#!/usr/bin/env python3
import argparse
import os
import sys
import pandas as pd
import numpy as np
from Bio import SeqIO, Align, pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess

def run_command(cmd, shell=False):
    """Run a shell command and check for errors."""
    try:
        subprocess.run(cmd, check=True, shell=shell, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
        print(f"Stderr: {e.stderr.decode()}")
        sys.exit(1)

def get_depth_per_position(bam_file, chrom, start, end, samtools_bin="samtools"):
    """
    Get depth for a specific region using samtools depth.
    Returns a dictionary {pos: depth}.
    """
    region = f"{chrom}:{start}-{end}"
    cmd = [samtools_bin, "depth", "-r", region, "-a", bam_file]
    
    # Run and capture output
    result = subprocess.run(cmd, stdout=subprocess.PIPE, text=True)
    depth_map = {}
    if result.returncode == 0:
        for line in result.stdout.splitlines():
            parts = line.split('\t')
            if len(parts) >= 3:
                 # samtools depth output: chrom, pos, depth
                 pos = int(parts[1])
                 d = int(parts[2])
                 depth_map[pos] = d
    return depth_map

def extract_cds_from_genbank(gb_file, gene_name):
    """
    Extracts the CDS sequence and its global coordinates from a GenBank file
    for a specific gene.
    """
    records = list(SeqIO.parse(gb_file, "genbank"))
    if not records:
        raise ValueError(f"No records found in {gb_file}")
    
    ref_record = records[0]
    
    for feature in ref_record.features:
        if feature.type == "CDS":
            # Check for gene name in qualifiers
            qualifiers = feature.qualifiers
            current_gene = qualifiers.get("gene", [""])[0] 
            
            # Sometimes gene name is in 'product' or 'note' if 'gene' is missing
            if gene_name.lower() in current_gene.lower():
                # Extract sequence
                seq = feature.extract(ref_record.seq)
                
                # Get global coordinates (1-based for user friendliness, but stored as logic)
                # Note: This simple logic assumes contiguous CDS. Splicing is handled by feature.extract() 
                # but coordinate mapping for coverage requires care.
                # For Influenza, CDS are usually contiguous per segment.
                start = feature.location.start + 1 # BioPython is 0-based
                end = feature.location.end
                
                return {
                    "seq": str(seq),
                    "translation": str(seq.translate()),
                    "start": start,
                    "end": end,
                    "location": feature.location,
                    "record": ref_record
                }
    
    raise ValueError(f"Gene {gene_name} not found in {gb_file}")

def align_and_map_cds(ref_seq_str, cons_seq_str):
    """
    Mapping helper that aligns two genomic sequences (Global)
    and returns a dictionary: ref_pos (1-based) -> cons_pos (1-based).
    Returns None if alignment fails.
    """
    try:
        from Bio import pairwise2
    except ImportError:
        print("Biopython pairwise2 not found.")
        return None

    # Global alignment with gap penalties
    # match=2, mismatch=-1, open=-5, extend=-0.5
    alns = pairwise2.align.globalms(ref_seq_str, cons_seq_str, 2, -1, -5, -0.5)
    if not alns:
        return None
    
    best_aln = alns[0]
    seqA = best_aln.seqA # Ref with gaps
    seqB = best_aln.seqB # Cons with gaps
    
    mapping = {}
    ref_idx = 0 # 0-based counter for ref
    cons_idx = 0 # 0-based counter for cons
    
    for a, b in zip(seqA, seqB):
        is_ref_base = (a != "-")
        is_cons_base = (b != "-")
        
        if is_ref_base:
            ref_pos = ref_idx + 1
            if is_cons_base:
                mapping[ref_pos] = cons_idx + 1
            else:
                mapping[ref_pos] = None # Deletion in Cons
            ref_idx += 1
        
        if is_cons_base:
            cons_idx += 1
            
    return mapping

def align_and_map_positions(ref_cds_seq, cons_seq_full):
    """
    Align ref CDS to full consensus sequence and build a mapping from
    ref CDS nucleotide positions (1-based) -> consensus nucleotide positions (1-based).
    Uses local alignment with gaps.
    """
    # Local alignment: ref CDS is a substring of full consensus segment
    aln = pairwise2.align.localms(
        ref_cds_seq, cons_seq_full, 2.0, -1.0, -3.0, -1.0, one_alignment_only=True
    )
    if not aln:
        return None, None, None

    aligned_ref, aligned_cons, _, _, _ = aln[0]
    ref_pos = 0
    cons_pos = 0
    mapping = {}

    for r, c in zip(aligned_ref, aligned_cons):
        if r != "-":
            ref_pos += 1
        if c != "-":
            cons_pos += 1
        if r != "-":
            mapping[ref_pos] = cons_pos if c != "-" else None

    return mapping, aligned_ref, aligned_cons


def cds_positions_from_location(location):
    parts = location.parts if hasattr(location, "parts") else [location]
    positions = []
    is_negative = getattr(location, "strand", 1) == -1
    if is_negative:
        parts = parts[::-1]
    for part in parts:
        start = int(part.start) + 1
        end = int(part.end)
        if getattr(part, "strand", 1) == -1:
            rng = range(end, start - 1, -1)
        else:
            rng = range(start, end + 1)
        positions.extend(list(rng))
    return positions

def _safe_name(name):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name)

def analyze_sample(sample_name, virus_type, bam_file, consensus_fasta, config, resistance_list, samtools_bin="samtools", out_prefix=None):
    """
    Analyze one sample for multiple segments/genes.
    """
    all_coverage_data = []
    all_mutation_data = []

    # Read Consensus Fasta
    cons_dict = SeqIO.to_dict(SeqIO.parse(consensus_fasta, "fasta"))

    for segment, info in config.items():
        gene = info['gene']
        ref_gb = info['ref_gb']
        cons_acc = info['cons_acc'] # Accession of this segment in the consensus file
        
        print(f"Processing {sample_name} - {gene}...")

        # 1. Get Reference Info
        try:
            ref_info = extract_cds_from_genbank(ref_gb, gene)
        except Exception as e:
            print(f"Skipping {gene}: {e}")
            continue

        # 2. Get Consensus Sequence for this segment
        if cons_acc not in cons_dict:
            print(f"Consensus accession {cons_acc} not found in {consensus_fasta}")
            continue
            
        cons_seq_full = str(cons_dict[cons_acc].seq)
        
        # 3. Build CDS mapping from GenBank (spliced CDS-aware)
        cds_positions = cds_positions_from_location(ref_info['location'])
        if not cds_positions:
            print(f"Skipping {gene}: CDS mapping failed")
            continue

        # Get Reference Full Seq for Length & Alignment
        ref_seq_full = str(ref_info['record'].seq)

        # 4. Get Depth Information (full segment, then map per codon)
        # Use Reference Length because BAM is aligned to Reference
        depth_map = get_depth_per_position(
            bam_file, cons_acc, 1, len(ref_seq_full), samtools_bin=samtools_bin
        )

        # Ensure CDS length is multiple of 3
        trim_len = (len(cds_positions) // 3) * 3
        cds_positions = cds_positions[:trim_len]

        ref_protein = ref_info['translation']
        if ref_protein.endswith("*"):
            ref_protein = ref_protein[:-1]

        # Consensus CDS from genomic coordinates (using Alignment)
        # ref_seq_full is already defined above
        
        # Map Ref Coordinates -> Cons Coordinates
        # (This handles INDELs which shift coordinates)
        mapping = align_and_map_cds(ref_seq_full, cons_seq_full)
        
        if not mapping:
            print(f"  Alignment failed for {gene}, skipping.")
            continue
            
        cons_cds_chars = []
        for p in cds_positions:
            # p is 1-based index in Reference
            mapped_p = mapping.get(p)
            if mapped_p is not None and mapped_p <= len(cons_seq_full):
                cons_cds_chars.append(cons_seq_full[mapped_p - 1])
            else:
                cons_cds_chars.append("-") # Deletion
                
        cons_cds_seq = "".join(cons_cds_chars)
        
        cds_strand = getattr(ref_info["location"], "strand", 1)
        if cds_strand == -1:
            cons_cds_seq = str(Seq(cons_cds_seq).reverse_complement())

        # 5. Build Analysis Tables per Amino Acid
        aa_len = len(ref_protein)
        cds_positions = cds_positions[:aa_len * 3]
        ref_cds_trim = ref_info['seq'][:aa_len * 3]
        cons_cds_trim = cons_cds_seq[:aa_len * 3]
        
        for i in range(aa_len):
            aa_pos = i + 1
            ref_aa = ref_protein[i]

            # Build ref/cons codons
            ref_codon = ref_cds_trim[i * 3:(i * 3) + 3]
            cons_codon = cons_cds_seq[i * 3:(i * 3) + 3]
            
            # Use Reference Genomic Positions for Depth Lookup
            # (cds_positions contains the Ref coordinates)
            ref_genomic_positions = cds_positions[i * 3:(i * 3) + 3]

            if len(ref_codon) < 3:
                cons_aa = "X"
            else:
                try:
                    # cons_codon is already strand-corrected and extracted from aligned Cons
                    if "-" in cons_codon:
                         cons_aa = "-" # Frameshift/Deletion representation
                    else:
                         cons_aa = str(Seq(cons_codon).translate())
                except Exception:
                    cons_aa = "X"
            
            # Calculate average depth for this codon
            # Use Reference positions to query depth_map (which matches BAM/Ref)
            codon_depths = [depth_map.get(p, 0) for p in ref_genomic_positions if p is not None]
            avg_depth = np.mean(codon_depths) if codon_depths else 0
            min_depth = min(codon_depths) if codon_depths else 0
            
            # Check Resistance
            # Filter where gene matches AND pos matches
            res_rows = resistance_list[
                (resistance_list['gene'] == gene) & 
                (resistance_list['aa_pos'] == aa_pos)
            ]
            
            is_resistance_site = not res_rows.empty
            
            # Check if current Mutation matches a resistance mutation
            resistance_mut_found = False
            if is_resistance_site:
                for _, rrow in res_rows.iterrows():
                    target_mut = rrow.get('mutation')
                    if target_mut and target_mut == cons_aa:
                        resistance_mut_found = True
                        break

            # Determine Status
            status = "Ref"
            if min_depth < 100:
                status = "Low Coverage"
            elif cons_aa != ref_aa:
                 if cons_aa == "X" or cons_aa == "*":
                     status = "Stop/Ambig"
                 elif resistance_mut_found:
                     status = "Resistance Mutation"
                 else:
                     status = "Other Mutation"
            
            row = {
                "sample": sample_name,
                "virus": virus_type,
                "gene": gene,
                "aa_pos": aa_pos,
                "ref_aa": ref_aa,
                "cons_aa": cons_aa,
                "depth": avg_depth,
                "status": status,
                "is_resistance": is_resistance_site,
                "res_mut_found": resistance_mut_found
            }
            all_mutation_data.append(row)

        # Write aligned nucleotide/protein FASTA for manual inspection
        if out_prefix:
            sample_tag = _safe_name(sample_name)
            gene_tag = _safe_name(gene)
            nt_path = f"{out_prefix}_{sample_tag}_{gene_tag}_aligned_nt.fasta"
            aa_path = f"{out_prefix}_{sample_tag}_{gene_tag}_aligned_aa.fasta"

            with open(nt_path, "w") as nt_out:
                nt_out.write(f">{sample_tag}|{gene_tag}|ref_cds\n{ref_cds_trim}\n")
                nt_out.write(f">{sample_tag}|{gene_tag}|cons_cds\n{cons_cds_trim}\n")

            # Protein alignment by reference positions (X for unmapped/ambig)
            cons_aa_list = [
                row["cons_aa"]
                for row in all_mutation_data
                if row["gene"] == gene and row["sample"] == sample_name
            ]
            ref_aa_list = list(ref_protein)
            with open(aa_path, "w") as aa_out:
                aa_out.write(f">{sample_tag}|{gene_tag}|ref_protein\n{''.join(ref_aa_list)}\n")
                aa_out.write(f">{sample_tag}|{gene_tag}|cons_protein\n{''.join(cons_aa_list)}\n")

    return pd.DataFrame(all_mutation_data)

import json

def main():
    parser = argparse.ArgumentParser(description="Analyze Viral Resistance")
    parser.add_argument("--config", required=True, help="Path to JSON configuration file")
    parser.add_argument("--resistance-csv", required=True, help="Path to resistance list CSV")
    parser.add_argument("--out-csv", required=True, help="Output CSV path")
    
    args = parser.parse_args()
    
    # Load Configuration
    with open(args.config, 'r') as f:
        config = json.load(f)
        
import re

def parse_resistance_csv(csv_path):
    """
    Parses the messy resistance CSV into a clean structured DataFrame.
    Target columns: virus_type, gene, aa_pos, ref_aa, mut_aa
    """
    df = pd.read_csv(csv_path)

    def norm_header(x):
        return re.sub(r"[\s_]+", " ", str(x).strip()).lower()

    col_map = {norm_header(c): c for c in df.columns}

    def find_col(aliases, required=True):
        for alias in aliases:
            key = norm_header(alias)
            if key in col_map:
                return col_map[key]
        if required:
            raise KeyError(f"Missing required column. Looked for: {aliases}. Available: {list(df.columns)}")
        return None

    col_strain = find_col(["Strain"])
    col_segment = find_col(["Segment"])
    col_aa_change = find_col(["Amino Acid Change", "Amino acid change"])
    col_resistance = find_col(["Resistance Level", "Resistance level"], required=False)
    col_drug = find_col(
        ["Drug Target", "Drug target", "Drug target or importance of the mutation"],
        required=False
    )
    
    clean_rows = []
    
    for _, row in df.iterrows():
        strains = str(row[col_strain])
        segment_raw = str(row[col_segment])
        aa_change = str(row[col_aa_change])
        resistance_level = str(row[col_resistance]) if col_resistance else ""
        drug_target = str(row[col_drug]) if col_drug else ""

        # Skip synergy/combination entries
        if re.search(r"synerg", resistance_level, flags=re.IGNORECASE) or re.search(r"synerg", drug_target, flags=re.IGNORECASE):
            continue
        
        # 1. Parse Strain (Handle "&")
        target_strains = []
        if "A(H1N1)pdm09" in strains:
            target_strains.append("A(H1N1)pdm09")
        if "A(H3N2)" in strains:
            target_strains.append("A(H3N2)")
            
        # 2. Parse Gene from Segment
        # e.g. "NA (N1)" -> "NA", "M2 ..." -> "M2"
        # Just take the first part before space or parenthesis
        gene = segment_raw.split()[0]
        
        # 3. Parse Amino Acid Change (e.g. H275Y, I223I/V)
        # Regex: Ref, Position, Mut(s)
        match = re.match(r"([A-Z])(\d+)([A-Z](?:/[A-Z])*)", aa_change)
        if not match:
            # Maybe ambiguous or complex like "I223I/V" or just "Global"
            continue
            
        ref_aa, pos_str, mut_aa = match.groups()
        aa_pos = int(pos_str)
        mut_list = mut_aa.split("/")
        
        for v in target_strains:
            for mut in mut_list:
                clean_rows.append({
                    "virus": v,
                    "gene": gene,
                    "aa_pos": aa_pos,
                    "ref_aa": ref_aa,
                    "mutation": mut,
                    "original_entry": aa_change
                })
            
    return pd.DataFrame(clean_rows)

def main():
    parser = argparse.ArgumentParser(description="Analyze Viral Resistance")
    parser.add_argument("--config", required=True, help="Path to JSON configuration file")
    parser.add_argument("--resistance-csv", required=True, help="Path to resistance list CSV")
    parser.add_argument("--out-csv", required=True, help="Output CSV path")
    parser.add_argument("--samtools", default="samtools", help="Path to samtools executable")
    
    args = parser.parse_args()
    
    # Load Configuration
    with open(args.config, 'r') as f:
        config = json.load(f)
        
    # Load and Parse Resistance List
    if not os.path.exists(args.resistance_csv):
         print(f"Error: Resistance CSV not found at {args.resistance_csv}")
         sys.exit(1)
         
    res_df = parse_resistance_csv(args.resistance_csv)
    print(f"Loaded {len(res_df)} resistance sites from CSV.")
    
    all_results = []
    
    # Iterate through samples defined in config
    for sample in config['samples']:
        sample_name = sample['name']
        bam_file = sample['bam']
        consensus_fasta = sample['consensus']
        virus_type = sample['virus_type'] # e.g. "A(H3N2)"
        
        print(f"Analyzing Sample: {sample_name} ({virus_type})")
        
        # Filter resistance list for this virus
        sample_res_list = res_df[res_df['virus'] == virus_type]
        
        # Analyze each segment defined for this sample
        # We pass the segments dict directly
        df = analyze_sample(
            sample_name=sample_name,
            virus_type=virus_type,
            bam_file=bam_file,
            consensus_fasta=consensus_fasta,
            config=sample['segments'],
            resistance_list=sample_res_list,
            samtools_bin=args.samtools,
            out_prefix=args.out_csv.replace(".csv", "")
        )
        all_results.append(df)

    if all_results:
        final_df = pd.concat(all_results, ignore_index=True)
        final_df.to_csv(args.out_csv, index=False)
        print(f"Analysis complete. Results saved to {args.out_csv}")

        # Write resistance-site report (positions + expected mutations)
        res_map = (
            res_df.groupby(["virus", "gene", "aa_pos"])["mutation"]
            .apply(lambda x: ";".join(sorted(set(x))))
            .reset_index()
            .rename(columns={"mutation": "expected_mutations"})
        )
        resistance_report = final_df.merge(
            res_map, on=["virus", "gene", "aa_pos"], how="inner"
        )
        resistance_csv = args.out_csv.replace(".csv", "_resistance_sites.csv")
        resistance_report.to_csv(resistance_csv, index=False)
        print(f"Resistance-site report saved to {resistance_csv}")

        # Full amino-acid check table with expected mutations
        aa_checks = final_df.merge(
            res_map, on=["virus", "gene", "aa_pos"], how="left"
        )
        aa_checks_csv = args.out_csv.replace(".csv", "_aa_checks.csv")
        aa_checks.to_csv(aa_checks_csv, index=False)
        print(f"AA check table saved to {aa_checks_csv}")

        # Alignment identity summary per gene
        alignment_summary = (
            final_df.assign(match=lambda d: d["ref_aa"] == d["cons_aa"])
            .groupby(["sample", "virus", "gene"])
            .agg(
                aa_len=("aa_pos", "max"),
                match_pct=("match", "mean")
            )
            .reset_index()
        )
        alignment_summary["match_pct"] = (alignment_summary["match_pct"] * 100).round(2)
        align_csv = args.out_csv.replace(".csv", "_alignment_summary.csv")
        alignment_summary.to_csv(align_csv, index=False)
        print(f"Alignment summary saved to {align_csv}")
    else:
        print("No results generated.")

if __name__ == "__main__":
    main()
