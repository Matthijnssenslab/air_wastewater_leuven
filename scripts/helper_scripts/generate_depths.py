#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys
import time


def resolve_samtools():
    samtools = os.environ.get("SAMTOOLS", "")
    if samtools:
        return samtools
    return "samtools"


def resolve_bcftools():
    bcftools = os.environ.get("BCFTOOLS", "")
    if bcftools:
        return bcftools
    return "bcftools"


def run(cmd):
    result = subprocess.run(cmd, shell=False)
    if result.returncode != 0:
        sys.exit(result.returncode)


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def run_with_heartbeat(proc, label, interval=30):
    last = time.time()
    while proc.poll() is None:
        time.sleep(1)
        now = time.time()
        if now - last >= interval:
            log(f"[heartbeat] {label} still running...")
            last = now


def ensure_index(samtools, bam):
    bai = bam + ".bai"
    if not os.path.exists(bai) or os.path.getmtime(bai) < os.path.getmtime(bam):
        run([samtools, "index", bam])


def trim_softclips(in_bam, out_bam):
    if os.path.exists(out_bam) and os.path.getmtime(out_bam) >= os.path.getmtime(in_bam):
        log(f"[trim] Using existing trimmed BAM: {out_bam}")
        return
    log(f"[trim] Trimming soft clips: {in_bam} -> {out_bam}")
    try:
        import pysam
    except ImportError:
        raise SystemExit("pysam is required for --trim-softclips. Install it in your environment.")

    with pysam.AlignmentFile(in_bam, "rb") as bam_in, \
            pysam.AlignmentFile(out_bam, "wb", template=bam_in) as bam_out:
        for read in bam_in.fetch(until_eof=True):
            if read.is_unmapped or read.cigartuples is None:
                bam_out.write(read)
                continue
            cigar = read.cigartuples
            if not cigar:
                bam_out.write(read)
                continue

            start_clip = cigar[0][1] if cigar[0][0] == 4 else 0
            end_clip = cigar[-1][1] if cigar[-1][0] == 4 else 0

            if start_clip == 0 and end_clip == 0:
                bam_out.write(read)
                continue

            seq = read.query_sequence
            qual = read.query_qualities
            if seq is None:
                bam_out.write(read)
                continue

            seq_end = len(seq) - end_clip if end_clip else len(seq)
            new_seq = seq[start_clip:seq_end]
            new_qual = qual[start_clip:seq_end] if qual is not None else None

            new_cigar = cigar[:]
            if start_clip:
                new_cigar = new_cigar[1:]
            if end_clip and len(new_cigar) > 0:
                new_cigar = new_cigar[:-1]

            if not new_cigar:
                continue

            read.query_sequence = new_seq
            if new_qual is not None:
                read.query_qualities = new_qual
            read.cigartuples = new_cigar
            bam_out.write(read)

    pysam.index(out_bam)
    log(f"[trim] Indexed: {out_bam}.bai")


def write_depth(samtools, bam, out_txt):
    if os.path.exists(out_txt) and os.path.getmtime(out_txt) >= os.path.getmtime(bam):
        log(f"[depth] Using existing depth: {out_txt}")
        return
    ensure_index(samtools, bam)
    log(f"[depth] {bam} -> {out_txt}")
    with open(out_txt, "w") as out_handle:
        proc = subprocess.run(
            [samtools, "depth", "-a", "-q", "30", "-Q", "30", bam],
            stdout=out_handle
        )
    if proc.returncode != 0:
        sys.exit(proc.returncode)


def compute_af_from_ad(ad_field):
    if not ad_field or ad_field == ".":
        return None
    parts = [p for p in ad_field.split(",") if p.isdigit()]
    if len(parts) < 2:
        return None
    counts = [int(x) for x in parts]
    total = sum(counts)
    if total == 0:
        return None
    alt_max = max(counts[1:])
    return alt_max / total


def write_variants_tsv(bcftools, bam, ref_fa, out_tsv):
    if os.path.exists(out_tsv) and os.path.getmtime(out_tsv) >= os.path.getmtime(bam):
        log(f"[vcf-tsv] Using existing TSV: {out_tsv}")
        return
    log(f"[vcf-tsv] {bam} -> {out_tsv}")
    # Use FORMAT/AD to compute AF since INFO/AF may be missing.
    cmd = [
        bcftools, "mpileup",
        "-f", ref_fa,
        "-d", "0",
        "-q", "30",
        "-Q", "30",
        "--redo-BAQ",
        "--adjust-MQ", "50",
        "--ff", "0x904",
        "-Ou", bam
    ]
    cmd2 = [
        bcftools, "call",
        "-mv", "-Ou"
    ]
    cmd3 = [
        bcftools, "query",
        "-f", "%CHROM\t%POS\t%REF\t%ALT\t[%DP]\t[%AD]\n"
    ]
    p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(cmd2, stdin=p1.stdout, stdout=subprocess.PIPE)
    proc = subprocess.run(cmd3, stdin=p2.stdout, stdout=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        sys.exit(proc.returncode)
    with open(out_tsv, "w") as out_handle:
        for line in proc.stdout.splitlines():
            chrom, pos, ref, alt, dp, ad = line.split("\t")
            af = compute_af_from_ad(ad)
            if af is None:
                continue
            out_handle.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{dp}\t{af:.6f}\n")


def write_vcf(bcftools, bam, ref_fa, out_vcf):
    if os.path.exists(out_vcf) and os.path.getmtime(out_vcf) >= os.path.getmtime(bam):
        log(f"[vcf] Using existing VCF: {out_vcf}")
        return
    log(f"[vcf] {bam} -> {out_vcf}")
    cmd = [
        bcftools, "mpileup",
        "-f", ref_fa,
        "-d", "0",
        "-q", "30",
        "-Q", "30",
        "--redo-BAQ",
        "--adjust-MQ", "50",
        "--ff", "0x904",
        "-Ou", bam
    ]
    cmd2 = [
        bcftools, "call",
        "-mv", "-Oz",
        "-o", out_vcf
    ]
    p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(cmd2, stdin=p1.stdout)
    run_with_heartbeat(p2, f"bcftools call ({os.path.basename(out_vcf)})")
    if p2.returncode != 0:
        sys.exit(p2.returncode)
    idx = subprocess.run([bcftools, "index", "-t", out_vcf])
    if idx.returncode != 0:
        sys.exit(idx.returncode)


def write_mask_from_depth(depth_txt, min_dp, out_bed):
    log(f"[mask-dp] {depth_txt} -> {out_bed}")
    with open(depth_txt, "r") as fh, open(out_bed, "w") as out_handle:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom, pos, depth = parts[0], int(parts[1]), int(parts[2])
            if depth < min_dp:
                start = pos - 1
                out_handle.write(f"{chrom}\t{start}\t{pos}\n")


def write_mask_from_vcf(bcftools, vcf_path, min_dp, min_af, out_bed):
    log(f"[mask-af] {vcf_path} -> {out_bed}")
    query = [
        bcftools, "query",
        "-f", "%CHROM\t%POS\t[%DP]\t[%AD]\n",
        vcf_path
    ]
    proc = subprocess.run(query, stdout=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        sys.exit(proc.returncode)
    with open(out_bed, "w") as out_handle:
        for line in proc.stdout.splitlines():
            chrom, pos, dp, ad = line.split("\t")
            dp_val = float(dp) if dp != "." else 0
            af = compute_af_from_ad(ad)
            if af is None:
                continue
            if dp_val >= min_dp and af < min_af:
                pos = int(pos)
                start = pos - 1
                out_handle.write(f"{chrom}\t{start}\t{pos}\n")


def merge_masks(mask_a, mask_b, out_bed):
    log(f"[mask-merge] {mask_a}, {mask_b} -> {out_bed}")
    seen = set()
    with open(out_bed, "w") as out_handle:
        for path in (mask_a, mask_b):
            if not os.path.exists(path):
                continue
            with open(path, "r") as fh:
                for line in fh:
                    if line in seen:
                        continue
                    seen.add(line)
                    out_handle.write(line)


def write_consensus(bcftools, ref_fa, vcf_path, mask_bed, out_fa):
    if os.path.exists(out_fa) and os.path.getmtime(out_fa) >= os.path.getmtime(vcf_path):
        log(f"[consensus] Using existing consensus: {out_fa}")
        return
    log(f"[consensus] {vcf_path} -> {out_fa}")
    cmd = [
        bcftools, "consensus",
        "-f", ref_fa,
        "-m", mask_bed,
        vcf_path
    ]
    with open(out_fa, "w") as out_handle:
        proc = subprocess.run(cmd, stdout=out_handle)
    if proc.returncode != 0:
        sys.exit(proc.returncode)


def main():
    parser = argparse.ArgumentParser(description="Generate samtools depth files and consensus.")
    parser.add_argument("--bam", action="append", required=True, help="BAM file path (repeatable).")
    parser.add_argument("--out", action="append", required=True, help="Output depth txt (repeatable).")
    parser.add_argument("--ref", help="Reference FASTA for bcftools mpileup/consensus.")
    parser.add_argument("--vcf", action="append", help="Output variant TSV (repeatable).")
    parser.add_argument("--vcf-raw", action="append", help="Output VCF (repeatable).")
    parser.add_argument("--consensus", action="append", help="Output consensus FASTA (repeatable).")
    parser.add_argument("--samtools", help="Path to samtools binary (optional).")
    parser.add_argument("--bcftools", help="Path to bcftools binary (optional).")
    parser.add_argument("--trim-softclips", action="store_true", help="Trim soft-clipped bases before depth/consensus.")
    parser.add_argument("--trim-bam", action="append", help="Output trimmed BAM (repeatable).")
    parser.add_argument("--min-dp", type=int, default=100, help="Minimum depth for consensus (default: 100).")
    parser.add_argument("--min-af", type=float, default=0.9, help="Minimum AF for consensus (default: 0.9).")
    args = parser.parse_args()
    if len(args.bam) != len(args.out):
        parser.error("--bam and --out must be the same length.")
    if args.vcf and len(args.vcf) != len(args.bam):
        parser.error("--vcf must be the same length as --bam when provided.")
    if args.vcf_raw and len(args.vcf_raw) != len(args.bam):
        parser.error("--vcf-raw must be the same length as --bam when provided.")
    if args.consensus and len(args.consensus) != len(args.bam):
        parser.error("--consensus must be the same length as --bam when provided.")

    samtools = args.samtools if args.samtools else resolve_samtools()
    bcftools = args.bcftools if args.bcftools else resolve_bcftools()

    bam_list = list(args.bam)
    if args.trim_softclips:
        if not args.trim_bam or len(args.trim_bam) != len(args.bam):
            parser.error("--trim-bam must be provided once per --bam when using --trim-softclips.")
        trimmed = []
        for in_bam, out_bam in zip(args.bam, args.trim_bam):
            trim_softclips(in_bam, out_bam)
            trimmed.append(out_bam)
        bam_list = trimmed

    for bam, out_txt in zip(bam_list, args.out):
        if not os.path.exists(bam):
            raise FileNotFoundError(bam)
        write_depth(samtools, bam, out_txt)

    if args.ref and args.vcf:
        if not os.path.exists(args.ref):
            raise FileNotFoundError(args.ref)
        for bam, out_tsv in zip(bam_list, args.vcf):
            write_variants_tsv(bcftools, bam, args.ref, out_tsv)

    if args.ref and args.vcf_raw and args.consensus:
        if not os.path.exists(args.ref):
            raise FileNotFoundError(args.ref)
        for bam, out_txt, out_vcf, out_fa in zip(bam_list, args.out, args.vcf_raw, args.consensus):
            write_vcf(bcftools, bam, args.ref, out_vcf)
            mask_depth = out_txt + ".lowdp.bed"
            mask_af = out_vcf + ".lowaf.bed"
            mask_all = out_vcf + ".mask.bed"
            write_mask_from_depth(out_txt, args.min_dp, mask_depth)
            write_mask_from_vcf(bcftools, out_vcf, args.min_dp, args.min_af, mask_af)
            merge_masks(mask_depth, mask_af, mask_all)
            write_consensus(bcftools, args.ref, out_vcf, mask_all, out_fa)


if __name__ == "__main__":
    main()
