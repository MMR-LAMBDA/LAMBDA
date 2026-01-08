#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
UMI-Tagged Amplicon Sequencing Analysis Pipeline for the att-IR Platform


Version: 1.1

DESCRIPTION:
This pipeline is designed to analyze amplicon sequencing data from an att-site mediated
Intracellular Recombination (att-IR) platform. The platform is used to physically
link a library of enzyme variants (e.g., MutS for DNA Mismatch Repair) to their
functional activity on a corresponding library of substrates (e.g., DNA mismatches
within a ccdB gene).

In the att-IR system, a ΦC31 integrase-mediated recombination event in a single E. coli
cell fuses the MutS-encoding plasmid with the repaired ccdB substrate plasmid. This
event reconstitutes a functional β-lactamase gene, enabling selection with ampicillin.
The resulting fused plasmid serves as a physical record linking a specific MutS
variant to the substrate it has successfully repaired.

This bioinformatics pipeline analyzes UMI-tagged amplicons that span both the MutS
variant region and the repaired site in the ccdB gene. Its ultimate goal is to
quantify which MutS variants are enriched for the repair of each specific type of
DNA mismatch.

PIPELINE LOGIC:
The analysis is divided into six core steps:
1.  Adapter Trimming: Preprocessing of raw FASTQ files.
2.  UMI Processing & SE Filtering: Establishes a high-confidence read list using a
    single-end (Read 1) workflow, UMI-based consensus collapsing, and majority-rule
    filtering on key MutS variant codons.
3.  Paired-End Alignment: Generates a complete PE alignment database.
4.  Information Merge & Filter: Links the high-confidence SE reads (with UMIs) to
    their corresponding full PE alignments.
5.  DNA Mutation Analysis (Substrate De-multiplexing): Segregates reads into
    different BAM files based on the identified repaired substrate in the ccdB gene.
6.  Amino Acid Mutation Calling (Quantification): Quantifies the frequency of each
    MutS variant for each substrate-specific BAM file.

USAGE:
python att-IR_pipeline.py --samples ./input ./selection --ref ./reference/integrated.fa --threads 15
"""

import argparse
import logging
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from collections import Counter, defaultdict
from itertools import groupby
from multiprocessing import Pool, cpu_count

# Configure logging for standardized output
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)


def _se_filter_by_majority(input_bam: Path, output_bam: Path, chromosome: str, positions: str):
    """
    Filters a BAM file based on the majority allele at key positions within each UMI cluster.
    """
    import pysam

    try:
        target_positions_0based = sorted([int(p.strip()) - 1 for p in positions.split(',')])
        target_positions_set = set(target_positions_0based)
    except ValueError:
        logging.error("FATAL: --positions argument is not formatted correctly.")
        sys.exit(1)

    logging.info(f"Starting majority-rule filtering on: {input_bam.name}")
    start_time = time.time()

    def get_all_target_alleles(read: pysam.AlignedSegment, t_pos_set: set):
        if read.reference_end is None or not t_pos_set: return {}
        first_pos, last_pos = min(t_pos_set), max(t_pos_set)
        if read.reference_end <= first_pos or read.reference_start > last_pos: return {}
        alleles = {}
        try:
            aligned_pairs = read.get_aligned_pairs(with_seq=True)
        except Exception:
            return {}
        for _, ref_pos, query_base in aligned_pairs:
            if ref_pos in t_pos_set:
                alleles[ref_pos] = '-' if query_base is None else query_base.upper()
        return alleles

    with pysam.AlignmentFile(str(input_bam), "rb") as bam_in, \
         pysam.AlignmentFile(str(output_bam), "wb", template=bam_in) as bam_out:
        keyfunc = lambda read: read.get_tag('MI')
        pcc, wrc, trp = 0, 0, 0
        for _, reads_in_cluster_iter in groupby(bam_in, key=keyfunc):
            cluster_reads, allele_matrix = [], []
            for read in reads_in_cluster_iter:
                trp += 1
                if read.reference_name == chromosome:
                    alleles = get_all_target_alleles(read, target_positions_set)
                    cluster_reads.append(read)
                    allele_matrix.append(alleles)

            if not cluster_reads: continue

            n_reads = len(cluster_reads)
            survivor_mask = [True] * n_reads
            for pos in target_positions_0based:
                counts = Counter(allele_matrix[i].get(pos) for i in range(n_reads) if survivor_mask[i] and allele_matrix[i].get(pos))
                if not counts: continue
                major_allele = counts.most_common(1)[0][0]
                for i in range(n_reads):
                    if survivor_mask[i] and allele_matrix[i].get(pos) is not None and allele_matrix[i].get(pos) != major_allele:
                        survivor_mask[i] = False

            for i in range(n_reads):
                if survivor_mask[i]:
                    bam_out.write(cluster_reads[i])
                    wrc += 1
            pcc += 1
            if pcc > 0 and pcc % 50000 == 0:
                logging.info(f"  ...Clusters processed: {pcc}, Reads scanned: {trp}, Time elapsed: {time.time() - start_time:.2f}s")

    logging.info(f"Filtering finished. Total clusters: {pcc}. Reads written: {wrc}. Total time: {time.time() - start_time:.2f}s.")


def _merge_tags(se_bam_path: Path, pe_bam_path: Path, output_bam_path: Path, tags_to_transfer: list):
    """
    Transfers specified tags from a single-end BAM to a paired-end BAM by matching read names.
    """
    import pysam
    logging.info("Starting tag transfer from SE to PE alignments...")
    reads_written, pe_groups_processed = 0, 0
    get_base_qname = lambda read: read.query_name.split('/')[0]

    with pysam.AlignmentFile(str(se_bam_path), "rb", check_sq=False) as se_bam, \
         pysam.AlignmentFile(str(pe_bam_path), "rb") as pe_bam, \
         pysam.AlignmentFile(str(output_bam_path), "wb", header=pe_bam.header) as out_bam:

        se_groups = groupby(se_bam, key=get_base_qname)
        pe_groups = groupby(pe_bam, key=get_base_qname)

        try:
            se_qname, se_reads = next(se_groups)
            pe_qname, pe_reads = next(pe_groups)
        except StopIteration:
            logging.warning("One of the input BAM files is empty.")
            return

        while True:
            if se_qname == pe_qname:
                try:
                    first_se_read = next(se_reads)
                    tags_data = {tag: first_se_read.get_tag(tag, with_value_type=True) for tag in tags_to_transfer if first_se_read.has_tag(tag)}
                    if tags_data:
                        for pe_read in pe_reads:
                            for tag, (val, vtype) in tags_data.items():
                                pe_read.set_tag(tag, val, value_type=vtype)
                            out_bam.write(pe_read)
                            reads_written += 1
                except StopIteration: pass
                pe_groups_processed += 1
                try:
                    se_qname, se_reads = next(se_groups)
                    pe_qname, pe_reads = next(pe_groups)
                except StopIteration: break
            elif pe_qname < se_qname:
                pe_groups_processed += 1
                try:
                    pe_qname, pe_reads = next(pe_groups)
                except StopIteration: break
            else:
                try:
                    se_qname, se_reads = next(se_groups)
                except StopIteration: break
            if pe_groups_processed > 0 and pe_groups_processed % 1000000 == 0:
                logging.info(f"  ...Processed {pe_groups_processed // 1000000}M PE read groups")
    logging.info(f"Tag transfer complete. Total reads written: {reads_written}.")


# --- Functions for Step 5: DNA Mutation Analysis ---
MUTATION_DEFINITIONS = {'mm_G_419', 'del_A_422', 'ins_C_434', 'ins_A_435', 'ins_T_435', 'ins_G_435'}

def _get_mutations_from_read(alignment):
    import pysam
    mutations = set()
    if not alignment.query_sequence: return mutations
    q_seq, last_ref_pos, ins_buf = alignment.query_sequence, -1, []
    try:
        pairs = alignment.get_aligned_pairs(with_seq=True)
    except Exception:
        return mutations
    
    for q_pos, r_pos, r_base in pairs:
        if r_pos is not None and ins_buf:
            ins, ins_pos = "".join(ins_buf).upper(), last_ref_pos + 1
            if ins_pos == 434 and ins == 'C': mutations.add('ins_C_434')
            elif ins_pos == 435:
                if ins == 'A': mutations.add('ins_A_435')
                elif ins == 'T': mutations.add('ins_T_435')
                elif ins == 'G': mutations.add('ins_G_435')
            ins_buf = []
        
        if q_pos is None:
            if r_pos + 1 == 422: mutations.add('del_A_422')
        elif r_pos is None:
            ins_buf.append(q_seq[q_pos])
        else:
            if r_pos + 1 == 419 and q_seq[q_pos].upper() == 'G' and (r_base is None or q_seq[q_pos].upper() != r_base.upper()):
                mutations.add('mm_G_419')
            last_ref_pos = r_pos
    
    if ins_buf:
        ins, ins_pos = "".join(ins_buf).upper(), last_ref_pos + 1
        if ins_pos == 434 and ins == 'C': mutations.add('ins_C_434')
        elif ins_pos == 435:
            if ins == 'A': mutations.add('ins_A_435')
            elif ins == 'T': mutations.add('ins_T_435')
            elif ins == 'G': mutations.add('ins_G_435')
            
    return mutations

def _process_cluster_worker(cluster_data):
    counts, reps = Counter(), {}
    for r1_str, r2_str, mutations in cluster_data:
        for mut in mutations:
            if mut in MUTATION_DEFINITIONS:
                counts[mut] += 1
                if mut not in reps:
                    reps[mut] = (r1_str, r2_str)
    if not counts:
        return None
    major_mut = counts.most_common(1)[0][0]
    return (major_mut, *reps.get(major_mut)) if reps.get(major_mut) else None

def _create_job_iterator(bam_iter, min_pairs):
    import pysam
    for _, reads_iter in groupby(bam_iter, key=lambda r: r.get_tag('MI')):
        pairs = defaultdict(list)
        for r in reads_iter:
            pairs[r.query_name].append(r)
        
        if len(pairs) < min_pairs:
            continue
        
        job_data = []
        for _, members in pairs.items():
            if len(members) != 2:
                continue
            r1 = next((r for r in members if r.is_read1), None)
            r2 = next((r for r in members if r.is_read2), None)
            if r1 and r2:
                mutations = _get_mutations_from_read(r2)
                job_data.append((r1.to_string(), r2.to_string(), mutations))
        if job_data:
            yield job_data

def _analyze_clusters_by_majority(input_bam: Path, output_dir: Path, threads: int, min_pairs: int):
    """
    Analyzes UMI clusters to identify the majority DNA mutation (repaired substrate)
    and splits the BAM file accordingly.
    """
    import pysam
    from tqdm import tqdm

    output_dir.mkdir(exist_ok=True)
    
    with pysam.AlignmentFile(str(input_bam), "rb") as bam_in:
        mut_map = {m: '_'.join(m.split('_')[:2]) for m in MUTATION_DEFINITIONS}
        s_names = set(mut_map.values())
        out_bams = {n: pysam.AlignmentFile(str(output_dir / f"{n}.bam"), "wb", template=bam_in) for n in s_names}
        
        job_iter = _create_job_iterator(bam_in, min_pairs)
        logging.info("De-multiplexing BAM by repaired substrate type...")
        
        with Pool(processes=threads) as pool:
            results = pool.imap_unordered(_process_cluster_worker, job_iter, chunksize=100)
            for res in tqdm(results, desc="  Analyzing UMI Clusters", unit=" cluster", ncols=100):
                if res:
                    full_mut, r1_str, r2_str = res
                    s_name = mut_map.get(full_mut)
                    if out_bams.get(s_name):
                        out_bams[s_name].write(pysam.AlignedSegment.fromstring(r1_str, out_bams[s_name].header))
                        out_bams[s_name].write(pysam.AlignedSegment.fromstring(r2_str, out_bams[s_name].header))

        for f in out_bams.values():
            f.close()
            
    logging.info("Sorting and indexing output BAM files...")
    for n in s_names:
        path = output_dir / f"{n}.bam"
        if path.exists() and path.stat().st_size > 0:
            try:
                pysam.sort("-o", str(path), str(path))
                pysam.index(str(path))
            except pysam.utils.SamtoolsError as e:
                logging.warning(f"Could not sort/index {path}: {e}")
    logging.info(f"Substrate de-multiplexing complete. Output is in '{output_dir}'.")

# --- Functions for Step 6: Amino Acid Mutation Calling ---
CODON_TABLE_3_LETTER = {'ATA':'Ile','ATC':'Ile','ATT':'Ile','ATG':'Met','ACA':'Thr','ACC':'Thr','ACG':'Thr','ACT':'Thr','AAC':'Asn','AAT':'Asn','AAA':'Lys','AAG':'Lys','AGC':'Ser','AGT':'Ser','AGA':'Arg','AGG':'Arg','CTA':'Leu','CTC':'Leu','CTG':'Leu','CTT':'Leu','CCA':'Pro','CCC':'Pro','CCG':'Pro','CCT':'Pro','CAC':'His','CAT':'His','CAA':'Gln','CAG':'Gln','CGA':'Arg','CGC':'Arg','CGG':'Arg','CGT':'Arg','GTA':'Val','GTC':'Val','GTG':'Val','GTT':'Val','GCA':'Ala','GCC':'Ala','GCG':'Ala','GCT':'Ala','GAC':'Asp','GAT':'Asp','GAA':'Glu','GAG':'Glu','GGA':'Gly','GGC':'Gly','GGG':'Gly','GGT':'Gly','TCA':'Ser','TCC':'Ser','TCG':'Ser','TCT':'Ser','TTC':'Phe','TTT':'Phe','TTA':'Leu','TTG':'Leu','TAC':'Tyr','TAT':'Tyr','TAA':'Ter','TAG':'Ter','TGC':'Cys','TGT':'Cys','TGA':'Ter','TGG':'Trp'}

def _reverse_complement(seq: str) -> str:
    return seq.upper().translate(str.maketrans('ATCGN', 'TAGCN'))[::-1]

def _get_sequence_for_ref_interval(read, start_1based: int, end_1based: int):
    ref_start_0based, ref_end_0based = start_1based - 1, end_1based - 1
    expected_len = (end_1based - start_1based) + 1
    try:
        ref_to_query_map, deletions = {}, set()
        for q_pos, r_pos in read.get_aligned_pairs():
            if r_pos is not None and ref_start_0based <= r_pos <= ref_end_0based:
                if q_pos is None:
                    deletions.add(r_pos)
                else:
                    ref_to_query_map[r_pos] = read.query_sequence[q_pos]
        if any(pos in deletions for pos in range(ref_start_0based, ref_end_0based + 1)):
            return None
        extracted_seq = "".join([ref_to_query_map.get(i, "") for i in range(ref_start_0based, ref_end_0based + 1)])
        return extracted_seq if len(extracted_seq) == expected_len else None
    except Exception:
        return None

def _extract_and_translate_codons(in_bam: Path, out_txt: Path, mut_count_file: Path):
    """
    Extracts, translates, and quantifies MutS variant codons from a BAM file.
    """
    import pysam
    from tqdm import tqdm
    
    logging.info(f"Reading from BAM file: {in_bam.name}")
    reads_processed, results_written, mutation_counter = 0, 0, Counter()
    
    with pysam.AlignmentFile(str(in_bam), "rb") as bamfile, open(out_txt, "w") as outfile:
        outfile.write("read_id\tmutation_info\n")
        
        for read in tqdm(bamfile, desc="    Quantifying MutS Variants", unit=" reads", leave=False, ncols=100):
            reads_processed += 1
            if read.is_unmapped or read.is_secondary or read.is_supplementary or not read.is_read1:
                continue
            
            seq_44_46 = _get_sequence_for_ref_interval(read, 44, 46)
            seq_50_52 = _get_sequence_for_ref_interval(read, 50, 52)
            
            if seq_44_46 is None or seq_50_52 is None:
                continue
                
            aa_44_46 = CODON_TABLE_3_LETTER.get(_reverse_complement(seq_44_46), "Unk") # Codon for aa 38
            aa_50_52 = CODON_TABLE_3_LETTER.get(_reverse_complement(seq_50_52), "Unk") # Codon for aa 36
            
            parts = []
            if aa_50_52 != 'Phe': parts.append(f"p.Phe36{aa_50_52}")
            if aa_44_46 != 'Glu': parts.append(f"p.Glu38{aa_44_46}")
            
            result_string = ", ".join(parts) if parts else "_wt"
            outfile.write(f"{read.query_name}\t{result_string}\n")
            mutation_counter[result_string] += 1
            results_written += 1
            
    logging.info(f"  Writing mutation counts to {mut_count_file.name}...")
    with open(mut_count_file, "w") as count_file:
        count_file.write("mutation_info\tcount\n")
        for mutation, count in mutation_counter.most_common():
            count_file.write(f"{mutation}\t{count}\n")
            
    logging.info(f"    - Total reads scanned: {reads_processed}")
    logging.info(f"    - Reads passing filters & quantified: {results_written}")

# =========================================================================================
# --- Main Pipeline Class ---
# =========================================================================================

class AttIRPipeline:
    def __init__(self, config):
        self.config = config
        self.sample_name = ""
        self.sample_dir = None
        self.main_output_dir = None
        self.dirs = {}

    def _run_cmd(self, cmd, **kwargs):
        """Helper to run an external command and log its execution."""
        logging.info(f"Executing command: {' '.join(map(str, cmd))}")
        try:
            subprocess.run(cmd, check=True, **kwargs)
        except subprocess.CalledProcessError as e:
            logging.error(f"Command failed with exit code {e.returncode}")
            if e.stderr:
                logging.error(f"Stderr: {e.stderr}")
            raise e

    def _setup_directories(self, sample_dir: Path):
        """Create the directory structure for a given sample."""
        self.sample_dir = sample_dir
        self.sample_name = sample_dir.name
        self.main_output_dir = self.sample_dir / "pipeline_results"
        
        if self.main_output_dir.exists():
            logging.warning(f"Output directory {self.main_output_dir} already exists. Removing it.")
            shutil.rmtree(self.main_output_dir)
        
        self.dirs = {
            "step1": self.main_output_dir / "01_trimmed_reads",
            "step2": self.main_output_dir / "02_umi_processed_se",
            "step3": self.main_output_dir / "03_pe_alignment",
            "step4": self.main_output_dir / "04_merged_bam",
            "step5": self.main_output_dir / "05_split_by_dna_mutation",
            "step6": self.main_output_dir / "06_final_aa_counts",
        }
        for d in self.dirs.values():
            d.mkdir(parents=True, exist_ok=True)

    def _run_step1_trimming(self, raw_r1, raw_r2):
        logging.info(f"--- [Step 1/6] Adapter Trimming & Quality Filtering for {self.sample_name} ---")
        trimmed_r1 = self.dirs['step1'] / f"{self.sample_name}.trimmed_1.fq.gz"
        trimmed_r2 = self.dirs['step1'] / f"{self.sample_name}.trimmed_2.fq.gz"
        report_path = self.dirs['step1'] / f"{self.sample_name}.cutadapt_report.txt"
        
        cmd = [
            'cutadapt', '--times', '2',
            '-a', self.config['ADAPTER_R1_3PRIME'],
            '-G', self.config['ADAPTER_R2_5PRIME'],
            '-A', self.config['ADAPTER_R2_3PRIME'],
            '-m', f"{self.config['MIN_LEN_R1']}:{self.config['MIN_LEN_R2']}",
            '-M', f"{self.config['MAX_LEN_R1']}:{self.config['MAX_LEN_R2']}",
            '-q', f"{self.config['QUALITY_CUTOFF']},{self.config['QUALITY_CUTOFF']}",
            '--cores', str(self.config['THREADS']),
            '-o', str(trimmed_r1),
            '-p', str(trimmed_r2),
            str(raw_r1), str(raw_r2)
        ]
        
        with open(report_path, 'w') as f:
            self._run_cmd(cmd, stdout=f, stderr=subprocess.STDOUT)
            
        logging.info("✅ Step 1 complete.")
        return trimmed_r1, trimmed_r2

    def _run_step2_umi_processing(self, trimmed_r1):
        logging.info(f"--- [Step 2/6] UMI Processing & High-Confidence SE Filtering for {self.sample_name} ---")
        
        s2_dir = self.dirs['step2']
        se_aligned_bam = s2_dir / f"{self.sample_name}.se.aligned.bam"
        se_collapsed_bam = s2_dir / f"{self.sample_name}.se.collapsed.bam"
        se_collapsed_sorted_mi_bam = s2_dir / f"{self.sample_name}.se.collapsed.mi_sorted.bam"
        se_filtered_bam = s2_dir / f"{self.sample_name}.se.filtered.bam"
        se_final_namesorted_bam = s2_dir / f"{self.sample_name}.se.filtered.namesorted.bam"
        cluster_dist_tsv = s2_dir / "cluster_size_dist.tsv"
        umi_log_file = s2_dir / "umi_extract.log"

        # (2a) Piped command: umi_tools extract -> bowtie2 -> samtools sort
        logging.info("  - (2a) Extracting UMI, Aligning R1, and Sorting via stream...")
        
        cmd_umi = [
            'umi_tools', 'extract',
            '--stdin', str(trimmed_r1),
            '--extract-method=regex',
            f"--bc-pattern={self.config['UMI_PATTERN']}",
            '--log', str(umi_log_file),
        ]
        cmd_bt2 = [
            'bowtie2', '-p', str(self.config['THREADS']), '--np', '0', 
            '--very-sensitive', '-x', self.config['REF_BT2_INDEX'], '-U', '-'
        ]
        cmd_sam_sort = [
            'samtools', 'sort', '--threads', str(self.config['THREADS']), 
            '-m', self.config['SORT_MEM_PER_THREAD'], 
            '-T', str(s2_dir / "sort_tmp"), 
            '-o', str(se_aligned_bam), '-'
        ]
        
        try:
            p_umi = subprocess.Popen(cmd_umi, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            p_bt2 = subprocess.Popen(cmd_bt2, stdin=p_umi.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if p_umi.stdout:
                p_umi.stdout.close()
            
            p_sam = subprocess.Popen(cmd_sam_sort, stdin=p_bt2.stdout, stderr=subprocess.PIPE)
            if p_bt2.stdout:
                p_bt2.stdout.close()

            _, p_sam_err = p_sam.communicate()
            _, p_bt2_err = p_bt2.communicate()
            _, p_umi_err = p_umi.communicate()

            if p_umi.returncode != 0:
                logging.error(f"umi_tools failed with exit code {p_umi.returncode}:\n{p_umi_err.decode()}")
                raise subprocess.CalledProcessError(p_umi.returncode, cmd_umi)
            if p_bt2.returncode != 0:
                logging.error(f"bowtie2 failed with exit code {p_bt2.returncode}:\n{p_bt2_err.decode()}")
                raise subprocess.CalledProcessError(p_bt2.returncode, cmd_bt2)
            if p_sam.returncode != 0:
                logging.error(f"samtools sort failed with exit code {p_sam.returncode}:\n{p_sam_err.decode()}")
                raise subprocess.CalledProcessError(p_sam.returncode, cmd_sam_sort)
            
            logging.info("  - (2a) Streaming alignment complete.")
            logging.info(f"  - Bowtie2 alignment summary:\n---\n{p_bt2_err.decode().strip()}\n---")

        except Exception as e:
            logging.error("An error occurred during the streaming alignment process.")
            raise e

        # (2b) Collapse reads
        logging.info("  - (2b) Collapsing reads by UMI...")
        self._run_cmd(['umicollapse', 'bam', '-i', str(se_aligned_bam), '-o', str(se_collapsed_bam), '--algo', 'adj', '-k', '4', '--tag'])

        # (2c) Sort by MI tag
        logging.info("  - (2c) Sorting collapsed BAM by MI tag...")
        self._run_cmd(['samtools', 'sort', '-t', 'MI', '-@', str(self.config['THREADS']), str(se_collapsed_bam), '-o', str(se_collapsed_sorted_mi_bam)])

        # (2d) Cluster size distribution
        logging.info("  - (2d) Generating UMI cluster size distribution...")
        dist_cmd = f"samtools view {se_collapsed_sorted_mi_bam} | awk -F '\\t' '{{for(i=1;i<=NF;i++){{if($i~/^MI:Z:/){{mi=substr($i,6); count[mi]++}}}}}} END{{for(mi in count) print mi, count[mi]}}' | awk '{{print $2}}' | sort -n | uniq -c | awk 'BEGIN{{OFS=\"\\t\"}}{{print $2, $1}}' > {cluster_dist_tsv}"
        subprocess.run(dist_cmd, shell=True, check=True, executable='/bin/bash')

        # (2e) Majority rule filtering
        logging.info("  - (2e) Filtering clusters by majority rule at MutS variant codons...")
        _se_filter_by_majority(
            se_collapsed_sorted_mi_bam, se_filtered_bam,
            self.config['SE_FILTER_CHROMOSOME'], self.config['SE_FILTER_POSITIONS']
        )
        
        # (2f) Sort by read name
        logging.info("  - (2f) Sorting filtered BAM by read name...")
        self._run_cmd(['samtools', 'sort', '-n', '-@', str(self.config['THREADS']), str(se_filtered_bam), '-o', str(se_final_namesorted_bam), '-T', str(s2_dir / "sort_name_tmp")])

        logging.info("✅ Step 2 complete.")
        return se_final_namesorted_bam

    def _run_step3_pe_alignment(self, trimmed_r1, trimmed_r2):
        logging.info(f"--- [Step 3/6] Paired-End Alignment for {self.sample_name} ---")
        pe_namesorted_bam = self.dirs['step3'] / f"{self.sample_name}.pe.namesorted.bam"
        temp_dir = self.dirs['step3'] / "temp"
        temp_dir.mkdir()
        
        r1_umi_trimmed = temp_dir / "r1.fq.gz"
        r2_umi_trimmed = temp_dir / "r2.fq.gz"

        logging.info("  - (3a) Trimming UMI from R1 for PE alignment...")
        self._run_cmd(['cutadapt', '-u', str(self.config['TRIM_BASES_R1_FOR_PE']), '-o', str(r1_umi_trimmed), '-p', str(r2_umi_trimmed), str(trimmed_r1), str(trimmed_r2)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        logging.info("  - (3b) Aligning PE reads and sorting by name...")
        cmd_bt2 = ['bowtie2', '-p', str(self.config['THREADS']), '--np', '0', '--very-sensitive', '-x', self.config['REF_BT2_INDEX'], '-1', str(r1_umi_trimmed), '-2', str(r2_umi_trimmed)]
        cmd_sam_sort = ['samtools', 'sort', '-n', '--threads', str(self.config['THREADS']), '-m', self.config['SORT_MEM_PER_THREAD'], '-T', str(temp_dir / "sort_tmp"), '-o', str(pe_namesorted_bam), '-']
        
        p_bt2 = subprocess.Popen(cmd_bt2, stdout=subprocess.PIPE)
        p_sam = subprocess.Popen(cmd_sam_sort, stdin=p_bt2.stdout)
        if p_bt2.stdout:
            p_bt2.stdout.close()
        
        ret_sam = p_sam.wait()
        ret_bt2 = p_bt2.wait()
        
        if ret_bt2 != 0: raise subprocess.CalledProcessError(ret_bt2, cmd_bt2)
        if ret_sam != 0: raise subprocess.CalledProcessError(ret_sam, cmd_sam_sort)
        
        shutil.rmtree(temp_dir)
        logging.info("✅ Step 3 complete.")
        return pe_namesorted_bam

    def _run_step4_merge_filter(self, se_namesorted_bam, pe_namesorted_bam):
        logging.info(f"--- [Step 4/6] Information Merge & Filtering for {self.sample_name} ---")
        final_merged_mi_sorted_bam = self.dirs['step4'] / f"{self.sample_name}.merged.tagged.misorted.bam"
        merged_unsorted_bam = self.dirs['step4'] / f"{self.sample_name}.merged.tmp.bam"

        logging.info("  - (4a) Transferring UMI tags to high-confidence PE reads...")
        _merge_tags(se_namesorted_bam, pe_namesorted_bam, merged_unsorted_bam, ['MI'])

        logging.info("  - (4b) Sorting final merged BAM by MI tag...")
        self._run_cmd(['samtools', 'sort', '-t', 'MI', '-@', str(self.config['THREADS']), str(merged_unsorted_bam), '-o', str(final_merged_mi_sorted_bam)])

        merged_unsorted_bam.unlink()
        logging.info("✅ Step 4 complete.")
        return final_merged_mi_sorted_bam

    def _run_step5_dna_mutation_analysis(self, merged_mi_sorted_bam):
        logging.info(f"--- [Step 5/6] DNA Mutation Analysis & Substrate De-multiplexing for {self.sample_name} ---")
        _analyze_clusters_by_majority(
            merged_mi_sorted_bam, self.dirs['step5'], self.config['THREADS'], self.config['MIN_PAIRS_PER_CLUSTER']
        )
        logging.info("✅ Step 5 complete.")

    def _run_step6_aa_mutation_calling(self):
        logging.info(f"--- [Step 6/6] Amino Acid Mutation Calling & Quantification for {self.sample_name} ---")
        
        split_bams = list(self.dirs['step5'].glob("*.bam"))
        if not split_bams:
             logging.warning("No split BAM files found in step 5 directory. Skipping step 6.")
             return

        for split_bam_file in split_bams:
            if not (split_bam_file.with_suffix(split_bam_file.suffix + '.bai')).exists():
                logging.warning(f"Index not found for {split_bam_file.name}. File may be empty. Skipping.")
                continue
                
            basename = split_bam_file.stem
            details_file = self.dirs['step6'] / f"{basename}.details.txt"
            counts_file = self.dirs['step6'] / f"{basename}.counts.tsv"
            
            _extract_and_translate_codons(split_bam_file, details_file, counts_file)
        
        logging.info("✅ Step 6 complete.")

    def run_pipeline_for_sample(self, sample_dir: Path):
        """Execute the full analysis pipeline for a single sample directory."""
        self._setup_directories(sample_dir)
        
        logging.info("=================================================================================")
        logging.info(f"▶️ Starting Full Pipeline for Sample: {self.sample_name}")
        logging.info("=================================================================================")
        
        try:
            raw_fastq1 = next(sample_dir.glob("*_1.fq.gz"))
            raw_fastq2 = next(sample_dir.glob("*_2.fq.gz"))
        except StopIteration:
            logging.warning(f"Paired FASTQ files not found for {sample_dir}. Skipping.")
            return

        try:
            trimmed_r1, trimmed_r2 = self._run_step1_trimming(raw_fastq1, raw_fastq2)
            se_namesorted_bam = self._run_step2_umi_processing(trimmed_r1)
            pe_namesorted_bam = self._run_step3_pe_alignment(trimmed_r1, trimmed_r2)
            merged_mi_sorted_bam = self._run_step4_merge_filter(se_namesorted_bam, pe_namesorted_bam)
            self._run_step5_dna_mutation_analysis(merged_mi_sorted_bam)
            self._run_step6_aa_mutation_calling()

            logging.info(f"\n🎉🎉🎉 Successfully completed the full pipeline for sample: {self.sample_name} 🎉🎉🎉")

        except Exception as e:
            logging.error(f"❌❌❌ Pipeline failed for sample {self.sample_name}. Error: {e} ❌❌❌")
            return


def check_dependencies():
    """Verify that all required command-line tools are in the system's PATH."""
    logging.info("Starting pre-run dependency checks...")
    dependencies = ['cutadapt', 'umi_tools', 'bowtie2', 'samtools', 'umicollapse']
    missing = []
    for tool in dependencies:
        if not shutil.which(tool):
            missing.append(tool)
    if missing:
        logging.error(f"FATAL: Required tools not found in PATH: {', '.join(missing)}")
        sys.exit(1)
    
    try:
        import pysam
        import tqdm
    except ImportError as e:
        logging.error(f"FATAL: Missing required Python package: {e.name}. Please install it (e.g., 'pip install {e.name}')")
        sys.exit(1)

    logging.info("All dependencies are satisfied.")


def main():
    parser = argparse.ArgumentParser(
        description="UMI-Tagged Amplicon Sequencing Analysis Pipeline for the att-IR Platform.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '--samples', 
        nargs='+', 
        required=True,
        help="One or more paths to sample directories."
    )
    parser.add_argument(
        '--ref',
        required=True,
        help="Path to the reference genome FASTA file."
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=cpu_count(),
        help=f"Number of threads to use for parallel tasks. (Default: {cpu_count()})"
    )
    args = parser.parse_args()

    # --- Pipeline Configuration ---
    CONFIG = {
        "THREADS": args.threads,
        "SORT_MEM_PER_THREAD": "1G",
        "REF_GENOME": Path(args.ref),
        "REF_BT2_INDEX": str(Path(args.ref).with_suffix('')),
        "ADAPTER_R1_3PRIME": 'CATGGGCGTATGGGCGTCGA',
        "ADAPTER_R2_5PRIME": '^AGACGATAAC',
        "ADAPTER_R2_3PRIME": 'TCAGCGCGCAAATACGCATA',
        "MIN_LEN_R1": 98,
        "MIN_LEN_R2": 96,
        "MAX_LEN_R1": 125,
        "MAX_LEN_R2": 102,
        "QUALITY_CUTOFF": 20,
        "UMI_PATTERN": '(?P<umi_1>.{22}).*',
        "SE_FILTER_CHROMOSOME": 'integrated',
        "SE_FILTER_POSITIONS": '44,45,46,50,51,52',
        "TRIM_BASES_R1_FOR_PE": 22,
        "MIN_PAIRS_PER_CLUSTER": 2,
    }

    check_dependencies()

    ref_bt2_index_path = Path(f"{CONFIG['REF_BT2_INDEX']}.1.bt2")
    if not ref_bt2_index_path.exists():
        logging.error(f"FATAL: Bowtie2 index not found! Expected '{ref_bt2_index_path}'.")
        logging.error(f"Please build it first: bowtie2-build {CONFIG['REF_GENOME']} {CONFIG['REF_BT2_INDEX']}")
        sys.exit(1)
        
    pipeline = AttIRPipeline(CONFIG)
    
    sample_dirs = [Path(p) for p in args.samples]
    for sample_dir in sample_dirs:
        if not sample_dir.is_dir():
            logging.warning(f"Provided sample path is not a directory: {sample_dir}. Skipping.")
            continue
        pipeline.run_pipeline_for_sample(sample_dir)

    logging.info("\n=================================================================================")
    logging.info("All specified samples have been processed.")
    logging.info("=================================================================================")


if __name__ == '__main__':
    main()