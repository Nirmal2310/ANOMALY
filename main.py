import yaml
import subprocess
import os
import argparse

def load_config(config_path):
	new_path=os.path.abspath(config_path)
	with open(new_path, 'r') as file:
		config = yaml.safe_load(file)
	return config

path=subprocess.run("conda info --envs | grep 'anomaly' | awk '{print $2}'", \
                    shell=True, stdout=subprocess.PIPE, encoding='utf-8')\
                    .stdout.strip()

current_path=os.path.abspath(os.getcwd())

def generate_alignment_shell_script(config):
    alignment_mode = config.get("read_type", "None")
    reference_file = config.get("reference_nuclear", "ref.fa")
    threads = config.get("threads_minimap2", 16)
    if alignment_mode == "ONT":
        mode = "map-ont"
    elif alignment_mode == "HiFi":
        mode = "map-hifi"
    elif alignment_mode == "pb":
        mode = "map-pb"
    elif alignment_mode == "None":
        mode = ''
    else:
        raise ValueError(f"Unknown alignment mode: {alignment_mode}")
    
    shell_script = f"""
    reference={reference_file}
    input_fastq={{input}}
    sorted_bam={{output}}
    threads={threads}
    mode={mode}

    minimap2 -ax "$mode" -t "$threads" -Y "$reference" "$input_fastq" | \
    samtools view -@ "$threads" -bS | \
    samtools sort -@ "$threads" -o "$sorted_bam" && samtools index -@ "$threads" "$sorted_bam"

    echo "Alignment, sorting, and indexing complete. Output file: $sorted_bam"
    """


    return shell_script
    
def generate_sniffle_shell_script(config):
    min_support = config.get("min_support", None)
    minsvlen = config.get("minimum_numt_length", None)
    minimum_mapq = config.get("minimum_mapq")
    long_ins_length = config.get("long_ins_length", 2500)
    genotype_ploidy = config.get("genotype_ploidy", 2)
    quiet = config.get("quiet", False)
    allow_overwrite = config.get("allow_overwrite", False)
    threads_sniffles = config.get("threads_sniffles", 16)
    threads = threads_sniffles if threads_sniffles else config.get("threads", 16)
    min_support = min_support if min_support else config.get("min_support", None)
    param_string = []

    if min_support != None:
        param_string.append(f"--minsupport {min_support}")

    if minimum_mapq is not None:
        param_string.append(f"--mapq {minimum_mapq}")
        
    if minsvlen != None:
        param_string.append(f"--minsvlen {minsvlen}")

    if long_ins_length:
        param_string.append(f"--long-ins-length {long_ins_length}")

    if genotype_ploidy:
        param_string.append(f"--genotype-ploidy {genotype_ploidy}")

    if quiet:
        param_string.append("--quiet")

    if allow_overwrite:
        param_string.append("--allow-overwrite")

    parameters = " ".join(param_string)

    basic_command = f"sniffles --threads {threads} {parameters} --input {{input}} --vcf {{output}} > {{log}} 2>&1"

    return basic_command



def generate_snakemake_rule_alignment(config):
    alignment_shell_script = generate_alignment_shell_script(config)
    sample_folder = config.get("sample_directory")
    working_directory = config.get("working_directory")
    rule_alignment = f"""
rule alignment:
    input:
        "{sample_folder}/{{sample}}.fastq.gz"
    output:
        "{working_directory}/bam_files/{{sample}}.bam"
    conda: "{path}"
    resources: 
        mem_mb=180000
    benchmark:
        "{working_directory}/bam_files/{{sample}}.mapping.benchmark.txt"
    shell:
        \"\"\"
        #!/bin/bash
        {alignment_shell_script}
        \"\"\"
"""
    return rule_alignment

def generate_snakemake_rule_sniffles(config):
    sniffles_shell_script=generate_sniffle_shell_script(config)
    bam_or_fastq = config.get("bam_or_fastq")
    sample_folder = config.get("sample_directory")
    working_directory = config.get("working_directory")
    threads = config.get("threads_sniffles", 16)

    if bam_or_fastq == 'b':
        rule_sniffles = f"""
rule sniffles:
    input:"{sample_folder}/{{sample}}.bam"
    output:"{working_directory}/vcf_files/{{sample}}_sniffles.vcf.gz"
    log:"{working_directory}/vcf_files/{{sample}}_sniffles.log"
    conda: "{path}"
    threads: {threads}
    resources: mem_mb=180000
    benchmark:
        "{working_directory}/vcf_files/{{sample}}.sniffles.benchmark.txt"
    shell:
        \"\"\"
        #!/bin/bash
        {sniffles_shell_script}
        \"\"\"
        """
    if bam_or_fastq == 'f':
        rule_sniffles = f"""
rule sniffles:
    input:"{working_directory}/bam_files/{{sample}}.bam"
    output:"{working_directory}/vcf_files/{{sample}}_sniffles.vcf.gz"
    log:"{working_directory}/vcf_files/{{sample}}_sniffles.log"
    conda: "{path}"
    threads: {threads}
    resources: mem_mb=180000
    benchmark:
        "{working_directory}/vcf_files/{{sample}}.sniffles.benchmark.txt"
    shell:
        \"\"\"
        #!/bin/bash
        {sniffles_shell_script}
        \"\"\"
        """
    return rule_sniffles

def generate_snakemake_rule_inserts(config):
    working_directory = config.get("working_directory")
    headers_file = config.get("headers_nuclear", "chr_headers.txt")
    rule_inserts = f""" 
rule inserts:
    input:
        vcf_file='{working_directory}/vcf_files/{{sample}}_sniffles.vcf.gz',
        ref_list='{headers_file}'
    output:
        '{working_directory}/inserts/{{sample}}_insertions.txt'
    benchmark:
        "{working_directory}/inserts/{{sample}}.get_inserts.benchmark.txt"
    shell:
        \"\"\"
        bash \"{current_path}/Scripts/get_insertion_calls.sh\" {{input.vcf_file}} {{output}} {{input.ref_list}}
        \"\"\"
    """
    return rule_inserts

def generate_snakemake_rule_get_inserts_fasta(config):
    working_directory = config.get("working_directory")
    rule_get_inserts_fasta = f"""
rule get_inserts_fasta:
    input: 
        '{working_directory}/inserts/{{sample}}_insertions.txt'
    output:
        '{working_directory}/fasta_files/{{sample}}_inserts.fasta'
    benchmark:
        "{working_directory}/fasta_files/{{sample}}.create_inserts_fasta.benchmark.txt"
    shell:
        '''
        bash \"{current_path}/Scripts/make_inserts_fasta.sh\" {{input}} {{output}}
        '''
    """
    return rule_get_inserts_fasta

def generate_snakemake_rule_inserts_blast(config):
    working_directory = config.get("working_directory")
    headers_file = config.get("headers_nuclear", "chr_headers.txt")
    evalue_cutoff = config.get("evalue_cutoff")
    threads = config.get("threads_sniffles", 16)
    rule_inserts_blast = f"""
rule inserts_blast:
    input: 
        fasta_file="{working_directory}/fasta_files/{{sample}}_inserts.fasta",
        ref_list="{headers_file}",
    
    output:'{working_directory}/filtered/{{sample}}_blast_result_filtered.txt'
    threads: {threads}
    params:
        evalue_cutoff={evalue_cutoff}
    conda: "{path}"
    benchmark:
        "{working_directory}/filtered/{{sample}}.blast.benchmark.txt"
    shell:
        '''
        bash \"{current_path}/Scripts/blast_run.sh\" {{input.fasta_file}} {{output}} {{threads}} {{input.ref_list}} {{params.evalue_cutoff}}
        '''
    """
    return rule_inserts_blast

def generate_snakemake_rule_inserts_numt_concat(config):
    working_directory = config.get("working_directory")
    coverage_cutoff = config.get("numt_coverage")
    rule_inserts_numt_concat = f"""
rule inserts_numt_concat:
    input: 
        filtered_file='{working_directory}/filtered/{{sample}}_blast_result_filtered.txt'

    output:'{working_directory}/NUMTs/{{sample}}_concatenated_numts.txt'

    conda: "{path}"

    params:
        coverage_cutoff={coverage_cutoff}
    
    benchmark:
        "{working_directory}/NUMTs/{{sample}}.numt_concatenation.benchmark.txt"

    shell:
        '''
        bash \"{current_path}/Scripts/numt_concat.sh\" {{input.filtered_file}} {{output}} {{params.coverage_cutoff}}
        '''
    """
    return rule_inserts_numt_concat

def generate_snakemake_rule_get_supplementary_alignments(config):
    working_directory = config.get("working_directory")
    bam_or_fastq = config.get("bam_or_fastq")
    sample_directory = config.get("sample_directory")
    threads = config.get("threads_sniffles",16)
    headers_file = config.get("headers_nuclear", "chr_headers.txt")

    if bam_or_fastq == 'b':
        rule_get_supplementary_alignments = f"""
rule get_supplementary_alignments:
    input: 
        main_bam="{sample_directory}/{{sample}}.bam",
        ref_headers="{headers_file}"
    
    output:'{working_directory}/SA_Data/{{sample}}_MT_SA.bam'
    threads: {threads}
    benchmark:
        "{working_directory}/SA_Data/{{sample}}.supplementary_alignment.benchmark.txt"
    shell:
        '''
        bash \"{current_path}/Scripts/get_supplementary_alignments.sh\" {{input.main_bam}} {{output}} {{threads}} {{input.ref_headers}}
        '''
    """
    if bam_or_fastq == 'f':
        rule_get_supplementary_alignments = f"""
rule get_supplementary_alignments:
    input: 
        main_bam="{working_directory}/bam_files/{{sample}}.bam",
        ref_headers="{headers_file}"
    
    output:"{working_directory}/SA_Data/{{sample}}_MT_SA.bam"
    threads: {threads}
    benchmark:
        "{working_directory}/SA_Data/{{sample}}.supplementary_alignment.benchmark.txt"
    shell:
        '''
        bash \"{current_path}/Scripts/get_supplementary_alignments.sh\" {{input.main_bam}} {{output}} {{threads}} {{input.ref_headers}}
        '''
    """
    return rule_get_supplementary_alignments

def generate_snakemake_rule_potential_numts_from_sa(config):
    working_directory = config.get("working_directory")
    threads = config.get("threads_sniffles",16)
    supporting_reads = config.get("numt_supporting_reads")
    rule_potential_numts_from_sa = f"""
rule potential_numts_from_sa:
    input: 
        data="{working_directory}/SA_Data/{{sample}}_MT_SA.bam",
    output: 
        sa_calls="{working_directory}/filtered/{{sample}}_MT_SA_calls.txt",
        potential_numts="{working_directory}/filtered/{{sample}}_potential_numts_from_sa.txt"
    threads: {threads}
    params:
        read_cutoff={supporting_reads}
    benchmark:
        "{working_directory}/filtered/{{sample}}.get_numts_from_sa.benchmark.txt"
    shell:
        '''
        bash \"{current_path}/Scripts/get_potential_numts_from_sa.sh\" {{input.data}} {{output.sa_calls}} {{output.potential_numts}} {{threads}} {{params.read_cutoff}}
        '''
    """
    return rule_potential_numts_from_sa

def generate_snakemake_rule_final_numts_from_sa(config):
    working_directory = config.get("working_directory")
    rule_final_numts_from_sa = f"""
rule final_numts_from_sa:
    input: 
        sa_calls="{working_directory}/filtered/{{sample}}_MT_SA_calls.txt",
        potential_numts="{working_directory}/filtered/{{sample}}_potential_numts_from_sa.txt"

    output: "{working_directory}/NUMTs/{{sample}}_final_numts_from_sa.txt"

    benchmark:
        "{working_directory}/NUMTs/{{sample}}.final_numts_from_sa.benchmark.txt"

    shell:
        '''
        bash \"{current_path}/Scripts/get_final_numts_from_sa.sh\" {{input.sa_calls}} {{output}} {{input.potential_numts}}
        '''
    """
    return rule_final_numts_from_sa

def generate_snakemake_rule_remove_duplicate_numts(config):
    working_directory = config.get("working_directory")
    rule_remove_duplicate_numts = f"""
rule remove_duplicate_numts:
    input: 
        sa_numt="{working_directory}/NUMTs/{{sample}}_final_numts_from_sa.txt",
        ins_numt="{working_directory}/NUMTs/{{sample}}_concatenated_numts.txt"

    output: "{working_directory}/Results/{{sample}}_final_numts.txt"

    benchmark:
        "{working_directory}/Results/{{sample}}.remove_duplicates.benchmark.txt"

    shell:
        '''
	bash \"{current_path}/Scripts/remove_duplicates.sh\" -i {{input.ins_numt}} -s {{input.sa_numt}} -o {{output}}
        '''
    """
    return rule_remove_duplicate_numts

def generate_snakemake_rule_visualize_numts(config):
    working_directory = config.get("working_directory")
    reference_file = config.get("reference_nuclear")
    headers_file = config.get("headers_nuclear", "chr_headers.txt")
    rule_visualize_numts = f"""
rule visualize_numts:
    input:
        input_file="{working_directory}/Results/{{sample}}_final_numts.txt",
        ref_genome="{reference_file}",
        ref_headers="{headers_file}"
    output:
        svg_out="{working_directory}/Results/{{sample}}_numt_circos_plot.svg",
        png_out="{working_directory}/Results/{{sample}}_numt_circos_plot.png"
    benchmark:
        "{working_directory}/Results/{{sample}}.ciros_plot.benchmark.txt"

    shell:
        '''
    bash \"{current_path}/Scripts/visualization.sh\" -i {{input.input_file}} -r {{input.ref_genome}} -v {{input.ref_headers}} -s {{output.svg_out}} -p {{output.png_out}}
        '''
    """
    return rule_visualize_numts

def generate_snakemake_rule_combined_circos_plot(config):
    working_directory = config.get("working_directory")
    reference_file = config.get("reference_nuclear")
    header_file = config.get("headers_nuclear", "chr_headers.txt")
    rule_combined_circos_numts = f"""
rule combined_circos_numts:
    input:
        sample_plots = expand("{working_directory}/Results/{{sample}}_numt_circos_plot.svg", sample=samples),
        sample_numts = expand("{working_directory}/Results/{{sample}}_final_numts.txt", sample=samples),
        ref_genome = "{reference_file}",
        ref_headers = "{header_file}"
    output:
        svg_out="{working_directory}/Results/cohort_numt_circos_plot.svg",
        png_out="{working_directory}/Results/cohort_numt_circos_plot.png"
    
    conda: "{path}"

    benchmark:
        "{working_directory}/Results/cohort.circos_plot.benchmark.txt"
    
    shell:
        '''
    bash \"{current_path}/Scripts/cohort_visualization.sh\" -r {{input.ref_genome}} -v {{input.ref_headers}} -s {{output.svg_out}} -p {{output.png_out}}
        '''
    """
    return rule_combined_circos_numts

def generate_snakemake_rule_all(config,config_path):
    file_type = config.get("bam_or_fastq", "b")
    sample_folder = config.get("sample_directory")
    working_directory = config.get("working_directory")
    new_path = os.path.abspath(config_path)
    if sample_folder is None:
        raise ValueError("The 'sample_directory' key must be specified in the config.")

    if file_type == "b":
        rule_all = f"""
import os
import glob

# Define the directory where samples are located
configfile: "{new_path}"
sample_dir = "{sample_folder}"
sample_files = glob.glob(os.path.join(sample_dir, "*.bam"))
samples = [os.path.basename(f).replace(".bam", "") for f in sample_files]

# Include the rest of the workflow from other files
include: "{working_directory}/sniffles.smk"
include: "{working_directory}/inserts.smk"
include: "{working_directory}/get_inserts_fasta.smk"
include: "{working_directory}/inserts_blast.smk"
include: "{working_directory}/inserts_numt_concat.smk"
include: "{working_directory}/get_supplementary_alignments.smk"
include: "{working_directory}/potential_numts_from_sa.smk"
include: "{working_directory}/final_numts_from_sa.smk"
include: "{working_directory}/remove_duplicate_numts.smk"
include: "{working_directory}/visualize_numts.smk"
include: "{working_directory}/combined_circos_plot.smk"


rule all:
    input:
        expand("{working_directory}/Results/{{sample}}_numt_circos_plot.svg", sample=samples),
        expand("{working_directory}/Results/{{sample}}_numt_circos_plot.png", sample=samples),
        "{working_directory}/Results/cohort_numt_circos_plot.svg",
        "{working_directory}/Results/cohort_numt_circos_plot.png"
"""
    else:
        # FASTQ file processing
        rule_all = f"""
import os
import glob
configfile: "{new_path}"
# Define the directory where samples are located
sample_dir = "{sample_folder}"
sample_files = glob.glob(os.path.join(sample_dir, "*.fastq.gz"))
samples = [os.path.basename(f).replace(".fastq.gz", "") for f in sample_files]

# Include the rest of the workflow from other files
include: "{working_directory}/alignment.smk" 
include: "{working_directory}/sniffles.smk"
include: "{working_directory}/inserts.smk"
include: "{working_directory}/get_inserts_fasta.smk"
include: "{working_directory}/inserts_blast.smk"
include: "{working_directory}/inserts_numt_concat.smk"
include: "{working_directory}/get_supplementary_alignments.smk"
include: "{working_directory}/potential_numts_from_sa.smk"
include: "{working_directory}/final_numts_from_sa.smk"
include: "{working_directory}/remove_duplicate_numts.smk"
include: "{working_directory}/visualize_numts.smk"
include: "{working_directory}/combined_circos_plot.smk"

rule all:
    input:
        expand("{working_directory}/Results/{{sample}}_numt_circos_plot.svg", sample=samples),
        expand("{working_directory}/Results/{{sample}}_numt_circos_plot.png", sample=samples),
        "{working_directory}/Results/cohort_numt_circos_plot.svg",
        "{working_directory}/Results/cohort_numt_circos_plot.png"
"""

    return rule_all
     

def write_snakefile(rule, snakefile_path):
    with open(snakefile_path, "w") as snakefile:
        snakefile.write(rule)


def run_snakemake(snakefile_path="Snakefile", config_path="config.yml"):
    try:
        subprocess.run(["snakemake", "-s", snakefile_path, "--configfile", config_path], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Snakemake failed: {e}")
        raise

def main():

    parser = argparse.ArgumentParser(description="config_file")
    parser.add_argument(
    '-c', '--config',
    type=str,
    required=True,
    help="Absolute Path to the configuration YAML file. (Working Directory + snake_config.yml)"
    )
    args = parser.parse_args()
    config_path = args.config
    config = load_config(config_path)
    
    working_directory = config.get("working_directory")
    
    sample_directory = config.get("sample_directory")
    
    directories = ['fasta_files','NUMTs','inserts','filtered','vcf_files','bam_files','SA_Data', 'Results']
    
    for directory in directories:
        folder = os.path.join(working_directory,directory)
        if not os.path.exists(folder):
            os.makedirs(folder)
            print(f"Directory '{folder}' created.")
    
    rule_alignment = generate_snakemake_rule_alignment(config)
    rule_alignment_path = os.path.join(working_directory,"alignment.smk")
    write_snakefile(rule_alignment,snakefile_path=rule_alignment_path)

    rule_sniffles = generate_snakemake_rule_sniffles(config)
    rule_sniffles_path = os.path.join(working_directory,"sniffles.smk")
    write_snakefile(rule_sniffles,snakefile_path=rule_sniffles_path)

    rule_inserts = generate_snakemake_rule_inserts(config)
    rule_inserts_path = os.path.join(working_directory,"inserts.smk")
    write_snakefile(rule_inserts,snakefile_path=rule_inserts_path)

    rule_fasta = generate_snakemake_rule_get_inserts_fasta(config)
    rule_fasta_path = os.path.join(working_directory,"get_inserts_fasta.smk")
    write_snakefile(rule_fasta,snakefile_path=rule_fasta_path)

    rule_blast = generate_snakemake_rule_inserts_blast(config)
    rule_blast_path = os.path.join(working_directory,"inserts_blast.smk")
    write_snakefile(rule_blast,snakefile_path=rule_blast_path)

    rule_numt = generate_snakemake_rule_inserts_numt_concat(config)
    rule_numt_path = os.path.join(working_directory,"inserts_numt_concat.smk")
    write_snakefile(rule_numt,snakefile_path=rule_numt_path)

    rule_get_alignments = generate_snakemake_rule_get_supplementary_alignments(config)
    rule_get_alignments_path = os.path.join(working_directory,"get_supplementary_alignments.smk")
    write_snakefile(rule_get_alignments,snakefile_path=rule_get_alignments_path)

    rule_potential_numts_from_sa = generate_snakemake_rule_potential_numts_from_sa(config)
    rule_potential_numts_from_sa_path = os.path.join(working_directory,"potential_numts_from_sa.smk")
    write_snakefile(rule_potential_numts_from_sa,snakefile_path=rule_potential_numts_from_sa_path)

    rule_final_numts_from_sa = generate_snakemake_rule_final_numts_from_sa(config)
    rule_final_numts_from_sa_path = os.path.join(working_directory,"final_numts_from_sa.smk")
    write_snakefile(rule_final_numts_from_sa,snakefile_path=rule_final_numts_from_sa_path)

    rule_remove_duplicate_numts = generate_snakemake_rule_remove_duplicate_numts(config)
    rule_remove_duplicate_numts_path = os.path.join(working_directory,"remove_duplicate_numts.smk")
    write_snakefile(rule_remove_duplicate_numts,snakefile_path=rule_remove_duplicate_numts_path)

    rule_visualize_numts = generate_snakemake_rule_visualize_numts(config)
    rule_visualize_numts_path = os.path.join(working_directory,"visualize_numts.smk")
    write_snakefile(rule_visualize_numts,snakefile_path=rule_visualize_numts_path)

    rule_combined_circos_plot = generate_snakemake_rule_combined_circos_plot(config)
    rule_combined_circos_plot_path = os.path.join(working_directory,"combined_circos_plot.smk")
    write_snakefile(rule_combined_circos_plot, rule_combined_circos_plot_path)

    rule_all = generate_snakemake_rule_all(config,config_path=config_path)
    rule_all_path = os.path.join(working_directory,"Snakefile")
    write_snakefile(rule_all,snakefile_path=rule_all_path)

if __name__ == "__main__":
    main()