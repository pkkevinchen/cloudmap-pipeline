import shlex, subprocess
import sys

def read_input_files( input_file ):
    file_dict = {}
    for in_line in input_file:
        if in_line[0] == "#":
            pass
        else:
            parsed_in_line = in_line.strip().split()
            file_dict[parsed_in_line[0]] = parsed_in_line[1:]
    return file_dict

def read_config( config_file ):
    config_dict = {}
    for i in config_file:
        if i[0] == "#":
            pass
        else:
            k = i.strip().split()
            config_dict[k[0]] = k[1]
    return config_dict

def popen_simple(cmd_line, log_file):
    cmd_args = shlex.split(cmd_line)
    cmd_process = subprocess.Popen(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    cmd_std = cmd_process.communicate()

    log_file.write(cmd_std[0])
    log_file.write(cmd_std[1])

def popen_write_to_file(cmd_line, write_file, log_file):
    cmd_args = shlex.split(cmd_line)
    cmd_process = subprocess.Popen(cmd_args, stdout=write_file, stderr=subprocess.PIPE)
    cmd_std = cmd_process.communicate()
    write_file.close()

    log_file.write(cmd_std[1])

def prinseq_filter(prinseq_env, file_to_filter, filtered_file, script_file, log_file):

    prinseq_cmd = "%s -fastq %s " % (prinseq_env, file_to_filter)
    prinseq_cmd += "-trim_qual_left 30 -trim_qual_right 30 -min_qual_mean 30 -min_len 50 -log prinseq_all.log "
    prinseq_cmd += "-out_good %s_good -out_bad null" % (filtered_file)

    script_file.write("# prinseq filter reads: %s\n" % (filtered_file))
    script_file.write(prinseq_cmd)
    script_file.write("\n\n")
    
    popen_simple(prinseq_cmd, log_file)

    print "prinseq filtered: " + file_to_filter + " to " + filtered_file

def bwa_align(bwa_env, ref_genome, out_dir, strain_name, script_file, log_file):

    # $bwa aln -I -t 6 c_elegans.genomic.fa test.fastq > test.sai
    bwa_aln_cmd = "%s aln -t 8 %s %s" % (bwa_env, ref_genome, out_dir+strain_name+"_good.fastq")

    # $bwa samse c_elegans.genomic.fa test.sai test.fastq > test.sam
    bwa_sam_cmd = "%s samse " % (bwa_env)
    bwa_sam_cmd += "-r '@RG\tID:%s\tSM:%s\tLB:%s\tPL:Illumina' " % (strain_name, strain_name, strain_name)
    bwa_sam_cmd += "%s %s.sai %s" % (ref_genome, out_dir+strain_name, out_dir+strain_name+"_good.fastq")

    script_file.write("# bwa align reads: %s\n" % (strain_name))
    script_file.write(bwa_aln_cmd + " > %s.sai \n" % (out_dir+strain_name))
    script_file.write(bwa_sam_cmd + " > %s.sam \n" % (out_dir+strain_name))
    script_file.write("\n")

    bwa_aln_out = open(out_dir+strain_name+".sai", 'w')
    popen_write_to_file(bwa_aln_cmd, bwa_aln_out, log_file)

    bwa_sam_out = open(out_dir+strain_name+".sam", 'w')
    popen_write_to_file(bwa_sam_cmd, bwa_sam_out, log_file)

    print "bwa aligned: " + strain_name

def process_sam(sam_env, out_dir, strain_name, script_file, log_file):

    # $samtools view -h -F 4 test.sam -o test_mapped.sam
    sam_view_cmd = "%s view -h -F 4 %s.sam -o %s_mapped.sam" % (sam_env, out_dir+strain_name, out_dir+strain_name)
    # $samtools sort --threads 8 -O bam -o test_mapped_sorted.bam test_mapped.sam
    sam_sort_cmd = "%s sort --threads 8 -O bam -o %s_mapped_sorted.bam %s_mapped.sam" % (sam_env, out_dir+strain_name, out_dir+strain_name)
    # $samtools index test_mapped_sorted.bam
    sam_index_cmd = "%s index %s_mapped_sorted.bam" % (sam_env, out_dir+strain_name)

    script_file.write("# samtools process reads: %s\n" % (strain_name))
    script_file.write(sam_view_cmd + "\n")
    script_file.write(sam_sort_cmd + "\n")
    script_file.write(sam_index_cmd + "\n")
    script_file.write("\n")

    popen_simple(sam_view_cmd, log_file)
    popen_simple(sam_sort_cmd, log_file)
    popen_simple(sam_index_cmd, log_file)

    print "SAM file processed: " + strain_name

def gatk_realign(gatk_env, ref_genome, out_dir, strain_name, script_file, log_file):

    #$GenomeAnalysisTK.jar -T RealignerTargetCreator -R c_elegans.genomic.fa -I test_mapped_sorted.bam -o test.intervals
    gatk_realign_target_cmd = "%s -T RealignerTargetCreator -R %s -I %s_mapped_sorted.bam -o %s.intervals" % (gatk_env, ref_genome, out_dir+strain_name, out_dir+strain_name)
    #$GenomeAnalysisTK.jar -T IndelRealigner -R c_elegans.genomic.fa -I test_mapped_sorted.bam -targetIntervals test.intervals -o test_realigned.bam
    gatk_realign_cmd = "%s -T IndelRealigner -R %s -I %s_mapped_sorted.bam -targetIntervals %s.intervals -o %s_realigned.bam" % (gatk_env, ref_genome, out_dir+strain_name, out_dir+strain_name, out_dir+strain_name)

    script_file.write("# GATK realign: %s\n" % (strain_name))
    script_file.write(gatk_realign_target_cmd + "\n")
    script_file.write(gatk_realign_cmd + "\n")
    script_file.write("\n")

    popen_simple(gatk_realign_target_cmd, log_file)
    popen_simple(gatk_realign_cmd, log_file)

    print "GATK realigned: " + strain_name

def picard_markdup(picard_env, sam_env, out_dir, strain_name, script_file, log_file):

    # $picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=test_1_realigned.bam I=test_2_realigned.bam O=test_all_marked_dup.bam M=test_all_marked_dup_metrics.txt
    picard_markdup_cmd = "%s MarkDuplicates REMOVE_DUPLICATES=true I=%s_realigned.bam" % (picard_env, out_dir+strain_name)
    picard_markdup_cmd += " O=%s_marked_dup.bam M=%s_marked_dup_metrics.txt " % (out_dir+strain_name, out_dir+strain_name)
    # $samtools index test_all_marked_dup.bam
    sam_index_cmd = "%s index %s_marked_dup.bam " % (sam_env, out_dir+strain_name)
    
    script_file.write("# picard mark duplicate reads: %s\n" % (strain_name))
    script_file.write(picard_markdup_cmd + "\n")
    script_file.write(sam_index_cmd + "\n")
    script_file.write("\n")

    popen_simple(picard_markdup_cmd, log_file)
    popen_simple(sam_index_cmd, log_file)

    print "picard duplicate marked: " + strain_name

def gatk_genotype(gatk_env, ref_genome, out_dir, strain_name, script_file, log_file):

    #$GenomeAnalysisTK.jar -T UnifiedGenotyper -R c_elegans.genomic.fa -I test_all_marked_dup.bam -o test_CM_variants.vcf
    gatk_genotype_cmd = "%s -T UnifiedGenotyper -R %s -I %s_marked_dup.bam -o %s_variants.vcf " % (gatk_env, ref_genome, out_dir+strain_name, out_dir+strain_name)
    gatk_genotype_cmd += "--num_threads 8 --standard_min_confidence_threshold_for_calling 30 --standard_min_confidence_threshold_for_emitting 30 "
    gatk_genotype_cmd += "--read_filter MappingQuality --interval_set_rule UNION --interval_merging ALL --downsampling_type NONE "
    gatk_genotype_cmd += "--baq OFF --baqGapOpenPenalty 40.0 --validation_strictness STRICT --pedigreeValidationType STRICT --logging_level INFO "
    gatk_genotype_cmd += "--genotype_likelihoods_model BOTH --p_nonref_model EXACT_ORIGINAL --heterozygosity 0.001 --pcr_error_rate 1.0E-4 "
    gatk_genotype_cmd += "--genotyping_mode DISCOVERY --output_mode EMIT_VARIANTS_ONLY --min_base_quality_score 17 --max_deletion_fraction 0.05 "
    gatk_genotype_cmd += "--max_alternate_alleles 5 --min_indel_count_for_genotyping 5 --indel_heterozygosity 1.25E-4 --indelGapContinuationPenalty 10 "
    gatk_genotype_cmd += "--indelGapOpenPenalty 45 --indelHaplotypeSize 80 --min_mapping_quality_score 10 "

    script_file.write("# GATK genotype: %s\n" % (strain_name))
    script_file.write(gatk_genotype_cmd + "\n")
    script_file.write("\n")

    popen_simple(gatk_genotype_cmd, log_file)

    print "GATK genotyped: " + strain_name

def subtract_variants(vcf_filter_env, out_dir, strain_name, parent_vcf, cross_vcf, filter_qual, script_file, log_file):

    subtract_parent_cmd = "%s %s %s_variants.vcf %s_noParent_q%s.vcf %s" % (vcf_filter_env, parent_vcf, out_dir+strain_name, out_dir+strain_name, filter_qual, filter_qual)
    subtract_cross_cmd = "%s %s %s_noParent_q%s.vcf %s_unique_q%s.vcf %s" % (vcf_filter_env, cross_vcf, out_dir+strain_name, filter_qual, out_dir+strain_name, filter_qual, filter_qual)

    script_file.write("# subtract variants: %s\n" % (strain_name))
    script_file.write(subtract_parent_cmd + "\n")
    script_file.write(subtract_cross_cmd + "\n")
    script_file.write("\n")

    popen_simple(subtract_parent_cmd, log_file)
    popen_simple(subtract_cross_cmd, log_file)

    print "Subtracted VCF from parental strain: " + strain_name

def snpeff_call(snpeff_env, out_dir, strain_name, script_file, log_file):

    snpeff_cmd = "%s WS250 %s_unique_q50.vcf " % (snpeff_env, out_dir+strain_name)

    script_file.write("# SnpEff calls: %s\n" % (strain_name))
    script_file.write(snpeff_cmd + " > " + out_dir+strain_name + "_snpeff.vcf \n")
    script_file.write("\n")

    snpeff_out = open(out_dir+strain_name+"_snpeff.vcf", 'w')
    popen_write_to_file(snpeff_cmd, snpeff_out, log_file)

    print "Snpeff called: " + strain_name

def VDM_mapping(vdm_env, out_dir, strain_name, filter_qual, script_file, log_file):

    vdm_map_cmd = "%s --sample_vcf %s_unique_q%s.vcf --output %s_vdm_variants_q%s.csv " % (vdm_env, out_dir+strain_name, filter_qual, out_dir+strain_name, filter_qual)
    vdm_map_cmd += "--location_plot_output %s_vdm_variants_q%s.pdf " % (out_dir+strain_name, filter_qual)
    vdm_map_cmd += "--loess_span \"0.4\" --d_yaxis \"1.0\" --h_yaxis \"0\" --points_color \"gray27\" "
    vdm_map_cmd += "--loess_color \"green\" --standardize \"true\" --normalize_bins \"true\" --break_file \"C.elegans\" "

    script_file.write("# VDM mapping: %s\n" % (strain_name))
    script_file.write(vdm_map_cmd + "\n")
    script_file.write("\n")

    popen_simple(vdm_map_cmd, log_file)

    print "VDM mapped: " + strain_name

def snpeff_cleanup(snpeff_clean_env, out_dir, strain_name, script_file, log_file):

    snpeff_clean_cmd = "%s %s_snpeff.vcf %s_snpeff_clean.csv" % (snpeff_clean_env, out_dir+strain_name, out_dir+strain_name)

    script_file.write("# snpeff clean: %s\n" % (strain_name))
    script_file.write(snpeff_clean_cmd + "\n")
    script_file.write("\n")

    popen_simple(snpeff_clean_cmd, log_file)

    print "Snpeff cleanup: " + strain_name

def main():
    try:
        input_file_list = open( sys.argv[1] ).readlines()
        config_file = open( sys.argv[2] ).readlines()
        script_file = open( sys.argv[3]+".sh", "w")
        log_file = open( sys.argv[3]+".log", "w")
    except IndexError:
        print "Please include: [input_file_list], [config_file], [run_name]"
        exit()

    all_files = read_input_files(input_file_list)
    strain_names = all_files.keys()

    config_dict = read_config(config_file)

    env_ngs_file_directory = config_dict["ngs_file_dicectory"]
    env_out_directory = config_dict["analysis_result_directory"]
    env_ref_genome_file = config_dict["ref_genome_file"]

    env_prinseq = "perl " + config_dict["prinseq"]
    env_fastqc = config_dict["fastqc"]
    env_bwa = config_dict["bwa"]
    env_samtools = config_dict["samtools"]
    env_gatk = config_dict["java"] + " -Xmx6g -jar " + config_dict["gatk_jar"]
    env_picard = config_dict["java"] + " -Xmx10g -jar " + config_dict["picard_jar"]
    env_snpeff = config_dict["java"] + " -Xmx4g -jar " + config_dict["snpeff_jar"]
    env_vdm_mapping = "python " + config_dict["vdm_mapping"]
    env_subtract_variants = "python " + config_dict["filter_varient"]
    env_snpeff_cleanup = "python " + config_dict["snpeff_to_csv"]
    

    #filter with prinseq and combine all files of one strain
    for strain in strain_names:
        strain_files = all_files[strain]
        if len(strain_files) < 1:
            print "No sequnce file found, please check input_file_table.csv"
            exit()
        else:
            #file_num = len(strain_files)
            combine_good_fq_cmd = []
            combine_good_fq_cmd.append('cat')
            file_counter = 1
            for file_name in strain_files:
                out_file = strain + "_" + str(file_counter)
                prinseq_filter(env_prinseq, env_ngs_file_directory+file_name, env_out_directory+out_file, script_file, log_file)
                file_counter += 1
                combine_good_fq_cmd.append(env_out_directory+out_file+"_good.fastq")

            #print combine_good_fq_cmd
            
            combine_fq_write = open(env_out_directory+strain+"_good.fastq", 'w')
            combine_process = subprocess.Popen(combine_good_fq_cmd, stdout=combine_fq_write)
            errcode = combine_process.wait()
            combine_fq_write.close()

    # bwa for each strain
    
    for strain in strain_names:
        bwa_align(env_bwa, env_ref_genome_file, env_out_directory, strain, script_file, log_file)
        process_sam(env_samtools, env_out_directory, strain, script_file, log_file)
        gatk_realign(env_gatk, env_ref_genome_file, env_out_directory, strain, script_file, log_file)
        picard_markdup(env_picard, env_samtools, env_out_directory, strain, script_file, log_file)
        gatk_genotype(env_gatk, env_ref_genome_file, env_out_directory, strain, script_file, log_file)
        subtract_variants(env_subtract_variants, env_out_directory, strain, config_dict["parent_variant"], config_dict["cross_variant"], "300", script_file, log_file)
        subtract_variants(env_subtract_variants, env_out_directory, strain, config_dict["parent_variant"], config_dict["cross_variant"], "200", script_file, log_file)
        subtract_variants(env_subtract_variants, env_out_directory, strain, config_dict["parent_variant"], config_dict["cross_variant"], "100", script_file, log_file)
        subtract_variants(env_subtract_variants, env_out_directory, strain, config_dict["parent_variant"], config_dict["cross_variant"], "50", script_file, log_file)
        VDM_mapping(env_vdm_mapping, env_out_directory, strain, "300", script_file, log_file)
        VDM_mapping(env_vdm_mapping, env_out_directory, strain, "200", script_file, log_file)
        VDM_mapping(env_vdm_mapping, env_out_directory, strain, "100", script_file, log_file)
        VDM_mapping(env_vdm_mapping, env_out_directory, strain, "50", script_file, log_file)
        snpeff_call(env_snpeff, env_out_directory, strain, script_file, log_file)
        snpeff_cleanup(env_snpeff_cleanup, env_out_directory, strain, script_file, log_file)

if __name__ == '__main__':
    main()

