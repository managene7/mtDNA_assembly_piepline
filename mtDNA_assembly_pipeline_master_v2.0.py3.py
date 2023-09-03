#Pre requisites
# 1. ncbi-blast+
# 2. Rust (https://www.rust-lang.org/learn/get-started)
# 3. msbwt2 (https://github.com/HudsonAlpha/rust-msbwt)
# 4. FMLRC2 (https://github.com/HudsonAlpha/fmlrc2)
# 5. NextDenovo (https://github.com/Nextomics/NextDenovo) with requirement softwares.
# 6. h5py (https://docs.h5py.org/en/latest/build.html)
# 7. rmblast (https://www.repeatmasker.org/rmblast/)
# 8. TRF (https://tandem.bu.edu/trf/downloads)
# 9. RepeatMasker (https://www.repeatmasker.org/RepeatMasker/)
# 10. MitoFinder (https://github.com/RemiAllio/MitoFinder) with requirement softwares. 
#________________ option parse _______________________________
import sys 
args = sys.argv[1:]
option_dict={'-out':"", '-out_type':"1", '-short_read2':"", '-gene_annotation':"yes", "-num_threads":"20","-gene_annotation":"no", '-min_ont_read_len':"10000", '-min_identity':'0.90','-estimated_mt_size':'600000'}
help_description="""
________________________________________________________________________________________________________________________________
Usage;
If default value exists, corresponding value can be omitted.
-long_read              (Required)  file name of ONT sequences for assembly 
                                    (compressed by 'gz' or uncompressed fastq format)
-min_ont_read_len       (option)    The length threshold of ONT-seq by base pair.
                                    (default:10000)
-short_read1            (Required)  file name of Illumina sequence for error correction
                                    (compressed by 'gz' or uncompressed fastq format)
-short_read2            (option)    the second file name of the paired Illumina sequence (option)
                                    (compressed by 'gz' or uncompressed fastq format)
-ref_mt_genome          (Required)  file name of reference mitochondrial genome sequence to isolate mtDNA containing ONT reads
                                    (FASTA format)
-ref_nucl_genome        (Required)  file name of reference nuclear genome sequence for filtering
                                    (FASTA format. mtDNA-masked sequence can be used with '.masked' extension)
-ref_mt_gene            (Required)  file name of reference mitochondrial gene for gene prediction
                                    (The data is available from NCBI. The file extension is '.gb')
-gene_annotation        (option)    the final folder name of the assembled mtDNA
                                    no: default
-num_threads            (option)    number of threads for the pipeline 
                                    (default: 20)
-min_identity           (option)    minimum sequence identity between ont-reads and reference mt genome sequence 
                                    (default: 0.90)
-estimated_mt_size      (option)    estimated mitochondrial genome size
                                    (default: 600000)
________________________________________________________________________________________________________________________________
"""
        if line1=="":
if args==[]:
    print (help_description)
    quit()

for i in range(len(args)):
    if args[i].startswith("-"):
        try:
            option_dict[args[i]]=args[i+1]
        except:
            if args[0]=="-help" or args[0]=="":
                print (help_description)
                quit()

import os

temp_file_list=os.listdir(".")

#uncompress the compressed long-read seq_________________________________________________________________________________________________
print ("Check long-read seq..\n")
if option_dict['-long_read'][-3:]==".gz":
    if option_dict['-long_read'][:-3] not in temp_file_list:
        print ("1/19. Uncompressing the long-read gz file..\n")
        gunzip=os.system("gunzip %s" % option_dict['-long_read'])
        if gunzip==0:
            option_dict['-long_read']=option_dict['-long_read'][:-3]
        else:
            print ("-Error!-: Uncompressing raised an error!\n")
            quit()
    else:
        option_dict['-long_read']=option_dict['-long_read'][:-3]
        print ("1/19. Uncompressed long-read file exists. The pipeline will use the existing file for the next step..\n")
else:
    print ("1/19. The long-read file is uncompressed format. The pipeline will use the existing file for the next step..\n")

# filtering min-len long-read sequence_________________________________________________________________________________________________
if option_dict['-long_read']+"_min%s.fasta" % option_dict['-min_ont_read_len'] not in temp_file_list and option_dict['-long_read']+"_min%s.fastq" % option_dict['-min_ont_read_len'] not in temp_file_list:
    if option_dict['-long_read'][-5:]=="fastq":
        print ("2/19. Filtering out long-read sequences under %s bp.. \n" % option_dict['-min_ont_read_len'])
        filtering_len=os.system("python /workdir/python_codes/ont_fastq_size_filtering_v1.0.py3.py -fastq %s -out_type 3 -min_len %s" %(option_dict['-long_read'], option_dict['-min_ont_read_len']))

        filtered_fasta_name=option_dict['-long_read']+"_min%s.fasta" % option_dict['-min_ont_read_len']
        filtered_fastq_name=option_dict['-long_read']+"_min%s.fastq" % option_dict['-min_ont_read_len']
    else:
        print ("-Error!-: the extension of the long-read file is not 'fastq'. Please check and try again..\n")
        quit()
else:
    filtered_fasta_name=option_dict['-long_read']+"_min%s.fasta" % option_dict['-min_ont_read_len']
    filtered_fastq_name=option_dict['-long_read']+"_min%s.fastq" % option_dict['-min_ont_read_len']
    filtering_len=0
    print ("2/19. Filtered min-len long-read sequence already exists. The pipeline will use the existing files..\n")

# make mt-ref blast db_________________________________________________________________________________________________
if option_dict['-ref_mt_genome']+".ndb" not in temp_file_list:
    if filtering_len==0:
        print ("3/19. Generating blastdb of reference mtDNA sequence..")
        makeblastdb=os.system("makeblastdb -dbtype nucl -in %s" % option_dict['-ref_mt_genome'])
    else:
        print ("-Error!-: Filtering min-len long-read sequence raised an error!\n")
        quit()
else:
    makeblastdb=0
    print ("3/19. blastdb of the %s already exists. The pipeline will use the existing blastdb for the next step..\n" % option_dict['-ref_mt_genome'] )

# run blastn to isolate mt-dna containing long-read seqs_________________________________________________________________________________________________
if filtered_fasta_name[:-6]+"-mtRef.xml" not in temp_file_list:
    if makeblastdb==0:
        print ("4/19. Run blast against reference mtDNA..\n")
        run_blast=os.system("blastn -outfmt 5 -num_threads %s -max_target_seqs 1 -db %s -query %s -out %s" % (option_dict['-num_threads'], option_dict['-ref_mt_genome'], filtered_fasta_name, filtered_fasta_name[:-6]+"-mtRef.xml"))
    else:
        print ("-Error!-: makeblastdb for reference mtDNA seq raised an error!\n")
        quit()
else:
    run_blast=0
    print ("4/19. The blast xml file for mtDNA isolation already exists. The pipeline will use the existing file for the next step..\n")

# parsing reference mtDNA blast result_________________________________________________________________________________________________
if filtered_fasta_name[:-6]+"-mtRef.csv" not in temp_file_list:
    if run_blast==0:
        print ("\n5/19. Parsing blast results..\n")
        pwd=os.getcwd()
        parsing_ref_blast=os.system("python /workdir/python_codes/0_blast_parsing___xml_format_v4.0_multi-threads.py3.py -xml %s -seq_parse 0 -num_hit 1 -num_hsp 1 -lim_score 200 -ov_len 1 -num_threads %s" % (pwd+"/"+filtered_fasta_name[:-6]+"-mtRef.xml", option_dict['-num_threads']))
    else:
        print ("-Error!-: blast search against reference mtDNA raised an error!\n")
        quit()
else:
    parsing_ref_blast=0
    print ("5/19. The parsed blast result for mtDNA isolation exists. The pipeline will use the existing file for the next step..\n")

# parsing mtDNA-containing long-reads_________________________________________________________________________________________________
if filtered_fasta_name[:-6]+"-mtReads.fasta" not in temp_file_list:
    if parsing_ref_blast==0:
        print ("6/19. Parsing mtDNA-containing long-reads..\n")
        parsing_mtDNA_long_reads=os.system("python /workdir/python_codes/csv_query_or_hit_seq_parsing_v1.0.py3.py -csv %s -seq %s -seq_type 2 -target 1 -method 1 -out %s " % (filtered_fasta_name[:-6]+"-mtRef.csv", filtered_fasta_name, filtered_fasta_name[:-6]+"-mtReads.fasta"))
    else:
        print ("-Error!-: Parsing mtDNA-containing long-read sequence raised an error!\n")
        quit()
else:
    parsing_mtDNA_long_reads=0
    print ("6/19. The parsed mtDNA-containing long-read file already exists. The pipeline will use the existing file for the next step..\n")

# masking reference genome sequence_________________________________________________________________________________________________
if option_dict['-ref_nucl_genome'][-7:]!=".masked":
    if option_dict['-ref_nucl_genome']+".masked" not in temp_file_list:
        print ("\n7/19. Masking reference genome sequence with mtDNA..\n")
        mtDNA_masking=os.system("/workdir/RepeatMasker/RepeatMasker -e rmblast -norna -no_is -pa %s -lib %s %s" % (option_dict['-num_threads'], option_dict['-ref_mt_genome'], option_dict['-ref_nucl_genome']))
        if mtDNA_masking==0:
            option_dict['-ref_nucl_genome']=option_dict['-ref_nucl_genome']+".masked"
        else:
            print ("-Error!-: Reference genome masking by mtDNA raised error!\n")
            quit()
    else:
        mtDNA_masking=0
        option_dict['-ref_nucl_genome']=option_dict['-ref_nucl_genome']+".masked"
        print ("7/19. Masked reference genome sequence already exists. The pipeline will use the existing file for the next step..\n")
else:
    print ("7/19. Masked reference genome sequence already exists. The pipeline will use the existing file for the next step..\n")

# make mtDNA-masked reference genome seq blastdb_________________________________________________________________________________________________
if option_dict['-ref_nucl_genome']+".ndb" not in temp_file_list:
    if mtDNA_masking==0:
        print ("8/19. Make blastdb of mt-masked reference genome sequence..\n")
        make_blastdb=os.system("makeblastdb -dbtype nucl -in %s" % option_dict['-ref_nucl_genome'])
    else:
        print ("-Errror!-: Masking reference genome by mtDNA raised an error.\n")
        quit()
else:
    make_blastdb=0
    print ("8/19. blastdb of mt-masked reference genome already exists. The pipeline will use the existing file for the next step..\n")

# run blast to filter nuclDNA containing long-reads_________________________________________________________________________________________________
if filtered_fasta_name[:-6]+"-mtReads_nuclRead-filtering.xml" not in temp_file_list:
    if parsing_mtDNA_long_reads==0:
        print ("9/19. Run blast to filter out nuclDNA-containing mtDNA long-reads..\n")
        run_blast_filtering_nuclDNA=os.system("blastn -outfmt 5 -num_threads %s -max_target_seqs 1 -db %s -query %s -out %s" % (option_dict['-num_threads'], option_dict['-ref_nucl_genome'], filtered_fasta_name[:-6]+"-mtReads.fasta", filtered_fasta_name[:-6]+"-mtReads_nuclRead-filtering.xml"))
    else:
        print ("-Error!-: Parsing mtDNA-containing long reads raised an error!")
        quit()
else:
    run_blast_filtering_nuclDNA=0
    print ("9/19. The blast xml file for filtering nuclDNA-containing long reads already exists. The pipeline will use the existing file for the next step..\n")

# parsing blast results_________________________________________________________________________________________________
if filtered_fasta_name[:-6]+"-mtReads_nuclRead-filtering.csv" not in temp_file_list:
    if run_blast_filtering_nuclDNA==0:
        print ("\n10/19. Parsing blast results to filter out nuclDNA containing long-reads.\n")
        pwd=os.getcwd()
        parsing_nuclDNA_blast=os.system("python /workdir/python_codes/0_blast_parsing___xml_format_v4.0_multi-threads.py3.py -xml %s -seq_parse 0 -num_hit 1 -num_hsp 1 -lim_score 600 -ov_len 1 -num_threads %s" % (pwd+"/"+filtered_fasta_name[:-6]+"-mtReads_nuclRead-filtering.xml", option_dict['-num_threads']))
    else:
        print ("-Error!-: Blast search to filter nuclDNA raised an error!\n")
        quit()
else:
    parsing_nuclDNA_blast=0
    print ("10/19. Parsed blast xml file for filtering nuclDNA-containing long reads already exists. The pipeline will use the existing file for the next step..\n")

# parsing nuclDNA filtered mtDNA long-reads_________________________________________________________________________________________________
if filtered_fasta_name[:-6]+"-mtReads_nuclRead-filtered.fasta" not in temp_file_list:
    if parsing_nuclDNA_blast==0:
        print ("11/19. Parsing nuclDNA-filtered mtDNA-containing long-reads..\n")
        pure_mtDNA_parsing=os.system("python /workdir/python_codes/csv_query_or_hit_seq_parsing_v1.0.py3.py -csv %s -seq %s -seq_type 1 -target 1 -method 2 -out %s" % (filtered_fasta_name[:-6]+"-mtReads_nuclRead-filtering.csv", filtered_fastq_name, filtered_fasta_name[:-6]+"-mtReads_nuclRead-filtered.fasta"))
    else:
        print ("-Error!-: Parsing blast results for filtering nuclDNA-containing long-reads raised an error!\n")
        quit()
else:
    pure_mtDNA_parsing=0
    print ("11/19. Parsed fastq file for nuclDNA-filtered mtDNA-containing long reads already exists. The pipeline will use the existing file for the next step..\n")

# gunzip illumina sequences_________________________________________________________________________________________________
if option_dict['-short_read1'][-3:]==".gz":
    print ("12-1/19. Uncompress short read sequence1..\n")
    os.system("gunzip %s" % option_dict['-short_read1'])
    option_dict['-short_read1']=option_dict['-short_read1'][:-3]
if option_dict['-short_read2']!="":
    if option_dict['-short_read2'][-3:]==".gz":
        print ("12-2/19. Uncompress short read sequence2..\n")
        os.system("gunzip %s" % option_dict['-short_read2'])
        option_dict['-short_read2']=option_dict['-short_read2'][:-3]

# run AdapterRemoval for the illumina sequences_________________________________________________________________________________________________
if option_dict['-short_read2']!="":
    if option_dict['-short_read1']+"_trimmed.fastq" not in temp_file_list or option_dict['-short_read2']+"_trimmed.fastq" not in temp_file_list:
        print ("13/19. AdapterRemoval is running for trimming low-quality base..\n")
        os.system("AdapterRemoval --threads %s --file1 %s --file2 %s --output1 %s --output2 %s" % (option_dict['-num_threads'], option_dict['-short_read1'], option_dict['-short_read2'], option_dict['-short_read1']+"_trimmed.fastq", option_dict['-short_read2']+"_trimmed.fastq"))
        option_dict['-short_read1']=option_dict['-short_read1']+"_trimmed.fastq"
        option_dict['-short_read2']=option_dict['-short_read2']+"_trimmed.fastq"
    else:
        option_dict['-short_read1']=option_dict['-short_read1']+"_trimmed.fastq"
        option_dict['-short_read2']=option_dict['-short_read2']+"_trimmed.fastq"
        print ("13/19. The trimmed short-read sequences already exist. The existing files will be used for the next step..\n")
else:
    if option_dict['-short_read1']+"_trimmed.fastq" not in temp_file_list:
        print ("13/19. AdapterRemoval is running for trimming low-quality base..\n")
        os.system("AdapterRemoval --threads %s --file1 %s --output1 %s" % (option_dict['-num_threads'], option_dict['-short_read1'], option_dict['-short_read1']+"_trimmed.fastq"))
        option_dict['-short_read1']=option_dict['-short_read1']+"_trimmed.fastq"
    else:
        option_dict['-short_read1']=option_dict['-short_read1']+"_trimmed.fastq"
        print ("13/19. The trimmed short-read sequence already exists. The existing files will be used for the next step..\n")

# make comp_msbwt.npy file for error correction of long-reads_________________________________________________________________________________________________
if option_dict['-short_read1']+"_comp_msbwt.npy" not in temp_file_list:
    if pure_mtDNA_parsing==0:
        print ("14/19. Generate 'comp_msbwt.npy' file for long-read error correction..\n")
        if option_dict['-short_read2']=="":
            comp_msbwt=os.system("msbwt2-build -o %s %s" %(option_dict['-short_read1']+"_comp_msbwt.npy", option_dict['-short_read1']))
        else:
            comp_msbwt=os.system("msbwt2-build -o %s %s %s" % (option_dict['-short_read1']+"_comp_msbwt.npy", option_dict['-short_read1'], option_dict['-short_read2']))
    else:
        print ("-Error!-: Parsing nuclDNA filtered long-reads raised an error!\n")
        quit()
else:
    comp_msbwt=0
    print ("14/19. 'comp_msbwt.npy' file for %s already exists. The pipeline will use the existing file for the next step..\n" % option_dict['-short_read1'] )

# Error correction of the long-reads_________________________________________________________________________________________________
if filtered_fasta_name[:-6]+"-mtReads_nuclRead-filtered_fmlrc2-corrected.fasta" not in temp_file_list:
    if comp_msbwt==0:
        print ("15/19. Error correction of the long-reads using fmlrc2..\n")
        err_corr=os.system("fmlrc2 -t %s %s %s %s" %(option_dict['-num_threads'],option_dict['-short_read1']+"_comp_msbwt.npy", filtered_fasta_name[:-6]+"-mtReads_nuclRead-filtered.fasta", filtered_fasta_name[:-6]+"-mtReads_nuclRead-filtered_fmlrc2-corrected.fasta"))
    else:
        print ("-Error!-: Generating 'comp_msbwt.npy' file raised an error!\n")
        quit()
else:
    err_corr=0
    print ("15/19. Error corrected long reads by FMLRC2 already exists. The pipeline will use the existing file for the next step..\n")

# Fasta filtering by identity threshold_________________________________________________________________________________________________
if filtered_fasta_name[:-6]+"-mtReads_nuclRead-filtered_fmlrc2-corrected_idty-%s.fasta" % option_dict['-min_identity'] not in temp_file_list:
    if err_corr==0:
        print ("16/19. Filtering corrected sequence by identity.. (threshold: %s)\n" % option_dict['-min_identity'])
        idty_filtering_blast=os.system("blastn -outfmt 5 -num_threads %s -max_target_seqs 5 -db %s -query %s -out %s" % (option_dict['-num_threads'], option_dict['-ref_mt_genome'], filtered_fasta_name[:-6]+"-mtReads_nuclRead-filtered_fmlrc2-corrected.fasta", filtered_fasta_name[:-6]+"-mtReads_nuclRead-filtered_fmlrc2-corrected.xml"))
        if idty_filtering_blast==0:
            idty_filtering_blast_parsing=os.system("python /workdir/python_codes/0_blast_parsing___xml_format_v4.0_multi-threads.py3.py -xml %s -seq_parse 0 -num_hit 1 -num_hsp 1 -lim_score 2000 -ov_len 1 -num_threads %s" % (filtered_fasta_name[:-6]+"-mtReads_nuclRead-filtered_fmlrc2-corrected.xml", option_dict['-num_threads']) )
        else:
            print ("-Error!-: Blast search for filtering by identity with mtDNA raised an error!\n")
            quit()
        if idty_filtering_blast_parsing==0:
            idty_filtering_seq_parsing=os.system("python /workdir/python_codes/csv_query_or_hit_seq_parsing_v1.0.py3.py -csv %s -seq %s -seq_type 2 -target 1 -method 1 -out %s -idty %s " % (filtered_fasta_name[:-6]+"-mtReads_nuclRead-filtered_fmlrc2-corrected.csv", filtered_fasta_name[:-6]+"-mtReads_nuclRead-filtered_fmlrc2-corrected.fasta", filtered_fasta_name[:-6]+"-mtReads_nuclRead-filtered_fmlrc2-corrected_idty-%s.fasta" % option_dict['-min_identity'], option_dict['-min_identity']))
        else:
            print ("-Error!-: Blast result parsing for filtering by identity with mtDNA raised an error!\n")
            quit()
    else:
        print ("-Error!-: Error correction by FMLRC2 raised an error!\n")
        quit()
else:
    idty_filtering_seq_parsing=0
    print ("16/19. The fasta file filtered by identity thresold already exists. The pipeline will use the existing file for the next step..\n")

# Fasta parsing by coverage from top-length reads_________________________________________________________________________________________________
coverage=10
assembly_info_list=[["Result folder", "Coverage", "Number of Contigs", "Total contig length", "Each contig length"]]
print ("17 and 18 steps will be iterated 15 times with 2x coverage interval from 10x to 40x coverage..\n\n")
for i in range(16):
    if filtered_fasta_name[:-6]+"-mtReads_nuclRead-filtered_fmlrc2-corrected_idty-%s_%sx.fasta" % (option_dict['-min_identity'], str(coverage)) not in temp_file_list:
        if idty_filtering_seq_parsing==0:
            print ("17/19. Isolate reads from top read length by coverage (%sx)..\n" % str(coverage) )
            tot_seq_size=int(option_dict['-estimated_mt_size'])*coverage
            coverage_parsing=os.system("python /workdir/python_codes/longest_fasta_seq_parser_by_tot_size_v1.0.py3.py -total_size %s -fasta_name %s -out %s" % (str(tot_seq_size), filtered_fasta_name[:-6]+"-mtReads_nuclRead-filtered_fmlrc2-corrected_idty-%s.fasta" % option_dict['-min_identity'], filtered_fasta_name[:-6]+"-mtReads_nuclRead-filtered_fmlrc2-corrected_idty-%s_%sx.fasta" % (option_dict['-min_identity'], str(coverage))))
        else:
            print ("-Error!-: Filtering sequence by identity threshold raised an error!\n")
            quit()
    else:
        coverage_parsing=0
        print ("17/19. The %sx coverage long-read sequence file already exists. The pipeline will use the existing file for the next step..\n" % str(coverage))


    # Assemble the error corrected long reads by NextDenovo_________________________________________________________________________________________________
    os.system("ls %s > %s" % (
        filtered_fasta_name[:-6]+"-mtReads_nuclRead-filtered_fmlrc2-corrected_idty-%s_%sx.fasta" % (option_dict['-min_identity'], str(coverage)), 
        filtered_fasta_name[:-6]+"_idty-%s_%sx_NextDenovo.fofn" % (option_dict['-min_identity'], str(coverage))
        ))
    run_cont="""
    [General]
    job_type = local # local, slurm, sge, pbs, lsf
    job_prefix = nextDenovo
    task = assemble # all, correct, assemble
    rewrite = yes # yes/no
    deltmp = yes 
    parallel_jobs = %s # number of tasks used to run in parallel
    input_type = corrected # raw, corrected
    read_type = ont # clr, ont, hifi
    input_fofn = %s 
    workdir = %s 

    [correct_option]
    read_cutoff = 1k
    #seed_cutoff=1000
    genome_size = %sk # estimated genome size
    sort_options = -m 20g -t 10
    minimap2_options_raw = -t 10
    pa_correction = 10 # number of corrected tasks used to run in parallel, each corrected task requires ~TOTAL_INPUT_BASES/4 bytes of memory usage.
    correction_options = -p 10

    [assemble_option]
    minimap2_options_cns = -t 10
    nextgraph_options = -a 1
    # see https://nextdenovo.readthedocs.io/en/latest/OPTION.html for a detailed introduction about all the parameters
    """ % (
        option_dict['-num_threads'], 
        filtered_fasta_name[:-6]+"_idty-%s_%sx_NextDenovo.fofn" % (option_dict['-min_identity'], str(coverage)), 
        filtered_fasta_name[:-6]+"_assembled_results_idty-%s_%sx" % (option_dict['-min_identity'], str(coverage)), 
        option_dict['-estimated_mt_size'][:-3]
        )

    run_cfg=open(filtered_fasta_name[:-6]+"-idty-%s_NextDenovo.cfg" % option_dict['-min_identity'],'w')
    run_cfg.write(run_cont)
    run_cfg.close()

    if filtered_fasta_name[:-6]+"_assembled_results_idty-%s_%sx" % (option_dict['-min_identity'], str(coverage)) not in temp_file_list:
        if coverage_parsing==0:
            print ("18/19. Assembly of the mtDNA-containing ONT-seq started..\n")
            assemble_mtDNA=os.system("/workdir/NextDenovo/nextDenovo %s" % (filtered_fasta_name[:-6]+"-idty-%s_NextDenovo.cfg" % option_dict['-min_identity']))
        else:
            print ("-Error!-: Parsing sequence by coverage raised an error!\n")
            quit()
    else:
        assemble_mtDNA=0
        print ("18/19. The %s folder is already exists." % (filtered_fasta_name[:-6]+"_assembled_results_idty-%s_%sx" % (option_dict['-min_identity'], str(coverage))))

    assembled_mtDNA=open(filtered_fasta_name[:-6]+"_assembled_results_idty-%s_%sx/03.ctg_graph/nd.asm.fasta" % (option_dict['-min_identity'], str(coverage)),'r')
    temp_assembly_info=[]
    while 1:
        line1=assembled_mtDNA.readline()
        line2=assembled_mtDNA.readline().strip()
        if line1=="":
            break
        temp_assembly_info.append(len(line2))
    cont=[filtered_fasta_name[:-6]+"_assembled_results_idty-%s_%sx" % (option_dict['-min_identity'], str(coverage)), str(len(temp_assembly_info)), str(sum(temp_assembly_info))]+[str(k) for k in temp_assembly_info]
    assembly_info_list.append(cont)
    coverage+=2

#write assembled contig info in a csv file
import csv
outcsv=csv.writer(open(filtered_fasta_name[:-6]+"_assembled_results_idty-%s_10-40x.csv" % (option_dict['-min_identity']),'w', newline=""))
for row in assembly_info_list:
    print ("\t".join(row))
    outcsv.writerow(row)


# Gene annotation of the assembled mitochondrial genome with MitoFinder _________________________________________________________________________________________________
if option_dict['-gene_annotation']!="no":
    if option_dict['-gene_annotation']+"__gene_annotation" not in temp_file_list:
        if assemble_mtDNA==0:
            print ("\n19/19. Gene annotation of the assembled mtDNA started..\n" )
            annotation=os.system("/workdir/MitoFinder/mitofinder -j %s -a %s -r %s -o 1 -p %s --new-genes --max-contig-size 700000" % (option_dict['-gene_annotation']+"__gene_annotation", option_dict['-gene_annotation']+"/03.ctg_graph/nd.asm.fasta", option_dict['-ref_mt_gene'], option_dict['-num_threads'] ))
            if annotation==0:
                print ("\nMitochondrial gene annotation has successfully completed. Check %s folder to see the results.\n" % (option_dict['-gene_annotation']+"__gene_annotation" ))
        else:
            print ("\n-Error!-: Assembly of mtDNA by NextDenovo raised an error!\n")
            quit()
    else:
        print ("\nThe folder for the gene annotation already exists. Remove the folder if you want to run gene annotation again.\n")