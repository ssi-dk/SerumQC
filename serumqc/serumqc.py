#!/usr/bin/env python3

import os
import re
import sys
import logging
import subprocess
import argparse
import tempfile
import configparser
import time
import shutil #apparently has a nice remove script for directory trees
import gzip
import zipfile
import inspect
import pathlib

sys.path.append(os.path.realpath(os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "lib")))
import serum
config = configparser.RawConfigParser()
config_path = os.path.realpath(os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "lib", "serum.config"))
#config_path = "/home/ssi.ad/kimn/git.repositories/ssi_scripts/serum/serum.config"
assert os.path.isfile(config_path)
config.read(config_path)

def program_initialization():   
    #Handle arguements
    parser = argparse.ArgumentParser(description='SerumQC, a batch pipeline for QC on reads by denovo assembling bacterial genomes')      
    parser.add_argument("-i", "--input_directory",
                        help="Give a working directory",
                        required=True)
    parser.add_argument("-o", "--output_directory",
                        help="Give a output directory",
                        required=True)
    parser.add_argument("-run", "--run_name",
                        help="Run name",
                        required=True)
    parser.add_argument("-g","--group",
                        help="group for accounting on computerome",
                        type=str,
                        default="")
    parser.add_argument("-ss", "--sample_sheet",
                        help="Sample sheeet in excel",
                        type=str,
                        default="")
    parser.add_argument("-fileids", "--identifier_file", 
                        help="Identifier file",
                        type=str,
                        default="")
    parser.add_argument("-t", "--threads",
                        help="Number of maximum concurrent threads",
                        type=int,
                        default=4)
    parser.add_argument("-mem", "--memory", 
                        help="Memory limit in Gb (where applicable)",
                        type=int,
                        default=12)
    parser.add_argument("-p", "--partition", 
                        help="Partition to use",
                        type=str,
                        default="standard")
    parser.add_argument("-qc", "--redo_qc_summary", 
                        help="redo_qc_summary",
                        action="store_true",
                        default=False)
    parser.add_argument("-clean", "--clean_temp_files", 
                        help="clean_temp_files",
                        action="store_true",
                        default=False)
    parser.add_argument("-keep", "--keep_temp_files", 
                        help="Keep temp files",
                        action="store_true",
                        default=False)
    parser.add_argument("-se","--send_email",
                        help="Send emails out",
                        action="store_true",
                        default=False),
    parser.add_argument("-e","--email_list",
                        help="Send emails out",
                        type=str,
                        default="")
    args = parser.parse_args()

    return args 

def run_qc_on_grid(args, argv):
    input_directory = os.path.realpath(os.path.expanduser(args.input_directory))
    output_directory = os.path.realpath(os.path.expanduser(args.output_directory))
    run_name = args.run_name
    sample_sheet_file = os.path.realpath(os.path.expanduser(args.sample_sheet))
    identifier_file = args.identifier_file
    threads = args.threads
    memory = args.memory
    partition = args.partition
    keep_temp_files = args.keep_temp_files
    group = args.group
    available_queues = config.get("categories","queues").split(",")
    email_default_to = config.get("email","default_to").split(";")
    send_email = args.send_email
    email_list = args.email_list.split(";")


    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    os.chdir(output_directory)
 
    serum.set_file_for_logger(run_name+"_SerumQC.log")
    serum.log_info("Current working directory: {}".format(os.getcwd()))
    serum.log_info("Run command: {}".format(" ".join(argv)))
    serum.script__create_run_summary_file([run_name+"_summary.tsv"])

    if identifier_file == "":
        if os.path.isdir(sample_sheet_file):
            serum.script__create_identifier_file(input_directory, run_name=run_name, output_directory=output_directory)
            identifier_file = "identifier_list.txt"
        else:
            shutil.copy(sample_sheet_file, os.path.join(output_directory,"sample_sheet.xlsx"))
            os.chmod(os.path.join(output_directory,"sample_sheet.xlsx"),0o664)
            sample_sheet_file = os.path.join(output_directory,"sample_sheet.xlsx")
            serum.script__create_identifier_file(input_directory, run_name=run_name, sample_sheet=sample_sheet_file, output_directory=output_directory)
            identifier_file = "identifier_list.txt"

    if os.path.isfile(identifier_file):
        identifier_dataframe = serum.get_from_file__pandas_dataframe([identifier_file])
        jobs_queued = []
        email_list.extend(email_default_to)

        for index, row in identifier_dataframe.iterrows():
            identifier = row['SampleID']
            experiment_name = row['ExperimentName']
            sample_output_directory = row['output_directory']
            emails = row["Initials"]
            if len(emails) > 0:
                email_list.extend([x.lower() for x in emails.split(";")])
                email_list = list(set(email_list))

            if experiment_name == run_name:
                priority = partition
                if "Priority" in identifier_dataframe:
                    if row["Priority"].lower() in available_queues:
                        priority = row["Priority"].lower()
                serum.log_info("Running sample {}".format(identifier))

                if not os.path.exists(sample_output_directory):
                    os.makedirs(sample_output_directory)

                sample_identifier_file = os.path.join(sample_output_directory, identifier+"_identifier.txt")
                identifier_dataframe.loc[[index]].to_csv(sample_identifier_file, sep='\t', index=False)

                lib_path = os.path.dirname(os.path.realpath(__file__))
                run_command = "serumqc.qc_sample('"+sample_identifier_file+"', '"+run_name+"', "+str(threads)+", "+str(memory)+", '"+priority+"', "+str(keep_temp_files)+")"
                job_command = (
                    "python3 -c \""
                        "import sys;"
                        "sys.path.append('{sys_path}');"
                        "import serumqc;"
                        "{run_command}\""
                    )
                variables = {"sys_path":lib_path,"run_command":run_command}
                run_command = job_command.format(**variables)

                job_id = serum.run_job(job_name = run_name+"_"+identifier, job_command = run_command, task = "", working_dir = os.getcwd(), log_file = run_name+"_SerumQC.log", partition = priority, cpu_requirement = threads, memory_requirement_in_Gb = memory, time_requirement = "03:00:00", group=group)
                jobs_queued.append(job_id)

        if len(jobs_queued) > 0:
            email_to = ""
            for email in email_list:
                if "@" in email:
                    email_to += email +";"
            email_to = email_to[:-1]
            pre_run_operations(email_to, run_name, input_directory, output_directory, send_email)

            run_command = "serumqc.post_run_operations('"+email_to+"', '"+run_name+"', '"+input_directory+"',"+str(send_email)+")"
            job_command = (
                "python3 -c \""
                    "import sys;"
                    "sys.path.append('{sys_path}');"
                    "import serumqc;"
                    "{run_command}\""
                )
            variables = {"sys_path":lib_path,"run_command":run_command}
            run_command = job_command.format(**variables)
            job_id = serum.run_job(job_name = run_name+"__post_operations", job_command = run_command, task = "", working_dir = os.getcwd(), log_file = run_name+"_SerumQC.log", partition = priority, cpu_requirement = threads, memory_requirement_in_Gb = memory, time_requirement = "03:00:00", dependent_job_ids = jobs_queued, dependency = "afterany", group=group)
    return 0

def pre_run_operations(email_to, run_name, email_read_directory, email_output_directory, send_email):
    if(send_email):
        email_read_directory = serum.convert_path_from_linux_to_msft(email_read_directory)
        email_output_directory = serum.convert_path_from_linux_to_msft(email_output_directory)
        serum.script__email_reads_available(email_to, run_name, email_read_directory, email_output_directory)
    return 0

def post_run_operations(email_to, run_name, email_read_directory, send_email):
    serum.script__convert_tsv_to_excel(run_name+"_summary.tsv",run_name+"_summary.xlsx")
    email_output_directory = os.path.realpath(run_name+"_summary.xlsx")
    if(send_email):
        email_read_directory = serum.convert_path_from_linux_to_msft(email_read_directory)
        email_output_directory = serum.convert_path_from_linux_to_msft(email_output_directory)
        serum.script__email_qc_complete(email_to, run_name, email_read_directory, email_output_directory)
    return 0

def qc_sample(identifier_file, run_name, threads=4, memory=12, partition='', keep_temp_files=False):
    #Want to adjust here to obtain ID list from excel file matching the run, perhaps a function which first makes identifier list

    identifier_dataframe = serum.get_from_file__pandas_dataframe([identifier_file])

    identifier = identifier_dataframe.loc[0]['SampleID'].strip()
    read_files = [identifier_dataframe.loc[0]['R1_location'], identifier_dataframe.loc[0]['R2_location']]

    output_directory = identifier_dataframe.loc[0]['output_directory']
    provided_species = identifier_dataframe.loc[0]['Organism']
    detected_species = ""
    priority = identifier_dataframe.loc[0]['Priority']
    email = identifier_dataframe.loc[0]['Initials']
    comments = identifier_dataframe.loc[0]['Comments']

    mlst_species = serum.get_from_file__mlst_species(provided_species)

    os.chdir(output_directory)

    print("Running sample: {}".format(identifier))

    log_file = "SerumQC.log"

    serum.set_file_for_logger(log_file)
    #done here so that each sample can rebuild from log if needed
    serum.log_info("Run name: {}".format(run_name))
    serum.log_version(os.path.realpath(os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", ".git")))
    serum.log_info("Working directory: {}".format(os.getcwd()))
    serum.log_info("Output directory: {}".format(output_directory))
    #serum.check_all_programs()

    to_remove = []

    #add check for if final file is here don't rerun so a check on qc file

    qc_file = "qc.txt"
    job_id_provided = False

    if(os.path.isfile(qc_file)):
        print("Sample {} already done, delete folder to rerun or delete {}_qc_files to continue from failure point".format(identifier, identifier))
        return 0

    try:
        #Prep reads for assembly
        step01 = serum.step__read_management(read_files, threads, memory)
        to_remove.extend(step01.get_temp_files())
        #Check reads for contaminants
        step02 = serum.step__kraken_on_reads(read_files, threads)
        to_remove.extend(step02.get_temp_files())
        #Grab detected species and use if needed
        detected_species = serum.get_from_file__first_occurance_of_pattern(["kraken_report_summary.txt"],"Species detected: (?P<detected_species>.*?)\n").strip()
        if mlst_species == "":
            serum.log_info("Species provided {} has has no mlst, attempting mlst with {}".format(provided_species, detected_species))
            mlst_species = serum.get_from_file__mlst_species(detected_species,"ncbi_species")
        #Assemble
        step03 = serum.step__assembler(["normalized_R1.fastq.gz","normalized_R2.fastq.gz"], threads, memory)
        to_remove.extend(step03.get_temp_files())
        #Map reads to contigs
        step04 = serum.step__denovo_mapping(["spades_contigs.fasta","trimmed_R1.fastq.gz","trimmed_R2.fastq.gz"], threads)
        to_remove.extend(step04.get_temp_files())
        #Read correction on contigs from mapping
        step05 = serum.step__contig_read_correction(["spades_contigs.fasta","spades_contigs_elprep.bam"], threads, memory)
        to_remove.extend(step05.get_temp_files())
        #Coverage and stats
        step06 = serum.step__clean_final_contigs(["pilon.fasta","coverage_summary.txt"],identifier,[identifier+"_contigs.fasta"])
        to_remove.extend(step06.get_temp_files())
        #Check contigs for contaminants
        step07 = serum.step__kraken_on_contigs([identifier+"_contigs.fasta"])
        to_remove.extend(step07.get_temp_files())

        #Analysis steps
        #Surveillance (species relative)
        step0A = serum.program__mlst(["spades_contigs.fasta"],mlst_species)
        to_remove.extend(step0A.get_temp_files())
        #Run finders on reads
        step0B = serum.step__ariba_finders(read_files, threads)
        to_remove.extend(step0B.get_temp_files())
        #Run additional scripts as described in DB
        step0C = serum.step__species_specific_analysis(detected_species)

    except Exception as e:
        print(str(e))

    else:
        print("completed successfully")
        if(not keep_temp_files):
            serum.remove_files_and_folders(to_remove)

    finally:
        stepQC = serum.step__generate_sample_QC([identifier_file, os.path.join(output_directory, "..", run_name+"_summary.tsv".format(identifier))])

    return 0

def redo_qc_summary(args,argv):
    run_name = args.run_name
    output_directory = os.path.realpath(args.output_directory)
    redo_qc_summary = os.path.join(output_directory,run_name+"_summary_redo.tsv")
    serum.script__create_run_summary_file([redo_qc_summary])
    for sample in os.listdir(output_directory):
        sample_directory = (os.path.join(output_directory,sample))
        if os.path.isdir(sample_directory):
            os.chdir(sample_directory)
            serum.step__generate_sample_QC([sample+"_identifier.txt",redo_qc_summary], output_files = ["qc_redo.txt"])
            os.chdir(output_directory)
    serum.script__convert_tsv_to_excel(redo_qc_summary,redo_qc_summary.replace(".tsv",".xlsx"))
    return 0

def clean_temp_files(args,argv):
    potential_temp_files = []
    potential_temp_files.extend(["trimmed_R1.fastq.gz","trimmed_R1_unpaired.fastq.gz","trimmed_R2.fastq.gz","trimmed_R2_unpaired.fastq.gz","normalized_R1.fastq.gz","normalized_R2.fastq.gz"])
    potential_temp_files.extend(["kraken.txt"])
    potential_temp_files.extend(["spades"])
    potential_temp_files.extend(["spades_contigs.fasta.amb","spades_contigs.fasta.bwt","spades_contigs.fasta.sa","spades_contigs.fasta.ann","spades_contigs.fasta.pac","spades_contigs.sam","spades_contigs_elprep.sam","spades_contigs_elprep.bam","spades_contigs_elprep.bam.bai","spades_contigs.cov"])
    potential_temp_files.extend(["pilon.vcf"])
    potential_temp_files.extend(["contigs_kraken.txt"])

    output_directory = os.path.realpath(args.output_directory)
    for sample in os.listdir(output_directory):
        sample_directory = (os.path.join(output_directory,sample))
        if os.path.isdir(sample_directory):
            os.chdir(sample_directory)
            serum.remove_files_and_folders(potential_temp_files)
            os.chdir(output_directory)
    return 0

def main(argv):
    args = program_initialization()
    if args.redo_qc_summary:
        redo_qc_summary(args,argv)
    elif args.clean_temp_files:
        clean_temp_files(args,argv)
    else:
        run_qc_on_grid(args, argv)
    return 0

if __name__ == "__main__":   
    main(sys.argv)