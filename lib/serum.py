#!/usr/bin/env python3
""" 
library of common functions used in the various pipelines
"""
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
import time
import collections
import pandas 
import inspect
import Bio.SeqIO
import statistics
import vcf #for script__summarize_vcf
import smtplib
import email.mime.text
import datetime
import collections #to avoid bug with statistics mode

#obtain program paths from config file
config = configparser.RawConfigParser()
lib_path = os.path.dirname(os.path.realpath(__file__))
config_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),"serum.config")
assert os.path.isfile(config_path)
config.read(config_path)

#requires use of a logger, simple one is created though a fancier one can be set and is set for QC pipeline
global logger #root logger
logger = logging.getLogger()
logFormatter = logging.Formatter(fmt = '%(asctime)s - %(levelname)-9s- %(message)s', datefmt = '%y-%m-%d %H:%M:%S')
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
logger.addHandler(consoleHandler)
logger.setLevel(logging.DEBUG)

global encoding
encoding = 'utf-8'

global SUCCESS, FAILURE, RUNNING
SUCCESS = 0
FAILURE = -1
RUNNING = 1

def resolve_config_path(category, path):
    stored_path = config.get(category,path)
    stored_path = stored_path.replace("<serum>",os.path.dirname(os.path.realpath(__file__))+"/..")
    stored_path = os.path.expanduser(stored_path)
    return stored_path

global GLOBAL_spades_program_path
GLOBAL_spades_program_path = resolve_config_path("programs","spades_path")
#global GLOBAL_quast_program_path
#GLOBAL_quast_program_path = resolve_config_path("programs","quast_path")
global GLOBAL_bbduk_program_path
GLOBAL_bbduk_program_path = resolve_config_path("programs","bbduk_path")
global GLOBAL_bbmerge_program_path
GLOBAL_bbmerge_program_path = resolve_config_path("programs","bbmerge_path")
global GLOBAL_bbnorm_program_path
GLOBAL_bbnorm_program_path = resolve_config_path("programs","bbnorm_path")
global GLOBAL_kraken_program_path
GLOBAL_kraken_program_path = resolve_config_path("programs","kraken_path")
global GLOBAL_kraken_report_program_path
GLOBAL_kraken_report_program_path = resolve_config_path("programs","kraken_report_path")
global GLOBAL_mlst_program_path
GLOBAL_mlst_program_path = resolve_config_path("programs","mlst_path")
global GLOBAL_trimmomatic_program_path
GLOBAL_trimmomatic_program_path = resolve_config_path("programs","trimmomatic")
global GLOBAL_ariba_program_path
GLOBAL_ariba_program_path = resolve_config_path("programs","ariba")
global GLOBAL_bwa_program_path
GLOBAL_bwa_program_path = resolve_config_path("programs","bwa")
global GLOBAL_elprep_program_path
GLOBAL_elprep_program_path = resolve_config_path("programs","elprep")
global GLOBAL_samtools_program_path
GLOBAL_samtools_program_path = resolve_config_path("programs","samtools")
global GLOBAL_pilon_jar_program_path
GLOBAL_pilon_jar_program_path = resolve_config_path("programs","pilon_jar")

global GLOBAL_kraken_db_path
GLOBAL_kraken_db_path = resolve_config_path("files","kraken_db")
global GLOBAL_species_db_path
GLOBAL_species_db_path = resolve_config_path("files","species_db")


class _function_helper(object):
    """
    Function helper used to normalize script and program functions to check for 
    inputs and whether the function needs to be run at all, also contains log
    functions.
    """
    def __init__(self, name, logger=logger):
        self.name = name
        self.input_files = []
        self.output_files = []
        self.temp_files = []
        self.program_command = ""
        self.logger = logger
        self.process_out = ""
        self.process_err = ""
    def start(self):
        self.logger.info("Function \'{}\' is starting".format(self.name))
        if self.output_files_exist() and len(self.output_files) > 0:
            self.logger.info("Function \'{}\' not ran as output_files {} already exist".format(self.name, self.output_files))
            return SUCCESS
        return RUNNING
    def completed(self):
        if self.output_files_exist():
            self.logger.info("Function \'{}\' has finished with output_files {}".format(self.name, self.output_files))
            return SUCCESS
        elif len(self.output_files) == 0:
            self.logger.info("Function \'{}\' has finished".format(self.name)) 
            return SUCCESS
        else:
            self.logger.warn("Function \'{}\' has finished but did not create expected output".format(self.name))
            return FAILURE
    def set_program_command(self, program_command):
        self.program_command = program_command
    def run_program(self):
        self.logger.info("Running command: {}".format(self.program_command))
        program_as_args = self.program_command.split()
        process = subprocess.Popen(program_as_args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        self.process_out, self.process_err = process.communicate()
        return (self.process_out, self.process_err)
    def get_stdout(self):
        return self.process_out.decode(encoding)
    def get_stderr(self):
        return self.process_err.decode(encoding)
    def set_output_files(self, output_files):
        self.output_files = output_files
    def set_temp_files(self, temp_files):
        self.temp_files = temp_files
    def get_temp_files(self):
        return self.temp_files
    def output_files_exist(self):
        for output_file in self.output_files:
            if not os.path.isfile(output_file):
                return False
        return True
    def exception_handler(self, e):
        self.logger.error("Function \"{}\" had an error. Error message:{}".format(self.name,str(e)))
        return FAILURE
class _Job_handler(object):
    """
    Job handler which helps manage executing jobs for the grid engine
    """
    def __init__(self, job_name = "", job_command = "", task = "", working_dir = "", log_file = "", dependent_job_arg = "", dependent_job_ids = [], dependency = "afterok", partition_arg = "", partition = "", cpu_requirement = 1, memory_requirement_in_Gb = 1, time_requirement = "02:00:00", group=""):
        self.job_name = job_name
        self.job_command = job_command
        self.task = task
        self.working_dir = working_dir
        self.log_file = log_file
        self.dependent_job_arg = dependent_job_arg
        self.dependent_job_ids = dependent_job_ids
        self.dependency = dependency
        self.partition_arg = partition_arg
        self.partition = partition
        self.cpu_requirement = cpu_requirement
        self.memory_requirement_in_Gb = memory_requirement_in_Gb
        self.time_requirement = time_requirement
        self.group = group
        self.job_id = ""
        self.grid_system = config.get("values","grid_system")

    #Done so I can create a basic instance then add other fields
    def set_working_dir(self, working_dir):
        self.working_dir = working_dir

    def set_log_file(self, log_file):
        self.log_file = log_file

    def set_cpu_requirement(self, cpu_requirement):
        self.cpu_requirement = cpu_requirement

    def set_memory_requirement_in_Gb(self, memory_requirement_in_Gb):
        self.memory_requirement_in_Gb = memory_requirement_in_Gb

    def set_partition(self, partition):
        self.partition = partition

    def set_dependency_info(self, dependent_job_ids, dependency = "afterok"):
        self.dependent_job_ids = dependent_job_ids
        self.dependency = dependency

    def set_time_requirements(self, time_requirement):
        self.time_requirement = time_requirement

    def set_job_command(self, job_command, task):
        self.job_command = job_command
        self.task = task

    def set_job_name(self, job_name):
        self.job_name = job_name

    def generate_dependent_job_arg(self):
        if self.grid_system == "slurm":
            self.dependent_job_arg = "--dependency={dependency}:{job_id}".format(dependency = self.dependency,job_id = ":".join(map(str,self.dependent_job_ids)))
        elif self.grid_system == "torque":
            self.dependent_job_arg = "-W depend={dependency}:{job_id}".format(dependency = self.dependency,job_id = ":".join(map(str,self.dependent_job_ids)))

    def submit_job(self):
        if self.job_name == "":
            # self.job_name = "j{run_name}_{identifier}__{task}".format(run_name = self.run_name, identifier = self.identifier,task = self.task)
            self.job_name = "jQC"
        if self.partition != "" and self.partition_arg == "":
            self.partition_arg = "-p {}".format(self.partition)
        if len(self.dependent_job_ids) > 0 and self.dependent_job_arg == "":
            if self.grid_system == "slurm":
                self.dependent_job_arg = "--dependency={dependency}:{job_id}".format(dependency = self.dependency,job_id = ":".join(map(str,self.dependent_job_ids)))
        self.job_id = submit_slurm_job(self)
        return self.job_id

    def check_job_submission(self):
        if self.job_id == "-1":
            return Exception("Job failed to submit")
        elif len(self.job_id) > 0:
            return self.job_id
        else:
            return Exception("Job not submitted")

    def set_function_as_job(self, function_name, function_params):
        run_command = "serum.{}({})".format(function_name,str(function_params)[1:len(str(function_params))-1])
        job_command = (
            "python3 -c \""
                "import sys;"
                "sys.path.append('{sys_path}');"
                "import serum;"
                "serum.set_file_for_logger('{log_file}');"
                "{run_command}\""
            )
        #should also be able to substitute more stuff
        variables = {"sys_path":lib_path,"log_file":self.log_file,"run_command":run_command}

        run_command = job_command.format(**variables)
        task = function_name

        self.set_job_command(run_command, task)

    def copy(self):
        return _Job_handler(job_name = self.job_name, job_command = self.job_command, task = self.task, working_dir = self.working_dir, log_file = self.log_file, dependent_job_arg = self.dependent_job_arg, dependent_job_ids = self.dependent_job_ids, dependency = self.dependency, partition_arg = self.partition_arg, partition = self.partition, cpu_requirement = self.cpu_requirement, memory_requirement_in_Gb = self.memory_requirement_in_Gb, time_requirement = self.time_requirement)

    def exception_handler(self,e):
        """
        Place to do exception handling
        """
        if self.job_name == "":
            self.job_name = "j{identifier}_{task}".format(identifier = self.identifier,task = self.task)
        logger.error("Job: \"{}\" had an error. Error message {}".format(self.job_name,str(e)))

#
# Small functions used by other functions or for initialization
#
def set_file_for_logger(log_file):
    logFormatter = logging.Formatter(fmt = '%(asctime)s - %(levelname)-9s- %(message)s', datefmt = '%y-%m-%d %H:%M:%S')
    fileHandler = logging.FileHandler(log_file)
    fileHandler.setFormatter(logFormatter)
    logger.addHandler(fileHandler)
    return 0
def log_info(message):
    logger.info("{}".format(message))
    return
def log_version(git_directory):
    try:
        args = ["git","--git-dir",git_directory,"rev-parse","origin/master"]
        process = subprocess.Popen(args,stdout = subprocess.PIPE,stderr = subprocess.STDOUT)
        process_out, process_err = process.communicate()
        logger.info("Git branch hash (version): {}".format(process_out.strip()))
    except:
        logger.error("pipeline version check failed")
        raise
    return True
def check_all_programs():
    logger.info("config: {}".format(config_path))
    check_program(GLOBAL_spades_program_path)
    logger.info("spades: {}".format(obtain_program_path(GLOBAL_spades_program_path)))
    # check_program(GLOBAL_quast_program_path)
    # logger.info("QUAST: {}".format(obtain_program_path(GLOBAL_quast_program_path)))
    check_program(GLOBAL_bbduk_program_path)
    logger.info("BBDuk: {}".format(obtain_program_path(GLOBAL_bbduk_program_path)))
    check_program(GLOBAL_bbmerge_program_path)
    logger.info("BBMerge: {}".format(obtain_program_path(GLOBAL_bbmerge_program_path)))
    check_program(GLOBAL_bbnorm_program_path)
    logger.info("BBNorm: {}".format(obtain_program_path(GLOBAL_bbnorm_program_path)))
    check_program(GLOBAL_kraken_program_path,["-h"])
    logger.info("kraken: {}".format(obtain_program_path(GLOBAL_kraken_program_path)))
    check_program(GLOBAL_kraken_report_program_path,["-h"])
    logger.info("kraken report: {}".format(obtain_program_path(GLOBAL_kraken_report_program_path)))
    #check_program(GLOBAL_samtools_program_path)
    #logger.info("samtools: {}".format(obtain_program_path(GLOBAL_samtools_program_path)))
    #check_program(GLOBAL_mlst_program_path)    #not done because mlst has non standard return codes
    logger.info("mlst: {}".format(obtain_program_path(GLOBAL_mlst_program_path)))

    return True
def obtain_program_path(program_path):
    #want the path as it'll contain version number
    try:
        path = os.getenv('PATH')
        for p in path.split(os.path.pathsep):
            p = os.path.join(p,program_path)
            if os.path.exists(p) and os.access(p, os.X_OK):
                return os.path.realpath(p)
    except Exception:
        logger.error("Error obtaining program path at {}".format(program_path))
        raise
    return ("NA")
def check_program(program_path, options = []):
    try:
        run_command = [program_path]
        run_command.extend(options)
        subprocess.check_call(run_command, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
    except Exception:
        logger.error("Program not found/working at {}".format(program_path))
        raise
    return True
def remove_files_and_folders(to_remove_list):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.start()

    try:
        for item_to_remove in to_remove_list:
            if os.path.isfile(item_to_remove):
                logger.info("Removing file:\t {}".format(item_to_remove))
                os.remove(item_to_remove)
            elif os.path.isdir(item_to_remove):
                logger.info("Removing directory:\t {}".format(item_to_remove))
                shutil.rmtree(item_to_remove)
            else:
                logger.info("Error item doesn't exist:\t {}".format(item_to_remove))

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def script__symlink_files(files_to_link):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.start()
    try:
        for input_file in files_to_link:
            #Assuming input file can be mix of directory and file
            file_symlink = os.path.join(os.path.basename(input_file))
            if os.path.isfile(input_file):
                if not os.path.isfile(file_symlink) and not os.path.islink(file_symlink):
                    os.symlink(os.path.realpath(input_file),os.path.basename(input_file))
            else:
                logger.error("Input file is not a file, not symlinking")
                raise Exception("no read files")

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def get_from_file__pandas_dataframe(input_files):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.start()

    dataframe = ""
    try:
        identifiers = []
        run_identifier_file = input_files[0]
        dataframe = pandas.DataFrame.from_csv(run_identifier_file, sep = '\t', index_col = None).fillna('NA')
        for column in dataframe:
            dataframe[column] = dataframe[column].astype(str)

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        return dataframe
def get_from_file__first_occurance_of_pattern(input_files, search_pattern):
    search_file = input_files[0]
    try:
        with open(search_file,"r") as input_file:
            lines = input_file.read()
            result = re.search(search_pattern, lines, re.MULTILINE)
            if result:
                return result.group(1)
            else:
                return "NA"
    except Exception as e:
        return "NA"
def get_from_file__mlst_species(provided_species, category = "organism"):
    try:
        species_database = pandas.DataFrame.from_csv(GLOBAL_species_db_path, sep = '\t', index_col = None).fillna('')
        #allow some wiggle room of non exact matches
        if category == "organism":
            species_database_species_index = species_database[species_database[category] == provided_species].index.tolist()
        else:
            if category in species_database:
                species_database_species_index = species_database[species_database[category] == provided_species].index.tolist()
        if(len(species_database_species_index) != 1):
            logger.warn("species_database issue (mlst species) {} on {}".format(provided_species, category))
            return ""
        else:
            logger.info("mlst species {} found".format(species_database.at[species_database_species_index[0],"mlst_species"]))
            return species_database.at[species_database_species_index[0],"mlst_species"]
    except Exception as e:
        logger.error("get_from_file__mlst_species {}".format(str(e)))
        return ""
def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False
def setIfNotNA(variable, value):
    if value == "NA":
        return variable
    else:
        return value
#
# Steps - groups of programs/scripts that represent a step in the QC pipeline
#
def step__read_management(input_files, threads, memory, output_files=["read_length.txt","trimmomatic.txt","normalized.log"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        read_files = input_files

        script__symlink_files(read_files)
        script__summarize_read_info(read_files)
        program__trimmomatic_read_trimming(read_files, threads)
        script__summarize_read_info(["trimmed_R1.fastq.gz","trimmed_R2.fastq.gz"])
        program__bbnorm_normalization(["trimmed_R1.fastq.gz","trimmed_R2.fastq.gz"], threads, memory)
        script__summarize_read_info(["normalized_R1.fastq.gz","normalized_R2.fastq.gz"])
        fh.set_temp_files(["trimmed_R1.fastq.gz","trimmed_R1_unpaired.fastq.gz","trimmed_R2.fastq.gz","trimmed_R2_unpaired.fastq.gz","normalized_R1.fastq.gz","normalized_R2.fastq.gz"])

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def step__kraken_on_reads(input_files, threads, output_files = ["kraken_report_summary.txt"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh
    
    try:
        read_files = input_files
        program__kraken_on_reads(read_files, threads)
        program__kraken_report(["kraken.txt"])
        script__summarize_kraken_report(["kraken_report.txt"])

        fh.set_temp_files(["kraken.txt"])

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def step__assembler(input_files, threads, memory, output_files = ["spades_contigs.fasta","spades.log"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        program__spades_assembler_only(input_files, threads, memory)
        fh.set_temp_files(["spades"])

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def step__denovo_mapping(input_files, threads, output_files = ["spades_contigs_elprep.bam"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        contigs_file = [input_files[0]]
        read_files = [input_files[1],input_files[2]]

        program__bwa_index(contigs_file)
        program__bwa_mem(contigs_file+read_files, threads)
        program__elprep(["spades_contigs.sam"], threads)
        program__samtools_sam_to_bam(["spades_contigs_elprep.sam"], threads)
        program__samtools_index_bam(["spades_contigs_elprep.bam"])
        program__samtools_calculate_coverage(["spades_contigs_elprep.bam"])
        script__get_coverage_and_length_from_cov(["spades_contigs.cov"])

        fh.set_temp_files(["spades_contigs.fasta.amb","spades_contigs.fasta.bwt","spades_contigs.fasta.sa","spades_contigs.fasta.ann","spades_contigs.fasta.pac","spades_contigs.sam","spades_contigs_elprep.sam","spades_contigs_elprep.bam","spades_contigs_elprep.bam.bai","spades_contigs.cov"])
    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def step__contig_read_correction(input_files, threads, memory, output_files = ["pilon_vcf_summary.txt"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        contigs_file = input_files[0]
        elprep_bam_file = input_files[1]
        program__pilon_with_vcf([contigs_file, elprep_bam_file], threads, memory)
        script__summarize_vcf(["pilon.vcf"])
        fh.set_temp_files(["pilon.vcf"])

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def step__kraken_on_contigs(input_files, output_files = ["contigs_kraken_report_summary.txt"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        contigs_file = input_files

        program__kraken_on_contigs(input_files)
        script__summarize_kraken_on_contigs(["contigs_kraken.txt"], output_files = ["contigs_kraken_report_summary.txt"])
        fh.set_temp_files(["contigs_kraken.txt"])

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def step__clean_final_contigs(input_files, identifier, output_files=["final_contigs.fasta"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        contigs_file = input_files[0]
        coverage_file = input_files[1]
        final_contigs_file = output_files[0]
        identifier = identifier

        script__remove_word_from_contigs("_pilon", contigs_file)
        script__append_contigs_with_bpcov([contigs_file, coverage_file], output_files = [final_contigs_file])
        script__remove_value_from_contigs("cov",final_contigs_file)
        script__update_contigs_length([final_contigs_file], output_files = [final_contigs_file])
        script__prepend_value_to_contigs(identifier+"_",final_contigs_file)
        script__contigs_check([final_contigs_file])
        script__filter_contigs([final_contigs_file],[final_contigs_file.replace("contigs","filtered_contigs")])

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def step__generate_sample_QC(input_files, output_files = ["qc.txt"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        identifier_file = input_files[0]
        run_summary_file = input_files[1]
        qc_file = output_files[0]
        script__qc_sample([identifier_file], output_files = [qc_file])
        script__append_run_summary_file([qc_file],[run_summary_file])
        master_summary_file = config.get("files","master_summary_file")
        if os.path.isfile(master_summary_file):
            script__append_run_summary_file([qc_file],[master_summary_file])

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def step__ariba_finders(input_files, threads, output_files = ["ariba_resfinder_report.tsv","ariba_virulencefinder_report.tsv","ariba_plasmidfinder_report.tsv"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        R1_files = input_files[0]
        R2_files = input_files[1]
        resfinder_report = output_files[0]
        virulencefinder_report = output_files[1]
        plasmidfinder_report = output_files[2]
        program__ariba_on_database([R1_files,R2_files], threads, "resfinder", output_files=[resfinder_report])
        program__ariba_on_database([R1_files,R2_files], threads, "virulencefinder", output_files=[virulencefinder_report])
        program__ariba_on_database([R1_files,R2_files], threads, "plasmidfinder", output_files=[plasmidfinder_report])

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def step__species_specific_analysis(provided_species, category="ncbi_species"):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.start()

    try:
        species_database = pandas.DataFrame.from_csv(GLOBAL_species_db_path, sep = '\t', index_col = None).fillna('')
        category = "ncbi_species"
        species_database_species_index = species_database[species_database[category] == provided_species].index.tolist()
        if(len(species_database_species_index) != 1):
            logger.warn("species_database issue (ncbi species) {} on {}".format(provided_species, "ncbi_species"))
        else:
            species_script = species_database.at[species_database_species_index[0],"additional_script"]
            if(len(species_script) > 0):
                logger.info("Will run script: {} for {}".format(species_script, provided_species))
                fh.set_program_command(species_script)
                fh.run_program()
            else:
                logger.info("No species specific script specified for {} ".format(provided_species))
    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def convert_path_from_linux_to_msft(path):
    if str(config.get("email","convert_paths_from_linux_to_msft"))=='True':
        replacement_path = config.get("email","replacement_path").split(";")
        path = path.replace(replacement_path[0],replacement_path[1]).replace("/","\\")
        return path
    else:
        return path
def script__email_reads_available(email_to, run_name, read_directory, output_directory):
    email_from = str(config.get("email","default_from"))
    email_subject = "SerumQC: {} reads available".format(run_name)
    email_template = (
        "Reads are available at:\n"
        "{read_directory} \n"
        "\n"
        "QC initiated at {current_time} and will be located at:\n"
        "{output_directory}\n"
        "\n"
        "{signature}"
        )
    context = {
        "read_directory":read_directory,
        "current_time":datetime.datetime.now().strftime('%Y-%m-%d'),
        "output_directory":output_directory,
        "signature":""
        }

    msg = email.mime.text.MIMEText(email_template.format(**context))
    msg['Subject'] = email_subject
    msg['From'] = email_from
    msg['To'] = email_to

    s = smtplib.SMTP('localhost')
    s.send_message(msg)
    s.quit()

def script__email_qc_complete(email_to, run_name, read_directory, qc_report_file):
    email_from = str(config.get("email","default_from"))
    email_subject = "SerumQC: {} QC complete".format(run_name)
    email_template = (
        "Reads are available at:\n"
        "{read_directory} \n"
        "\n"
        "QC finished with qc report:\n"
        "{qc_report_file}\n"
        "\n"
        "{signature}"
        )
    context = {
        "read_directory":read_directory,
        "qc_report_file":qc_report_file,
        "signature":""
        }

    msg = email.mime.text.MIMEText(email_template.format(**context))
    msg['Subject'] = email_subject
    msg['From'] = email_from
    msg['To'] = email_to

    s = smtplib.SMTP('localhost')
    s.send_message(msg)
    s.quit()

def script__convert_tsv_to_excel(tsv_file, excel_file):
    if os.path.isfile(tsv_file):
        df = get_from_file__pandas_dataframe([tsv_file])
        writer = pandas.ExcelWriter(excel_file)
        df.to_excel(writer,'Sheet1')
        writer.save()
    else:
        return FAILURE

#
# Sets of scripts and programs run through the pipeline
#
def script__create_identifier_file(input_directory, run_name="", sample_sheet="", output_files = ["identifier_list.txt"], output_directory = os.getcwd()):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh
    try:
        # io.update_status()

        illumina_run_directory = input_directory
        identifier_file = output_files[0]
        read_file_pattern = config.get("values","read_file_pattern")
        samples_to_ignore = config.get("values","samples_to_ignore").split(",")
        run_folder = output_directory
        sample_sheet = sample_sheet

        sample_sheet_data = None
        if sample_sheet == "":
            sample_sheet_data = pandas.DataFrame()
        if sample_sheet != "":
            try:
                sample_sheet_data = pandas.read_excel(sample_sheet)
            except FileNotFoundError as e:
                print("sample_sheet not provided")

        required_columns = config.get("categories","identifier_file").split(",")

        for column in required_columns:
            if column not in sample_sheet_data:
                sample_sheet_data[column] = ''
        column_names = sample_sheet_data.columns.values.tolist()

        identifier_list = []
        for read_file in sorted(os.listdir(illumina_run_directory)):
            result = re.match(read_file_pattern,read_file)
            if result:
                identifier = result.group("identifier")
                if identifier not in samples_to_ignore:
                    if identifier in sample_sheet_data["SampleID"].astype(str).values:
                        sample_sheet_index = sample_sheet_data[sample_sheet_data["SampleID"].astype(str) == identifier].index
                        if len(sample_sheet_index) != 1:
                            logger.error("{} not unique in sample_sheet for SampleID".format(identifier))
                            raise
                        else:
                            sample_sheet_index = sample_sheet_index[0]
                        if result.group("read_number") == "_R1":
                            sample_sheet_data.loc[sample_sheet_index,"R1_location"] = os.path.join(illumina_run_directory,result.group(0))
                        elif result.group("read_number") == "_R2":
                            sample_sheet_data.loc[sample_sheet_index,"R2_location"] = os.path.join(illumina_run_directory,result.group(0))
                    else:
                        column_values = [[""]*len(column_names)]
                        next_index = [len(sample_sheet_data.index)]
                        new_row = pandas.DataFrame(column_values,columns = column_names,index = next_index)
                        new_row["SampleID"] = identifier
                        new_row["ExperimentName"] = run_name
                        if result.group("read_number") == "_R1":
                            new_row["R1_location"] = os.path.join(illumina_run_directory,result.group(0))
                        elif result.group("read_number") == "_R2":
                            new_row["R2_location"] = os.path.join(illumina_run_directory,result.group(0))
                        sample_sheet_data = sample_sheet_data.append(new_row)
        for index in sample_sheet_data.index:
            sample_sheet_data.loc[index,"output_directory"] = os.path.join(run_folder,str(sample_sheet_data.loc[index,"SampleID"]))

        sample_sheet_data.to_csv(identifier_file,sep = "\t",index = False)
    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def script__summarize_read_info(input_files, output_files = ["read_length.txt"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.start()

    try:
        """
        Run program/function here
        """
        R1_read_file = input_files[0]
        R2_read_file = input_files[1]
        read_length_file = output_files[0]
        minimum_reads_for_assembly = int(config.get("values","minimum_reads_for_assembly"))

        output_exists = False
        if os.path.isfile(read_length_file):
            with open(read_length_file,"r") as output_file:
                for line in output_file:
                    if re.search(R1_read_file+"|"+R2_read_file,line):
                        fh.set_output_files(output_files)
                        if(fh.start()==SUCCESS):
                            fh.completed()
                            return fh

        R1_file_handler = gzip.open(R1_read_file,"rt")
        R1_reads = Bio.SeqIO.parse(R1_file_handler,"fastq")
        R1_sizes = [len(read) for read in R1_reads] 
        R2_file_handler = gzip.open(R2_read_file,"rt")
        R2_reads = Bio.SeqIO.parse(R2_file_handler,"fastq")
        R2_sizes = [len(read) for read in R2_reads]
        Bt_sizes = R1_sizes+R2_sizes

        R1_reads = len(R1_sizes)
        R2_reads = len(R2_sizes)
        Bt_reads = min(len(R1_sizes),len(R2_sizes))
        if R1_reads == 0:
            R1_sizes = [0]
        if R2_reads == 0:
            R2_sizes = [0]
        if Bt_reads == 0:
            Bt_sizes = [0]

        output_template = (
            "Input files: {R1_read_file},{R2_read_file}\n"
            "R1-\tReads: {R1_reads}\t Mean: {R1_mean}\t Mode: {R1_mode}\t Max: {R1_max}\t Min: {R1_min}\n" 
            "R2-\tReads: {R2_reads}\t Mean: {R2_mean}\t Mode: {R2_mode}\t Max: {R2_max}\t Min: {R2_min}\n" 
            "Both-\tReads: {Bt_reads}\t Mean: {Bt_mean}\t Mode: {Bt_mode}\t Max: {Bt_max}\t Min: {Bt_min}\n" 
            "\n"
            )

        context = {
            "R1_read_file":R1_read_file,"R2_read_file":R2_read_file,
            "R1_reads":R1_reads,"R1_mean":statistics.mean(R1_sizes),"R1_mode":collections.Counter(R1_sizes).most_common()[0][0],"R1_max":max(R1_sizes),"R1_min":min(R1_sizes),
            "R2_reads":R2_reads,"R2_mean":statistics.mean(R2_sizes),"R2_mode":collections.Counter(R2_sizes).most_common()[0][0],"R2_max":max(R2_sizes),"R2_min":min(R2_sizes),
            "Bt_reads":Bt_reads,"Bt_mean":statistics.mean(Bt_sizes),"Bt_mode":collections.Counter(Bt_sizes).most_common()[0][0],"Bt_max":max(Bt_sizes),"Bt_min":min(Bt_sizes)
            }

        with open(read_length_file,"a") as output_file:
            output_file.write(output_template.format(**context))

        if Bt_reads <= minimum_reads_for_assembly:
            raise Exception("not enough reads")

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def program__kraken_on_reads(input_files, threads, output_files = ["kraken.txt"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh
    try:
        db_path = GLOBAL_kraken_db_path
        program_path = GLOBAL_kraken_program_path
        R1_read_file = input_files[0]
        R2_read_file = input_files[1]
        kraken_file = output_files[0]
        threads = threads

        command_template = (
            "{program_path} --paired --threads {threads} --db {db_path} --output {kraken_file} {R1_read_file} {R2_read_file}"
            )
        context = {
            "program_path":program_path,
            "threads":threads,
            "db_path":db_path,
            "kraken_file":kraken_file,
            "R1_read_file":R1_read_file,
            "R2_read_file":R2_read_file
            }

        fh.set_program_command(command_template.format(**context))
        fh.run_program()

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def program__kraken_report(input_files, output_files = ["kraken_report.txt"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        db_path = GLOBAL_kraken_db_path
        program_path = GLOBAL_kraken_report_program_path
        kraken_file = input_files[0]
        kraken_report_file = output_files[0]

        command_template = (
            "{program_path} --db {db_path} {kraken_file}"
            )
        context = {
            "program_path":program_path,
            "db_path":db_path,
            "kraken_file":kraken_file
            }

        fh.set_program_command(command_template.format(**context))
        fh.run_program()
        with open(kraken_report_file,"w") as output_file:
            output_file.write(fh.get_stdout())

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def script__summarize_kraken_report(input_files, output_files = ["kraken_report_summary.txt"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        kraken_report_file = input_files[0]
        kraken_summary_report_file = output_files[0]
        taxonomic_level_space_shift = int(config.get("values","summarize_kraken_report__taxonomic_level_space_shift"))
        cutoff_threshold = float(config.get("values","summarize_kraken_report__cutoff_threshold"))
        kraken_report_pattern = re.compile("(?P<proportion>[0-9]+\.[0-9]+)\t[0-9]+\t[0-9]+\t(?P<classification_level>[A-Z])\t[0-9]+\t(?P<spaces>\s*)(?P<classification_name>[a-zA-Z]+.*)")
        number_of_spaces_for_genus = 0
        percent_unclassified = 0.0

        with open(kraken_summary_report_file,"w") as output_file:
            with open(kraken_report_file,"r") as input_file:
                for line in input_file:
                    result = re.search(kraken_report_pattern,line)
                    if(result):
                        if result.group("classification_level") == "U":
                            output_file.write(line)
                            percent_unclassified = float(result.group("proportion"))
                        if float(result.group("proportion")) > cutoff_threshold:
                            output_file.write(line)

            with open(kraken_report_file,"r") as kraken_report:
                going_down = True
                classifications_found = []
                previous_num_of_spaces = -1
                classification_level = ""
                saved_line = ""
                number_of_genuses = 0
                current_domain = "NA"
                for line in kraken_report:
                    result = re.search(kraken_report_pattern,line)
                    if(result):
                        if result.group("classification_level") == "D":
                            current_domain = result.group("classification_name")
                            number_of_spaces_for_genus = len(result.group("spaces")) + taxonomic_level_space_shift
                        proportion = float(result.group("proportion"))
                        if(proportion >= cutoff_threshold):
                            current_num_of_spaces = len(result.group("spaces"))
                            if current_num_of_spaces <= number_of_spaces_for_genus:
                                current_classification_level = result.group("classification_level")
                                if previous_num_of_spaces >= current_num_of_spaces: #has gone up a level or multiple levels
                                    going_down = False
                                    #print previous classification level
                                else:
                                    going_down = True
                                if going_down != True:
                                    # print going_down, line.strip(), saved_line.strip()
                                    classifications_found.append(re.search(kraken_report_pattern,saved_line).group("classification_level")+" "+re.search(kraken_report_pattern,saved_line).group("classification_name"))
                                    number_of_genuses += 1
                                    # print classification_level, saved_line
                                if current_num_of_spaces <= number_of_spaces_for_genus:
                                    classification_level = current_classification_level
                                    saved_line = line
                                previous_num_of_spaces = current_num_of_spaces
                number_of_genuses += 1
                result = re.search(kraken_report_pattern,saved_line)
                if(result):
                    classifications_found.append(result.group("classification_level")+" "+result.group("classification_name"))

                output_file.write("Genuses found: {}\n".format(",".join(classifications_found)))
                output_file.write("Number of genuses in sample: {}\n".format(len(classifications_found)))
                output_file.write("Percent unclassified: {}\n".format(percent_unclassified))

                if(len(classifications_found) != 1):
                    logger.warning("{} genuses found >{}%, possible comtamination".format(len(classifications_found),cutoff_threshold))

            with open(kraken_report_file,"r") as input_file:
                highest_proportion = 0
                most_common_species = "None"
                for line in input_file:
                    result = re.search(kraken_report_pattern,line)
                    if result:
                        if result.group("classification_level") == "S":
                            if float(result.group("proportion")) > highest_proportion:
                                highest_proportion = float(result.group("proportion"))
                                most_common_species = result.group("classification_name")

                    if most_common_species == "":
                        logger.error("No species found in sample")
                output_file.write("Species detected: {}\n".format(most_common_species))

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def program__kraken_on_contigs(input_files, output_files = ["contigs_kraken.txt"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        db_path = GLOBAL_kraken_db_path
        program_path = GLOBAL_kraken_program_path
        contigs_file = input_files[0]
        kraken_contigs_file = output_files[0]

        command_template = (
            "{program_path} --db {db_path} --output {kraken_contigs_file} {contigs_file}"
            )
        context = {
            "program_path":program_path,
            "db_path":db_path,
            "kraken_contigs_file":kraken_contigs_file,
            "contigs_file":contigs_file
            }

        fh.set_program_command(command_template.format(**context))
        fh.run_program()

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def script__summarize_kraken_on_contigs(input_files, output_files = ["contigs_kraken_report_summary.txt"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        db_path = GLOBAL_kraken_db_path
        kraken_contigs_file = input_files[0]
        kraken_name_file = os.path.join(db_path,"taxonomy/names.dmp")
        name_dict = {'0':"unclassified"}
        with open(kraken_name_file,"r") as name_file:
            for line in name_file:
                info = line.split("|")
                for i in range(0,len(info)):
                    info[i] = info[i].strip()
                if info[0] not in name_dict:
                    name_dict[info[0]] = info[1]
                if info[3] == "scientific name":
                    name_dict[info[0]] = info[1]

        with open(kraken_contigs_file,"r") as kraken_file, open(output_files[0], "w") as output_file:
            for line in kraken_file:
                info = line.split("\t")
                output_file.write("{}\t{}\n".format(info[1],name_dict[info[2]]))

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def program__bbnorm_normalization(input_files, threads, memory, output_files = ["normalized.log"], temp_files = ["normalized_R1.fastq.gz", "normalized_R2.fastq.gz"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        fh.set_temp_files(temp_files)

        program_path = GLOBAL_bbnorm_program_path
        R1_reads = input_files[0]
        R2_reads = input_files[1]
        R1_normalized_reads = temp_files[0]
        R2_normalized_reads = temp_files[1]
        normalized_log = output_files[0]
        target_coverage = int(config.get("values","bbnorm_normalization__target_coverage"))
        threads = threads
        memory = memory #in Gb

        command_template = (
            "{program_path} in={R1_reads} in2={R2_reads} out={R1_normalized_reads} out2={R2_normalized_reads} target={target_coverage} threads={threads} -Xmx{memory}g"
            )
        context = {
            "program_path":program_path,
            "R1_reads":R1_reads,
            "R2_reads":R2_reads,
            "R1_normalized_reads":R1_normalized_reads,
            "R2_normalized_reads":R2_normalized_reads,
            "target_coverage":target_coverage,
            "threads":threads,
            "memory":memory
            }

        fh.set_program_command(command_template.format(**context))
        fh.run_program()

        with open(normalized_log,"w") as output_file:
            output_file.write(fh.get_stderr())

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def program__trimmomatic_read_trimming(input_files, threads, output_files = ["trimmomatic.txt"], temp_files = ["trimmed_R1.fastq.gz","trimmed_R1_unpaired.fastq.gz","trimmed_R2.fastq.gz","trimmed_R2_unpaired.fastq.gz"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        fh.set_temp_files(temp_files)

        program_path = GLOBAL_trimmomatic_program_path
        R1_reads = input_files[0]
        R2_reads = input_files[1]
        R1_trimmed_reads = temp_files[0]
        R1_trimmed_reads_unpaired = temp_files[1]
        R2_trimmed_reads = temp_files[2]
        R2_trimmed_reads_unpaired = temp_files[3]
        trimmed_log = output_files[0]
        trim_quality = int(config.get("values","trimmomatic_trimming__trim_quality"))
        window_size = int(config.get("values","trimmomatic_trimming__window_size"))
        minimum_length = int(config.get("values","trimmomatic_trimming__minimum_length"))
        threads = threads

        command_template = (
            "{program_path} PE -threads {threads} {R1_reads} {R2_reads} {R1_trimmed_reads} {R1_trimmed_reads_unpaired} {R2_trimmed_reads} {R2_trimmed_reads_unpaired} LEADING:{trim_quality} TRAILING:{trim_quality} SLIDINGWINDOW:{window_size}:{trim_quality} MINLEN:{minimum_length}"
            )
        context = {
            "program_path":program_path,
            "threads":threads,
            "R1_reads":R1_reads,
            "R2_reads":R2_reads,
            "R1_trimmed_reads":R1_trimmed_reads,
            "R1_trimmed_reads_unpaired":R1_trimmed_reads_unpaired,
            "R2_trimmed_reads":R2_trimmed_reads,
            "R2_trimmed_reads_unpaired":R2_trimmed_reads_unpaired,
            "trim_quality":trim_quality,
            "window_size":window_size,
            "minimum_length":minimum_length
            }

        fh.set_program_command(command_template.format(**context))
        fh.run_program()
        with open(trimmed_log,"w") as output_file:
            output_file.write(fh.get_stderr())

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def program__bbduk_trimming(input_files, threads, output_files = ["trimmed.log"], temp_files = ["trimmed_R1.fastq.gz", "trimmed_R2.fastq.gz"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        fh.set_temp_files(temp_files)

        program_path = GLOBAL_bbduk_program_path
        R1_reads = input_files[0]
        R2_reads = input_files[1]
        R1_trimmed_reads = temp_files[0]
        R2_trimmed_reads = temp_files[1]
        trimmed_log = output_files[0]
        trim_sides = config.get("values","bbduk_trimming__trim_sides") 
        trim_quality = int(config.get("values","bbduk_trimming__trim_quality"))
        minimum_average_quality = int(config.get("values","bbduk_trimming__minimum_average_quality"))
        minimum_length = int(config.get("values","bbduk_trimming__minimum_length"))
        memory = 1 #in Gb
        threads = threads

        command_template = (
            "{program_path} in={R1_reads} in2={R2_reads} out={R1_trimmed_reads} out2={R2_trimmed_reads} qtrim={trim_sides} trimq={trim_quality} minavgquality={minimum_average_quality} minlength={minimum_length} threads={threads} -Xmx{memory}g "
            )
        context = {
            "program_path":program_path,
            "R1_reads":R1_reads,
            "R2_reads":R2_reads,
            "R1_trimmed_reads":R1_trimmed_reads,
            "R2_trimmed_reads":R2_trimmed_reads,
            "trim_sides":trim_sides,
            "trim_quality":trim_quality,
            "minimum_average_quality":minimum_average_quality,
            "minimum_length":minimum_length,
            "threads":threads,
            "memory":memory
            }

        fh.set_program_command(command_template.format(**context))
        fh.run_program()
        with open(trimmed_log,"w") as output_file:
            output_file.write(fh.get_stdout())

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def program__spades_assembler_only(input_files, threads, memory, output_files = ["spades_contigs.fasta","spades.log"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        program_path = GLOBAL_spades_program_path
        R1_reads = input_files[0]
        R2_reads = input_files[1]
        assembler_contigs = output_files[0]
        assembler_log = output_files[1]
        spades_output_directory = "spades"
        threads = threads #redundant I know but keeping temporary for formatting
        memory = memory #in Gb

        command_template = (
            "{program_path} -1 {R1_reads} -2 {R2_reads} --only-assembler -o {spades_output_directory} -t {threads} -m {memory}"
            )
        context = {
            "program_path":program_path,
            "R1_reads":R1_reads,
            "R2_reads":R2_reads,
            "spades_output_directory":spades_output_directory,
            "threads":threads,
            "memory":memory
            }

        fh.set_program_command(command_template.format(**context))
        fh.run_program()

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        shutil.copy(os.path.join(spades_output_directory,"contigs.fasta"),assembler_contigs)
        shutil.copy(os.path.join(spades_output_directory,"spades.log"),assembler_log)
        fh.completed()
        return fh
# def program__quast(input_files, threads, output_files = ["quast_report.txt","quast.cov"]):
#     fh = _function_helper(inspect.currentframe().f_code.co_name)
#     fh.set_output_files(output_files)
#     if(fh.start()==SUCCESS):
#         fh.completed()
#         return fh

#     try:
#         program_path = GLOBAL_quast_program_path
#         contigs_file = input_files[0]
#         reference_file = input_files[1]
#         bam_file = input_files[2]
#         quast_output_directory = "quast"
#         quast_report = output_files[0]
#         quast_coverage  = output_files[1]
#         threads = threads

#         command_template = (
#             "{program_path} -s --bam {bam_file} --fragmented -R {reference_file} -o {quast_output_directory} -t {threads} {contigs_file}"
#             )
#         context = {
#             "program_path":program_path,
#             "bam_file":bam_file,
#             "reference_file":reference_file,
#             "contigs_file":contigs_file,
#             "quast_output_directory":quast_output_directory,
#             "threads":threads
#             }

#         fh.set_program_command(command_template.format(**context))
#         fh.run_program()

#     except Exception as e:
#         fh.exception_handler(e)
#         raise
#     else:
#         shutil.copy(os.path.join(quast_output_directory,"report.txt"),quast_report)
#         directory, filename = os.path.split(reference_file)
#         shutil.copy(os.path.join(quast_output_directory,"structural_variations/{}".format(filename.replace(".fasta",".cov"))),quast_coverage)
#         fh.completed()
#         return fh
def program__nucmer(input_files):
    """
    TODO: code for function
    """
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        print("hi")
    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def script__nucmer_duplicate_parser(input_files):
    #IN PROGRESS
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        for line in nucmer_file:
            nucmer_contig_line_pattern = re.compile(config.get("values","contig_length_and_coverage_pattern"))
            nucmer_contig_line_pattern = re.compile(config.get("values","contig_length_and_coverage_pattern"))

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def program__bwa_index(input_files, output_files = ["spades_contigs.fasta.amb","spades_contigs.fasta.ann","spades_contigs.fasta.bwt","spades_contigs.fasta.pac","spades_contigs.fasta.sa"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        program_path = GLOBAL_bwa_program_path
        contigs_file = input_files[0]

        command_template = (
            "{program_path} index {contigs_file}"
            )
        context = {
            "program_path":program_path,
            "contigs_file":contigs_file
            }

        fh.set_program_command(command_template.format(**context))
        fh.run_program()

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def program__bwa_mem(input_files, threads, output_files = ["spades_contigs.sam"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        program_path = GLOBAL_bwa_program_path
        contigs_file = input_files[0]
        R1_trimmed_reads = input_files[1]
        R2_trimmed_reads = input_files[2]
        sam_file = output_files[0]
        threads = threads

        command_template = (
            "{program_path} mem -t {threads} {contigs_file} {R1_trimmed_reads} {R2_trimmed_reads}"
            )
        context = {
            "program_path":program_path,
            "threads":threads,
            "contigs_file":contigs_file,
            "R1_trimmed_reads":R1_trimmed_reads,
            "R2_trimmed_reads":R2_trimmed_reads
            }

        fh.set_program_command(command_template.format(**context))
        fh.run_program()
        with open(sam_file,"w") as output_file:
            output_file.write(fh.get_stdout())

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def program__elprep(input_files, threads, output_files = ["spades_contigs_elprep.sam","elprep.log"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        program_path = GLOBAL_elprep_program_path
        sam_file = input_files[0]
        elprep_sam_file = output_files[0]
        log_file = output_files[1]
        threads = threads

        command_template = (
            "{program_path} {sam_file} {elprep_sam_file} --nr-of-threads {threads} --sorting-order coordinate --mark-duplicates --clean-sam"
            )
        context = {
            "program_path":program_path,
            "sam_file":sam_file,
            "elprep_sam_file":elprep_sam_file,
            "threads":threads
            }

        fh.set_program_command(command_template.format(**context))
        fh.run_program()
        with open(log_file,"w") as output_file:
            output_file.write(fh.get_stderr())

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def program__samtools_sam_to_bam(input_files, threads, output_files = ["spades_contigs_elprep.bam"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        program_path = GLOBAL_samtools_program_path
        sam_file = input_files[0]
        bam_file = output_files[0]
        threads = threads

        command_template = (
            "{program_path} view -@ {threads} -b -o {bam_file} {sam_file}"
            )
        context = {
            "program_path":program_path,
            "threads":threads,
            "bam_file":bam_file,
            "sam_file":sam_file
            }

        fh.set_program_command(command_template.format(**context))
        fh.run_program()

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def program__samtools_index_bam(input_files, output_files = ["spades_contigs_elprep.bam.bai"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        program_path = GLOBAL_samtools_program_path
        bam_file = input_files[0]

        command_template = (
            "{program_path} index {bam_file}"
            )
        context = {
            "program_path":program_path,
            "bam_file":bam_file
            }

        fh.set_program_command(command_template.format(**context))
        fh.run_program()

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def program__samtools_calculate_coverage(input_files, output_files = ["spades_contigs.cov"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        program_path = GLOBAL_samtools_program_path
        bam_file = input_files[0]
        coverage_file = output_files[0]

        command_template = (
            "{program_path} depth {bam_file}"
            )
        context = {
            "program_path":program_path,
            "bam_file":bam_file
            }

        fh.set_program_command(command_template.format(**context))
        fh.run_program()
        with open(coverage_file,"w") as output_file:
             output_file.write(fh.get_stdout())

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def script__get_coverage_and_length_from_cov(input_files, output_files = ["coverage_summary.txt", "coverage_qa.txt"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        coverage_file = input_files[0]
        coverage_summary_file = output_files[0]
        coverage_qa_file = output_files[1]

        contigs_dict = {}
        coverage_range = [0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100]
        binned_coverage = {}
        for coverage_size in coverage_range:
            binned_coverage[coverage_size]=[0,0]
        contig_list = []
        total_length = 0 #not actually total length with quast output as it's already binned
        total_coverage = 0
        with open(coverage_file,"r") as sample_coverage_file:
            for line in sample_coverage_file:
                coverage_per_position = line.split("\t")
                if len(coverage_per_position) == 3:
                    current_key = coverage_per_position[0]
                    total_length = total_length+1
                    total_coverage = total_coverage+int(coverage_per_position[2])
                    if current_key not in contig_list:
                        contig_list.append(current_key)
                        contigs_dict[current_key] = [1,int(coverage_per_position[2])]
                    else:
                        contigs_dict[current_key][0] = contigs_dict[current_key][0] + 1
                        contigs_dict[current_key][1] = contigs_dict[current_key][1] + int(coverage_per_position[2])
                    #binned_coverage
                    for coverage_size in coverage_range:
                        if int(coverage_per_position[2]) > coverage_size:
                            binned_coverage[coverage_size][0] = binned_coverage[coverage_size][0] + 1
                            binned_coverage[coverage_size][1] = binned_coverage[coverage_size][1] + int(coverage_per_position[2])
                        else:
                            break

        with open(coverage_summary_file,"w") as coverage_summary:
            coverage_summary.write("contig\taverage_coverage\tlength\ttotal_coverage\n")
            for contig in contig_list: 
                coverage_summary.write("{identifier!s}\t{average_coverage:0.1f}\t{length}\t{total_coverage}\n".format(
                    identifier = contig, 
                    average_coverage = float(contigs_dict[contig][1]/contigs_dict[contig][0]), 
                    length = contigs_dict[contig][0], 
                    total_coverage = contigs_dict[contig][1]))
            coverage_summary.write("Total\t{average_coverage:0.1f}\t{length}\t{total_coverage}\n".format(
                average_coverage = float(total_coverage/total_length), 
                length = total_length, 
                total_coverage = total_coverage))

        with open(coverage_qa_file,"w") as coverage_qa:
            coverage_qa.write("coverage_bin\taverage_coverage\tlength\ttotal_coverage\n")
            for coverage in coverage_range:
                coverage_qa.write("{coverage_bin}\t{average_coverage}\t{length}\t{total_coverage}\n".format(
                    coverage_bin = coverage, 
                    average_coverage = float(binned_coverage[coverage][1]/binned_coverage[coverage][0]),
                    length = binned_coverage[coverage][0], 
                    total_coverage = binned_coverage[coverage][1]))

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def program__pilon_with_vcf(input_files, threads, memory, output_files = ["pilon.fasta","pilon.vcf"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        program_path = GLOBAL_pilon_jar_program_path
        contigs_file = input_files[0]
        bam_file = input_files[1]

        command_template = (
            "java -Xmx{memory}G -jar {program_path} --changes --genome {contigs_file} --frags {bam_file} --vcf --output pilon --threads {threads}"
            )
        context = {
            "memory":memory,
            "program_path":program_path,
            "contigs_file":contigs_file,
            "bam_file":bam_file,
            "threads":threads
            }

        fh.set_program_command(command_template.format(**context))
        fh.run_program()

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def script__summarize_vcf(input_files, output_files = ["pilon_vcf_summary.txt","pilon_summary.vcf"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        vcf_file = input_files[0]
        vcf_summary_file = output_files[0]
        summary_vcf = output_files[1]

        number_of_records = 0
        number_of_single_site_records = 0
        number_of_single_site_records_with_depth = 0
        number_of_single_site_records_with_depth_and_ambiguous = 0
        number_of_single_site_records_with_depth_and_snp = 0
        number_of_single_site_records_low_depth = 0
        number_of_indels_and_deletions = 0

        ambiguous_cutoff = float(config.get("values","summarize_vcf__ambiguous_cutoff"))
        min_depth = int(config.get("values","summarize_vcf__min_depth"))

        vcf_reader = vcf.Reader(open(vcf_file, 'r'))
        vcf_writer = vcf.Writer(open(summary_vcf, 'w'), vcf_reader)
        for record in vcf_reader:
            number_of_records += 1
            if(record.affected_end - record.affected_start == 1):
                number_of_single_site_records += 1
                if "DP" in record.INFO and "AF" in record.INFO:
                    if record.INFO["DP"] > min_depth:
                        number_of_single_site_records_with_depth += 1
                        if sum(record.INFO["AF"]) < 1-ambiguous_cutoff and sum(record.INFO["AF"]) > ambiguous_cutoff:
                            number_of_single_site_records_with_depth_and_ambiguous += 1
                            vcf_writer.write_record(record)
                        if sum(record.INFO["AF"]) >= 1-ambiguous_cutoff:
                            number_of_single_site_records_with_depth_and_snp += 1
                            vcf_writer.write_record(record)
                    else:
                        number_of_single_site_records_low_depth += 1
                else:
                    vcf_writer.write_record(record)
            else:
                if "DP" in record.INFO:
                    if record.INFO["DP"] > min_depth:
                        number_of_indels_and_deletions += 1
                        vcf_writer.write_record(record)
                else:
                    vcf_writer.write_record(record)

        with open(vcf_summary_file, "w") as output_file:
            template = (
                "number_of_records\t{number_of_records}\n"
                "number_of_single_site_records\t{number_of_single_site_records}\n"
                "number_of_single_site_records_with_depth\t{number_of_single_site_records_with_depth}\n"
                "number_of_single_site_records_with_depth_and_ambiguous\t{number_of_single_site_records_with_depth_and_ambiguous}\n"
                "number_of_single_site_records_with_depth_and_snp\t{number_of_single_site_records_with_depth_and_snp}\n"
                "number_of_single_site_records_low_depth\t{number_of_single_site_records_low_depth}\n"
                "number_of_indels_and_deletions\t{number_of_indels_and_deletions}\n"
                )
            context = {
                "number_of_records":number_of_records,
                "number_of_single_site_records":number_of_single_site_records,
                "number_of_single_site_records_with_depth":number_of_single_site_records_with_depth,
                "number_of_single_site_records_with_depth_and_ambiguous":number_of_single_site_records_with_depth_and_ambiguous,
                "number_of_single_site_records_with_depth_and_snp":number_of_single_site_records_with_depth_and_snp,
                "number_of_single_site_records_low_depth":number_of_single_site_records_low_depth,
                "number_of_indels_and_deletions":number_of_indels_and_deletions
                }
            output_file.write(template.format(**context))

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def program__ariba_on_database(input_files, threads, database, output_files = ["ariba_report.tsv"]):
    #removed threads for now due to issues with ariba and threads
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        R1_reads = input_files[0]
        R2_reads = input_files[1]
        ariba_output_directory = "ariba_"+database
        ariba_report = output_files[0]
        program_path = GLOBAL_ariba_program_path
        ariba_available_dbs = {}
        ariba_dbs = config.get("files","ariba_db").split(",")
        for item in ariba_dbs:
            if len(item.split(":")) == 2:
                ariba_available_dbs[item.split(":")[0]]=item.split(":")[1].replace("<serum>",os.path.dirname(os.path.realpath(__file__))+"/..")

        if database not in ariba_available_dbs:
            logger.error("{} not found in ariba_db for config")
            raise
            
        command_template = (
           "{program_path} run {db_path} {R1_reads} {R2_reads} {output_directory}"
            )
        context = {
            "program_path":program_path,
            "db_path":ariba_available_dbs[database],
            "R1_reads":R1_reads,
            "R2_reads":R2_reads,
            "output_directory":ariba_output_directory
            }
        fh.set_program_command(command_template.format(**context))
        fh.run_program()

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        shutil.copy(os.path.join(ariba_output_directory,"report.tsv"),ariba_report)
        fh.completed()
        return fh
def program__mlst(input_files, mlst_species = "", output_files = ["mlst.txt"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        program_path = GLOBAL_mlst_program_path
        contigs_file = input_files[0]
        mlst_file = output_files[0]

        mlst_arg = ""
        if mlst_species != "" and mlst_species != "NA":
            mlst_arg = "--scheme="+mlst_species.strip()

        command_template = (
            "{program_path} {contigs_file} {mlst_arg}"
            )
        context = {
            "program_path":program_path,
            "contigs_file":contigs_file,
            "mlst_arg":mlst_arg
            }

        fh.set_program_command(command_template.format(**context))
        fh.run_program()
        with open(mlst_file,"w") as output_file:
            output_file.write(fh.get_stdout())

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def script__contigs_check(input_files, output_files = ["contigs_qa.txt"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        coverage_summary_file = input_files[0]
        contigs_check_file = output_files[0]
        min_contig_length = int(config.get("values","contigs_check__min_contig_length"))

        min_contig_coverage = [0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100]
        total_length_over_min = [0]*len(min_contig_coverage)
        total_contigs_over_min = [0]*len(min_contig_coverage)

        contig_length_and_coverage_pattern = re.compile(config.get("values","contig_length_and_coverage_pattern"))

        with open(coverage_summary_file, "r") as coverage_summary_data:
            for line in coverage_summary_data:
                result = re.match(contig_length_and_coverage_pattern,line)
                if result:
                    length = int(result.group("length"))
                    coverage = float(result.group("bpcov"))

                    for i in range(len(min_contig_coverage)):
                        if length > min_contig_length and coverage > min_contig_coverage[i]:
                            total_length_over_min[i] += length
                            total_contigs_over_min[i] += 1

        with open(contigs_check_file, "w") as output_file:
            template = (
                "min coverage\ttotal length\ttotal contigs\n"
                "{min_contig_coverage_0}\t{total_length_over_min_0}\t{total_contigs_over_min_0}\n"
                "{min_contig_coverage_1}\t{total_length_over_min_1}\t{total_contigs_over_min_1}\n"
                "{min_contig_coverage_2}\t{total_length_over_min_2}\t{total_contigs_over_min_2}\n"
                "{min_contig_coverage_3}\t{total_length_over_min_3}\t{total_contigs_over_min_3}\n"
                "{min_contig_coverage_4}\t{total_length_over_min_4}\t{total_contigs_over_min_4}\n"
                "{min_contig_coverage_5}\t{total_length_over_min_5}\t{total_contigs_over_min_5}\n"
                "{min_contig_coverage_6}\t{total_length_over_min_6}\t{total_contigs_over_min_6}\n"
                "{min_contig_coverage_7}\t{total_length_over_min_7}\t{total_contigs_over_min_7}\n"
                "{min_contig_coverage_8}\t{total_length_over_min_8}\t{total_contigs_over_min_8}\n"
                "{min_contig_coverage_9}\t{total_length_over_min_9}\t{total_contigs_over_min_9}\n"
                "{min_contig_coverage_10}\t{total_length_over_min_10}\t{total_contigs_over_min_10}\n"
                "{min_contig_coverage_11}\t{total_length_over_min_11}\t{total_contigs_over_min_11}\n"
                "{min_contig_coverage_12}\t{total_length_over_min_12}\t{total_contigs_over_min_12}\n"
                "{min_contig_coverage_13}\t{total_length_over_min_13}\t{total_contigs_over_min_13}\n"
                "{min_contig_coverage_14}\t{total_length_over_min_14}\t{total_contigs_over_min_14}\n"
                "{min_contig_coverage_15}\t{total_length_over_min_15}\t{total_contigs_over_min_15}\n"
                "{min_contig_coverage_16}\t{total_length_over_min_16}\t{total_contigs_over_min_16}\n"
                "{min_contig_coverage_17}\t{total_length_over_min_17}\t{total_contigs_over_min_17}\n"
                "{min_contig_coverage_18}\t{total_length_over_min_18}\t{total_contigs_over_min_18}\n"
                "{min_contig_coverage_19}\t{total_length_over_min_19}\t{total_contigs_over_min_19}\n"
                "{min_contig_coverage_20}\t{total_length_over_min_20}\t{total_contigs_over_min_20}\n"
                "* using a minimum contig length filter of {min_contig_length}bp"
                )
            context = {
                "min_contig_coverage_0":min_contig_coverage[0],"total_length_over_min_0":total_length_over_min[0],"total_contigs_over_min_0":total_contigs_over_min[0],
                "min_contig_coverage_1":min_contig_coverage[1],"total_length_over_min_1":total_length_over_min[1],"total_contigs_over_min_1":total_contigs_over_min[1],
                "min_contig_coverage_2":min_contig_coverage[2],"total_length_over_min_2":total_length_over_min[2],"total_contigs_over_min_2":total_contigs_over_min[2],
                "min_contig_coverage_3":min_contig_coverage[3],"total_length_over_min_3":total_length_over_min[3],"total_contigs_over_min_3":total_contigs_over_min[3],
                "min_contig_coverage_4":min_contig_coverage[4],"total_length_over_min_4":total_length_over_min[4],"total_contigs_over_min_4":total_contigs_over_min[4],
                "min_contig_coverage_5":min_contig_coverage[5],"total_length_over_min_5":total_length_over_min[5],"total_contigs_over_min_5":total_contigs_over_min[5],
                "min_contig_coverage_6":min_contig_coverage[6],"total_length_over_min_6":total_length_over_min[6],"total_contigs_over_min_6":total_contigs_over_min[6],
                "min_contig_coverage_7":min_contig_coverage[7],"total_length_over_min_7":total_length_over_min[7],"total_contigs_over_min_7":total_contigs_over_min[7],
                "min_contig_coverage_8":min_contig_coverage[8],"total_length_over_min_8":total_length_over_min[8],"total_contigs_over_min_8":total_contigs_over_min[8],
                "min_contig_coverage_9":min_contig_coverage[9],"total_length_over_min_9":total_length_over_min[9],"total_contigs_over_min_9":total_contigs_over_min[9],
                "min_contig_coverage_10":min_contig_coverage[10],"total_length_over_min_10":total_length_over_min[10],"total_contigs_over_min_10":total_contigs_over_min[10],
                "min_contig_coverage_11":min_contig_coverage[11],"total_length_over_min_11":total_length_over_min[11],"total_contigs_over_min_11":total_contigs_over_min[11],
                "min_contig_coverage_12":min_contig_coverage[12],"total_length_over_min_12":total_length_over_min[12],"total_contigs_over_min_12":total_contigs_over_min[12],
                "min_contig_coverage_13":min_contig_coverage[13],"total_length_over_min_13":total_length_over_min[13],"total_contigs_over_min_13":total_contigs_over_min[13],
                "min_contig_coverage_14":min_contig_coverage[14],"total_length_over_min_14":total_length_over_min[14],"total_contigs_over_min_14":total_contigs_over_min[14],
                "min_contig_coverage_15":min_contig_coverage[15],"total_length_over_min_15":total_length_over_min[15],"total_contigs_over_min_15":total_contigs_over_min[15],
                "min_contig_coverage_16":min_contig_coverage[16],"total_length_over_min_16":total_length_over_min[16],"total_contigs_over_min_16":total_contigs_over_min[16],
                "min_contig_coverage_17":min_contig_coverage[17],"total_length_over_min_17":total_length_over_min[17],"total_contigs_over_min_17":total_contigs_over_min[17],
                "min_contig_coverage_18":min_contig_coverage[18],"total_length_over_min_18":total_length_over_min[18],"total_contigs_over_min_18":total_contigs_over_min[18],
                "min_contig_coverage_19":min_contig_coverage[19],"total_length_over_min_19":total_length_over_min[19],"total_contigs_over_min_19":total_contigs_over_min[19],
                "min_contig_coverage_20":min_contig_coverage[20],"total_length_over_min_20":total_length_over_min[20],"total_contigs_over_min_20":total_contigs_over_min[20],
                "min_contig_length":min_contig_length
                }
            output_file.write(template.format(**context))

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def script__filter_contigs(input_files, output_files = ["filtered_contigs.fasta"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        contigs_file = input_files[0]
        filtered_contigs_file = output_files[0]
        min_length = int(config.get("values","filter_contigs__min_length"))
        min_depth = float(config.get("values","filter_contigs__min_depth"))
        contig_length_and_coverage_pattern = re.compile(config.get("values","contig_length_and_coverage_pattern"))

        with open(contigs_file,"r")  as fasta_input:
            records = list(Bio.SeqIO.parse(fasta_input,"fasta"))
            filtered_sequences = []
            for record in records:
                result = re.search(contig_length_and_coverage_pattern,record.id)
                if result:
                    length = int(result.group("length"))
                    coverage = float(result.group("bpcov"))
                    if length >= min_length and coverage >= min_depth:
                        filtered_sequences.append(record)

        with open(filtered_contigs_file,"w") as filtered_contigs:
            Bio.SeqIO.write(filtered_sequences, filtered_contigs,"fasta")

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def script__append_run_summary_file(input_files, output_files = ["run_summary_file.txt"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.start()

    try:
        sample_summary_file = input_files[0]
        run_summary_file = output_files[0]
        sample_dict = {}
        qc_categories = config.get("categories","qc").split(",")

        for item in qc_categories:
            sample_dict[item] = "NA"

        with open(sample_summary_file,"r") as sample_summary:
            for line in sample_summary:
                key_value = line.split(":",1)
                if len(key_value) == 2:
                    key = key_value[0]
                    value = key_value[1].replace("\t"," ").strip()
                    sample_dict[key] = value
                else:
                    logger.warning("run_script__append_run_summary_file failed to handle {}".format(line))

        sample_line = ""
        with open(run_summary_file,"a") as output_file:
            for item in qc_categories:
                sample_line = sample_line+sample_dict[item]+"\t"
            sample_line = sample_line.strip()+"\n"
            output_file.write(sample_line)

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def script__create_run_summary_file(output_files = ["run_summary_file.txt"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        run_summary_file = output_files[0]
        header = ""
        qc_categories = config.get("categories","qc").split(",")
        for item in qc_categories:
            header = header+item+"\t"
        header = header.strip()+"\n"

        with open(run_summary_file,"w") as output_file:
            output_file.write(header)

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def script__qc_sample(input_files, output_files = ["qc.txt"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        #sample_name,supplying_lab,initials,run_name,qc_action,R1_location,R2_location,
        #output_directory,num_of_reads,trimmed_num_of_reads,normalized_num_of_reads,
        #read_length,mean_read_length,provided_species,detected_species,genuses_detected,
        #number_of_genuses_detected,mlst_species,mlst_type,mlst_alleles,read_insert_size,
        #read_deviation,total_length_at_0,total_length_at_25,total_length_difference,
        #total_contigs_at_0,total_contigs_at_25,total_contigs_difference,
        #theoretical_coverage,total_coverage,failed_proportion_filter_count,called_snps,
        #qc_organism,qc_genuses_count,qc_length_at_0,qc_length_at_25,qc_length_difference,
        #qc_contigs_at_0,qc_contigs_at_25,qc_coverage,comments
        qc_categories = config.get("categories","qc").split(",")
        identifier_file_categories = config.get("categories","identifier_file").split(",")
        qc_labels = config.get("categories","qc_labels").split(",")
        qc_labels_dict = {}

        for item in qc_labels:
            qc_labels_dict[item.split(":")[0]] = item.split(":")[1]
        #Output stored in dict and all values as strings
        qc_dict = {}
        for item in qc_categories:
            qc_dict[item] = "NA"

        identifier_file = input_files[0]
        qc_file = output_files[0]

        qc_dict["identifier_file"] = identifier_file
        qc_dict["qc_file"] = qc_file

        identifier_dataframe = get_from_file__pandas_dataframe([identifier_file])

        for category in identifier_file_categories:
            if category in identifier_dataframe:
                qc_dict[category] = identifier_dataframe.loc[0][category]
            else:
                qc_dict[category] = "NA"

        #step__read_management QC
        read_lengths_file = ["read_length.txt"]
        if os.path.isfile(read_lengths_file[0]):
            qc_dict["num_of_reads"] = get_from_file__first_occurance_of_pattern(read_lengths_file,"Both-\tReads: (?P<reads>[0-9]+)\t Mean: [0-9]+[\.]?[0-9]*\t Mode: [0-9]+\t Max: [0-9]+\t Min: [0-9]+")
            qc_dict["mean_read_length"] = get_from_file__first_occurance_of_pattern(read_lengths_file,".*trimmed.*\n.*\n.*\n.*Both-\tReads: [0-9]+\t Mean: (?P<mean_read_length>[0-9]+[\.]?[0-9]*)\t Mode: [0-9]+\t Max: [0-9]+\t Min: [0-9]+")
            qc_dict["read_length"] = get_from_file__first_occurance_of_pattern(read_lengths_file,"Both-\tReads: [0-9]+\t Mean: [0-9]+[\.]?[0-9]*\t Mode: (?P<read_length>[0-9]+)\t Max: [0-9]+\t Min: [0-9]+")
            qc_dict["trimmed_num_of_reads"] = get_from_file__first_occurance_of_pattern(read_lengths_file,".*trimmed.*\n.*\n.*\n.*Both-\tReads: (?P<reads>[0-9]+)\t Mean: [0-9]+[\.]?[0-9]*\t Mode: [0-9]+\t Max: [0-9]+\t Min: [0-9]+")
            qc_dict["normalized_num_of_reads"] = get_from_file__first_occurance_of_pattern(read_lengths_file,".*normalized.*\n.*\n.*\n.*Both-\tReads: (?P<reads>[0-9]+)\t Mean: [0-9]+[\.]?[0-9]*\t Mode: [0-9]+\t Max: [0-9]+\t Min: [0-9]+")

        #step__kraken_on_reads QC
        kraken_summary_report_file = ["kraken_report_summary.txt"]
        if os.path.isfile(kraken_summary_report_file[0]):
            qc_dict["detected_species_ncbi"] = get_from_file__first_occurance_of_pattern(kraken_summary_report_file,"Species detected: (?P<detected_species>.*?)\n")
            qc_dict["genuses_detected"] = get_from_file__first_occurance_of_pattern(kraken_summary_report_file,"Genuses found:(?P<genuses>.*?)\n")
            qc_dict["number_of_genuses_detected"] = get_from_file__first_occurance_of_pattern(kraken_summary_report_file,"Number of genuses in sample: (?P<number_of_genuses_detected>[0-9]+)\n")
            qc_dict["percent_unclassified"] = get_from_file__first_occurance_of_pattern(kraken_summary_report_file,"Percent unclassified: (?P<percent_unclassified>[0-9]+[\.]?[0-9]*)\n")

        #program__mlst
        mlst_file = ["mlst.txt"]
        if os.path.isfile(mlst_file[0]):
            qc_dict["mlst_species"] = get_from_file__first_occurance_of_pattern(mlst_file,".*_contigs\.fasta\s(?P<mlst_species>.+?)\s.+?\s.*")
            qc_dict["mlst_type"] = get_from_file__first_occurance_of_pattern(mlst_file,".*_contigs\.fasta\s.+?\s(?P<mlst_type>.+?)\s.*")
            qc_dict["mlst_alleles"] = get_from_file__first_occurance_of_pattern(mlst_file,".*_contigs\.fasta\s.+?\s.+?\s(?P<mlst_alleles>.*)")

        #step__assembler QC
        spades_log_file = ["spades.log"]
        if os.path.isfile(spades_log_file[0]):
            qc_dict["read_insert_size"] = get_from_file__first_occurance_of_pattern(spades_log_file,"Insert size = (?P<insert_size>.*?), deviation = .*?,")
            qc_dict["read_deviation"] = get_from_file__first_occurance_of_pattern(spades_log_file,"Insert size = .*?, deviation = (?P<deviation>.*?),")

        #step__contig_read_correction QC
        vcf_summary_file = ["pilon_vcf_summary.txt"]
        if os.path.isfile(vcf_summary_file[0]):
            qc_dict["failed_proportion_filter_count"] = get_from_file__first_occurance_of_pattern(vcf_summary_file,"number_of_single_site_records_with_depth_and_ambiguous\t(?P<failed_proportion_filter_count>[0-9]+)")
            qc_dict["called_snps"] = get_from_file__first_occurance_of_pattern(vcf_summary_file,"number_of_single_site_records_with_depth_and_snp\t(?P<number_of_single_site_records_with_depth_and_snp>[0-9]+)")

        #QC actions and warnings
        qc_action = ""
        qc_coverage = "NA"
        qc_organism = "OK"
        any_warnings = False

        qc_db = get_from_file__pandas_dataframe([GLOBAL_species_db_path])
        qc_db__organism = ""
        qc_db__mlst_species = ""
        qc_db__ncbi_species = ""
        qc_db__min_length = 0
        qc_db__max_length = 0
        qc_db__max_contigs = 0
        qc_db__coverage_base = 0
        qc_db__coverage_compare = 0
        qc_db__lack_of_reads_cov = 0.0
        qc_db__not_enough_reads_cov = 0.0
        qc_db__low_reads_cov = 0.0
        qc_db__max_genuses = 0
        qc_db__allowable_percent_unclassified = 0.0
        qc_db__length_difference = 0
        qc_db__contig_difference = 0

        qc_db_default_index = qc_db[qc_db["organism"] == "default"].index.tolist()
        if len(qc_db_default_index) == 1:
            qc_db_default_index = qc_db_default_index[0]
            qc_db__organism = str(setIfNotNA(qc_db__organism, qc_db.at[qc_db_default_index,"organism"]))
            qc_db__mlst_species = str(setIfNotNA(qc_db__mlst_species, qc_db.at[qc_db_default_index,"mlst_species"]))
            qc_db__ncbi_species = str(setIfNotNA(qc_db__ncbi_species, qc_db.at[qc_db_default_index,"ncbi_species"]))
            qc_db__min_length = float(setIfNotNA(qc_db__min_length, qc_db.at[qc_db_default_index,"min_length"]))
            qc_db__max_length = float(setIfNotNA(qc_db__max_length, qc_db.at[qc_db_default_index,"max_length"]))
            qc_db__max_contigs = float(setIfNotNA(qc_db__max_contigs, qc_db.at[qc_db_default_index,"max_contigs"]))
            qc_db__coverage_base = int(float(setIfNotNA(qc_db__coverage_base, qc_db.at[qc_db_default_index,"coverage_base"])))
            qc_db__coverage_compare  = int(float(setIfNotNA(qc_db__coverage_compare, qc_db.at[qc_db_default_index,"coverage_compare"])))
            qc_db__lack_of_reads_cov = float(setIfNotNA(qc_db__lack_of_reads_cov, qc_db.at[qc_db_default_index,"lack_of_reads_cov"]))
            qc_db__not_enough_reads_cov = float(setIfNotNA(qc_db__not_enough_reads_cov, qc_db.at[qc_db_default_index,"not_enough_reads_cov"]))
            qc_db__low_reads_cov = float(setIfNotNA(qc_db__low_reads_cov, qc_db.at[qc_db_default_index,"low_reads_cov"]))
            qc_db__max_genuses = float(setIfNotNA(qc_db__max_genuses, qc_db.at[qc_db_default_index,"max_genuses"]))
            qc_db__allowable_percent_unclassified = float(setIfNotNA(qc_db__allowable_percent_unclassified, qc_db.at[qc_db_default_index,"allowable_percent_unclassified"]))
            qc_db__length_difference = float(setIfNotNA(qc_db__length_difference, qc_db.at[qc_db_default_index,"length_difference"]))
            qc_db__contig_difference = float(setIfNotNA(qc_db__contig_difference, qc_db.at[qc_db_default_index,"contig_difference"]))
        else:
            logger.error("Default entry for species database not found or duplicated, this will likely cause issues with QC")
        
        provided_species = qc_dict["Organism"]
        detected_species_ncbi = qc_dict["detected_species_ncbi"]
        qc_db_species_index = qc_db[qc_db["organism"] == provided_species].index.tolist()
        if len(qc_db_species_index) == 0:
            if provided_species == "NA":
                qc_organism = "no species provided, using {} for QC values".format(detected_species_ncbi)
            else:
                qc_organism = "WARNING: {} provided and not found in DB, using {} for QC values".format(provided_species, detected_species_ncbi)
            qc_db_species_index = qc_db[qc_db["ncbi_species"] == detected_species_ncbi].index.tolist()
        if len(qc_db_species_index) == 0:
            qc_organism = "WARNING: Dominant species {} not found in database, using default values".format(detected_species_ncbi)
        elif len(qc_db_species_index) == 1:
            qc_db_species_index = qc_db_species_index[0]
            qc_db__organism = str(setIfNotNA(qc_db__organism, qc_db.at[qc_db_species_index,"organism"]))
            qc_db__mlst_species = str(setIfNotNA(qc_db__mlst_species, qc_db.at[qc_db_species_index,"mlst_species"]))
            qc_db__ncbi_species = str(setIfNotNA(qc_db__ncbi_species, qc_db.at[qc_db_species_index,"ncbi_species"]))
            qc_db__min_length = float(setIfNotNA(qc_db__min_length, qc_db.at[qc_db_species_index,"min_length"]))
            qc_db__max_length = float(setIfNotNA(qc_db__max_length, qc_db.at[qc_db_species_index,"max_length"]))
            qc_db__max_contigs = float(setIfNotNA(qc_db__max_contigs, qc_db.at[qc_db_species_index,"max_contigs"]))
            qc_db__coverage_base = int(float(setIfNotNA(qc_db__coverage_base, qc_db.at[qc_db_species_index,"coverage_base"])))
            qc_db__coverage_compare  = int(float(setIfNotNA(qc_db__coverage_compare, qc_db.at[qc_db_species_index,"coverage_compare"])))
            qc_db__lack_of_reads_cov = float(setIfNotNA(qc_db__lack_of_reads_cov, qc_db.at[qc_db_species_index,"lack_of_reads_cov"]))
            qc_db__not_enough_reads_cov = float(setIfNotNA(qc_db__not_enough_reads_cov, qc_db.at[qc_db_species_index,"not_enough_reads_cov"]))
            qc_db__low_reads_cov = float(setIfNotNA(qc_db__low_reads_cov, qc_db.at[qc_db_species_index,"low_reads_cov"]))
            qc_db__max_genuses = float(setIfNotNA(qc_db__max_genuses, qc_db.at[qc_db_species_index,"max_genuses"]))
            qc_db__allowable_percent_unclassified = float(setIfNotNA(qc_db__allowable_percent_unclassified, qc_db.at[qc_db_default_index,"allowable_percent_unclassified"]))
            qc_db__length_difference = float(setIfNotNA(qc_db__length_difference, qc_db.at[qc_db_species_index,"length_difference"]))
            qc_db__contig_difference = float(setIfNotNA(qc_db__contig_difference, qc_db.at[qc_db_species_index,"contig_difference"]))

        if(qc_db__ncbi_species != "" and qc_db__ncbi_species != detected_species_ncbi):
            qc_organism = "WARNING: {} provided but dominant species detected was {}".format(qc_db__ncbi_species, detected_species_ncbi)

        R1_location = qc_dict["R1_location"]
        R2_location = qc_dict["R2_location"]
        if R1_location == "NA" or R2_location == "NA":
            qc_action = "Action: Core facility"
            qc_coverage = "REDO: not found"

        number_of_genuses = qc_dict["number_of_genuses_detected"]
        if(number_of_genuses.isdigit()):
            qc_dict["qc_genuses_count"] = "OK"
            if(int(number_of_genuses) > qc_db__max_genuses):
                qc_dict["qc_genuses_count"] = "WARNING: {} genuses detected".format(number_of_genuses)

        percent_unclassified = qc_dict["percent_unclassified"]
        if(isfloat(percent_unclassified)):
            qc_dict["qc_allowable_percent_unclassified"] = "OK"
            if(float(percent_unclassified) > qc_db__allowable_percent_unclassified):
                qc_dict["qc_allowable_percent_unclassified"] = "WARNING: {} unclassified reads".format(percent_unclassified)

        trimmed_num_of_reads = qc_dict["trimmed_num_of_reads"]
        mean_read_length = qc_dict["mean_read_length"]
        if trimmed_num_of_reads.isdigit() and isfloat(mean_read_length):
            number_of_bases = float(trimmed_num_of_reads) * float(mean_read_length) * 2 #for paired reads derp
            qc_dict["theoretical_coverage"] = "{coverage_min:0.1f} - {coverage_max:0.1f}".format(coverage_min = number_of_bases/qc_db__max_length, coverage_max = number_of_bases/qc_db__min_length)

        #step__denovo_mapping QC
        coverage_base = qc_db__coverage_base
        coverage_compare = qc_db__coverage_compare
        coverage_qa_file = ["coverage_qa.txt"]
        contigs_qa_file = ["contigs_qa.txt"]
        coverage_summary_file = ["coverage_summary.txt"]
        if(os.path.isfile(coverage_qa_file[0]) and os.path.isfile(contigs_qa_file[0]) and os.path.isfile(coverage_summary_file[0])):
            qc_dict["coverage_base"] = coverage_base
            qc_dict["coverage_compare"] = coverage_compare
            qc_dict["bp_length_at_coverage_base"] = get_from_file__first_occurance_of_pattern(coverage_qa_file,"{}\t[0-9]+[\.]?[0-9]*\t(?P<length>[0-9]+)\t[0-9]+\n".format(coverage_base))
            qc_dict["bp_length_at_coverage_compare"] = get_from_file__first_occurance_of_pattern(coverage_qa_file,"{}\t[0-9]+[\.]?[0-9]*\t(?P<length>[0-9]+)\t[0-9]+\n".format(coverage_compare))
            qc_dict["length_at_coverage_base"] = get_from_file__first_occurance_of_pattern(contigs_qa_file,"{}\t(?P<length>[0-9]+)\t[0-9]+\n".format(coverage_base))
            qc_dict["length_at_coverage_compare"] = get_from_file__first_occurance_of_pattern(contigs_qa_file,"{}\t(?P<length>[0-9]+)\t[0-9]+\n".format(coverage_compare))
            qc_dict["contigs_at_coverage_base"] = get_from_file__first_occurance_of_pattern(contigs_qa_file,"{}\t[0-9]+\t(?P<contigs>[0-9]+)\n".format(coverage_base))
            qc_dict["contigs_at_coverage_compare"] = get_from_file__first_occurance_of_pattern(contigs_qa_file,"{}\t[0-9]+\t(?P<contigs>[0-9]+)\n".format(coverage_compare))
            qc_dict["total_coverage"] = get_from_file__first_occurance_of_pattern(coverage_summary_file,"Total\t(?P<total_coverage>[0-9]+[\.]?[0-9]*)\t[0-9]*\t[0-9]*\n")
            qc_dict["bp_length_difference"] = "NA"
            qc_dict["length_difference"] = "NA"
            qc_dict["contig_difference"] = "NA"
        coverage_comparisons = [["bp_length_at_coverage_base", "bp_length_at_coverage_compare", "bp_length_difference", qc_db__min_length, qc_db__max_length, qc_db__length_difference, "testqc_bp_length_at_coverage_base", "testqc_bp_length_at_coverage_compare", "testqc_bp_length_difference"],
                                ["length_at_coverage_base", "length_at_coverage_compare", "length_difference", qc_db__min_length, qc_db__max_length, qc_db__length_difference, "qc_length_at_coverage_base", "qc_length_at_coverage_compare", "qc_length_difference"],
                                ["contigs_at_coverage_base", "contigs_at_coverage_compare", "contig_difference", 1, qc_db__max_contigs, qc_db__contig_difference, "qc_contigs_at_coverage_base", "qc_contigs_at_coverage_compare", "qc_contigs_difference"],]
        for item in coverage_comparisons:
            if qc_dict[item[0]].isdigit() and qc_dict[item[1]].isdigit():
                qc_dict[item[2]] = str(int(qc_dict[item[0]]) - int(qc_dict[item[1]]))
            if qc_dict[item[0]].isdigit() and (int(qc_dict[item[0]]) >= item[3] and int(qc_dict[item[0]]) <= item[4]):
                qc_dict[item[6]] = "OK"
            else:
                qc_dict[item[6]] = "WARNING: size mismatch between {}-{} expected, {} found".format(item[3], item[4], qc_dict[item[0]])
            if qc_dict[item[1]].isdigit() and (int(qc_dict[item[1]]) >= item[3] and int(qc_dict[item[1]]) <= item[4]):
                qc_dict[item[7]] = "OK"
            else:
                qc_dict[item[7]] = "WARNING: size mismatch between {}-{} expected, {} found".format(item[3], item[4], qc_dict[item[1]])
            if qc_dict[item[2]].isdigit() and int(qc_dict[item[2]]) <= item[5]:
                qc_dict[item[8]] = "OK"
            else:
                qc_dict[item[8]] = "WARNING: {} size difference detected, {} is maximum allowed".format(qc_dict[item[2]], item[5])

        total_coverage = qc_dict["total_coverage"]
        num_of_reads = qc_dict["num_of_reads"]
        if(not num_of_reads.isdigit()):
            num_of_reads = 0
        if(int(num_of_reads) < int(config.get("values","minimum_reads_for_assembly"))):
                qc_action = "Action: Core facility"
                qc_coverage = "REDO: not enough reads to preform assembly"
        elif(isfloat(total_coverage)):
            if(float(total_coverage) < qc_db__lack_of_reads_cov):
                qc_action = "Action: Core facility"
                qc_coverage = "REDO: lack of reads"
            elif(float(total_coverage) < qc_db__not_enough_reads_cov):
                qc_action = "Action: Core facility"
                qc_coverage = "REDO: not enough reads"
            elif(float(total_coverage) < qc_db__low_reads_cov):
                qc_coverage = "WARNING: low reads"
            else:
                qc_coverage = "OK"
        else:
            qc_action = "Action: Core facility"
            qc_coverage = "REDO: not found"

        for item in qc_dict:
            if item.startswith("qc_"):
                if "WARNING" in qc_dict[item]:
                    any_warnings = True

        #Special case
        if(qc_action == "" and any_warnings == True):
            qc_action = "Action: Supplying lab"
            if( "WARNING" in qc_coverage and 
                "WARNING" in qc_dict["qc_length_difference"] and
                not "WARNING" in qc_dict["qc_length_at_coverage_base"] and 
                not "WARNING" in qc_dict["qc_genuses_count"] and 
                not "WARNING" in qc_organism
            ):
                qc_action = "Action: Core facility"
            if( "WARNING" in qc_organism or 
                "WARNING" in qc_dict["qc_genuses_count"]
            ):
                qc_action = "Action: Supplying lab"
        elif(qc_action == "" and any_warnings == False):
            qc_action = "Sample passed all QC"

        qc_dict["qc_action"] = qc_action
        qc_dict["qc_coverage"] = qc_coverage
        qc_dict["qc_organism"] = qc_organism

        for item in qc_labels_dict:
            if item in qc_dict:
                qc_dict[qc_labels_dict[item]] = qc_dict[item]

        with open(qc_file,"w") as output_file:
            for item in qc_categories:
                output_file.write("{}:{}\n".format(item,qc_dict[item]))

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def script__generate_nasp_config(input_files, identifier, partition, output_files = ["nasp.xml"]):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.set_output_files(output_files)
    if(fh.start()==SUCCESS):
        fh.completed()
        return fh

    try:
        reference_name = os.path.splitext(os.path.basename(input_files[0]))[0]
        reference_location = os.path.realpath(input_files[0])
        R1_read_file = os.path.basename(input_files[1])
        R2_read_file = os.path.basename(input_files[2])
        read_location =  os.path.dirname(os.path.realpath(input_files[1]))
        output_folder = os.path.join(os.getcwd(),"nasp")
        nasp_config = output_files[0]

        command_template = (
            "<?xml version = \"1\" ?>\n"
            "<NaspInputData>\n"
            "   <Options>\n"
            "       <RunName>nasp</RunName>\n"
            "       <OutputFolder>{output_folder}</OutputFolder>\n"
            "       <Reference name = \"{reference_name}\" path = \"{reference_location}\">\n"
            "           <FindDups>True</FindDups>\n"
            "       </Reference>\n"
            "       <Filters>\n"
            "           <ProportionFilter>0.9</ProportionFilter>\n"
            "           <CoverageFilter>10</CoverageFilter>\n"
            "       </Filters>\n"
            "       <JobSubmitter>SLURM</JobSubmitter>\n"
            "   </Options>\n"
            "   <Files>\n"
            "       <ReadFolder path = \"{read_location}\">\n"
            "           <ReadPair sample = \"{identifier}\">\n"
            "               <Read1Filename>{R1_read_file}</Read1Filename>\n"
            "               <Read2Filename>{R2_read_file}</Read2Filename>\n"
            "           </ReadPair>\n"
            "       </ReadFolder>\n"
            "   </Files>\n"
            "   <ExternalApplications>\n"
            "       <Index name = \"Index\" path = \"/tools/linuxbrew/bin\">\n"
            "           <AdditionalArguments/>\n"
            "           <JobParameters name = \"nasp_index\">\n"
            "               <MemRequested>2</MemRequested>\n"
            "               <NumCPUs>1</NumCPUs>\n"
            "               <Walltime>4</Walltime>\n"
            "               <Queue>{partition}</Queue>\n"
            "               <JobSubmitterArgs/>\n"
            "           </JobParameters>\n"
            "       </Index>\n"
            "       <MatrixGenerator name = \"MatrixGenerator\" path = \"/tools/linuxbrew/opt/python3/lib/python3.6/site-packages/nasp/nasptool_linux_64\">\n"
            "           <AdditionalArguments/>\n"
            "           <JobParameters name = \"nasp_matrix\">\n"
            "               <MemRequested>8</MemRequested>\n"
            "               <NumCPUs>4</NumCPUs>\n"
            "               <Walltime>90</Walltime>\n"
            "               <Queue>{partition}</Queue>\n"
            "               <JobSubmitterArgs/>\n"
            "           </JobParameters>\n"
            "       </MatrixGenerator>\n"
            "       <Picard name = \"Picard\" path = \"/tools/bin/picard.jar\">\n"
            "           <AdditionalArguments/>\n"
            "       </Picard>\n"
            "       <Samtools name = \"Samtools\" path = \"/tools/linuxbrew/Cellar/samtools/0.1.19/bin/samtools\">\n"
            "           <AdditionalArguments/>\n"
            "       </Samtools>\n"
            "       <DupFinder name = \"DupFinder\" path = \"/tools/linuxbrew/bin/nucmer\">\n"
            "           <AdditionalArguments/>\n"
            "           <JobParameters>\n"
            "               <MemRequested>4</MemRequested>\n"
            "               <NumCPUs>1</NumCPUs>\n"
            "               <Walltime>4</Walltime>\n"
            "               <Queue>{partition}</Queue>\n"
            "               <JobSubmitterArgs/>\n"
            "           </JobParameters>\n"
            "       </DupFinder>\n"
            "       <Aligner name = \"BWA-mem\" path = \"/tools/linuxbrew/bin/bwa\">\n"
            "           <AdditionalArguments/>\n"
            "           <JobParameters>\n"
            "               <MemRequested>10</MemRequested>\n"
            "               <NumCPUs>4</NumCPUs>\n"
            "               <Walltime>36</Walltime>\n"
            "               <Queue>{partition}</Queue>\n"
            "               <JobSubmitterArgs/>\n"
            "           </JobParameters>\n"
            "       </Aligner>\n"
            "       <SNPCaller name = \"GATK\" path = \"/tools/bin/GenomeAnalysisTK.jar\">\n"
            "           <AdditionalArguments>-stand_call_conf 100 -stand_emit_conf 100 -ploidy 1</AdditionalArguments>\n"
            "           <JobParameters>\n"
            "               <MemRequested>10</MemRequested>\n"
            "               <NumCPUs>4</NumCPUs>\n"
            "               <Walltime>36</Walltime>\n"
            "               <Queue>{partition}</Queue>\n"
            "               <JobSubmitterArgs/>\n"
            "           </JobParameters>\n"
            "       </SNPCaller>\n"
            "   </ExternalApplications>\n"
            "</NaspInputData>\n"
            )
        context = {
            "output_folder":output_folder,
            "reference_name":reference_name,
            "reference_location":reference_location,
            "read_location":read_location,
            "identifier":identifier,
            "R1_read_file":R1_read_file,
            "R2_read_file":R2_read_file,
            "partition":partition
        }

        with open(output_files[0],"w") as output_file:
            output_file.write(command_template.format(**context))

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def program__nasp_on_config(input_files):
    nasp_config = input_files[0]
    logger.info("Running nasp on {}".format(nasp_config))
    print("Running nasp:{}".format(" ".join(["nasp","--config",nasp_config])))
    process = subprocess.Popen(["nasp","--config",nasp_config],stdout = subprocess.PIPE,stderr = subprocess.STDOUT)
    process_out, process_err = process.communicate()

    last_job_id = "-1"
    with open(os.path.join("nasp","runlog.txt")) as runlog:
        for line in runlog:
            jobid = re.compile("^jobid = (?P<jobid>.+$)")
            result = re.search(jobid, line)
            if result:
                if result.group("jobid") > last_job_id:
                    last_job_id = result.group("jobid")

    logger.info("NASP last job id is {}".format(last_job_id))

    return last_job_id
def script__remove_word_from_contigs(word, contigs_file):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.start()

    try:
        with open(contigs_file,"r")  as fasta_input:
            records = list(Bio.SeqIO.parse(fasta_input,"fasta"))

        for record in records:
            record.id = record.id.replace(word,"")
            record.name = record.id
            record.description = record.id
        with open(contigs_file,"w") as output_file:
            Bio.SeqIO.write(records, output_file,"fasta")

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def script__remove_value_from_contigs(value, contigs_file):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.start()

    try:
        with open(contigs_file,"r")  as fasta_input:
            records = list(Bio.SeqIO.parse(fasta_input,"fasta"))

        pattern = re.compile(".*?(^|_)(?P<value>{}_.*?(_|$))".format(value))
        for record in records:
            result = re.search(pattern, record.id)
            if result:
                record.id = record.id.replace(result.group("value"),"")
                if record.id[-1] == "_":
                    record.id = record.id[:-1]
                record.name = record.id
                record.description = record.id
        with open(contigs_file,"w") as output_file:
            Bio.SeqIO.write(records, output_file,"fasta")

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def script__prepend_value_to_contigs(value, contigs_file):
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.start()

    try:
        with open(contigs_file,"r")  as fasta_input:
            records = list(Bio.SeqIO.parse(fasta_input,"fasta"))
        
        for record in records:
            record.id = value+record.id
            record.name = record.id
            record.description = record.id

        with open(contigs_file,"w") as output_file:
            Bio.SeqIO.write(records, output_file,"fasta")

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def script__append_contigs_with_bpcov(input_files, output_files = ["final_contigs.fasta"]):
    """
    Starts with assumption of spades contigs and ends with one formatted for QC needs
    """
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.start()

    try:
        contigs_file = input_files[0]
        coverage_file = input_files[1]
        final_contigs = output_files[0]

        coverage_dict = {}
        with open(coverage_file,"r") as coverage_input:
            for line in coverage_input:
                identifier = line.split("\t")[0]
                average_coverage = line.split("\t")[1]
                coverage_dict[identifier] = average_coverage

        with open(contigs_file,"r")  as fasta_input:
            records = list(Bio.SeqIO.parse(fasta_input,"fasta"))

        for record in records:
            if record.id in coverage_dict:
                record.id = "{identifier}_bpcov_{coverage}".format(identifier = record.id,coverage = coverage_dict[record.id])
                record.name = record.id
                record.description = record.id
            else:
                record.id = record.id+"_bpcov_ERROR"
                record.name = record.id
                record.description = record.id

        with open(final_contigs,"w") as output_file:
            Bio.SeqIO.write(records, output_file,"fasta")

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh
def script__update_contigs_length(input_files, output_files = ["final_contigs.fasta"]):
    """
    Small potential of length changes due to pilon, which can also affect contig numbers
    """
    fh = _function_helper(inspect.currentframe().f_code.co_name)
    fh.start()

    try:
        contigs_file = input_files[0]
        final_contigs = output_files[0]

        with open(contigs_file,"r")  as fasta_input:
            records = list(Bio.SeqIO.parse(fasta_input,"fasta"))

        pattern = re.compile(".*?(^|_)(?P<value>length_.*?(_|$))")

        for record in records:
            result = re.search(pattern, record.id)
            if result:
                record.id = record.id.replace(result.group("value"),"length_{}_".format(len(record.seq)))
                record.name = record.id
                record.description = record.id

        with open(final_contigs,"w") as output_file:
            Bio.SeqIO.write(records, output_file,"fasta")

    except Exception as e:
        fh.exception_handler(e)
        raise
    else:
        fh.completed()
        return fh

def run_job(job_name = "", job_command = "", task = "", working_dir = "", log_file = "", dependent_job_arg = "", dependent_job_ids = [], dependency = "afterok", partition_arg = "", partition = "", cpu_requirement = 1, memory_requirement_in_Gb = 1, time_requirement = "02:00:00", group = ""):
    #Wrapper to run jobs and be able to change queue system
    job = _Job_handler(
        job_name = job_name, 
        job_command = job_command, 
        task = task, 
        working_dir = working_dir, 
        log_file = log_file, 
        dependent_job_arg = dependent_job_arg, 
        dependent_job_ids = dependent_job_ids, 
        dependency = dependency, 
        partition_arg = partition_arg, 
        partition = partition, 
        cpu_requirement = cpu_requirement, 
        memory_requirement_in_Gb = memory_requirement_in_Gb, 
        time_requirement = time_requirement,
        group = group)
    if job.partition != "" and job.partition_arg == "":
        job.partition_arg = "-p {}".format(job.partition)
    grid_system = config.get("values","grid_system")

    if len(job.dependent_job_ids) > 0 and job.dependent_job_arg == "":
        job.generate_dependent_job_arg()

    if grid_system == "slurm":
        job_id = run_slurm_job(job)
    elif grid_system == "torque":
        job_id = run_torque_job(job)
    else:
        logger.error("Improper grid_system provided {}".format(grid_system))
        sys.exit(1)
    return job_id
def run_slurm_job(Job_handler):
    try:
        job_submit_command = "sbatch {dependant_job} -D {working_dir} -c {cpus} --mem={mem_in_Gb}G --time={time} -J \'{job_name}\' {partition}"
        job_submit_variables = {
            'dependant_job':Job_handler.dependent_job_arg,
            'working_dir':Job_handler.working_dir,
            'cpus':Job_handler.cpu_requirement,
            'mem_in_Gb':Job_handler.memory_requirement_in_Gb,
            'time':Job_handler.time_requirement, #D-HH:MM
            'job_name':Job_handler.job_name,
            'partition':Job_handler.partition_arg
            }
        run_command = Job_handler.job_command.replace('"','\\"')

        logger.info("Job submit command:{} --wrap=\"{}\"".format(job_submit_command.format(**job_submit_variables), run_command))
        process = subprocess.Popen("{} --wrap=\"{}\"".format(job_submit_command.format(**job_submit_variables), run_command), stdout = subprocess.PIPE, stderr = subprocess.STDOUT, shell = True)
        process_out, process_err = process.communicate()
        job_id_pattern = re.compile("^Submitted\sbatch\sjob\s(?P<job_id>.+)$")
        result = re.match(job_id_pattern,process_out.decode())
        if result:
            logger.info("{} last job id is {}".format(Job_handler.job_name,result.group("job_id")))
            job_id = result.group("job_id")
        else:
            job_id = "-1"
    except Exception as e:
        logger.error("Task: \"{}\" had an error. Error message {}".format(Job_handler.job_name,str(e)))
        raise
    else:
        pass
    finally:
        pass
    return job_id
def run_torque_job(Job_handler):
    try:
        job_submit_command = "qsub -V -d {working_dir} -w {working_dir} -l ncpus={cpus},mem={mem_in_Gb}gb,walltime={time} -N \'{job_name}\' -W group_list={group} -A {group} {dependant_job} "
        job_submit_variables = {
            'dependant_job':Job_handler.dependent_job_arg,
            'working_dir':Job_handler.working_dir,
            'cpus':Job_handler.cpu_requirement,
            'mem_in_Gb':Job_handler.memory_requirement_in_Gb,
            'time':Job_handler.time_requirement, #D-HH:MM
            'job_name':Job_handler.job_name,
            'group':Job_handler.group,
            'partition':Job_handler.partition_arg
            }
        run_command = Job_handler.job_command.replace('"','\\"')

        logger.info("Job submit command:echo \"{}\" | {}".format(run_command, job_submit_command.format(**job_submit_variables)))
        process = subprocess.Popen("echo \"{}\" | {}".format(run_command, job_submit_command.format(**job_submit_variables)), stdout = subprocess.PIPE, stderr = subprocess.STDOUT, shell = True)
        process_out, process_err = process.communicate()
        job_id_pattern = re.compile("^(?P<job_id>.+)$")
        result = re.match(job_id_pattern,process_out.decode())
        if result:
            logger.info("{} last job id is {}".format(Job_handler.job_name,result.group("job_id")))
            job_id = result.group("job_id")
        else:
            job_id = "-1"
    except Exception as e:
        logger.error("Task: \"{}\" had an error. Error message {}".format(Job_handler.job_name,str(e)))
        raise
    else:
        pass
    finally:
        pass
    return job_id
def submit_slurm_job(Job_handler):
    try:
        job_submit_command = "sbatch --hold {dependant_job} -D {working_dir} -c {cpus} --mem={mem_in_Gb}G --time={time} -J \'{job_name}\' {partition}"
        job_submit_variables = {
            'dependant_job':Job_handler.dependent_job_arg,
            'working_dir':Job_handler.working_dir,
            'cpus':Job_handler.cpu_requirement,
            'mem_in_Gb':Job_handler.memory_requirement_in_Gb,
            'time':Job_handler.time_requirement, #D-HH:MM
            'job_name':Job_handler.job_name,
            'partition':Job_handler.partition_arg
            }
        run_command = Job_handler.job_command.replace('"','\\"')

        logger.info("Job submit command:{} --wrap=\"{}\"".format(job_submit_command.format(**job_submit_variables), run_command))
        process = subprocess.Popen("{} --wrap=\"{}\"".format(job_submit_command.format(**job_submit_variables), run_command), stdout = subprocess.PIPE, stderr = subprocess.STDOUT, shell = True)
        process_out, process_err = process.communicate()
        job_id_pattern = re.compile("^Submitted\sbatch\sjob\s(?P<job_id>.+)$")
        result = re.match(job_id_pattern,process_out.decode())
        if result:
            logger.info("{} last job id is {}".format(Job_handler.job_name,result.group("job_id")))
            job_id = result.group("job_id")
        else:
            job_id = "-1"
    except Exception as e:
        logger.error("Task: \"{}\" had an error. Error message {}".format(Job_handler.job_name,str(e)))
        raise
    else:
        pass
    finally:
        pass
    return job_id
def release_slurm_job(job_ids):
    try:
        job_submit_command = "scontrol release JobId={job_ids}"
        job_submit_variables = {
            'job_ids':",".join([str(x) for x in job_ids])
            }
        run_command = job_submit_command.format(**job_submit_variables)

        process = subprocess.Popen(run_command, stdout = subprocess.PIPE, stderr = subprocess.STDOUT, shell = True)
        process_out, process_err = process.communicate()

    except Exception as e:
        logger.error("Task: release_slurm_job had an error. Error message {}".format(str(e)))
        raise

    return 0

def main():
    print("Class definition to be included, what does the candy shop provide?")

if __name__ == "__main__": main()
