#!/usr/bin/python
# coding: utf-8

 
# File: apytram.py
# Created by: Carine Rey
# Created on: Nov 2016
# 
# 
# Copyright or © or Copr. Carine Rey
# This software is a computer program whose purpose is to assembly 
# sequences from RNA-Seq data (paired-end or single-end) using one or
# more reference homologous sequences. 
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
# 


import os
import re
import sys
import time
import tempfile
import shutil
import logging
import argparse
import subprocess

from lib import ApytramNeeds
from lib import BlastPlus
from lib import Trinity
from lib import Aligner

start_time = time.time()

### Option defining
parser = argparse.ArgumentParser(prog = "apytram.py",
                                 description='''
    Run apytram.py on a fastq file to retrieve
    homologous sequences of bait sequences.''')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

##############
requiredOptions = parser.add_argument_group('Required arguments')
requiredOptions.add_argument('-d', '--database', type=str,
                             help='Database prefix name. If a database with the same name already exists, the existing database will be kept and the database will NOT be rebuilt.', required=True)
requiredOptions.add_argument('-dt', '--database_type', type=str, choices=["single","paired","FR","RF","F","R"],
                             help="""
                                  single: single unstranded data ______________________
                                  paired: paired unstranded data ______________________
                                  RF: paired stranded data (/1 = reverse ; /2 = forward)
                                  FR: paired stranded data (/2 = reverse ; /1 = forward)
                                  F: single stranded data (reads = forward) ____________
                                  R: single stranded data (reads = reverse) ____________
                                  WARNING: Paired read names must finished by 1 or 2""",
                            required=True)       
##############


##############
InUSOptions = parser.add_argument_group('Input Files')
InUSOptions.add_argument('-fa', '--fasta',  type=str, nargs='*',
                   help = "Fasta formated RNA-seq data to build the database of reads (only one file).")
InUSOptions.add_argument('-fq', '--fastq',  type=str, nargs='*',
                   help = "Fastq formated RNA-seq data to build the database of reads (several space delimited fastq file names are allowed). WARNING: Paired read names must finished by 1 or 2. (fastq files will be first converted to a fasta file. This process can require some time.)")
##############

##############
QueryOptions = parser.add_argument_group('Query File')
QueryOptions.add_argument('-q', '--query',  type=str,
                    help = "Fasta file (nucl) with homologous bait sequences which will be treated together for the apytram run. If no query is submitted, the program will just build the database. WARNING: Sequences must not contain other characters that a t g c n (eg. - * . )." )
QueryOptions.add_argument('-pep', '--query_pep',  type=str,
                   default = "",       
                   help = "Fasta file containing the query in the peptide format. It will be used at the first iteration as bait sequences to fish reads. It is compulsory to include also the query in nucleotide format (-q option)")
##############


##############
IterationOptions = parser.add_argument_group('Number of iterations')
IterationOptions.add_argument('-i', '--iteration_max',  type=int,
                    help = "Maximum number of iterations. (Default 5)",
                    default = 5 )
IterationOptions.add_argument('-i_start','--iteration_start',  type=int,
                    help = "Number of the first iteration. If different of 1, the tmp option must be used. (Default: 1)",
                    default = 1 )
##############


##############
OutOptions = parser.add_argument_group('Output Files')
OutOptions.add_argument('-out', '--output_prefix',  type=str, default = "./apytram",
                   help = "Output prefix (Default ./apytram)")
OutOptions.add_argument('-log', type=str, default="apytram.log",
                   help = "a log file to report avancement (default: apytram.log)")
OutOptions.add_argument('-tmp',  type=str,
                    help = "Directory to stock all intermediary files for the apytram run. (default: a directory in /tmp which will be removed at the end)",
                    default = "" )
OutOptions.add_argument('--keep_iterations',  action='store_true',
                    help = "A fasta file containing reconstructed sequences will be created at each iteration. (default: False)")
OutOptions.add_argument('--no_best_file',  action='store_true',
                        default = False,
                        help = "By default, a fasta file (Outprefix.best.fasta) containing only the best sequence is created. If this option is used, it will NOT be created.")

OutOptions.add_argument('--only_best_file',  action='store_true',
                        default = False,
                        help = "By default, a fasta file (Outprefix.fasta) containing all sequences from the last iteration is created. If this option is used, it will NOT be created.")

OutOptions.add_argument('--stats', action='store_true',
                             help='Create files with statistics on each iteration. (default: False)')
# comm Marie : tu expliques les détails des stats quelque part?                             
OutOptions.add_argument('--plot', action='store_true',
                             help='Create plots to represent the statistics on each iteration. (default: False)')
OutOptions.add_argument('--plot_ali', action='store_true',
                             help='Create file with a plot representing the alignement of all sequences from the last iteration on the query sequence. Take some seconds. (default: False)')
##############


##############
SearchOptions = parser.add_argument_group('Thresholds for EACH ITERATION')
SearchOptions.add_argument('-e', '--evalue',  type=float,
                    help = "Evalue threshold of the blastn of the bait queries on the database of reads. (Default 1e-3)",
                    default = 1e-3 )

SearchOptions.add_argument('-id', '--min_id',  type=int,
                    help = "Minimum identity percentage of a sequence with a query on the length of their alignment so that the sequence is kept at the end of a iteration (Default 60)",
                    default = 60 )
SearchOptions.add_argument('-mal', '--min_ali_len',  type=int,
                    help = "Minimum alignment length of a sequence on a query to be kept at the end of a iteration (Default 180)",
                    default = 180 )
SearchOptions.add_argument('-len', '--min_len',  type=int,
                    help = "Minimum length to keep a sequence at the end of a iteration. (Default 200)",
                    default = 200 )
##############


##############
StopOptions = parser.add_argument_group('Criteria to stop iteration')
StopOptions.add_argument('-required_coverage',  type=float,
                    help = "Required coverage of a bait sequence to stop iteration (Default: No threshold)",
                    default = 200 )
StopOptions.add_argument('--finish_all_iter', action='store_true',
                    help = "By default, iterations are stop if there is no improvment, if this option is used apytram will finish all iteration (-i).",
                    default = False)              
                    
##############


##############
FinalFilterOptions = parser.add_argument_group('Thresholds for Final output files')
FinalFilterOptions.add_argument('-flen', '--final_min_len',  type=int,
                    help = "Minimum PERCENTAGE of the query length to keep a sequence at the end of the run. (Default: 0)",
                    default = 0 )
FinalFilterOptions.add_argument('-fid', '--final_min_id',  type=int,
                    help = "Minimum identity PERCENTAGE of a sequence with a query on the length of their alignment so that the sequence is kept at the end of the run (Default 0)",
                    default = 0 )
FinalFilterOptions.add_argument('-fmal', '--final_min_ali_len',  type=int,
                     help = "Alignment length between a sequence and a query must be at least this PERCENTAGE of the query length to keep this sequence at the end of the run. (Default: 0)",
                    default = 0 )
##############


##############
MiscellaneousOptions = parser.add_argument_group('Miscellaneous options')
MiscellaneousOptions.add_argument('-threads',  type=int,
                    help = "Number of available threads. (Default 1)",
                    default = 1 )
MiscellaneousOptions.add_argument('-memory',  type=int,
                    help = "Memory available for the assembly in Giga. (Default 1)",
                    default = 1 )
MiscellaneousOptions.add_argument('-time_max',  type=int,
                    help = "Do not begin a new iteration if the job duration (in seconds) has exceed this threshold. (Default 7200)",
                    default = 7200 )
##############


### Option parsing
args = parser.parse_args()

os.system("echo '[Running in process...]\n[Warning messages may print, but they are no error. If real errors appear, the process will stop]'")


### Read the arguments
StartIteration = args.iteration_start
MaxIteration = args.iteration_max

Evalue = args.evalue

MinIdentityPercentage = args.min_id
MinAliLength = args.min_ali_len
MinLength = args.min_len
RequiredCoverage = args.required_coverage

FinalMinLength = args.final_min_len
FinalMinIdentityPercentage = args.final_min_id
FinalMinAliLength = args.final_min_ali_len

KeepIterations = args.keep_iterations
FinishAllIter = args.finish_all_iter

Threads = args.threads
Memory = args.memory
MaxTime = args.time_max

if args.database_type in ["RF","FR","R","F"]:
    StrandedData = True
else:
    StrandedData = False
    
if args.database_type in ["paired","RF","FR"]:
     PairedData = True
else:
     PairedData = False
    
    
if args.plot:
    args.stats = True

### Set up the log directory
if args.log:
    LogDirName = os.path.dirname(args.log)
    if not os.path.isdir(LogDirName) and LogDirName:
        os.makedirs(LogDirName)

### Set up the logger
LogFile = args.log
# create logger with 'spam_application'
logger = logging.getLogger('apytram')
logger.setLevel(logging.DEBUG)
# create file handler which logs even debug messages
fh = logging.FileHandler(LogFile)
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.WARN)
# create formatter and add it to the handlers
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(fh)
logger.addHandler(ch)

### If iteration begin not from 1, the temporary directory must be given by the user
if StartIteration != 1 and not args.tmp:
    logger.error("If you want to restart a previous job, the previous temporary directory must be given.")
    sys.exit(1)    

### Set up the working directory
if args.tmp:
    if os.path.isdir(args.tmp):
        logger.info("The temporary directory %s exists" %(args.tmp) )
    else:
        logger.info("The temporary directory %s does not exist, it will be created" % (args.tmp))
        os.makedirs(args.tmp)
    TmpDirName = args.tmp
else:
    TmpDirName = tempfile.mkdtemp(prefix='tmp_apytram')

### Set up the output directory
if args.output_prefix:
    OutDirName = os.path.dirname(args.output_prefix)
    OutPrefixName = args.output_prefix
    if os.path.isdir(OutDirName):
        logger.info("The output directory %s exists" %(os.path.dirname(args.output_prefix)) )
    elif OutDirName: # if OutDirName is not a empty string we create the directory
        logger.info("The temporary directory %s does not exist, it will be created" % (os.path.dirname(args.output_prefix)))
        os.makedirs(os.path.dirname(args.output_prefix))
else:
    logger.error("The output prefix must be defined")
    sys.exit(1)

### Check that query files exist
if args.query:
    QueryFile = args.query
    AliQueryFile = args.query
    if not os.path.isfile(QueryFile):
        logger.error(QueryFile+" (-q) is not a file.")
        sys.exit(1)
    elif not os.stat(QueryFile).st_size:
        logger.error(QueryFile+" (-q) is empty.")
        sys.exit(1)
    elif ApytramNeeds.count_sequences(QueryFile) !=1:
        logger.warning("%s (-q) contains more than one query. They are %s sequences." %(QueryFile,ApytramNeeds.count_sequences(QueryFile)))
        # If there are multiple probes, align them for the future coverage counter
        # Use Mafft
        start_mafft_time = time.time()
        MafftProcess = Aligner.Mafft(QueryFile)
        MafftProcess.QuietOption = True
        MafftProcess.AutoOption = True
        (MafftResult, err) = MafftProcess.get_output()
        AliQueryFile = "%s/References.ali.fasta" %TmpDirName
        ApytramNeeds.write_in_file(MafftResult,AliQueryFile)
        logger.debug("mafft --- %s seconds ---" % (time.time() - start_mafft_time))


# If the -pep option is used, the -q option must be precised
if args.query_pep:
    if not os.path.isfile(args.query_pep):
        logger.error(args.query_pep+" (-pep) is not a file.")
        sys.exit(1)
    
    if not args.query:
        logger.error("-pep option must be accompanied of the query in nucleotide format (-q option)")
        sys.exit(1)

### Check that there is a database, otherwise build it
DatabaseName = args.database
CheckDatabase_BlastdbcmdProcess = BlastPlus.Blastdbcmd(DatabaseName, "", "")
if not CheckDatabase_BlastdbcmdProcess.is_database():
    logger.info("Database %s does not exist" % DatabaseName)
    #Concatenate input files
    if args.fastq or args.fasta:
        if args.fastq:
            for fastq in args.fastq:
                if not os.path.isfile(fastq):
                    logger.error("The fastq file (%s) does not exist." %fastq)
                    sys.exit(1)
                # Format the fastq file in fasta
            InputFasta = "%s/input_fastq.fasta" %(TmpDirName)
            logger.info("Convert the fastq file in fasta format")
            start_convert = time.time()
            ExitCode = ApytramNeeds.fastq2fasta(" ".join(args.fastq),InputFasta)
            logger.info("Convertion takes %s seconds" %(time.time() - start_convert))
        elif args.fasta:
            for fasta in args.fasta:
                if not os.path.isfile(fasta):
                    logger.error("The fasta file (%s) does not exist." %fasta)
                    sys.exit(1)
                # Concatenate fasta files
            if len(args.fasta) > 1:
                InputFasta = "%s/input_fasta.fasta" %(TmpDirName)
                logger.info("Concatenate fasta files")
                start_convert = time.time()
                ExitCode = ApytramNeeds.cat_fasta(" ".join(args.fasta),InputFasta)
                logger.info("Concatenation takes %s seconds" %(time.time() - start_convert))
            else:
                InputFasta = args.fasta[0]
                
    else :
        logger.error("The database could not be formatted because fasta files (-fa) or fastq files (-fq) are required!")
        sys.exit(1)
    
    #Check if the end of sequence name of paired data are 1 or 2
    if PairedData:
        BadReadName = ApytramNeeds.check_paired_data(InputFasta)
        if BadReadName:
            logger.error("Paired read names must finished by 1 or 2. %s is uncorrect" %BadReadName)
            sys.exit(1)
    #Build blast formated database from a fasta file
    if not os.path.isfile(InputFasta):
        logger.error("Error during concatenation or conversion of input files.")
        sys.exit(1)
    if os.path.isdir(os.path.dirname(DatabaseName)):
        logger.info("Database directory exists")
    else:
        logger.info("Database directory does not exist, we create it")
        os.makedirs(os.path.dirname(DatabaseName))
    
    # Database building
    logger.info(DatabaseName + " database building")
    MakeblastdbProcess = BlastPlus.Makeblastdb(InputFasta,DatabaseName)
    out,err = MakeblastdbProcess.launch()

CheckDatabase_BlastdbcmdProcess = BlastPlus.Blastdbcmd(DatabaseName, "", "")
if not CheckDatabase_BlastdbcmdProcess.is_database():
    logger.error("Problem in the database building.\nAre you sure of your input format?\nAre all read names unique?")
    logger.info("Database %s does not exist" % DatabaseName)
    sys.exit(1)
else:
    logger.info("Database %s exists" % DatabaseName)

### If there is a query continue, else stop
if not args.query:
    logger.info("There is no query (-q), apytram has finished.")
    quit()
elif not os.path.isfile(args.query):
    logger.error(args.query+" (-q) is not a file.")
    sys.exit(1)
else:
    logger.info("DB: \"%s\"\tQuery: \"%s\"" %(DatabaseName,QueryFile))

### Make iterations
# Initialisation
i = 0
Stop = False
BaitSequences = QueryFile
IterationNotFinished = False

StatsDict = {0:{"IterationTime": 0,
                "CumulTime": time.time() - start_time,
                "LargeCoverage": 0,
                "StrictCoverage": 0,
                "NbContigs": 0,
                "AverageLength": 0,
                "TotalLength": 0,
                "BestLength":0,
                "AverageScore": 0,
                "TotalScore": 0,
                "BestScore":0,
                "AverageIdentity": 0,
                "TotalIdentity": 0,
                "BestIdentity":0,
                "ReadsNumber":0,
                "BlastTime": 0,
                "TrinityTime": 0,
                "Exonerate1Time":0,
                "Exonerate2Time":0,
                "MafftTime":0,
                "PythonTime":time.time() - start_time
                }}
    
logger.info("Iterations begin")
start_iter = time.time()

### If apytram restart a job
if StartIteration != 1 :
    for i in range(1,StartIteration):
        logger.debug("Iteration %s has already been executed" %i)
        StatsDict[i] = StatsDict[(i-1)].copy()
    Reali = i
    ### Last Trinity filtered file is the baitfile
    BaitSequences = "%s/Trinity_iter_%d.filtered.fasta" % (TmpDirName, i)
    ### Chech that the bait file exists and is not empty
    if not os.path.isfile(BaitSequences):
        logger.error("%s does not exit, apytram can not restart this job at the iteration %s because the iteration %s is not present in the temporary directory %s" %(BaitSequences,i+1,i,TmpDirName))
        exit(1)
    elif not os.stat(BaitSequences).st_size:
        logger.error("%s is empty, apytram can not restart this job at the iteration %sbecause the iteration %s is not present in the temporary directory %s" %(BaitSequences,i+1,i,TmpDirName))
        exit(1)
        
while (i < MaxIteration) and (Stop == False):
    start_iter_i = time.time()
    i += 1
    StatsDict[i] = StatsDict[(i-1)].copy()
    StatsDict[i].update({"IterationTime": 0,
                "CumulTime": 0,
                "BlastTime": 0,
                "TrinityTime": 0,
                "Exonerate1Time":0,
                "Exonerate2Time":0,
                "MafftTime":0,
                "PythonTime":0,
                })
    
    logger.info("Iteration %d/%d" %(i,MaxIteration))
    
    ### Blast bait sequences on database of reads
    
    logger.info("Blast bait sequences on reads database")
    start_blast_time = time.time()
    ReadNamesFile = "%s/ReadNames.%d.txt" % (TmpDirName,i)

    if args.query_pep and i == 1:
        BlastnProcess = BlastPlus.Blast("tblastn", DatabaseName, args.query_pep)
        BlastnProcess.Evalue = Evalue
    else:
        BlastnProcess = BlastPlus.Blast("blastn", DatabaseName, BaitSequences)
        BlastnProcess.Evalue = Evalue
        BlastnProcess.Task = "blastn"

    BlastnProcess.Threads = Threads
    BlastnProcess.OutFormat = "6 sacc"

    # Write read names in ReadNamesFile if the file does not exist
    if not os.path.isfile(ReadNamesFile):
        (out,err) = BlastnProcess.launch(ReadNamesFile)
    else:
        logger.warn("%s has already been created, it will be used" %ReadNamesFile )
        
    StatsDict[i]["BlastTime"] = time.time() - start_blast_time
    logger.debug("blast --- %s seconds ---" % (StatsDict[i]["BlastTime"]))
    if PairedData:
        # Get paired reads names and remove duplicated names
        logger.info("Get paired reads names and remove duplicated names")
        ExitCode = ApytramNeeds.add_paired_read_names(ReadNamesFile)
    else:
        # Remove duplicated names
        logger.info("Remove duplicated names")
        ExitCode = ApytramNeeds.remove_duplicated_read_names(ReadNamesFile)
    
    # Count the number of reads which will be used in the Trinity assembly
    logger.info("Count the number of reads")
    StatsDict[i]["ReadsNumber"] = ApytramNeeds.count_lines(ReadNamesFile)
    
    if not StatsDict[i]["ReadsNumber"]:
        logger.error("No read recruted by Blast at the iteration %s" %i)
        Stop = True
        IterationNotFinished = True
        i -= 1
    else:    
        # Compare the read list names with the list of the previous iteration:
        Identical = ApytramNeeds.are_identical(ReadNamesFile,"%s/ReadNames.%d.txt" % (TmpDirName,i-1))
        if Identical and not FinishAllIter:
            logger.info("Reads from the current iteration are identical to reads from the previous iteration")
            Stop = True
            IterationNotFinished = True
            i -= 1
        else:

            ### Retrieve sequences

            logger.info("Split read names depending on 1/ or 2/")
            if args.database_type in ["RF","FR"]:
                ReadNamesFile_Right = "%s/ReadNames.%d.1.txt" % (TmpDirName,i)
                ReadNamesFile_Left = "%s/ReadNames.%d.2.txt" % (TmpDirName,i)
                ReadFasta_Right = "%s/Reads.%d.1.fasta" % (TmpDirName,i)
                ReadFasta_Left = "%s/Reads.%d.2.fasta" % (TmpDirName,i)                
                ApytramNeeds.split_readnames_in_right_left(ReadNamesFile,ReadNamesFile_Right,ReadNamesFile_Left)
                StrandList = [".1",".2"]              
            else:
                StrandList = [""] 
                
            logger.info("Retrieve reads sequences")
            start_blastdbcmd_time = time.time()                
            for strand in StrandList:
                ReadFasta = "%s/Reads.%d%s.fasta" % (TmpDirName,i,strand)
                ReadNamesFile = "%s/ReadNames.%d%s.txt" % (TmpDirName,i,strand)
                BlastdbcmdProcess = BlastPlus.Blastdbcmd(DatabaseName, ReadNamesFile, ReadFasta)
                if not os.path.isfile(ReadFasta):
                    (out,err) = BlastdbcmdProcess.launch()
                else:
                    logger.warn("%s has already been created, it will be used" %(ReadFasta) ) 
                        
            
            StatsDict[i]["BlastTime"] += time.time() - start_blastdbcmd_time
            logger.debug("blastdbcmd --- %s seconds ---" %(time.time() - start_blastdbcmd_time))

            ### Launch Trinity

            start_trinity_time = time.time()
            logger.info("Launch Trinity")
            ExitCode = 0
            TrinityFasta = "%s/Trinity_iter_%d" %(TmpDirName, i)
            if StrandedData:
                if args.database_type in ["RF","FR"]:
                    TrinityProcess = Trinity.Trinity(TrinityFasta, right = ReadFasta_Right,
                                                    left = ReadFasta_Left)
                else:
                    TrinityProcess = Trinity.Trinity(TrinityFasta, single = ReadFasta)
                TrinityProcess.SS_lib_type = args.database_type
            else:
                TrinityProcess = Trinity.Trinity(TrinityFasta, single = ReadFasta)
                if PairedData:
                    TrinityProcess.RunAsPaired = True
            # If there is a huge number of reads, remove duplicated reads
            if StatsDict[i]["ReadsNumber"] > 1000:
                TrinityProcess.NormalizeReads
                
            TrinityProcess.CPU = Threads
            TrinityProcess.max_memory = Memory
            # Keep only contig with a length superior to MinLength
            TrinityProcess.MinLength = MinLength

            # Use the  --full_cleanup Trinity option to keep only the contig file
            TrinityProcess.FullCleanup = True
            if not os.path.isfile(TrinityFasta+".Trinity.fasta"):
                (out,err,ExitCode) = TrinityProcess.launch()
            else:
                logger.warn("%s has already been created, it will be used" %(TrinityFasta+".Trinity.fasta") ) 

            TrinityFasta = TrinityFasta + ".Trinity.fasta"
            StatsDict[i]["TrinityTime"] = time.time() - start_trinity_time
            logger.debug("trinity --- %s seconds ---" %(StatsDict[i]["TrinityTime"]))
            if not os.path.isfile(TrinityFasta): # Trinity found nothing
                logger.debug("Trinity found nothing...\n[...]\n"+"\n".join(out.strip().split("\n")[-15:]))
                if ExitCode == 2 or ExitCode == 0 : # Trinity exit 0 if "No butterfly assemblies to report"
                    logger.error("Trinity has assembled no contigs at the end of the iteration %s (ExitCode: %d)" %(i,ExitCode) )
                elif ExitCode != 0:
                    logger.error("Trinity has crashed (ExitCode: %d). Are all dependencies satisfied?" %ExitCode)
                Stop = True
                IterationNotFinished = True
                i -=1
            else:

                ### Filter Trinity contigs to keep only homologous sequences of the reference genes

                logger.info("Compare Trinity results with query sequences")
                # Use Exonerate 
                TrinityExonerate = "%s/Trinity_iter_%d.exonerate_cdna2g" % (TmpDirName, i)
                start_exo_time = time.time()
                TrinityExonerateProcess = Aligner.Exonerate(QueryFile,TrinityFasta)
                # Keep only the best hit for each contig from Trinity 
                TrinityExonerateProcess.Bestn = 1 
                TrinityExonerateProcess.Model = "cdna2genome"
                # Customize the output format
                TrinityExonerateProcess.Ryo = "%ti\t%qi\t%ql\t%tal\t%tl\t%tab\t%tae\t%s\t%pi\t%qab\t%qae\n"
                (out,err,TrinityExonerateResult) = TrinityExonerateProcess.get_output()
                # Write the result in a file
                TrinityExonerateFile = open(TrinityExonerate,"w")
                TrinityExonerateFile.write(TrinityExonerateResult)
                TrinityExonerateFile.close()
                if not TrinityExonerateResult:
                    logger.info("Reconstructed sequences but no homologous with references")
                    logger.info("Try to get homologies with a more sensible model")
                    ### Try to get homologies with a more sensible model
                    TrinityExonerate = "%s/Trinity_iter_%d.exonerate_coding2g" % (TmpDirName, i)
                    start_exo_time = time.time()
                    TrinityExonerateProcess = Aligner.Exonerate(QueryFile,TrinityFasta)
                    # Keep only the best hit for each contig from Trinity 
                    TrinityExonerateProcess.Bestn = 1 
                    TrinityExonerateProcess.Model = "coding2genome"
                    # Customize the output format
                    TrinityExonerateProcess.Ryo = "%ti\t%qi\t%ql\t%tal\t%tl\t%tab\t%tae\t%s\t%pi\t%qab\t%qae\n"
                    (out,err,TrinityExonerateResult) = TrinityExonerateProcess.get_output()
                    # Write the result in a file
                    TrinityExonerateFile = open(TrinityExonerate,"w")
                    TrinityExonerateFile.write(TrinityExonerateResult)
                    TrinityExonerateFile.close()

                if not TrinityExonerateResult:
                    logger.info("Reconstructed sequences but no homologous with references (even with the more sensible model)")
                    Stop = True
                    IterationNotFinished = True
                    i -= 1
                else:
                    # Keep only sequence with a identity percentage > MinIdentitypercentage on the whole hit
                    BestScoreNames, ReverseNames, TrinityExonerateResultsDict, StatsIter = ApytramNeeds.parse_exonerate_results(TrinityExonerateResult, MinIdentityPercentage,
                                             minalilength = MinAliLength)
                    StatsDict[i].update(StatsIter)
                    FilteredSequenceNames = TrinityExonerateResultsDict.keys()
                    StatsDict[i]["Exonerate1Time"] = time.time() - start_exo_time
                    logger.debug("exonerate on trinity --- %s seconds ---" % (StatsDict[i]["Exonerate1Time"]))

                    # Write filtered sequences in a file
                    logger.info("Filter sequence with a identity percentage superior to %d and a alignment len %d" %(MinIdentityPercentage, MinAliLength)) 
                    FileteredTrinityFasta =  "%s/Trinity_iter_%d.filtered.fasta" % (TmpDirName, i)
                    ExitCode = ApytramNeeds.filter_fasta(TrinityFasta, FilteredSequenceNames, FileteredTrinityFasta,
                                                        ReverseNames = ReverseNames)

                    ### Validated sequences become bait sequences

                    BaitSequences = FileteredTrinityFasta
                    if not os.stat(BaitSequences).st_size:
                        logger.error("No sequence has passed the iteration filter at the iteration %s" %(i))
                        Stop = True
                        IterationNotFinished = True
                        i -=1
                    else:
                        ### Compare sequences of the current iteration to those of the previous iteration

                        logger.info("Compare results with the previous iteration")

                        #Check if the number of contigs has changed

                        logger.info("Check if the number of contigs has changed")
                        StatsDict[i]["NbContigs"] = len(FilteredSequenceNames)

                        if StatsDict[i]["NbContigs"] != StatsDict[i-1]["NbContigs"]:
                             logger.info("The number of contigs has changed")          
                        elif i >= 2:
                            logger.info("Refind the \"parent\" contig from the previous contig for each contig and check they are different")
                            # Use Exonerate 
            # Comm Marie : pas compris le coup des 2 exonerate (par rapport à celui d'avant) -- en vrai j'ai compris et indiqué dans readme mais peux tu etre un peu plus explicite ici?

                            start_exo_time = time.time()
                            ExonerateProcess = Aligner.Exonerate(FileteredTrinityFasta, "%s/Trinity_iter_%d.filtered.fasta" % (TmpDirName,i-1) )
                            # Keep only the best hit for each contigs
                            Exonerate = "%s/iter_%d_%d.exonerate" % (TmpDirName, i-1, i)
                            ExonerateProcess.Bestn = 1
                            ExonerateProcess.Model =  "est2genome"
                            # Customize the output format
                            ExonerateProcess.Ryo = "%ti\t%qi\t%ql\t%qal\t%tal\t%tl\t%pi\n"
                            (out,err,ExonerateResult) = ExonerateProcess.get_output()
                            ExonerateFile = open(Exonerate,"w")
                            ExonerateFile.write(ExonerateResult)
                            ExonerateFile.close()
                            AlmostIdenticalResults = ApytramNeeds.check_almost_identical_exonerate_results(ExonerateResult)
                            if AlmostIdenticalResults and not FinishAllIter:
                                logger.info("Contigs are almost identical than the previous iteration (Same size (~98%), > 99% identity)")
                                Stop =True
                            StatsDict[i]["Exonerate2Time"] = time.time() - start_exo_time
                            logger.debug("exonerate on previous iter --- %s seconds ---" % (StatsDict[i]["Exonerate2Time"]))

                        # Check that the coverage has increased compared to the previous iteration

                        logger.info("Check that the coverage has inscreased compared to the previous iteration")
                        # Use Mafft
                        start_mafft_time = time.time()
                        MafftProcess = Aligner.Mafft(AliQueryFile)
                        MafftProcess.QuietOption = True
                        MafftProcess.AutoOption = True
                        #MafftProcess.AdjustdirectionOption = True
                        MafftProcess.AddOption = FileteredTrinityFasta
                        (MafftResult,err) = MafftProcess.get_output()
                        StatsDict[i]["StrictCoverage"], StatsDict[i]["LargeCoverage"], DicPlotCov = ApytramNeeds.calculate_coverage(MafftResult)
                        logger.info("Strict Coverage: %s\tLarge Coverage: %s" %(StatsDict[i]["StrictCoverage"], StatsDict[i]["LargeCoverage"]))
                        StatsDict[i]["MafftTime"] = time.time() - start_mafft_time
                        logger.debug("mafft --- %s seconds ---" % (StatsDict[i]["MafftTime"]))

                        if not FinishAllIter:
                            # Stop iteration if both Largecoverage and Total length are not improved
                            if StatsDict[i]["AverageLength"] > StatsDict[i-1]["AverageLength"]:
                                pass
                            elif StatsDict[i]["AverageScore"] > StatsDict[i-1]["AverageScore"]:
                                pass
                            elif StatsDict[i]["TotalLength"] > StatsDict[i-1]["TotalLength"]:
                                pass
                            elif StatsDict[i]["TotalScore"] > StatsDict[i-1]["TotalScore"]:
                                pass
                            elif StatsDict[i]["BestScore"] > StatsDict[i-1]["BestScore"]:
                                pass
                            elif StatsDict[i]["LargeCoverage"] <= StatsDict[i-1]["LargeCoverage"]:
                                logger.info("This iteration have a large coverage inferior (or equal) to the previous iteration")
                                Stop = True

                            # Stop iteration if the RequiredCoverage is reached
                            if StatsDict[i]["StrictCoverage"] >= RequiredCoverage:
                                logger.info("This iteration attains the required bait sequence coverage (%d >= %d)" % (StatsDict[i]["StrictCoverage"],RequiredCoverage))
                                Stop = True

                        ### Write a fasta file for this iteration if the option --keep_iterations was selected

                        if KeepIterations:
                            if not args.no_best_file:
                            # Best sequences of the iteration
                                ExitCode = ApytramNeeds.write_apytram_output(FileteredTrinityFasta, TrinityExonerateResultsDict,
                                                                 "%s.iter_%d.best.fasta" %(OutPrefixName,i), 
                                                                 Header = TrinityExonerateProcess.Ryo.replace('%',"").replace("\n","").split(),
                                                                 Names = BestScoreNames.values(),
                                                                 Message = "iter_%d.best." %i)
                            # All sequences of the iteration
                            ExitCode = ApytramNeeds.write_apytram_output(FileteredTrinityFasta,
                                                             TrinityExonerateResultsDict,
                                                             "%s.iter_%d.fasta" %(OutPrefixName,i),
                                                             Header = TrinityExonerateProcess.Ryo.replace('%',"").replace("\n","").split(),
                                                             Message = "iter_%d." %i)
                            # Mafft alignment
                            ApytramNeeds.write_in_file(MafftResult,"%s.iter_%s.ali.fasta" %(OutPrefixName,i))

    if IterationNotFinished:
            logger.debug("Iteration stop before end")
            Reali = i + 1
    else: 
        Reali = i
    
    NoPythonTime = StatsDict[Reali]["BlastTime"] + StatsDict[Reali]["TrinityTime"] +\
                 StatsDict[Reali]["MafftTime"] + StatsDict[Reali]["Exonerate1Time"] +\
                 StatsDict[Reali]["Exonerate2Time"]
    StatsDict[Reali].update({"IterationTime": time.time() - start_iter_i,
                             "CumulTime": time.time() - start_time,
                             "PythonTime": time.time() - start_iter_i - NoPythonTime })
    logger.debug("iteration %d --- %s seconds ---" % (Reali, time.time() - start_iter_i))
    
    if (time.time() - start_time) > MaxTime and Stop == False:
        logger.warn("No new iteration will begin because the maximum duration (%s seconds) of the job is attained. (%s seconds)" %(MaxTime, (time.time() - start_time)))
        Stop = True



logger.info("End of Iterations. Iterative process takes %s seconds." %(time.time() - start_iter))
start_output = time.time()
if i: #We check that there is at least one iteration with a result
    if FinalMinLength or FinalMinIdentityPercentage or FinalMinAliLength:
        start_iter_i = time.time()
        Reali +=1
        #### Final filter which is equivalent at a new iteration
        StatsDict[Reali] = StatsDict[Reali-1].copy()
        StatsDict[Reali].update({"IterationTime": 0,
                    "CumulTime": 0,
                    "BlastTime": 0,
                    "TrinityTime": 0,
                    "Exonerate1Time":0,
                    "Exonerate2Time":0,
                    "MafftTime":0,
                    "PythonTime":0,
                    })
 
        # Keep only sequence with a identity percentage > FinalMinIdentitypercentage on the whole hit
        BestScoreNames, ReverseNames, TrinityExonerateResultsDict, StatsIter = ApytramNeeds.parse_exonerate_results(TrinityExonerateResult,
                                         FinalMinIdentityPercentage,
                                         minalilengthpercentage = FinalMinAliLength,
                                         minlengthpercentage = FinalMinLength)
        StatsDict[Reali].update(StatsIter)
        FilteredSequenceNames = TrinityExonerateResultsDict.keys()
        logger.info("Filter sequence with a identity percentage superior to %d and a percentage alignment len %d" %(FinalMinIdentityPercentage, FinalMinAliLength))
        
        if FilteredSequenceNames: # If sequences pass the last filter
            # Write Filter hit
            FileteredTrinityFasta =  "%s/Trinity_iter_%d.filtered.fasta" % (TmpDirName, Reali+1)
            ExitCode = ApytramNeeds.filter_fasta(TrinityFasta, FilteredSequenceNames,
                                                 FileteredTrinityFasta, ReverseNames = ReverseNames)

    start_output = time.time()
    if FilteredSequenceNames: # If sequences pass the last filter
        StatsDict[Reali]["NbContigs"] = len(FilteredSequenceNames)
        #### Write output files
        logger.info("Write outputfiles")
        if not args.no_best_file:
            # Best sequences
            ExitCode = ApytramNeeds.write_apytram_output(FileteredTrinityFasta, TrinityExonerateResultsDict,
                                                         OutPrefixName+".best.fasta", 
                                                         Header = TrinityExonerateProcess.Ryo.replace('%',"").replace("\n","").split(),
                                                         Names = BestScoreNames.values(),
                                                         Message = "best_")
        if not args.only_best_file:
            # Last iteration seqeunces
            ExitCode = ApytramNeeds.write_apytram_output(FileteredTrinityFasta,
                                                         TrinityExonerateResultsDict,
                                                         OutPrefixName+".fasta",
                                                         Header = TrinityExonerateProcess.Ryo.replace('%',"").replace("\n","").split(),
                                                         Names = FilteredSequenceNames)
        
        if args.plot_ali or args.stats:
            ### Calculate the coverage
            logger.info("Calculate the final coverage")
            # Use Mafft
            start_mafft_time = time.time()
            MafftProcess = Aligner.Mafft(AliQueryFile)
            MafftProcess.QuietOption = True
            MafftProcess.AutoOption = True
            #MafftProcess.AdjustdirectionOption = True
            MafftProcess.AddOption = OutPrefixName+".fasta"
            (MafftResult,err) = MafftProcess.get_output()
            StatsDict[Reali]["StrictCoverage"], StatsDict[Reali]["LargeCoverage"], DicPlotCov = ApytramNeeds.calculate_coverage(MafftResult)
            logger.info("Strict Coverage: %s\tLarge Coverage: %s" %(StatsDict[Reali]["StrictCoverage"], StatsDict[Reali]["LargeCoverage"]))
            StatsDict[Reali]["MafftTime"] += time.time() - start_mafft_time
            logger.debug("mafft --- %s seconds ---" % (time.time() - start_mafft_time))

            NoPythonTime = StatsDict[Reali]["BlastTime"] + StatsDict[Reali]["TrinityTime"] +\
                 StatsDict[Reali]["MafftTime"] + StatsDict[Reali]["Exonerate1Time"] +\
                 StatsDict[Reali]["Exonerate2Time"]
            StatsDict[Reali].update({"IterationTime": time.time() - start_iter_i,
                                         "CumulTime": time.time() - start_time,
                                         "PythonTime": time.time() - start_iter_i - NoPythonTime })

    # Stats files
    
    if args.plot_ali:
        start_output_ali = time.time()
        LengthAlignment = len(DicPlotCov[DicPlotCov.keys()[0]])
        if LengthAlignment <= 3100:
            logger.info("Create plot of the final alignment (OutPrefix.ali.png)")         
            ApytramNeeds.create_plot_ali(DicPlotCov, OutPrefixName)
        else:
            logger.warn("Final alignment is longger than 3100 pb, the plot of the final alignment (OutPrefix.ali.png) can NOT be created. See the final alignement (OutPrefix.ali.fasta).")         
        logger.info("Write the final alignment in OutPrefix.ali.fasta")
        ApytramNeeds.write_in_file(MafftResult,"%s.ali.fasta" %OutPrefixName)
        logger.debug("Writing alignment plot and fasta --- %s seconds ---" % (time.time() - start_output_ali))
        
else:
    logger.warn("No results")

if args.stats:
    start_output_stat = time.time()
    logger.info("Write statistics file (OutPrefix.stats.csv)")
    ApytramNeeds.write_stats(StatsDict,OutPrefixName)

    if args.plot:
        logger.info("Create plot from the statistics file (OutPrefix.stats.pdf)")
        ApytramNeeds.create_plot(StatsDict, OutPrefixName)
        logger.debug("Writing stats file --- %s seconds ---" % (time.time() - start_output_stat))


logger.debug("Writing outputs --- %s seconds ---" % (time.time() - start_output))
        
### Remove tempdir if the option --tmp have not been use
if not args.tmp:
    logger.debug("Remove the temporary directory")
    #Remove the temporary directory :
    if "tmp_apytram" in TmpDirName:
        shutil.rmtree(TmpDirName)

logger.info("--- %s seconds ---" % (time.time() - start_time))

