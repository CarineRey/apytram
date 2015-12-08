#!/usr/bin/python
# coding: utf-8
import os
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


requiredOptions = parser.add_argument_group('Required arguments')
requiredOptions.add_argument('-d', '--database', nargs='?', type=str,
                             help='Database preffix name. If a database with the same name already exist, the existing database will be kept and the database will NOT be rebuilt.', required=True)
requiredOptions.add_argument('-dt', '--database_type', type=str, choices=["single","paired"],
                             help='single or paired end RNA-seq data', required=True)



InOutOptions = parser.add_argument_group('Input and Output Files')

InOutOptions.add_argument('-fa', '--fasta',  type=str,
                   help = "Fasta formatted RNA-seq data to build the database of reads")
InOutOptions.add_argument('-fq', '--fastq',  type=str,
                   help = "Fastq formatted RNA-seq data to build the database of reads. (The fastq will be first convert in fasta file")
InOutOptions.add_argument('-q', '--query',  type=str,
                    help = "Fasta file (nucl) with bait sequence for the apytram run. If no query is submitted, the database will be only build." )
InOutOptions.add_argument('-pep', '--query_pep',  type=str,
                   default = "",       
                   help = "Fasta file containing the query in the peptide format. It will be used at the first iteration as bait sequences to fish reads. It must be accompany of the query in nucleotide format (-q option)")
InOutOptions.add_argument('-out', '--output_preffix',  type=str, default = "./apytram",
                   help = "Output preffix (Default ./apytram)")
InOutOptions.add_argument('-tmp',  type=str,
                    help = "Directory to stock all intermediary files for the apytram run. (default: a directory in /tmp which will be removed at the end)",
                    default = "" )
InOutOptions.add_argument('--keep_iterations',  action='store_true',
                    help = "A fasta file will be created at each iteration.")
InOutOptions.add_argument('--stats', action='store_true',
                             help='Create files with statistics on each iteration')
InOutOptions.add_argument('--plot', action='store_true',
                             help='Create file with plot on statistics on each iteration')
InOutOptions.add_argument('-log', type=str, default="apytram.log",
                   help = "a log file to report avancement (default: apytram.log)")


SearchOptions = parser.add_argument_group('Arguments for search thresholds')

SearchOptions.add_argument('-i', '--iteration_max',  type=int,
                    help = "Maximum number of iteration. (Default 5)",
                    default = 5 )
InOutOptions.add_argument('--finish_all_iter',  action='store_true',
                    default = False,
                    help = "apytram will finish all iteration (-i) even if there is no improvment) (default: False)")
SearchOptions.add_argument('-e', '--evalue',  type=float,
                    help = "Evalue threshold of the blastn of the bait queries on the database of reads . (Default 1e-3)",
                    default = 1e-3 )
SearchOptions.add_argument('--required_coverage',  type=float,
                    help = "Required coverage of a bait sequence to stop iteration (Default: No threshold)",
                    default = 200 )
SearchOptions.add_argument('-id', '--min_id',  type=int,
                    help = "Minimum identity percentage with a query to keep a sequence at the end of a iteration (Default 20)",
                    default = 20 )
SearchOptions.add_argument('-l', '--min_len',  type=int,
                    help = "Minimum length to keep a sequence at the end of a iteration (Default 200)",
                    default = 200 )

MiscellaneousOptions = parser.add_argument_group('Miscellaneous options')
MiscellaneousOptions.add_argument('--threads',  type=int,
                    help = "Number of available threads. (Default 1)",
                    default = 1 )

### Option parsing
args = parser.parse_args()

## Example for the dev
#args = parser.parse_args(''' -d example_exec/db/examplefq
#                             -out Out_test/apytram
#                             -fq example/example_db.fastq
#                            -q example/ref_gene.fasta 
#                             --database_type paired
#                             --keep_iterations
#                             --plot
#                             -i 5'''.split())


### Arguments reading
MaxIteration = args.iteration_max
Threads = args.threads
Evalue = args.evalue
MinIdentityPercentage = args.min_id
MinLength = args.min_len
KeepIterations = args.keep_iterations
RequiredCoverage = args.required_coverage
FinishAllIter = args.finish_all_iter

if args.database_type == "paired":
    PairedData = True
else:
    PairedData = False
    
if args.plot:
    args.stats = True
    


### Set up the output directory
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
ch.setLevel(logging.ERROR)
# create formatter and add it to the handlers
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(fh)
logger.addHandler(ch)

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
if args.output_preffix:
    OutDirName = os.path.dirname(args.output_preffix)
    OutPreffixName = args.output_preffix
    if os.path.isdir(OutDirName):
        logger.info("The output directory %s exists" %(os.path.dirname(args.output_preffix)) )
    elif OutDirName: # if OutDirName is not a empty string we create the directory
        logger.info("The temporary directory %s does not exist, it will be created" % (os.path.dirname(args.output_preffix)))
        os.makedirs(os.path.dirname(args.output_preffix))
else:
    logger.error("The output preffix must be defined")
    sys.exit(1)

### Check that query files exist

if args.query:
    if not os.path.isfile(args.query):
        logger.error(args.query+" (-q) is not a file.")
        sys.exit(1)

# If the -pep option is used, the -q option must be precise
if args.query_pep:
    if not os.path.isfile(args.query_pep):
        logger.error(args.query_pep+" (-pep) is not a file.")
        sys.exit(1)
    
    if not args.query:
        logger.error("-pep option must be accompanied of the query in nucleotide format (-q option)")
        sys.exit(1)

### Check that there is a database else built it
DatabaseName = args.database
CheckDatabase_BlastdbcmdProcess = BlastPlus.Blastdbcmd(DatabaseName, "", "")
if not CheckDatabase_BlastdbcmdProcess.is_database():
    logger.info("Database %s does not exist" % DatabaseName)
    #Build blast formated database from a fasta file
    if args.fastq or args.fasta:
        if args.fastq:
            if not os.path.isfile(args.fastq):
                logger.error("The fastq file (-fq) does not exist.")
                sys.exit(1)
            else:
                # Format the fastq file in fasta
                InputFasta = TmpDirName + "/" + os.path.basename(args.fastq) + ".fasta"
                logger.info("Convert the fastq file in fasta format")
                ExitCode = ApytramNeeds.fastq2fasta(args.fastq,InputFasta)
        elif args.fasta:
            InputFasta = args.fasta
            
        if not os.path.isfile(InputFasta):
            logger.error("The fasta file (-fa) does not exist.")
            sys.exit(1)
        if os.path.isdir(os.path.dirname(DatabaseName)):
            logger.info("Database directory exists")
        else:
            logger.info("Database directory does not exist, we create it")
            os.makedirs(os.path.dirname(DatabaseName))
        # database building
        logger.info(DatabaseName + " database building")
        MakeblastdbProcess = BlastPlus.Makeblastdb(InputFasta,DatabaseName)
        ExitCode = MakeblastdbProcess.launch()
    else :
        logger.error("The database is not formatted ! A fasta file (-fa) or a fastq file (-fq) is required !")
        sys.exit(1)

CheckDatabase_BlastdbcmdProcess = BlastPlus.Blastdbcmd(DatabaseName, "", "")
if not CheckDatabase_BlastdbcmdProcess.is_database():
    logger.error("Problem in the database building")
    logger.info("Database %s does not exist" % DatabaseName)
    sys.exit(1)
else:
    logger.info("Database %s exists" % DatabaseName)

### If there is a query continue, else stop
if not args.query:
    logger.info("There is no query (-q), apytram have finished.")
    quit()
elif not os.path.isfile(args.query):
    logger.error(args.query+" (-q) is not a file.")
    sys.exit(1)
else:
    QueryFile = args.query
    logger.info("DB: \"%s\"\tQuery: \"%s\"" %(DatabaseName,QueryFile))
    
    

### Make iterations
# Initialisation
i = 0
Stop = False
BaitSequences = QueryFile
IterationNotFinished = False

StatsDict = {"iter_0":{"IterationTime": 0,
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
                "BlastTime": 0,
                "TrinityTime": 0,
                "Exonerate1Time":0,
                "Exonerate2Time":0,
                "MafftTime":0,
                "PythonTime":0
                }}
    
logger.info("Iterations begin")
start_iter = time.time()
while (i < MaxIteration) and (Stop == False):
    start_iter_i = time.time()
    i+=1
    StatsDict["iter_%d" %(i)] = StatsDict["iter_%d" %(i-1)].copy()
    StatsDict["iter_%d" %(i)].update({"IterationTime": 0,
                "CumulTime": 0,
                "BlastTime": 0,
                "TrinityTime": 0,
                "Exonerate1Time":0,
                "Exonerate2Time":0,
                "PythonTime":0
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
    # Write read names in ReadNamesFile
    ExitCode = BlastnProcess.launch(ReadNamesFile)
    StatsDict["iter_%d" %(i)]["BlastTime"] = time.time() - start_blast_time
    logger.debug("blast --- %s seconds ---" % (StatsDict["iter_%d" %(i)]["BlastTime"]))
    if PairedData:
        # Get paired reads names and remove duplicated names
        ExitCode = ApytramNeeds.add_paired_read_names(ReadNamesFile)
    else:
        # Remove duplicated names
        ExitCode = ApytramNeeds.remove_duplicated_read_names(ReadNamesFile)
        
    # Compare the read list names with the list of the previous iteration:
    Identical = ApytramNeeds.are_identical(ReadNamesFile,"%s/ReadNames.%d.txt" % (TmpDirName,i-1))
    if Identical and not FinishAllIter:
        logger.info("Reads from the current iteration are identical than the previous")
        Stop = True
        IterationNotFinished = True
        i -= 1
    else:
        
        ### Retrieve sequences
        
        logger.info("Retrieve sequences")
        ReadFasta = TmpDirName + "/Reads.%d.fasta" % (i)
        BlastdbcmdProcess = BlastPlus.Blastdbcmd(DatabaseName, ReadNamesFile, ReadFasta)
        BlastdbcmdProcess.launch()
        
        ### Launch Trinity
        
        start_trinity_time = time.time()
        logger.info("Launch Trinity")
        TrinityFasta = "%s/Trinity_iter_%d" %(TmpDirName, i)
        TrinityProcess = Trinity.Trinity(ReadFasta,TrinityFasta)
        TrinityProcess.CPU = Threads
        # Keep only contig with a length superior to MinLength
        TrinityProcess.MinLength = MinLength
        if PairedData:
            TrinityProcess.Paired = True
        # Use the  --full_cleanup Trinity option to keep only the contig file
        ExitCode = 0
        TrinityProcess.FullCleanup = True
        ExitCode = TrinityProcess.launch()
        TrinityFasta = TrinityFasta + ".Trinity.fasta"
        StatsDict["iter_%d" %(i)]["TrinityTime"] = time.time() - start_trinity_time
        logger.debug("trinity --- %s seconds ---" %(StatsDict["iter_%d" %(i)]["TrinityTime"]))
        if ExitCode != 0: # Trinity found nothing
            logger.error("Trinity found nothing (ExitCode: %d)" %ExitCode)
            Stop = True
            IterationNotFinished = True
            i -=1
        else:
            
            ### Filter Trinity contigs to keep only homologous sequences of the reference genes
            
            logger.info("Compare Trinity results with query sequences")
            # We use Exonerate 
            TrinityExonerate = "%s/Trinity_iter_%d.exonerate" % (TmpDirName, i)
            start_exo_time = time.time()
            TrinityExonerateProcess = Aligner.Exonerate(QueryFile,TrinityFasta)
            # We want to keep only the best hit for each Trinity sequences
            TrinityExonerateProcess.Bestn = 1
            TrinityExonerateProcess.Model = "cdna2genome"
            # We customize our output format
            TrinityExonerateProcess.Ryo = "%ti\t%qi\t%ql\t%tal\t%tl\t%tab\t%tae\t%s\t%pi\t%qab\t%qae\n"
            TrinityExonerateResult = TrinityExonerateProcess.get_output()
            # We write the result in a file
            TrinityExonerateFile = open(TrinityExonerate,"w")
            TrinityExonerateFile.write(TrinityExonerateResult)
            TrinityExonerateFile.close()
            # Keep only sequence with a identity percentage > MinIdentitypercentage on the whole hit
            BestScoreNames, TrinityExonerateResultsDict, StatsIter = ApytramNeeds.parse_exonerate_results(TrinityExonerateResult, MinIdentityPercentage)
            StatsDict["iter_%d" %(i)].update(StatsIter)
            FilteredSequenceNames = TrinityExonerateResultsDict.keys()
            StatsDict["iter_%d" %(i)]["Exonerate1Time"] = time.time() - start_exo_time
            logger.debug("exonerate on trinity --- %s seconds ---" % (StatsDict["iter_%d" %(i)]["Exonerate1Time"]))
            
            # Filter hit
            logger.info("Filter sequence with a identity percentage superior to %d" %(MinIdentityPercentage)) 
            FileteredTrinityFasta =  "%s/Trinity_iter_%d.filtered.fasta" % (TmpDirName, i)
            ExitCode = ApytramNeeds.filter_fasta(TrinityFasta, FilteredSequenceNames, FileteredTrinityFasta)
            
            ### Validated sequences become bait sequences
            
            BaitSequences = FileteredTrinityFasta

            ### Compare to the previous iteration
            
            logger.info("Compare results with the previous iteration")
            
            #Check if the number of contigs has changed
            
            logger.info("Check if the number of contigs has changed")
            StatsDict["iter_%d" %(i)]["NbContigs"] = len(FilteredSequenceNames)
            
            if StatsDict["iter_%d" %(i)]["NbContigs"] != StatsDict["iter_%d" %(i-1)]["NbContigs"]:
                 logger.info("The number of contigs has changed")          
            elif i >= 2:
                logger.info("Refind the \"brother\" contig from the previous contig for each contig and check they are different")
                # We use Exonerate 
                start_exo_time = time.time()
                ExonerateProcess = Aligner.Exonerate(FileteredTrinityFasta, "%s/Trinity_iter_%d.filtered.fasta" % (TmpDirName,i-1) )
                # We want to keep only the best hit for each contigs
                ExonerateProcess.Bestn = 1
                ExonerateProcess.Model = "est2genome"
                # We customize our output format
                ExonerateProcess.Ryo = "%ti\t%qi\t%ql\t%qal\t%tal\t%tl\t%pi\n"
                ExonerateResult = ExonerateProcess.get_output()
                AlmostIdenticalResults = ApytramNeeds.check_almost_identical_exonerate_results(ExonerateResult)
                if AlmostIdenticalResults and not FinishAllIter:
                    logger.info("Contigs are almost identical than the previous iteration (Same size (~98%), > 99% identity)")
                    Stop =True
                StatsDict["iter_%d" %(i)]["Exonerate2Time"] = time.time() - start_exo_time
                logger.debug("exonerate on previous iter --- %s seconds ---" % (StatsDict["iter_%d" %(i)]["Exonerate2Time"]))
     
            # Check that the coverage has inscreased compared to the previous iteration
        
            logger.info("Check that the coverage has inscreased compared to the previous iteration")
            # We use Mafft
            start_mafft_time = time.time()
            MafftProcess = Aligner.Mafft(QueryFile)
            MafftProcess.QuietOption = True
            MafftProcess.AdjustdirectionOption = True
            MafftProcess.AddOption = FileteredTrinityFasta
            MafftResult = MafftProcess.get_output()
            StatsDict["iter_%d" %(i)]["StrictCoverage"], StatsDict["iter_%d" %(i)]["LargeCoverage"] = ApytramNeeds.calculate_coverage(MafftResult)
            logger.info("Strict Coverage: %s\tLarge Coverage: %s" %(StatsDict["iter_%d" %(i)]["StrictCoverage"], StatsDict["iter_%d" %(i)]["LargeCoverage"]))
            StatsDict["iter_%d" %(i)]["MafftTime"] = time.time() - start_mafft_time
            logger.debug("mafft --- %s seconds ---" % (StatsDict["iter_%d" %(i)]["MafftTime"]))
            
            if not FinishAllIter:
                # Stop iteration if the Largecoverage is not improved and also the Total length
                if StatsDict["iter_%d" %(i)]["TotalLength"] > StatsDict["iter_%d" %(i-1)]["TotalLength"]:
                    pass
                elif StatsDict["iter_%d" %(i)]["TotalScore"] > StatsDict["iter_%d" %(i-1)]["TotalScore"]:
                    pass
                elif StatsDict["iter_%d" %(i)]["LargeCoverage"] <= StatsDict["iter_%d" %(i-1)]["LargeCoverage"]:
                    logger.info("This iteration have a large coverage inferior (or equal) to the previous iteration")
                    Stop = True

                # Stop iteration if the RequiredCoverage is attained
                if StatsDict["iter_%d" %(i)]["LargeCoverage"] >= RequiredCoverage:
                    logger.info("This iteration attains the required bait sequence coverage (%d >= %d)" % (StrictCoverage,RequiredCoverage))
                    Stop = True

            ### Write a fasta file for this iteration if tere is the option --keep_iterations
            
            if KeepIterations:
                # Best sequences of the iteration
                ExitCode = ApytramNeeds.write_apytram_output(FileteredTrinityFasta, TrinityExonerateResultsDict,
                                                 "%s.iter_%d.best.fasta" %(OutPreffixName,i), 
                                                 Header = TrinityExonerateProcess.Ryo.replace('%',"").replace("\n","").split(),
                                                 Names = BestScoreNames.values(),
                                                 Message = "iter_%d.best." %i)
                # All sequences of the iteration
                ExitCode = ApytramNeeds.write_apytram_output(FileteredTrinityFasta,
                                                 TrinityExonerateResultsDict,
                                                 "%s.iter_%d.fasta" %(OutPreffixName,i),
                                                 Header = TrinityExonerateProcess.Ryo.replace('%',"").replace("\n","").split(),
                                                 Message = "iter_%d." %i)

    if IterationNotFinished:
            Reali = i + 1
    else: 
        Reali = i
    
    NoPythonTime = StatsDict["iter_%d" %(Reali)]["BlastTime"] + StatsDict["iter_%d" %(Reali)]["TrinityTime"] +\
                 StatsDict["iter_%d" %(Reali)]["MafftTime"] + StatsDict["iter_%d" %(Reali)]["Exonerate1Time"] +\
                 StatsDict["iter_%d" %(Reali)]["Exonerate2Time"]
    StatsDict["iter_%d" %(Reali)].update({"IterationTime": time.time() - start_iter_i,
                                          "CumulTime": time.time() - start_iter,
                                          "PythonTime": time.time() - start_iter_i - NoPythonTime })
    
    logger.debug("iteration %d --- %s seconds ---" % (Reali, time.time() - start_iter_i))


#### Write output


logger.info("End of Iterations")
if i: #We check that there is at least one iteration with a result
    logger.info("Write outputfiles")
    # Best sequences
    ExitCode = ApytramNeeds.write_apytram_output(FileteredTrinityFasta, TrinityExonerateResultsDict,
                                                 OutPreffixName+".best.fasta", 
                                                 Header = TrinityExonerateProcess.Ryo.replace('%',"").replace("\n","").split(),
                                                 Names = BestScoreNames.values(),
                                                 Message = ".best." %i)
    # Last iteration
    ExitCode = ApytramNeeds.write_apytram_output(FileteredTrinityFasta,
                                                 TrinityExonerateResultsDict,
                                                 OutPreffixName+".fasta",
                                                 Header = TrinityExonerateProcess.Ryo.replace('%',"").replace("\n","").split(),)
    if args.stats:
        logger.info("Write statistics file (OutPreffix.stats.csv)")
        ApytramNeeds.write_stats(StatsDict,OutPreffixName)
        if args.plot:
            logger.info("Create plot from the statistics file (OutPreffix.stats.pdf)")
            ApytramNeeds.create_plot(StatsDict,OutPreffixName)
        
else:
    logger.warn("No results")
    
### Remove tempdir if the option --tmp have not been use
if not args.tmp:
    logger.debug("Remove the temporary directory")
    #Remove the temporary directory :
    if "tmp_apytram" in TmpDirName:
        shutil.rmtree(TmpDirName)
        
logger.debug("--- %s seconds ---" % (time.time() - start_time))