
# coding: utf-8

# In[ ]:

#!/usr/bin/python
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


# In[ ]:

### Option defining
parser = argparse.ArgumentParser(prog = "apytram.py",
                                 description='''
    Run apytram.py on a fastq file to retrieve
    homologous sequences of bait sequences.''')

requiredOptions = parser.add_argument_group('required arguments')
requiredOptions.add_argument('-d', '--database', nargs='?', type=str,
                             help='Database preffix name', required=True)
requiredOptions.add_argument('-t', '--database_type', type=str, choices=["single","paired"],
                             help='single or paired end database', required=True)

parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.add_argument('-log', nargs='?', type=str, default="apytram.log",
                   help = "a log file to report avancement (default: apytram.log)"
                   )
parser.add_argument('--threads',  type=int,
                    help = "Available threads. (Default 1)",
                    default = 1 )



parser.add_argument('-fa', '--fasta',  type=str,
                   help = "Fasta formatted RNA-seq data to build the database of reads")
parser.add_argument('-fq', '--fastq',  type=str,
                   help = "Fastq formatted RNA-seq data to build the database of reads")
parser.add_argument('-out', '--output_preffix',  type=str, default = "./apytram",
                   help = "Output preffix (Default ./apytram)")

parser.add_argument('--tmp',  type=str,
                    help = "Directory to stock all intermediary files for the apytram run. (default: a directory in /tmp which will be removed at the end)",
                    default = "" )
parser.add_argument('--keep_iterations',  action='store_true',
                    help = "A fasta file will be created at each iteration."
                   )

parser.add_argument('-q', '--query',  type=str,
                    help = "Fasta file (nt) with bait sequence for the apytram run." )
parser.add_argument('-i', '--iteration_max',  type=int,
                    help = "Maximum number of iteration. (Default 5)",
                    default = 5 )
parser.add_argument('-e', '--evalue',  type=float,
                    help = "Evalue. (Default 1e-3)",
                    default = 1e-3 )
parser.add_argument('--required_coverage',  type=float,
                    help = "Required coverage of a bait sequence to stop iteration (Default: No threshold)",
                    default = 200 )

parser.add_argument('-id', '--min_id',  type=int,
                    help = "Minimum identity percentage with a query to keep a sequence  (Default 50)",
                    default = 50 )
parser.add_argument('-l', '--min_len',  type=int,
                    help = "Minimum length to keep a sequence  (Default 200)",
                    default = 200 )

StatOptions = parser.add_argument_group('Arguments for statistics and plots')
StatOptions.add_argument('--stats', action='store_true',
                             help='Create files with statistics on each iteration')
StatOptions.add_argument('--plot', action='store_true',
                             help='Create file with plot on statistics on each iteration')



### Option parsing
args = parser.parse_args()

## Example for the dev
#args = parser.parse_args(''' -d example_exec/db/examplefq
#                             -out Out_test/apytram
#                             -fq example/example_db.fastq
#                             -q example/ref_gene.fasta 
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

if args.database_type == "paired":
    PairedData = True
else:
    PairedData = False

if args.plot:
    args.stats = True


# In[ ]:

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


# In[ ]:

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


# In[ ]:

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


# In[ ]:

### Check that there is a database else built it
DatabaseName = args.database
CheckDatabase_BlastdbcmdProcess = BlastPlus.Blastdbcmd(DatabaseName, "", "")

if not CheckDatabase_BlastdbcmdProcess.is_database():
    logger.info(DatabaseName+".nhr does not exist")
    #Build blast formated database from a fasta file
    if args.fastq or args.fasta:
        if args.fastq:
            if not os.path.isfile(args.fastq):
                logger.error("The fastq file (-fq) does not exist.")
                sys.exit(1)
            else:
                # Format the fastq file in fasta
                InputFasta = TmpDirName + "/" + os.path.basename(args.fastq) + ".fasta"
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
if not os.path.isfile(DatabaseName+".nhr"):
    logger.error("Problem in the database building")
    logger.info(DatabaseName+".nhr does not exist")
    sys.exit(1)
else:
    logger.info(DatabaseName+".nhr exists")


# In[ ]:

### If there is a query continue, else stop
if not args.query:
    logger.info("There is no query (-q), apytram have finished.")
    quit()
elif not os.path.isfile(args.query):
    logger.error(args.query+" (-q) is not a file.")
    sys.exit(1)
else:
    queryFile = args.query
    logger.info("apytram will run with \"%s\" as reads database and \"%s\" as bait sequences" %(DatabaseName,queryFile))
    


# In[ ]:

### Make iterations
# Initialisation
i = 0
Stop = False
LastIterCov = [0, 0, 0]
BaitSequences = queryFile
IterationNotFinished = False
if args.stats:  
    StatsDict = {"iter_0":{"iteration_time": 0,
                "cumul_time": time.time() - start_time,
                "large_coverage": 0,
                "strict_coverage": 0,
                "Nb_Contigs": 0}}
    
    

logger.info("Iterations begin")
start_iter = time.time()
while (i < MaxIteration) and (Stop == False):
    i+=1
    start_iter_i = time.time()
    logger.info("Iteration %d/%d" %(i,MaxIteration))
    # Blast bait seqeunce on database of reads
    logger.info("Blast bait sequences on reads database")
    ReadNamesFile = "%s/ReadNames.%d.txt" % (TmpDirName,i)
    BlastnProcess = BlastPlus.Blast("blastn", DatabaseName, BaitSequences)
    BlastnProcess.Evalue = Evalue
    BlastnProcess.Threads = Threads
    BlastnProcess.OutFormat = "6 sacc"
    # Write read names in ReadNamesFile
    ExitCode = BlastnProcess.launch(ReadNamesFile)
    if PairedData:
        # Get paired reads names
        ExitCode = ApytramNeeds.add_paired_read_names(ReadNamesFile)
    # Compare the read list names with the list of the previous iteration:
    Identical = ApytramNeeds.are_identical(ReadNamesFile,"%s/ReadNames.%d.txt" % (TmpDirName,i-1))
    if Identical:
        logger.info("Reads from the current iteration are identical than the previous")
        Stop = True
        IterationNotFinished = True
        i -= 1
    else:
        # Retrieve sequences
        logger.info("Retrieve sequences")
        ReadFasta = TmpDirName + "/Reads.%d.fasta" % (i)
        BlastdbcmdProcess = BlastPlus.Blastdbcmd(DatabaseName, ReadNamesFile, ReadFasta)
        BlastdbcmdProcess.launch()
        # Launch Trinity
        start_trinity_time = time.time()
        logger.info("Launch Trinity")
        TrinityFasta = TmpDirName + "/Trinity_iter_%d" % (i)
        TrinityProcess = Trinity.Trinity(ReadFasta,TrinityFasta)
        # Keep only contig with a length superior to MinLength
        TrinityProcess.MinLength = MinLength
        if PairedData:
            TrinityProcess.Paired = True
        # Use the  --full_cleanup Trinity option to keep only the contig file
        ExitCode = 0
        TrinityProcess.FullCleanup = True
        ExitCode = TrinityProcess.launch()
        TrinityFasta = TrinityFasta + ".Trinity.fasta"
        print("trinity --- %s seconds ---" % (time.time() - start_trinity_time))
        if ExitCode != 0: # Trinity found nothing
            logger.info("Trinity found nothing (ExitCode: %d)" %ExitCode)
            Stop = True
            IterationNotFinished = True
            i -=1
        else:
            ## Filter Trinity contigs to keep only homologous sequences of the reference genes
            logger.info("Compare Trinity results with query sequences")

            # We use Exonerate 
            TrinityExonerate = TmpDirName + "/Trinity_iter_%d.exonerate" % (i)
            start_exo_time = time.time()
            TrinityExonerateProcess = Aligner.Exonerate(queryFile,TrinityFasta)
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
            FilteredSequenceNames, BestScoreNames, TrinityExonerateResultsDict = ApytramNeeds.parse_exonerate_results(TrinityExonerateResult, MinIdentityPercentage)
            print("exonerate on trinity --- %s seconds ---" % (time.time() - start_exo_time))

            # Filter hit
            logger.info("Filter sequence with a identity percentage superior to %d" %(MinIdentityPercentage)) 
            FileteredTrinityFasta =  "%s/Trinity_iter_%d.filtered.fasta" % (TmpDirName,i)
            ExitCode = ApytramNeeds.filter_fasta(TrinityFasta, FilteredSequenceNames, FileteredTrinityFasta)
            
            # Validated sequences become bait sequences
            BaitSequences = FileteredTrinityFasta

            ## Compare to the previous iteration
            logger.info("Compare results with the previous iteration")
            
            logger.info("Check if the number of contigs has changed")
            NbContigs = len(FilteredSequenceNames)
            if NbContigs != LastIterCov[2]:
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
                if AlmostIdenticalResults:
                    logger.info("Contigs are almost identical than the previous iteration (Same size (~98%), > 99% identity)")
                    Stop =True
                print("exonerate on previous iter --- %s seconds ---" % (time.time() - start_exo_time))
     
            # Check that the coverage has inscreased compared to the previous iteration
            logger.info("Check that the coverage has inscreased compared to the previous iteration")
            start_mafft_time = time.time()
            MafftProcess = Aligner.Mafft(queryFile)
            MafftProcess.QuietOption = True
            MafftProcess.AdjustdirectionOption = True
            MafftProcess.AddOption = FileteredTrinityFasta
            MafftResult = MafftProcess.get_output()
            StrictCoverage, LargeCoverage = ApytramNeeds.calculate_coverage(MafftResult)
            logger.info("Strict Coverage: %d\tLarge Coverage: %d" %(StrictCoverage, LargeCoverage))
            print("mafft --- %s seconds ---" % (time.time() - start_mafft_time))

            if LargeCoverage <= LastIterCov[1]:
                logger.info("This iteration have a large coverage inferior (or equal) to the previous iteration")
                Stop = True
            elif StrictCoverage >= RequiredCoverage:
                logger.info("This iteration attains the required bait sequence coverage (%d >= %d)" % (StrictCoverage,RequiredCoverage))
                Stop = True

            # Write a fasta file for this iteration if tere is the option --keep_iterations
            if KeepIterations:
                # Best sequences
                ExitCode = ApytramNeeds.write_apytram_output(FileteredTrinityFasta, TrinityExonerateResultsDict,
                                                 "%s.iter_%d.best.fasta" %(OutPreffixName,i), 
                                                 Header = TrinityExonerateProcess.Ryo.replace('%',"").replace("\n","").split(),
                                                 Names = BestScoreNames,
                                                 Message = "iter_%d." %i)
                # Iteration
                ExitCode = ApytramNeeds.write_apytram_output(FileteredTrinityFasta,
                                                 TrinityExonerateResultsDict,
                                                 "%s.iter_%d.fasta" %(OutPreffixName,i),
                                                 Header = TrinityExonerateProcess.Ryo.replace('%',"").replace("\n","").split(),
                                                 Message = "iter_%d." %i)

            LastIterCov = [StrictCoverage, LargeCoverage, NbContigs]
    if IterationNotFinished:
            reali = i + 1
    else: 
        reali = i
    
    if args.stats:
        StatsDict["iter_%d" %(reali)] = {"iteration_time": time.time() - start_iter_i,
                                "cumul_time": time.time() - start_iter,
                                "large_coverage": LargeCoverage,
                                "strict_coverage": StrictCoverage,
                                "Nb_Contigs": NbContigs
                
                                            }
    print("iteration %d --- %s seconds ---" % (reali,time.time() - start_iter_i))


# In[ ]:

#### Write output
logger.info("End of Iterations")
if i: #We check that there is at least one iteration with a result
    logger.info("Write outputfiles")
    # Best sequences
    ExitCode = ApytramNeeds.write_apytram_output(FileteredTrinityFasta, TrinityExonerateResultsDict,
                                                 OutPreffixName+".best.fasta", 
                                                 Header = TrinityExonerateProcess.Ryo.replace('%',"").replace("\n","").split(),
                                                 Names = BestScoreNames)
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
    logger.info("No results")

    
### Remove tempdir if the option --tmp have not been use
if not args.tmp:
    logger.info("Remove the temporary directory")
    #Remove the temporary directory :
    if "tmp_apytram" in TmpDirName:
        shutil.rmtree(TmpDirName)


        

        
    


# In[ ]:

print("--- %s seconds ---" % (time.time() - start_time))


# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:



