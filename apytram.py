#!/usr/bin/python
# coding: utf-8

# File: apytram.py
# Created by: Carine Rey
# Created on: Nov 2015
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
import logging
import argparse

import ApytramLib


def strict_positive_integer(x):
    x = int(x)
    if type(x) != type(1) or x <= 0:
        raise argparse.ArgumentTypeError("Must be an integer superior to 0")
    return x

def positive_integer(x):
    x = int(x)
    if type(x) != type(1) or x < 0:
        raise argparse.ArgumentTypeError("Must be a positive integer")
    return x

def strict_positive_float(x):
    x = float(x)
    if type(x) != type(0.1) or x <= 0:
        raise argparse.ArgumentTypeError("Must be an float superior to 0")
    return x

def positive_float(x):
    x = float(x)
    if type(x) != type(0.1) or x < 0:
        raise argparse.ArgumentTypeError("Must be an float superior to 0")
    return x




start_time = time.time()

### Option defining
parser = argparse.ArgumentParser(prog="apytram.py",
                                 description='''
    Run apytram.py on a fastq file to retrieve
    homologous sequences of bait sequences.''')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

##############
requiredOptions = parser.add_argument_group('Required arguments')
requiredOptions.add_argument('-d', '--database', type=str,
                             help='Database prefix name. If a database with the same name already exists, the existing database will be kept and the database will NOT be rebuilt.', required=True)
requiredOptions.add_argument('-dt', '--database_type', type=str,
                             help="""
                                  single: single unstranded data ______________________
                                  paired: paired unstranded data ______________________
                                  RF: paired stranded data (/1 = reverse ; /2 = forward)
                                  FR: paired stranded data (/1 = forward ; /2 = reverse)
                                  F: single stranded data (reads = forward) ____________
                                  R: single stranded data (reads = reverse) ____________
                                  WARNING: Paired read names must finished by 1 or 2""",
                            required=True)
requiredOptions.add_argument('-out', '--output_prefix', type=str,
                   help="Output prefix", required=True)
##############


##############
InUSOptions = parser.add_argument_group('Input Files')
InUSOptions.add_argument('-fa', '--fasta', type=str,
                   help="Fasta formated RNA-seq data to build the database of reads (only one file).")
InUSOptions.add_argument('-fq', '--fastq', type=str,
                   help="Fastq formated RNA-seq data to build the database of reads (several space delimited fastq file names are allowed). For paired data, fq must be previously concatenated. WARNING: Paired read names must finished by 1 or 2. (fastq files will be first converted to a fasta file. This process can require some time.)")
##############

##############
QueryOptions = parser.add_argument_group('Query File')
QueryOptions.add_argument('-q', '--query', type=str,
                    help="""
                            Fasta file (nucl) with homologous bait sequences which will be treated together for the apytram run.
                            If no query is submitted, the program will just build the database.
                            WARNING: Sequences must not contain "- * . "
                           """,
                           )
#QueryOptions.add_argument('-pep', '--query_pep',  type=str,
#                   default="",
#                   help="Fasta file containing the query in the peptide format. It will be used at the first iteration as bait sequences to fish reads. It is compulsory to include also the query in nucleotide format (-q option)")
##############


##############
IterationOptions = parser.add_argument_group('Number of iterations')
IterationOptions.add_argument('-i', '--iteration_max', type=strict_positive_integer,
                    help="Maximum number of iterations. (Default 5)",
                    default=5)
IterationOptions.add_argument('-i_start', '--iteration_start', type=positive_integer,
                    help="Number of the first iteration. If different of 1, the tmp option must be used. (Default: 1)",
                    default=1)
##############


##############
OutOptions = parser.add_argument_group('Output Files')
OutOptions.add_argument('-log', type=str, default="apytram.log",
                   help="a log file to report avancement (default: apytram.log)")
OutOptions.add_argument('-tmp', type=str,
                    help="Directory to stock all intermediary files for the apytram run. (default: a directory in /tmp which will be removed at the end)",
                    default="")
OutOptions.add_argument('--keep_tmp', action='store_true',
                        default=False,
                        help="By default, the temporary directory will be remove.")

#OutOptions.add_argument('--keep_iterations',  action='store_true',
#                    help="A fasta file containing reconstructed sequences will be created at each iteration. (default: False)")

OutOptions.add_argument('--no_best_file', action='store_true',
                        default=False,
                        help="By default, a fasta file (Outprefix.best.fasta) containing only the best sequence is created. If this option is used, it will NOT be created.")

OutOptions.add_argument('--only_best_file', action='store_true',
                        default=False,
                        help="By default, a fasta file (Outprefix.fasta) containing all sequences from the last iteration is created. If this option is used, it will NOT be created.")

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
SearchOptions.add_argument('-e', '--evalue', type=positive_float,
                    help="Evalue threshold of the blastn of the bait queries on the database of reads. (Default 1e-5)",
                    default=1e-5)

SearchOptions.add_argument('-id', '--min_id', type=positive_integer,
                    help="Minimum identity percentage of a sequence with a query on the length of their alignment so that the sequence is kept at the end of a iteration (Default 70)",
                    default=70)
SearchOptions.add_argument('-mal', '--min_ali_len', type=positive_integer,
                    help="Minimum alignment length of a sequence on a query to be kept at the end of a iteration (Default 180)",
                    default=180)
SearchOptions.add_argument('-len', '--min_len', type=positive_integer,
                    help="Minimum length to keep a sequence at the end of a iteration. (Default 200)",
                    default=200)
##############


##############
StopOptions = parser.add_argument_group('Criteria to stop iteration')
StopOptions.add_argument('-required_coverage', type=positive_integer,
                    help="Required coverage of a bait sequence to stop iteration (Default: No threshold)",
                    default=200)
StopOptions.add_argument('--finish_all_iter', action='store_true',
                    help="By default, iterations are stop if there is no improvment, if this option is used apytram will finish all iteration (-i).",
                    default=False)

##############


##############
FinalFilterOptions = parser.add_argument_group('Thresholds for Final output files')
FinalFilterOptions.add_argument('-flen', '--final_min_len', type=positive_integer,
                    help="Minimum PERCENTAGE of the query length to keep a sequence at the end of the run. (Default: 0)",
                    default=0)
FinalFilterOptions.add_argument('-fid', '--final_min_id', type=positive_integer,
                    help="Minimum identity PERCENTAGE of a sequence with a query on the length of their alignment so that the sequence is kept at the end of the run (Default 0)",
                    default=0)
FinalFilterOptions.add_argument('-fmal', '--final_min_ali_len', type=positive_integer,
                     help="Alignment length between a sequence and a query must be at least this PERCENTAGE of the query length to keep this sequence at the end of the run. (Default: 0)",
                    default=0)
##############


##############
MiscellaneousOptions = parser.add_argument_group('Miscellaneous options')
MiscellaneousOptions.add_argument('-threads', type=positive_integer,
                    help="Number of available threads. (Default 1)",
                    default=1)
MiscellaneousOptions.add_argument('-memory', type=positive_integer,
                    help="Memory available for the assembly in Giga. (Default 1)",
                    default=1)
MiscellaneousOptions.add_argument('-time_max', type=positive_integer,
                    help="Do not begin a new iteration if the job duration (in seconds) has exceed this threshold. (Default 7200)",
                    default=7200)
MiscellaneousOptions.add_argument('--UseMapper', action='store_true',
                    help="Use NextGenMapper instead of balstn to fish reads",
                    default=False)
MiscellaneousOptions.add_argument('--write_even_empty', action='store_true',
                        default=False,
                        help="Write output fasta files, even if they must be empty. (Default: False)")
MiscellaneousOptions.add_argument('--out_by_species', action='store_true',
                        default=False,
                        help="Write output fasta files for each species. (Default: False)")
MiscellaneousOptions.add_argument('--debug', action='store_true', default=False,
                   help="debug mode, default False")
##############


### Option parsing
args = parser.parse_args()

### Set up the log directory
if args.log:
    LogDirName = os.path.dirname(args.log)
    if not os.path.isdir(LogDirName) and LogDirName:
        os.makedirs(LogDirName)

### Set up the logger
LogFile = args.log
# create logger with 'spam_application'
logger = logging.getLogger('apytram')
logger.setLevel(logging.INFO)
# create file handler which logs even debug messages
fh = logging.FileHandler(LogFile)
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
if args.debug:
    ch.setLevel(logging.DEBUG)
    fh.setLevel(logging.DEBUG)
    logger.setLevel(logging.DEBUG)
else:
    ch.setLevel(logging.WARN)

#ch.setLevel(logging.DEBUG)
# create formatter and add it to the handlers
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
formatter = logging.Formatter('%(message)s')
ch.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(fh)
logger.addHandler(ch)

logger.warning("[Running in process...]\n[Warning messages may print, but they are no error. If real errors appear, the process will stop]")
logger.info(" ".join(sys.argv))


error = 0
empty_queries = 0


### Set up the temporary directory
if args.tmp:
    if not "apytram" in args.tmp:
        logger.error("""ERROR: Temporary directory name (-tmp) must contain "apytram" to safety reasons because it will be completly remove.""")
        error += 1
    elif os.path.isdir(args.tmp):
        logger.info("The temporary directory %s exists", args.tmp)
    else:
        logger.info("The temporary directory %s does not exist, it will be created", args.tmp)
        os.makedirs(args.tmp)
    TmpDirName = args.tmp
else:
    TmpDirName = tempfile.mkdtemp(prefix='tmp_apytram_')


### Read the arguments


# Define global parameters
StartIteration = args.iteration_start
MaxIteration = args.iteration_max

UseMapper = args.UseMapper

if MaxIteration < 1:
    logger.error("The number of iteration (-i) must be superior to 0")
    error += 1

Evalue = args.evalue

MinIdentityPercentage = args.min_id
MinAliLength = args.min_ali_len
MinLength = args.min_len
RequiredCoverage = args.required_coverage

FinalMinLength = args.final_min_len
FinalMinIdentityPercentage = args.final_min_id
FinalMinAliLength = args.final_min_ali_len

#KeepIterations = args.keep_iterations
FinishAllIter = args.finish_all_iter

Threads = args.threads
Memory = args.memory
MaxTime = args.time_max


if args.plot:
    args.stats = True

# Define query files
if args.query:
    Queries = args.query.split(",")
else:
    Queries = []

logger.warning("Query:")

QueriesNamesList = []
QueriesList = []

for query in Queries:
    if re.search(":", query):
        (query, name) = query.split(":")
    else:
        name = ""
    if not os.path.isfile(query):
        logger.error("\t-%s ... ERROR (Don't exist)", query)
        error += 1
    elif not os.stat(query).st_size:
        logger.error("\t-%s ... ERROR (empty)", query)
        error += 1
        empty_queries += 1
        if not name:
            name = os.path.basename(os.path.splitext(query)[0])
        new_query = ApytramLib.ApytramClasses.Query(name, query, logger)
        new_query.TmpDirName = "%s/%s" %(TmptrinityDirName, name)
        QueriesList.append(new_query)
        QueriesNamesList.append(name)
    else:
        if not name:
            name = os.path.basename(os.path.splitext(query)[0])
        new_query = ApytramLib.ApytramClasses.Query(name, query, logger)
        new_query.TmpDirName = "%s/%s/" %(TmpDirName, name)
        ApytramLib.ApytramNeeds.set_directory_from_prefix(new_query.TmpDirName, "temporary", logger)
        logger.warning("\t-%s %s ... ok (%s sequences)", name, new_query.RawQuery, new_query.SequenceNb)
        if not name in QueriesNamesList:
            new_query.initialization()
            QueriesList.append(new_query)
            QueriesNamesList.append(name)
        else:
            logger.error("""The name "%s" must have only one associated query file.You must chose between:\n\t%s\n\t%s""",
                        new_query.Name,
                        new_query.RawQuery,
                        QueriesList[QueriesNamesList.index(new_query.Name)].RawQuery
                        )
            error += 1

if not len(QueriesList):
    logger.warning("No query")


# Define species
SpeciesList = []
SpeciesNamesList = []

DBs = args.database.split(",")
UniqSpecies_flag = True

if len(DBs) > 1:
    # There are several species:
    UniqSpecies_flag = False
elif len(DBs) == 1:
    if not re.search(":", DBs[0]):
        DBs[0] += ":SP"

DB_types_dict = {}
for item in args.database_type.split(","):
    db_type_sp = item.split(":")
    logger.debug(db_type_sp)
    if len(db_type_sp) == 2:
        db_type = db_type_sp[1]
        DB_types_dict[db_type] = db_type_sp[0]
    elif len(db_type_sp) == 1 and len(DBs) == 1:
        sp = DBs[0].split(":")[1]
        DB_types_dict[sp] = db_type_sp[0]


FAs = []
FQs = []
if args.fasta:
    FAs += args.fasta.split(",")
if args.fastq:
    FQs += args.fastq.split(",")



# Get species names and check if a database exist
for item in DBs:
    new_species = ApytramLib.ApytramClasses.RNA_species(start_time, logger)
    new_species.Evalue = Evalue
    new_species.MinLength = MinLength
    new_species.MinIdentityPercentage = MinIdentityPercentage
    new_species.MinAliLength = MinAliLength
    new_species.FinalMinLength = FinalMinLength
    new_species.FinalMinIdentityPercentage = FinalMinIdentityPercentage
    new_species.FinalMinAliLength = FinalMinAliLength
    new_species.keep_tmp = args.keep_tmp

    s = item.split(":")

    if len(s) == 2:
        (new_species.DatabaseName, new_species.Species) = s
        new_species.DatabaseType = DB_types_dict.get(new_species.Species ,"")
        # Check db_type
        if new_species.DatabaseType in ["single", "paired", "FR", "RF", "F", "R"]:
            if new_species.DatabaseType in ["RF", "FR", "R", "F"]:
                 new_species.StrandedData = True
            else:
                 new_species.StrandedData = False

            if new_species.DatabaseType in ["paired", "RF", "FR"]:
                new_species.PairedData = True
            else:
                new_species.PairedData = False
        else:
             logger.error("""The species "%s" must have data type in ["single", "paired", "FR", "RF", "F", "R"] and not %s""",
                        new_species.Species,
                        new_species.DatabaseType
                        )

        if not new_species.Species in SpeciesList:
            new_species.FormatedDatabase = new_species.has_a_formated_database()
            SpeciesList.append(new_species)
            SpeciesNamesList.append(new_species.Species)
        else:
            logger.error("""The species "%s" must have only one associated database.You must chose between:\n\t%s\n\t%s""",
                        new_species.Species,
                        new_species.DatabaseName,
                        SpeciesList[SpeciesNamesList.index(new_species.Species)].DatabaseName
                        )
    else:
        logger.error(""" "%s" must be formatted as "db:species" if you want to use the multispecies option" """, item)

# associate fasta files with a species
for item in FAs:
    f = item.split(":")
    if len(f) == 1 and UniqSpecies_flag:
        f.append(":SP")
    if len(f) == 2:
        (fa, species) = f
        if os.path.isfile(fa):
            SpeciesList[SpeciesNamesList.index(species)].Fasta.append(fa)
        else:
            logger.error("%s (-fa) is not a file.", fa)
            ApytramLib.ApytramNeeds.end(1, TmpDirName, keep_tmp=args.keep_tmp)
    else:
        logger.error(""" "%s" must be formatted as "fa_path:species" if you want to use the multispecies option" """, item)

# associate fastq files with a species
for item in FQs:
    f = item.split(":")
    if len(f) == 1 and UniqSpecies_flag:
        f.append(":SP")
    if len(f) == 2:
        (fq, species) = f
        if os.path.isfile(fq):
            if species in SpeciesNamesList:
                SpeciesList[SpeciesNamesList.index(species)].Fastq.append(fq)
            else:
                logger.error("The species associted with %s is %s. But there is no database associated with the species %s.", fq, species, species)
                ApytramLib.ApytramNeeds.end(1, TmpDirName, keep_tmp=args.keep_tmp)
        else:
            logger.error("%s (-fq) is not a file.", fq)
            ApytramLib.ApytramNeeds.end(1, TmpDirName, keep_tmp=args.keep_tmp)
    else:
        logger.error(""" "%s" must be formatted as "fq_path:species" if you want to use the multispecies option" """, item)


#Check all species has a formated database or input fasta/fastq files:

logger.warning("Species:")
for Species in SpeciesList:
    if not Species.DatabaseType in ["single", "paired", "FR", "RF", "F", "R"]:
        logger.warning("\t-%s ... Unknown data type (%s) -> ERROR", Species.Species, Species.DatabaseType)
        error += 1
    elif Species.FormatedDatabase:
        logger.warning("\t-%s ... Formated database (%s)", Species.Species, Species.DatabaseType)
    elif Species.Fasta and Species.Fastq:
        logger.warning("\t-%s ... NO formated database. fasta AND fastq input files -> ERROR", Species.Species)
        error += 1
    elif Species.Fasta or Species.Fastq:
        logger.warning("\t-%s ... NO formated database. (A database will be built)", Species.Species)
    else:
        logger.warning("\t-%s ... NO formated database and NO input fasta/fastq files -> ERROR", Species.Species)
        error += 1

logger.debug("Time to parse input command line: %s", time.time() - start_time)

if error > 0:
    if args.write_even_empty and error == empty_queries:
        logger.warning("Some query files are empty, but you use --write_even_empty option -> empty files will be create")
    else:
        logger.error("Error(s) occured, see above")
        ApytramLib.ApytramNeeds.end(1, TmpDirName, keep_tmp=args.keep_tmp)

### If iteration begin not from 1, the temporary directory must be given by the user
if StartIteration != 1 and not args.tmp:
    logger.error("If you want to restart a previous job, the previous temporary directory must be given.")
    ApytramLib.ApytramNeeds.end(1, TmpDirName, keep_tmp=args.keep_tmp)


### Get the available free space of the tmp dir
FreeSpaceTmpDir = ApytramLib.ApytramNeeds.get_free_space(TmpDirName)
logger.debug("%s free space in %s", FreeSpaceTmpDir, TmpDirName)

if UseMapper:
    logger.warn("Use NextGenMapper instead of Blastn to fish reads.\nRequire raw reads:")

### Check that there is a database for each species, otherwise build it
for Species in SpeciesList:
    if not Species.FormatedDatabase:
        Species.set_TmpDir(TmpDirName + "/db/" + Species.Species)
        Species.prepare_database(FreeSpaceTmpDir, TmpDirName)
        Species.build_database(FreeSpaceTmpDir, TmpDirName)
    ### If Use mapper, apytram needs raw reads
    if UseMapper:
        if Species.InputFastaFilename:
            logger.warn("\t-%s ... Raw reads available (%s)",Species.Species, Species.InputFastaFilename)
            pass
        elif Species.Fasta or Species.Fastq:
            logger.warn("\t-%s ... Raw reads NOT available. Get it from Input Fasta or Fastq",Species.Species)
            Species.set_TmpDir(TmpDirName + "/db/" + Species.Species)
            Species.prepare_database(FreeSpaceTmpDir, TmpDirName)
        elif Species.FormatedDatabase:
            logger.warn("\t-%s ... Raw reads NOT available. Get it from the database",Species.Species)
            Species.set_TmpDir(TmpDirName + "/db/" + Species.Species)
            Species.get_all_reads()

        if not Species.InputFastaFilename:
            logger.error("No raw reads available for %s.", Species.Species)
            ApytramLib.ApytramNeeds.end(1, TmpDirName, keep_tmp=args.keep_tmp)

### If there is a query continue, else stop
if not args.query:
    logger.info("There is no query (-q), apytram has finished.")
    ApytramLib.ApytramNeeds.end(0, TmpDirName, keep_tmp=args.keep_tmp)
else:
    ApytramLib.ApytramNeeds.set_directory_from_prefix(args.output_prefix, "output", logger)

### Set up the output directory
if args.output_prefix:
    OutDirName = os.path.dirname(args.output_prefix)
    OutPrefixName = args.output_prefix
    if os.path.isdir(OutDirName):
        logger.info("The output directory %s exists", os.path.dirname(args.output_prefix))
    elif OutDirName: # if OutDirName is not a empty string we create the directory
        logger.info("The output directory %s does not exist, it will be created", os.path.dirname(args.output_prefix))
        os.makedirs(os.path.dirname(args.output_prefix))
else:
    logger.error("The output prefix must be defined")
    ApytramLib.ApytramNeeds.end(1, TmpDirName, keep_tmp=args.keep_tmp)


for Query in QueriesList:
    logger.warning("NEW QUERY: %s", Query.Name)

    Query.OutPrefixName = "%s.%s" %(OutPrefixName, Query.Name)
    Query.NbSpecies = len(SpeciesNamesList)
    Query.StartTime = time.time()

    for Species in SpeciesList:
        Species.new_query(Query)

    #Iterative process
    while (Query.AbsIteration < MaxIteration) and (Query.continue_iter()) and Query.SequenceNb:
        Query.AbsIteration += 1
        Query.SpeciesWithoutImprovment[Query.AbsIteration] = []
        logger.warning("\tIteration %d/%d", Query.AbsIteration, MaxIteration)

        for Species in SpeciesList:
            if not Species.Finished:
                # Build new baitsequences file
                Query.new_species_iteration(SpeciesList)

            if Species.Improvment:
                Species.new_iteration()
                logger.warning("\t\t Start iteration %d/%d for %s", Species.CurrentIteration, MaxIteration, Species.Species)

                ### Blast bait sequences on database of reads
                # Write read names in ReadNamesFile if the file does not exist
                if not os.path.isfile(Species.ReadNamesFilename):
                    Species.fish_reads(Query.BaitSequences, Threads, mapper=UseMapper)
                else:
                    logger.warn("%s has already been created, it will be used", Species.ReadNamesFilename)

                if Species.PairedData:
                    # Get paired reads names and remove duplicated names
                    logger.info("Get paired reads names and remove duplicated names")
                    ApytramLib.ApytramNeeds.add_paired_read_names(Species.ReadNamesFilename, Species.ParsedReadNamesFilename, logger)
                else:
                    # Remove duplicated names
                    logger.info("Remove duplicated names")
                    ApytramLib.ApytramNeeds.remove_duplicated_read_names(Species.ReadNamesFilename, Species.ParsedReadNamesFilename, logger)

                # Count the number of reads which will be used in the Trinity assembly
                logger.info("Count the number of reads")
                Species.ReadsNumber = ApytramLib.ApytramNeeds.count_lines(Species.ParsedReadNamesFilename)
                Species.add_iter_statistic("ReadsNumber", Species.ReadsNumber)

                if not Species.ReadsNumber:
                    logger.warning("No read recruted by Blast at the iteration %s", Species.CurrentIteration)
                    Species.Improvment = False
                    Species.CompletedIteration = False

            if Species.Improvment:
                # Compare the read list names with the list of the previous iteration:
                NbNewReads = ApytramLib.ApytramNeeds.number_new_reads(Species.PreviousReadNamesFilename, Species.ParsedReadNamesFilename, nb_intial=Species.ReadsNumber)
                logger.warning("Iteration: %s - Species: %s - Number of new reads: %s", Species.CurrentIteration, Species.Species, NbNewReads)

                if (NbNewReads == 0) and not FinishAllIter:
                    logger.info("Reads from the current iteration are identical to reads from the previous iteration")
                    Species.Improvment = False
                    Species.CompletedIteration = False

            if Species.Improvment:
                ### Retrieve reads sequences
                if Species.InputFastaFilename:
                    Species.get_read_sequences(Threads, Memory, meth="seqtk")
                else:
                    Species.get_read_sequences(Threads, Memory, meth="blastdbcmd")

                ### Launch Trinity
                Species.launch_Trinity(Threads, Memory)

                if not os.path.isfile(Species.TrinityFastaFilename): # Trinity found nothing
                    Species.Improvment = False
                    Species.CompletedIteration = False

                if Species.Improvment:
                    ### Filter Trinity contigs to keep only homologous sequences of the reference genes
                    logger.info("Compare Trinity results with query sequences")
                    Species.get_homology_between_trinity_results_and_references(Query)

                    if not Species.TrinityExonerateResult:
                        logger.info("Reconstructed sequences but no homologous with references (even with the more sensible model)")
                        Species.Improvment = False
                        Species.CompletedIteration = False

                    if Species.Improvment:
                        # Keep only sequence with a identity percentage > MinIdentitypercentage on the whole hit
                        # and write filtered sequences in a file
                        Species.filter_trinity_results_according_homology_results()

                        ### Validated sequences (Species.FilteredTrinityFasta)  become bait sequences

                        if not Species.FilteredTrinityFasta.Sequences:
                            logger.warning("No sequence has passed the iteration filter at the iteration %s for %s", Species.CurrentIteration, Species.Species)
                            Species.Improvment = False
                            Species.CompletedIteration = False

                        else:
                            ### Compare sequences of the current iteration to those of the previous iteration
                            logger.info("Compare results with the previous iteration")

                            #Check if the number of contigs has changed
                            logger.info("Check if the number of contigs has changed")

                            if Species.get_iter_statistic("NbContigs") != Species.get_iter_statistic("NbContigs", RelIter=-1):
                                logger.info("The number of contigs has changed")
                            elif Query.AbsIteration >= 2:
                                # Use Exonerate to compare the current iteration with the previous
                                Species.compare_current_and_previous_iterations()

                            # Check that the coverage has increased compared to the previous iteration

                            logger.info("Check that the coverage has inscreased compared to the previous iteration")
                            Species.measure_coverage(Query)

                            if Species.CompletedIteration:
                                # Stop iteration if both Largecoverage and Total length are not improved
                                if Species.get_iter_statistic("AverageLength") != Species.get_iter_statistic("AverageLength", RelIter=-1):
                                    pass
                                elif Species.get_iter_statistic("AverageScore") != Species.get_iter_statistic("AverageScore", RelIter=-1):
                                    pass
                                elif Species.get_iter_statistic("TotalLength") != Species.get_iter_statistic("TotalLength", RelIter=-1):
                                    pass
                                elif Species.get_iter_statistic("TotalScore") != Species.get_iter_statistic("TotalScore", RelIter=-1):
                                    pass
                                elif Species.get_iter_statistic("BestScore") != Species.get_iter_statistic("BestScore", RelIter=-1):
                                    pass
                                elif Species.get_iter_statistic("LargeCoverage") != Species.get_iter_statistic("LargeCoverage", RelIter=-1):
                                    logger.info("This iteration have a large coverage inferior (or equal) to the previous iteration")
                                    Species.Improvment = False

                                # Stop iteration if the RequiredCoverage is reached
                                if Species.get_iter_statistic("StrictCoverage") >= RequiredCoverage:
                                    logger.info("This iteration attains the required bait sequence coverage (%d >= %d)", Species.get_iter_statistic("StrictCoverage"), RequiredCoverage)
                                    Species.Improvment = False


                            ### Write a fasta file for this iteration if the option --keep_iterations was selected
                            #if KeepIterations:
                            #    if not args.no_best_file:
                            #    # Best sequences of the iteration
                            #        ExitCode = ApytramLib.ApytramNeeds.write_apytram_output(FilteredTrinityFasta, TrinityExonerateResultsDict,
                            #                                        "%s.iter_%d.best.fasta" %(OutPrefixName,i),
                            #                                        Header = TrinityExonerateProcess.Ryo.replace('%',"").replace("\n","").split(),
                            #                                        Names = BestScoreNames.values(),
                            #                                        Message = "iter_%d.best." %i)
                            #    # All sequences of the iteration
                            #    ExitCode = ApytramLib.ApytramNeeds.write_apytram_output(FilteredTrinityFasta,
                            #                                    TrinityExonerateResultsDict,
                            #                                    "%s.iter_%d.fasta" %(OutPrefixName,i),
                            #                                    Header = TrinityExonerateProcess.Ryo.replace('%',"").replace("\n","").split(),
                            #                                    Message = "iter_%d." %i)
                            #    # Mafft alignment
                            #    ApytramLib.ApytramNeeds.write_in_file(MafftResult,"%s.iter_%s.ali.fasta" %(OutPrefixName,i))


            # End iteration
            Species.FinalIteration = Species.CurrentIteration
            if not Species.CompletedIteration:
                logger.debug("Iteration stop before end")
                Species.FinalIteration -= 1
            if not Species.Improvment:
                Query.SpeciesWithoutImprovment[Query.AbsIteration].append(Species.Species)

            Species.end_iteration() # just stop timer


            logger.info("End of the iteration %s for %s : --- %s seconds ---", Species.CurrentIteration, Species.Species, Species.get_iter_statistic("IterationTime"))


            if (time.time() - Query.StartTime) > MaxTime:
                logger.warn("No new iteration for this query and this species will begin because the maximum duration (%s seconds) of the job is attained. (%s seconds)", MaxTime, str(time.time() - Query.StartTime))
                Query.Stop = True


    logger.info("End of Iterations for %s. Iterative process takes %s seconds.", Query.Name, str(time.time() - Query.StartTime))

    ### Final filter
    start_output = time.time()
    for Species in SpeciesList:
        if Species.FinalIteration: #We check that there is at least one iteration with a result
            if Species.FinalMinLength or Species.FinalMinIdentityPercentage or Species.FinalMinAliLength: # A final filter is required
                start_iter_i = time.time()
                Species.new_iteration()
                logger.info("Start final filter for %s", Species.Species)
                # Keep only sequence with a identity percentage > FinalMinIdentitypercentage on the whole hit
                # and write filtered sequences in a file
                Species.filter_trinity_results_according_homology_results(final_iteration=True)
        else:
            logger.warn("No results for %s", Species.Species)

        if Species.FilteredTrinityFasta.Sequences or args.write_even_empty: # If sequences pass the last filter
            Species.rename_sequences()
            logger.debug("output: %s", Species.Species)
            # Prepare fasta output files by species
            if not args.no_best_file:
                QuerySpeciesBestFastaOutput = Species.get_output_fasta(fasta="best")
                Query.BestOutFileContent.extend(QuerySpeciesBestFastaOutput)
                if args.out_by_species:
                    QuerySpeciesBestFinalFastaFileName = Query.OutPrefixName + "." + Species.Species + ".best.fasta"
                    if QuerySpeciesBestFastaOutput or args.write_even_empty:
                         ApytramLib.ApytramNeeds.write_in_file("".join(QuerySpeciesBestFastaOutput), QuerySpeciesBestFinalFastaFileName)

            if not args.only_best_file:
                QuerySpeciesFastaOutput = Species.get_output_fasta(fasta="all")
                Query.OutFileContent.extend(QuerySpeciesFastaOutput)
                if args.out_by_species:
                    QuerySpeciesFinalFastaFileName = Query.OutPrefixName + "." + Species.Species + ".fasta"
                    if QuerySpeciesFastaOutput or args.write_even_empty:
                        ApytramLib.ApytramNeeds.write_in_file("".join(QuerySpeciesFastaOutput), QuerySpeciesFinalFastaFileName)

            if args.plot_ali or args.stats:
                 ### Calculate the coverage
                 logger.info("Calculate the final coverage")
                 Species.measure_coverage(Query)

            Species.end_iteration()

            # Stats files

        if args.stats:
            logger.info("Write statistics file (OutPrefix.stats.csv)")
            Query.StatsFileContent.extend(Species.get_stats())
            Query.TimeStatsDictList.append(Species.ExecutionStats.TimeStatsDict)
            Query.IterStatsDictList.append(Species.ExecutionStats.IterStatsDict)

    ### Write fasta outputfiles
    if Query.BestOutFileContent or args.write_even_empty:
        if not args.no_best_file:
            Query.FinalFastaFileName = Query.OutPrefixName + ".best.fasta"
            ApytramLib.ApytramNeeds.write_in_file("".join(Query.BestOutFileContent), Query.FinalFastaFileName)
    if Query.OutFileContent or args.write_even_empty:
        if not args.only_best_file:
            Query.FinalFastaFileName = Query.OutPrefixName + ".fasta"
            ApytramLib.ApytramNeeds.write_in_file("".join(Query.OutFileContent), Query.FinalFastaFileName)
    if Query.StatsFileContent:
        ApytramLib.ApytramNeeds.write_stats(Query.StatsFileContent, Query.OutPrefixName + ".stats.csv")
        if args.plot:
            logger.info("Create plot from the statistics file (OutPrefix.stats.pdf)")
            ApytramLib.ApytramNeeds.create_plot(Query.TimeStatsDictList,
                                     Query.IterStatsDictList,
                                     SpeciesNamesList,
                                     Query.OutPrefixName)

    if args.plot_ali and Query.FinalFastaFileName and (Query.OutFileContent or Query.BestOutFileContent):
        start_output_ali = time.time()
        Query.measure_final_coverage()
        LengthAlignment = len(Species.DicPlotCov[Species.DicPlotCov.keys()[0]])
        if LengthAlignment <= 3100:
            logger.info("Create plot of the final alignment (OutPrefix.ali.png)")
            ApytramLib.ApytramNeeds.create_plot_ali(Query.DicPlotCov, Query.OutPrefixName)
        else:
            logger.warn("Final alignment is longger than 3100 pb, the plot of the final alignment (OutPrefix.ali.png) can NOT be created. See the final alignement (OutPrefix.ali.fasta).")
        logger.info("Write the final alignment in OutPrefix.ali.fasta")
        ApytramLib.ApytramNeeds.write_in_file(Query.MafftResult, "%s.ali.fasta" %(Query.OutPrefixName))
        logger.debug("Writing alignment plot and fasta --- %s seconds ---", str(time.time() - start_output_ali))

    logger.debug("Writing outputs --- %s seconds ---", str(time.time() - start_output))


end_time = str(time.time() - start_time)
logger.info("--- %s seconds ---", end_time )
logger.warning("END (%s seconds)", end_time)
ApytramLib.ApytramNeeds.end(0, TmpDirName, keep_tmp=args.keep_tmp)
