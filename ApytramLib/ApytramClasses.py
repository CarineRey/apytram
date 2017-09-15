#!/usr/bin/python
# coding: utf-8

# File: ApytramNeeds.py
# Created by: Carine Rey
# Created on: June 2016
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


import os
import re
import string
import numpy as np
import sys
import subprocess
import time


import ApytramNeeds
import Aligner
import BlastPlus
import Trinity
import Seqtk

import pandas
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

#### Execution statistics class
class Exec_stats(object):
    def __init__(self,start_time):
        New_Time_stat_dic = {"SPDatabaseBuilding": 0,
                             "QDatabaseBuilding": 0,
                             "SPDatabasePreparation": 0,
                             "Blast_fish": 0,
                             "index_db": 0,
                             "Seqtk": 0,
                             "Blastdbcmd": 0,
                             "Trinity": 0,
                             "Exonerate_on_ref":0,
                             "Exonerate_between_iter":0,
                             "Blast_on_ref":0,
                             "Mafft":0,
                             "Prep_fasta_index":0,
                             "Prep_clstr_index":0,
                             "Python":0
                             }

        New_Iter_stat_dic = {"IterationTime": 0,
                             "CumulTime": 0,
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
                             }

        self.TimeStatsDict =  {0:New_Time_stat_dic}
        self.IterStatsDict =  {0:New_Iter_stat_dic}
        self.TimeAtributes = New_Time_stat_dic.keys()
        self.IterAtributes = New_Iter_stat_dic.keys()
        self.StartTime = start_time
        self.StartIterTime = start_time

    def new_iteration(self, i):
        self.TimeStatsDict[i] = {}
        self.IterStatsDict[i] = {"CumulTime": self.IterStatsDict[i-1]["CumulTime"]}



#### RNA sample class

class RNA_species(object):
    def __init__(self,start_time,logger):
        self.logger = logger
        self.Species = ""
        self.OutputPrefix = ""
        self.TmpDirName = ""

        self.Fasta = []
        self.Fastq = []
        self.InputFastaFilename = ""
        self.IndexFilename = ""
        self.IndexDB = None
        self.ClstrFilename = ""
        self.ClstrIndexFilename = ""
        self.ClstrIndexDB = None
        self.ClstrRepFilename = ""
        self.DatabaseName = ""
        self.DatabaseDirName = ""
        self.FormatedDatabase = False

        self.DatabaseType = "" # ["single","paired","FR","RF","F","R"]
        self.StrandedData = False # True/False
        self.PairedData = False # True/False
        self.Evalue = 0
        self.MinLength = 0
        self.MinIdentityPercentage = 0
        self.MinAliLength = 0
        self.FinalMinLength = 0
        self.FinalMinIdentityPercentage = 0
        self.FinalMinAliLength = 0

        self.ExecutionStats = Exec_stats(start_time)
        self.CurrentIteration = 0
        self.FinalIteration = 0
        self.Improvment = True
        self.CompletedIteration = True
        self.Finished = False

        # Intermediary files

        self.ReadNamesFilename = ""
        self.ParsedReadNamesFilename = ""
        self.TrinityFastaFilename = ""
        self.FilteredTrinityFastaFilename = ""
        self.HomologyOnRefResult = ""

        # Constante

        self.TrinityExonerateRyo = "%ti\t%qi\t%ql\t%tal\t%tl\t%tab\t%tae\t%s\t%pi\t%qab\t%qae\n"

    def add_time_statistic(self, attribute, start = 0, inter = 0, end = 0):
        if not start:
            start = self.ExecutionStats.StartIterTime

        if not end:
            end = time.time()

        if attribute in self.ExecutionStats.TimeAtributes:
            self.ExecutionStats.TimeStatsDict[self.CurrentIteration][attribute] = end - inter - start

    def get_time_statistic(self,attribute):
        if attribute in self.ExecutionStats.TimeStatsDict[self.CurrentIteration].keys():
            return(self.ExecutionStats.TimeStatsDict[self.CurrentIteration][attribute])
        else:
            return (0)

    def add_iter_statistic(self, attribute, value, mode_add = False):
        if attribute in self.ExecutionStats.IterAtributes:
            if mode_add:
                self.ExecutionStats.IterStatsDict[self.CurrentIteration].setdefault(attribute,0)
                self.ExecutionStats.IterStatsDict[self.CurrentIteration][attribute] += value
            else:
                self.ExecutionStats.IterStatsDict[self.CurrentIteration][attribute] = value

    def get_iter_statistic(self,attribute, IterNb = None, RelIter = 0):
        if not IterNb:
            IterNb = self.CurrentIteration
        if RelIter:
            IterNb += RelIter

        if IterNb in self.ExecutionStats.IterStatsDict.keys():
            if attribute in self.ExecutionStats.IterStatsDict[IterNb].keys():
                return(self.ExecutionStats.IterStatsDict[IterNb][attribute])
        else:
            return (0)

    def has_a_formated_database(self):
        CheckDatabase_BlastdbcmdProcess = BlastPlus.Blastdbcmd(self.DatabaseName,"","")
        return CheckDatabase_BlastdbcmdProcess.is_database()

    def set_TmpDir(self,TmpDirName):
        self.TmpDirName = TmpDirName
        if not os.path.isdir(self.TmpDirName):
            self.logger.info("The temporary directory %s does not exist, it will be created" % (self.TmpDirName))
            os.makedirs(self.TmpDirName)

    def set_OutputDir(self,OutPrefixName):
        ### Set up the output directory
        if OutPrefixName:
            OutDirName = os.path.dirname(OutPrefixName)
            if os.path.isdir(OutDirName):
                self.logger.info("The output directory %s exists" %(OutDirName) )
            elif OutDirName: # if OutDirName is not a empty string we create the directory
                self.logger.info("The output directory %s does not exist, it will be created" % (OutDirName))
                os.makedirs(OutDirName)
        self.OutDirPrefix = OutPrefixName

    def prepare_database(self, FreeSpaceTmpDir, TmpDirName):
        start = time.time()
        self.logger.info("Prepare data to build for the species: %s" %(self.Species))
        self.DatabaseDirName = os.path.dirname(self.DatabaseName)
        if os.path.isdir(self.DatabaseDirName) or not self.DatabaseDirName :
            self.logger.info("Database directory exists")
        else:
            self.logger.info("Database directory does not exist, we create it")
            os.makedirs(self.DatabaseDirName)

        #Get the available free space of the DB dir
        FreeSpaceDBDir = ApytramNeeds.get_free_space(self.DatabaseDirName)
        self.logger.debug("%s free space in %s" %(FreeSpaceDBDir,self.DatabaseDirName))

        #Concatenate input files
        if self.Fastq:
            for fastq in self.Fastq:
                if not os.path.isfile(fastq):
                    logger.error("The fastq file (%s) does not exist." %(fastq))
                    ApytramNeeds.end(1,self.TmpDirName,keep_tmp = self.keep_tmp)
                # Format the fastq file in fasta
            self.InputFastaFilename = "%s/input_fastq.fasta" %(self.TmpDirName)
            self.logger.info("Convert the fastq file in fasta format")
            start_convert = time.time()
            input_size = ApytramNeeds.get_size(self.Fastq)
            self.logger.debug("size of input files: %s " %input_size)
            if input_size*1.1 > FreeSpaceTmpDir:
                self.logger.error("Not enough available free space in %s to convert input files" %TmpDirName)
                ApytramNeeds.end(1,self.TmpDirName,keep_tmp = self.keep_tmp)
            if input_size*2.5 > FreeSpaceDBDir:
                self.logger.error("Not enough available free space in %s to build the database" %DatabaseDirName)
                ApytramNeeds.end(1,self.TmpDirName,keep_tmp = self.keep_tmp)
            (out,err) = ApytramNeeds.fastq2fasta(" ".join(self.Fastq),self.InputFastaFilename)
            if err:
                self.logger.error(err)
                ApytramNeeds.end(1, self.TmpDirName, keep_tmp=self.keep_tmp)
            self.logger.info("Convertion takes %s seconds" %(time.time() - start_convert))
        elif self.Fasta:
            for fasta in self.Fasta:
                if not os.path.isfile(fasta):
                    self.logger.error("The fasta file (%s) does not exist." %fasta)
                    ApytramNeeds.end(1,self.TmpDirName,keep_tmp = self.keep_tmp)
            # Check enough available free space
            input_size = ApytramNeeds.get_size(self.Fasta)
            self.logger.debug("size of input files: %s " %input_size)
            if input_size*2.5 > FreeSpaceDBDir:
                self.logger.error("Not enough available free space in %s to build the database" %DatabaseDirName)
                ApytramNeeds.end(1,self.TmpDirName,keep_tmp = self.keep_tmp)
            # Concatenate fasta files
            if len(self.Fasta) > 1:
                if input_size*1.1 > FreeSpaceTmpDir:
                    self.logger.error("Not enough available free space in %s to concatenate input files" %TmpDirName)
                    ApytramNeeds.end(1,self.TmpDirName,keep_tmp = self.keep_tmp)
                self.InputFastaFilename = "%s/input_fasta.fasta" %(self.TmpDirName)
                self.logger.info("Concatenate fasta files")
                start_convert = time.time()
                (out, err) = ApytramNeeds.cat_fasta(" ".join(args.fasta),InputFastaFilename)
                if err:
                    self.logger.error(err)
                    ApytramNeeds.end(1,self.TmpDirName,keep_tmp = self.keep_tmp)
                self.logger.info("Concatenation takes %s seconds" %(time.time() - start_convert))
            else:
                self.InputFastaFilename = self.Fasta[0]

        #Check if the end of sequence name of paired data are 1 or 2
        if self.PairedData:
            self.logger.debug("Fasta file: %s",self.InputFastaFilename)
            BadReadName = ApytramNeeds.check_paired_data(self.InputFastaFilename)
            if BadReadName:
                self.logger.error("Paired read names must finished by 1 or 2. %s is uncorrect" %(BadReadName))
                ApytramNeeds.end(1,self.TmpDirName,keep_tmp = self.keep_tmp)

        self.add_time_statistic("DatabasePreparation", start = start)
        #Build blast formated database from a fasta file
        if not os.path.isfile(self.InputFastaFilename):
            self.logger.error("Error during concatenation or conversion of input files. %s is not a file" %(self.InputFastaFilename))
            ApytramNeeds.end(1,self.TmpDirName,keep_tmp = self.keep_tmp)


    def build_database(self, FreeSpaceTmpDir, TmpDirName):
        start = time.time()

        if os.path.isfile(self.ClstrRepFilename):
            # Database building
            self.logger.info(self.DatabaseName + " database building")
            MakeblastdbProcess = BlastPlus.Makeblastdb(self.ClstrRepFilename, self.DatabaseName)
            (out,err) = MakeblastdbProcess.launch()
        elif os.path.isfile(self.InputFastaFilename):
            self.logger.info(self.DatabaseName + " database building")
            MakeblastdbProcess = BlastPlus.Makeblastdb(self.InputFastaFilename, self.DatabaseName)
            (out,err) = MakeblastdbProcess.launch()
        else:
            self.logger.error("Error during concatenation or conversion of input files. %s is not a file" %(self.InputFastaFilename))
            ApytramNeeds.end(1,self.TmpDirName,keep_tmp = self.keep_tmp)


        self.FormatedDatabase = self.has_a_formated_database()

        if not self.FormatedDatabase:
            self.logger.error("Problem in the database building.\nAre you sure of your input format?\nAre all read names unique?")
            self.logger.info("Database %s does not exist" % self.DatabaseName)
            ApytramNeeds.end(1,self.TmpDirName,keep_tmp = self.keep_tmp)
        else:
            self.add_time_statistic("SPDatabaseBuilding", start = start)
            self.logger.info("Database %s build in %s" %(self.DatabaseName,self.get_time_statistic("SPDatabaseBuilding")))

    def get_all_reads(self):
         self.InputFastaFilename = "%s/input_fastq.fasta" %(self.TmpDirName)
         BlastdbcmdProcess = BlastPlus.Blastdbcmd(self.DatabaseName, "", "")
         (out,err) = BlastdbcmdProcess.get_all_seq(self.InputFastaFilename)

    def new_iteration(self):

        self.ExecutionStats.StartIterTime = time.time()

        # Cleaning

        # Previous Iteration
        self.PreviousReadNamesFilename = self.ParsedReadNamesFilename
        self.PreviousFilteredTrinityFastaFilename = self.FilteredTrinityFastaFilename

        # New names files

        self.CurrentIteration +=1
        self.Improvment = True
        self.CompletedIteration = True

        self.ExecutionStats.new_iteration(self.CurrentIteration)

        self.ReadNamesFilename = "%s/ReadNames.%d.txt" %(self.TmpDirName,self.CurrentIteration)
        self.ParsedReadNamesFilename = "%s/ParsedReadNames.%d.txt" %(self.TmpDirName,self.CurrentIteration)
        self.ReadsNumber = 0

        self.ReadNamesFilename_Right = "%s/ReadNames.%d.1.txt" % (self.TmpDirName,self.CurrentIteration)
        self.ReadNamesFilename_Left  = "%s/ReadNames.%d.2.txt" % (self.TmpDirName,self.CurrentIteration)
        self.ReadFastaFilename_Right = "%s/Reads.%d.1.fasta"   % (self.TmpDirName,self.CurrentIteration)
        self.ReadFastaFilename_Left  = "%s/Reads.%d.2.fasta"   % (self.TmpDirName,self.CurrentIteration)

        self.ReadFastaFilename      = "%s/Reads.%d.fasta"   % (self.TmpDirName,self.CurrentIteration)

        self.TrinityFastaFilename = "%s/Trinity_iter_%d" %(self.TmpDirName,self.CurrentIteration)

        self.HomologyOnRefFilename = "%s/Trinity_iter_%d.blast_on_ref" %(self.TmpDirName,self.CurrentIteration)
        self.HomologyOnRefResult = ""

        self.HomologyBetweenIterFilename = "%s/iter_%d_%d" %(self.TmpDirName,self.CurrentIteration -1, self.CurrentIteration)

    def fish_reads(self, BaitSequencesFilename,Threads):

        self.launch_Blastn(BaitSequencesFilename.replace(".s2","").replace(".s",""),Threads, ext="test_without_reduce")
        self.launch_Blastn(BaitSequencesFilename,Threads)

    def launch_Blastn(self,BaitSequencesFilename,Threads,ext=""):
        start = time.time()
        self.logger.info("Blast bait sequences on reads database")
        BlastnProcess = BlastPlus.Blast("blastn", self.DatabaseName, BaitSequencesFilename)
        BlastnProcess.Evalue = self.Evalue
        BlastnProcess.Task = "blastn"
        BlastnProcess.Threads = Threads
        BlastnProcess.OutFormat = "6 sacc" #qseq
        if  self.CurrentIteration == 1: # first iteration bait sequence can be divergente
            BlastnProcess.perc_identity = 85
        else:
            BlastnProcess.perc_identity = 97 # other iteration -> less authorized divergente because same species

        (out,err) = BlastnProcess.launch(self.ReadNamesFilename+ext)
        self.add_time_statistic("Blast_fish", start = start)
        self.logger.info("End Blast (%s seconds) %s",self.get_time_statistic("Blast_fish"),ext)

    def get_read_sequences(self, Threads, Memory, meth="seqtk"):
        if self.PairedData:
            self.logger.info("Split read names depending on 1/ or 2/")
            (out, err) = ApytramNeeds.split_readnames_in_right_left(self.ParsedReadNamesFilename,self.ReadNamesFilename_Right,self.ReadNamesFilename_Left)
            if err:
                self.logger.error(err)
            StrandList = [".1",".2"]
        else:
            StrandList = [""]

        self.logger.info("Retrieve reads sequences")
        start_blastdbcmd_time = time.time()
        for strand in StrandList:
            ReadFastaFilename = "%s/Reads.%d%s.fasta" %(self.TmpDirName,self.CurrentIteration,strand)
            ReadNamesFilename = "%s/ReadNames.%d%s.txt" % (self.TmpDirName,self.CurrentIteration,strand)
            start = time.time()
            if (not os.path.isfile(ReadFastaFilename)) or (os.stat(ReadFastaFilename).st_size == 0):
                if meth == "index":
                    nb_reads = ApytramNeeds.retrieve_reads_from_index(self.IndexDB, ReadNamesFilename, ReadFastaFilename)
                    self.logger.debug("Add %s reads (get from index DB)" %(nb_reads))
                    self.add_time_statistic("index_db", start = start)
                    self.logger.debug("index_db %s --- %s seconds ---" %(strand, self.get_time_statistic("index_db")))

                elif meth == "blastdbcmd":
                    BlastdbcmdProcess = BlastPlus.Blastdbcmd(self.DatabaseName, ReadNamesFilename, ReadFastaFilename)
                    (out,err) = BlastdbcmdProcess.launch()
                    self.add_time_statistic("Blastdbcmd", start = start)
                    self.logger.debug("Blastdbcmd %s --- %s seconds ---" %(strand, self.get_time_statistic("Blastdbcmd")))

                elif meth == "index_seqtk":
                    NewReadNamesFilename = "%s/NewReads.%d%s.txt" %(self.TmpDirName,self.CurrentIteration, strand)
                    CommonReadNamesFilename = "%s/CommonReads.%d%s.txt" %(self.TmpDirName,self.CurrentIteration, strand)
                    PreviousReadNamesFilename = "%s/ReadNames.%d%s.txt" % (self.TmpDirName,self.CurrentIteration-1,strand)
                    PreviousReadFastaFilename = "%s/Reads.%d%s.fasta" %(self.TmpDirName,self.CurrentIteration-1,strand)

                    #CommonReadNamesFilename = Get read in ReadNamesFilename and in PreviousReadNamesFilename
                    Nb_common_reads = ApytramNeeds.common_reads(PreviousReadNamesFilename, ReadNamesFilename, CommonReadNamesFilename)

                    #NewReadNamesFilename = Get read in ReadNamesFilename and not in CommonReadNamesFilename
                    (Nb_new_reads, NewReadNamesFilename) = ApytramNeeds.new_reads(CommonReadNamesFilename, ReadNamesFilename, NewReadNamesFilename)

                    if Nb_new_reads:
                        nb_reads = ApytramNeeds.retrieve_reads_from_index(self.IndexDB, NewReadNamesFilename, ReadFastaFilename)
                        self.logger.debug("Add %s new reads (get from index DB)" %(nb_reads))
                        self.add_time_statistic("index_db", start = start)
                        self.logger.debug("index_db %s --- %s seconds ---" %(strand, self.get_time_statistic("index_db")))

                    if Nb_common_reads:
                        SeqtkProcess = Seqtk.Seqtk(fasta=PreviousReadFastaFilename)
                        self.logger.debug("Add %s old reads (get from previous iteration)" %(Nb_common_reads))
                        (out,err) = SeqtkProcess.launch_fasta_subseq(CommonReadNamesFilename, ReadFastaFilename, mode="a")

                    self.add_time_statistic("index_db", start = start)
                    self.logger.debug("index_db Seqtk (total) %s --- %s seconds ---" %(strand, self.get_time_statistic("index_db")))

                elif meth == "blastdbcmd_seqtk":

                    NewReadNamesFilename = "%s/NewReads.%d%s.txt" %(self.TmpDirName,self.CurrentIteration, strand)
                    CommonReadNamesFilename = "%s/CommonReads.%d%s.txt" %(self.TmpDirName,self.CurrentIteration, strand)
                    PreviousReadNamesFilename = "%s/ReadNames.%d%s.txt" % (self.TmpDirName,self.CurrentIteration-1,strand)
                    PreviousReadFastaFilename = "%s/Reads.%d%s.fasta" %(self.TmpDirName,self.CurrentIteration-1,strand)

                    #CommonReadNamesFilename = Get read in ReadNamesFilename and in PreviousReadNamesFilename
                    Nb_common_reads = ApytramNeeds.common_reads(PreviousReadNamesFilename, ReadNamesFilename, CommonReadNamesFilename)

                    #NewReadNamesFilename = Get read in ReadNamesFilename and not in CommonReadNamesFilename
                    (Nb_new_reads, NewReadNamesFilename) = ApytramNeeds.new_reads(CommonReadNamesFilename, ReadNamesFilename, NewReadNamesFilename)

                    if Nb_new_reads:
                        BlastdbcmdProcess = BlastPlus.Blastdbcmd(self.DatabaseName, NewReadNamesFilename, ReadFastaFilename)
                        self.logger.debug("Add %s new reads (get from DB)" %(Nb_new_reads))
                        (out,err) = BlastdbcmdProcess.launch()

                    if Nb_common_reads:
                        SeqtkProcess = Seqtk.Seqtk(fasta=PreviousReadFastaFilename)
                        self.logger.debug("Add %s old reads (get from previous iteration)" %(Nb_common_reads))
                        (out,err) = SeqtkProcess.launch_fasta_subseq(CommonReadNamesFilename, ReadFastaFilename, mode="a")

                    self.add_time_statistic("Blastdbcmd", start = start)
                    self.logger.debug("Blastdbcmd Seqtk (total) %s --- %s seconds ---" %(strand, self.get_time_statistic("Blastdbcmd")))

                elif meth == "seqtk":
                    SeqtkProcess = Seqtk.Seqtk(fasta=self.InputFastaFilename)
                    (out,err) = SeqtkProcess.launch_fasta_subseq(ReadNamesFilename, ReadFastaFilename)
                    self.add_time_statistic("Seqtk", start = start)
                    self.logger.debug("Seqtk %s --- %s seconds ---" %(strand, self.get_time_statistic("Seqtk")))
                else:
                    self.logger.error("No method to retrieve read")
            else:
                self.logger.warn("%s has already been created, it will be used" %(ReadFastaFilename) )

    def launch_Trinity(self, Threads, Memory, long_read=False):
        start = time.time()
        self.logger.info("Launch Trinity")
        ExitCode = 0
        if self.PairedData:
            TrinityProcess = Trinity.Trinity(self.TrinityFastaFilename,
                                             right = self.ReadFastaFilename_Right,
                                             left = self.ReadFastaFilename_Left)
        else:
            TrinityProcess = Trinity.Trinity(self.TrinityFastaFilename, single = self.ReadFastaFilename)
        if self.StrandedData:
            TrinityProcess.SS_lib_type = self.DatabaseType
        # If there is a huge number of reads, remove duplicated reads
        if self.ReadsNumber < 500:
            TrinityProcess.NoNormalizeReads = True
        if long_read and os.path.isfile(self.PreviousFilteredTrinityFastaFilename):
            TrinityProcess.LongReads = self.PreviousFilteredTrinityFastaFilename

        TrinityProcess.CPU = Threads
        TrinityProcess.max_memory = Memory
        # Keep only contig with a length superior to MinLength
        TrinityProcess.MinLength = self.MinLength

        # Use the  --full_cleanup Trinity option to keep only the contig file
        TrinityProcess.FullCleanup = True
        if not os.path.isfile(self.TrinityFastaFilename+".Trinity.fasta"):
            (out,err,ExitCode) = TrinityProcess.launch()
        else:
            self.logger.warn("%s has already been created, it will be used" %(self.TrinityFastaFilename+".Trinity.fasta") )

        self.TrinityFastaFilename = self.TrinityFastaFilename + ".Trinity.fasta"

        if not os.path.isfile(self.TrinityFastaFilename):
            if ExitCode == 2 or ExitCode == 0 : # Trinity exit 0 if "No butterfly assemblies to report"
               self.logger.debug("Trinity found nothing...\n[...]\n"+"\n".join(out.strip().split("\n")[-15:]))
               self.logger.warning("Trinity has assembled no contigs at the end of the iteration %s (ExitCode: %d)" %(self.CurrentIteration,ExitCode) )
            elif ExitCode == 255:
               self.logger.debug("Trinity found nothing...\n[...]\n"+"\n".join(out.strip().split("\n")[-15:]))
               self.logger.error("Trinity has crashed (ExitCode: %d). Do you use the last version of Trinity (>= 2.3)?" %(ExitCode))
            elif ExitCode != 0:
               self.logger.debug("Trinity found nothing...\n[...]\n"+"\n".join(out.strip().split("\n")[-15:]))
               self.logger.error("Trinity has crashed (ExitCode: %d). Are all dependencies satisfied?" %(ExitCode))


        self.add_time_statistic("Trinity", start = start)
        self.logger.info("End Trinity (%s seconds)" %(self.get_time_statistic("Trinity")))

    def get_homology_between_trinity_results_and_references(self,Query, method="blastn"):
        if method == "exonerate":
            # Use Exonerate
            start = time.time()
            TrinityExonerateProcess = Aligner.Exonerate(Query.RawQuery, self.TrinityFastaFilename)
            # Keep only the best hit for each contig from Trinity
            TrinityExonerateProcess.Bestn = 1
            TrinityExonerateProcess.Model = "ungapped" #aligned:affine #"est2genome" #cdna2genome ##too long !!!
            self.HomologyOnRefFilename = "%s/Trinity_iter_%d.exonerate_un" %(self.TmpDirName,self.CurrentIteration)
            # Customize the output format
            TrinityExonerateProcess.Ryo = self.TrinityExonerateRyo
            (out,err,self.HomologyOnRefResult) = TrinityExonerateProcess.get_output()
            # Write the result in a file
            if self.HomologyOnRefResult:
                ApytramNeeds.write_in_file(self.HomologyOnRefResult, self.HomologyBetweenIterFilename)
            else:
                self.logger.info("Reconstructed sequences but no homologous with references")
                self.logger.info("Try to get homologies with a more sensible model")
                ### Try to get homologies with a more sensible model
                self.HomologyOnRefFilename =   "%s/Trinity_iter_%d.exonerate_al" %(self.TmpDirName,self.CurrentIteration)
                TrinityExonerateProcess = Aligner.Exonerate(Query.RawQuery, self.TrinityFastaFilename)
                # Keep only the best hit for each contig from Trinity
                TrinityExonerateProcess.Bestn = 1
                TrinityExonerateProcess.Model = "affine:local" #"coding2genome"
                # Customize the output format
                TrinityExonerateProcess.Ryo = self.TrinityExonerateRyo
                (out,err,self.HomologyOnRefResult) = TrinityExonerateProcess.get_output()
                # Write the result in a file
                if self.TrinityExonerateResult:
                    ApytramNeeds.write_in_file(self.HomologyOnRefResult,self.HomologyBetweenIterFilename)

            self.add_time_statistic("Exonerate_on_ref", start = start)
            self.logger.debug("Exonerate_on_ref --- %s seconds ---" %(self.get_time_statistic("Exonerate_on_ref")))
        else:
            # Use Blast
            start = time.time()
            BlastnProcess = BlastPlus.Blast( "blastn", Query.QueryDB, self.TrinityFastaFilename)
            BlastnProcess.Evalue = 1e-5
            BlastnProcess.Task = "blastn"
            BlastnProcess.Threads = 1
            BlastnProcess.max_target_seqs = 1
            BlastnProcess.OutFormat = "6 sacc qacc qlen length slen sstart send score pident qstart qend "
            self.HomologyOnRefFilename = "%s/Trinity_iter_%d.blast_on_ref" %(self.TmpDirName,self.CurrentIteration)

            (out,err) = BlastnProcess.launch(self.HomologyOnRefFilename)
            self.add_time_statistic("Blast_on_ref", start = start)
            self.logger.info("End Blast (%s seconds)" %(self.get_time_statistic("Blast_on_ref")))

        if os.path.isfile(self.HomologyOnRefFilename) and (os.stat(self.HomologyOnRefFilename).st_size != 0):
            with open(self.HomologyOnRefFilename, "r") as FILEONREF:
                self.HomologyOnRefResult = FILEONREF.read()

    def read_and_parse_homology_results(self, final_iteration = False):
        "Return a list of Sequences if the identity percentage is superior to MinIdentityPercentage and the alignment length is superior to MinAliLen "

        if final_iteration:
            self.logger.info("Filter sequence with a identity percentage superior to %d and a percentage alignment len %d" %(self.FinalMinIdentityPercentage, self.FinalMinAliLength))
            minidentypercentage =  self.FinalMinIdentityPercentage
            minalilength = 0
            minlengthpercentage = self.FinalMinLength
            minalilengthpercentage = self.FinalMinAliLength
        else:
            self.logger.info("Filter sequence with a identity percentage superior to %d and a alignment len %d" %(self.MinIdentityPercentage, self.MinAliLength))
            minidentypercentage = self.MinIdentityPercentage
            minalilength = self.MinAliLength
            minlengthpercentage = 0
            minalilengthpercentage = 0

        self.FilteredTrinityFasta = ApytramNeeds.Fasta()
        BestScoreNames = {}

        HomologyScoreList = self.HomologyOnRefResult.strip().split("\n")

        for line in HomologyScoreList:
            ListLine = line.split("\t")
            (ti,qi) = ListLine[:2]
            (ql,tal,tl,tab,tae,score,pi,qab,qae) = [float(x) for x in ListLine[2:]]

            if (pi >=  minidentypercentage) and (tal >= minalilength) and (ql >= minlengthpercentage*tl/100) and (tal >= minalilengthpercentage*tl/100) :
                # We keep this sequence
                # A same sequence can be present 2 time if the hit scores are identical.
                # We keep only the first
                ever_seen = False
                for old_sequence in self.FilteredTrinityFasta.Sequences:
                    if old_sequence.Name == qi:
                        ever_seen = True
                        break
                if not ever_seen:
                    new_sequence = ApytramNeeds.Sequence()
                    new_sequence.add_attribute_from_exonerate(ListLine)

                    self.add_iter_statistic("TotalIdentity", pi, mode_add = True)
                    if self.get_iter_statistic("BestIdentity") <= pi:
                        self.add_iter_statistic("BestIdentity",pi)

                    self.add_iter_statistic("TotalLength", ql, mode_add = True)
                    if self.get_iter_statistic("BestLength") <= ql:
                        self.add_iter_statistic("BestLength", ql)

                    self.add_iter_statistic("TotalScore", score, mode_add = True)
                    if self.get_iter_statistic("BestScore") <= score :
                        self.add_iter_statistic("BestScore", score)
                        BestScoreNames[ti] = qi

                    # Check if the sequence is complement
                    if ((qae-qab)*(tae-tab) < 0):
                        new_sequence.Complement = True

                    self.FilteredTrinityFasta.append(new_sequence)

        for (Target,Query) in BestScoreNames.items():
            for i in range(len(self.FilteredTrinityFasta.Sequences)):
                if Query == self.FilteredTrinityFasta.Names[i]:
                    self.FilteredTrinityFasta.Sequences[i].BestSequence = Target

        NbContigs = len(self.FilteredTrinityFasta.Sequences)
        AverageIdentity = 0
        AverageLength = 0
        AverageScore = 0

        self.add_iter_statistic("NbContigs",NbContigs)

        if NbContigs:
            AverageIdentity = self.get_iter_statistic("TotalIdentity") / NbContigs
            AverageLength   = self.get_iter_statistic("TotalLength") / NbContigs
            AverageScore    = self.get_iter_statistic("TotalScore") / NbContigs

        self.add_iter_statistic("AverageIdentity", AverageIdentity)
        self.add_iter_statistic("AverageLength", AverageLength)
        self.add_iter_statistic("AverageScore", AverageScore)

    def filter_trinity_results_according_homology_results(self, final_iteration = False):

        if final_iteration:
            self.TrinityFastaFilename = self.FilteredTrinityFastaFilename


        self.FilteredTrinityFastaFilename = "%s/Trinity_iter_%d.filtered.fasta" %(self.TmpDirName,self.CurrentIteration)

        # Get and save filtered sequence names and their homology score in self.FilteredTrinityFasta
        self.read_and_parse_homology_results(final_iteration = final_iteration)

        if self.FilteredTrinityFasta.Sequences: # If sequences pass the filter
            # Read fasta
            TrinityFasta = ApytramNeeds.Fasta()
            TrinityFasta.read_fasta(FastaFilename = self.TrinityFastaFilename)

            # get sequence for filtered sequences in the trinityfasta
            self.FilteredTrinityFasta.complete_fasta(TrinityFasta)

            # rename seq with species name
            ## buil dic
            new_name_dic = { n: self.Species+"_"+n for n in self.FilteredTrinityFasta.Names}
            self.FilteredTrinityFasta = self.FilteredTrinityFasta.rename_fasta(new_name_dic)
            # Write fasta
            self.FilteredTrinityFasta.write_fasta(self.FilteredTrinityFastaFilename)

    def rename_sequences(self):
        # build dictionnary to rename sequences
        NewnameDict = {}
        i = 0
        for Sequence in self.FilteredTrinityFasta.Sequences:
            Message = self.Species + "_"
            i +=1
            OldName = Sequence.Name
            if Sequence.BestSequence:
                Message += "Best_"
            NewName = "APYTRAM_%s%d.len=%d.[%s.id=%d.len=%d]" %(Message,i,Sequence.ql,Sequence.ti,Sequence.pi,Sequence.tl)
            NewnameDict[OldName] = NewName

        self.FilteredTrinityFasta = self.FilteredTrinityFasta.rename_fasta(NewnameDict)

    def compare_current_and_previous_iterations(self, method="blastn"):
        self.logger.info("Refind the \"parent\" contig from the previous contig for each contig and check they are different")
        start = time.time()
        HomologyResult = ""
        if method == "exonerate":
            ExonerateProcess = Aligner.Exonerate(self.FilteredTrinityFastaFilename, self.PreviousFilteredTrinityFastaFilename)
            # Keep only the best hit for each contigs
            ExonerateProcess.Bestn = 1
            ExonerateProcess.Model =  "ungapped"
            # Customize the output format
            self.HomologyBetweenIterFilename = "%s/iter_%d_%d.exonerate" %(self.TmpDirName,self.CurrentIteration -1, self.CurrentIteration)

            ExonerateProcess.Ryo = "%ti\t%qi\t%ql\t%qal\t%tal\t%tl\t%pi\n"
            (out,err,HomologyResult) = ExonerateProcess.get_output()

            ApytramNeeds.write_in_file(HomologyResult,self.HomologyBetweenIterFilename)

            self.add_time_statistic("Exonerate_between_iter", start = start)
            self.logger.debug("Exonerate_between_iter --- %s seconds ---" %(self.get_time_statistic("Exonerate_between_iter")))

        else:

            if os.path.isfile(self.TrinityFastaFilename):
                tmp_DB_name = "%s/Trinity_iter_%d.db" %(self.TmpDirName,(self.CurrentIteration - 1))
                self.logger.info(tmp_DB_name + " database building")
                # Build DB
                MakeblastdbProcess = BlastPlus.Makeblastdb(self.TrinityFastaFilename, tmp_DB_name)
                (out,err) = MakeblastdbProcess.launch()
                # Blast ister i on iter i-1
                BlastnProcess = BlastPlus.Blast( "blastn", tmp_DB_name, self.PreviousFilteredTrinityFastaFilename)
                BlastnProcess.Evalue = 1e-5
                BlastnProcess.Task = "blastn"
                BlastnProcess.Threads = 1
                BlastnProcess.max_target_seqs = 1
                BlastnProcess.OutFormat = "6 sacc qacc qlen length length slen pident"
                self.HomologyBetweenIterFilename = "%s/iter_%d_%d.blast" %(self.TmpDirName,self.CurrentIteration -1, self.CurrentIteration)

                (out,err) = BlastnProcess.launch(self.HomologyBetweenIterFilename)
                self.add_time_statistic("Blast_between_iter", start = start)
                self.logger.info("End Blast (%s seconds)" %(self.get_time_statistic("Blast_between_iter")))

                with open(self.HomologyBetweenIterFilename, "r") as BLASTFILE:
                    HomologyResult = BLASTFILE.read()

        AlmostIdenticalResults = ApytramNeeds.check_almost_identical_exonerate_results(HomologyResult)
        if AlmostIdenticalResults:
            self.logger.info("Contigs are almost identical than the previous iteration (Same size (~98%), > 99% identity)")
            self.Improvment = False



    def measure_coverage(self,Query):
        if not Query.Aligned:
            Query.align_query()
        # Use Mafft
        start = time.time()
        MafftProcess = Aligner.Mafft(Query.AlignedQuery)
        MafftProcess.QuietOption = True
        MafftProcess.AutoOption = True
        #MafftProcess.AdjustdirectionOption = True
        MafftProcess.AddOption = self.FilteredTrinityFastaFilename
        (self.MafftResult,err) = MafftProcess.get_output()
        self.add_time_statistic("Mafft", start = start)

        self.logger.debug("Mafft --- %s seconds ---" %(self.get_time_statistic("Mafft")))

        (StrictCoverage, LargeCoverage, self.DicPlotCov) = ApytramNeeds.calculate_coverage(self.MafftResult,Query.ReferenceNames)

        self.add_iter_statistic("StrictCoverage", StrictCoverage)
        self.add_iter_statistic("LargeCoverage", LargeCoverage)
        self.logger.info("Strict Coverage: %s\tLarge Coverage: %s" %(StrictCoverage, LargeCoverage))

    def tmp_dir_clean_up(TmpDirName,i):
        if i == 1 :
            FilesToRemoves = ["%s/input_fastq.fasta" %TmpDirName]
        else:
            FilesToRemoves = ["%s/ReadNames.%d.txt" %(TmpDirName,i),
                            "%s/ReadNames.%d.1.txt" %(TmpDirName,i),
                            "%s/ReadNames.%d.2.txt" %(TmpDirName,i),
                            "%s/Reads.%d.fasta" %(TmpDirName,i),
                            "%s/Reads.%d.1.fasta" %(TmpDirName,i),
                            "%s/Reads.%d.2.fasta" %(TmpDirName,i),
                            "%s/Trinity_iter_%d.blast_on_ref" %(TmpDirName,i),
                            "%s/Trinity_iter_%d.exonerate_coding2g" % (TmpDirName, i),
                            "%s/Trinity_iter_%d.filtered.fasta" %(TmpDirName,i),
                            "%s/Trinity_iter_%d.Trinity.fasta" %(TmpDirName,i)]

    def end_iteration(self):

        iter_time = time.time() - self.ExecutionStats.StartIterTime
        self.ExecutionStats.IterStatsDict[self.CurrentIteration]["IterationTime"] = iter_time
        self.ExecutionStats.IterStatsDict[self.CurrentIteration]["CumulTime"] += iter_time

        NoPythonTime = 0

        for key in self.ExecutionStats.TimeStatsDict[self.CurrentIteration]:
            if key != "Python":
                NoPythonTime += self.get_time_statistic(key)

        self.HomologyOnRefResult = ""

        self.add_time_statistic("Python", inter = NoPythonTime )
        self.logger.debug("Python --- %s seconds ---" %(self.get_time_statistic("Python")))

    def new_query(self, Query):

        # Cleaning
        FilesToRemoves = [ self.TrinityFastaFilename,
                           self.FilteredTrinityFastaFilename]

        #Build a new tmp dir
        self.set_TmpDir("%s/%s" %(Query.TmpDirName, self.Species))

        self.ExecutionStats = Exec_stats(time.time())
        self.CurrentIteration = 0
        self.FinalIteration = 0
        self.Improvment = True
        self.CompletedIteration = True

        # Intermediary files

        self.ReadNamesFilename = ""
        self.TrinityFastaFilename = ""
        self.FilteredTrinityFastaFilename = ""
        self.HomologyOnRefResult = ""
        self.HomologyBetweenIterFilename = ""
        self.HomologyOnRefFilename = ""

        self.FilteredTrinityFasta = ApytramNeeds.Fasta()

    def get_output_fasta(self, fasta = "all"):
        assert fasta in ["all","best"], "fasta must be all or best"
        if fasta == "all":
            filename = "OutPrefix.fasta"
        else:
            filename = "OutPrefix.best.fasta"
        "Return fasta files with new names depending homology scores and options"
        self.logger.info("Prepare %s outputfile for %s" %(filename,self.Species))
        OutputFastaExtensions = []
        if fasta == "best":
            # Best sequences
            BestSequenceNames = []
            for Sequence in self.FilteredTrinityFasta.Sequences:
                if Sequence.BestSequence:
                    BestSequenceNames.append(Sequence.Name)
            return([str(self.FilteredTrinityFasta.filter_fasta(BestSequenceNames))])
        else:
            # Last iteration sequences
            return([str(self.FilteredTrinityFasta)])

    def get_stats(self):
        df_time = pandas.DataFrame(self.ExecutionStats.TimeStatsDict)
        df_iter = pandas.DataFrame(self.ExecutionStats.IterStatsDict)
        df = pandas.concat([df_time,df_iter]).T
        df.insert(0, "Iteration", df.index)
        df.insert(0, "Species", self.Species)
        return([df])



#### query class

class Query(object):
    def __init__(self,Name,QueryPath,logger):
        self.logger = logger
        self.Name = Name
        self.RawQuery = QueryPath
        self.SequenceNb = ApytramNeeds.count_sequences(QueryPath)
        self.FinalFastaFileName = ""

        #Read fasta
        QueryFasta = ApytramNeeds.Fasta()
        QueryFasta.read_fasta(FastaFilename = QueryPath)

        #Get Sequences names
        self.ReferenceNames = QueryFasta.Names

        self.AlignedQuery = ""
        self.QueryDB = ""
        self.Aligned = False

        self.CumulIteration = 0
        self.AbsIteration = 0

        self.Stop = False
        self.SpeciesWithoutImprovment = {0:[]}
        self.BaitSequences = ""
        self.PreviousBaitSequences = ""

        self.BestOutFileContent = []
        self.OutFileContent = []
        self.StatsFileContent = []
        self.TimeStatsDictList = []
        self.IterStatsDictList = []

    def initialization(self):
        #remove - in sequences
        self.logger.debug("Remove - in %s" ,self.RawQuery)
        query = ApytramNeeds.Fasta()
        query.read_fasta(FastaFilename=self.RawQuery)
        query.dealign_fasta()
        self.RawQuery = "%s/References.nogap.fasta" %(self.TmpDirName)
        self.QueryDB = "%s/References_DB" %(self.TmpDirName)
        query.write_fasta(self.RawQuery)
            # # If the -pep option is used, the -q option must be precised
            # if args.query_pep:
            #   if not os.path.isfile(args.query_pep):
            #        logger.error(args.query_pep+" (-pep) is not a file.")
            #        ApytramNeeds.end(1,self.TmpDirName,keep_tmp = self.keep_tmp)
            #
            #    if not args.query:
            #        logger.error("-pep option must be accompanied of the query in nucleotide format (-q option)")
            #        ApytramNeeds.end(1,self.TmpDirName,keep_tmp = self.keep_tmp)
                    # Database building
        self.logger.info(self.QueryDB + " database building")
        MakeblastdbProcess = BlastPlus.Makeblastdb(self.RawQuery, self.QueryDB)
        (out,err) = MakeblastdbProcess.launch()

    def align_query(self):
        self.AlignedQuery = self.RawQuery
        self.logger.warn("Align query %s", self.Name)
        if self.SequenceNb != 1:
            # If there are multiple probes, align them for the future coverage counter
            # Use Mafft
            aligned = ApytramNeeds.Fasta()
            aligned.read_fasta(FastaFilename=self.AlignedQuery)
            if not aligned.isalign():
                self.logger.debug("%s is not aligned, align it", self.RawQuery)
                start_mafft_time = time.time()
                MafftProcess = Aligner.Mafft(self.RawQuery)
                MafftProcess.QuietOption = True
                MafftProcess.AutoOption = True
                (MafftResult, err) = MafftProcess.get_output()
                self.AlignedQuery = "%s/References.ali.fasta" %(self.TmpDirName)
                ApytramNeeds.write_in_file(MafftResult, self.AlignedQuery)
                self.logger.debug("mafft --- %s seconds ---", str(time.time() - start_mafft_time))
            else:
                self.logger.warn("%s is already aligned, use it", self.RawQuery)
        self.Aligned = True

    def continue_iter(self):
        NbSpeciesWithoutImprovment = len(self.SpeciesWithoutImprovment[self.AbsIteration])
        if NbSpeciesWithoutImprovment == self.NbSpecies:
            return(False)
        elif self.Stop:
            return(False)
        else:
            return(True)

    def new_species_iteration(self,SpeciesList):
        self.CumulIteration += 1
        self.PreviousBaitSequences = self.BaitSequences
        self.BaitSequences = "%s/BaitSequences.%d.fasta" %(self.TmpDirName, self.CumulIteration)
        SpeciesCurrentReconstructedSequencesFileList = [Species.FilteredTrinityFastaFilename for Species in SpeciesList if Species.FilteredTrinityFastaFilename]

        if self.AbsIteration == 1 :
            SpeciesCurrentReconstructedSequencesFileList.append(self.RawQuery)

        if SpeciesCurrentReconstructedSequencesFileList:
            ApytramNeeds.cat_fasta(" ".join(SpeciesCurrentReconstructedSequencesFileList), self.BaitSequences)

        # Simplify baitsequences
        self.reduce_complexity_baitSequences()


    def reduce_complexity_baitSequences(self):
        start = time.time()
        Complex_BaitSequences = self.BaitSequences
        Simplified_BaitSequences = "%s/BaitSequences.%d.s.fasta" %(self.TmpDirName, self.CumulIteration)
        RepCluster_BaitSequences = "%s/BaitSequences.%d.s.fasta.rep" %(self.TmpDirName, self.CumulIteration)
        Cluster_BaitSequences = "%s/BaitSequences.%d.s.fasta.rep.clstr" %(self.TmpDirName, self.CumulIteration)

        # Run cd-hit-est to get cluster
        CdHitEstProcess = Aligner.CdHitEst(Complex_BaitSequences,RepCluster_BaitSequences, c=0.85)
        CdHitEstProcess.run()

        self.logger.info("reduce_complexity_baitSequences clustering --- %s seconds ---" %(time.time() - start))

        # Read clustr
        clstr_dic = ApytramNeeds.read_clstr(Cluster_BaitSequences)

        seq_ok = []
        nb_subseq = 0
        for (cluster,seq_l) in clstr_dic.items():
            rep_seq_name = ""
            rep_seq = ""
            if len(seq_l) == 1:
                rep_seq_name = seq_l[0][0]
                #Only 1 seq = rep
                file_with_seq_names = "%s/BaitSequences.%d.%s.tmp" %(self.TmpDirName, self.CumulIteration, cluster)
                ApytramNeeds.write_in_file("%s" %(rep_seq_name), file_with_seq_names)
                SeqtkProcess = Seqtk.Seqtk(fasta=self.BaitSequences)
                (f_rep_string,err) = SeqtkProcess.launch_fasta_subseq_stdout(file_with_seq_names)
                f_rep = ApytramNeeds.Fasta()
                f_rep.read_fasta(String = f_rep_string)
                rep_seq = f_rep.get(rep_seq_name)
                seq_ok.append((rep_seq_name+"_rep", rep_seq))
            else:
                #align all seqs
                # write all seq names in a file
                file_with_seq_names = "%s/BaitSequences.%d.%s.tmp" %(self.TmpDirName, self.CumulIteration, cluster)
                file_with_seq_fasta = "%s/BaitSequences.%d.%s.tmp.fasta" %(self.TmpDirName, self.CumulIteration, cluster)
                ApytramNeeds.write_in_file("\n".join([s for (s, _, _) in seq_l])+"\n", file_with_seq_names)
                # get sequence with seqtk
                SeqtkProcess = Seqtk.Seqtk(fasta=self.BaitSequences)
                (out,err) = SeqtkProcess.launch_fasta_subseq(file_with_seq_names, file_with_seq_fasta)
                # align file
                MafftProcess = Aligner.Mafft(file_with_seq_fasta)
                MafftProcess.QuietOption = True
                MafftProcess.AutoOption = True
                (ali_string,err) = MafftProcess.get_output()
                # read alignment
                ali = ApytramNeeds.Fasta()
                ali.read_fasta(String = ali_string)

                # get rep:
                for (seq, rep, rev) in seq_l:
                    if rep:
                        rep_seq_name = seq
                        #rep_seq = str(record_dict[rep_seq_name].seq)
                        break
                starti = time.time()
                for (seq_name, rep, rev) in seq_l:
                    if rep:
                        file_with_seq_names = "%s/BaitSequences.%d.%s.tmp1" %(self.TmpDirName, self.CumulIteration, cluster)
                        ApytramNeeds.write_in_file("%s" %(rep_seq_name), file_with_seq_names)
                        SeqtkProcess = Seqtk.Seqtk(fasta=self.BaitSequences)
                        (f_rep_string,err) = SeqtkProcess.launch_fasta_subseq_stdout(file_with_seq_names)
                        f_rep = ApytramNeeds.Fasta()
                        f_rep.read_fasta(String = f_rep_string)
                        rep_seq = f_rep.get(rep_seq_name)
                        seq_ok.append((rep_seq_name+"_rep", rep_seq))
                    else:
                        # compare seq with rep_seq:
                        aln_rep = ali.get(rep_seq_name)
                        aln_seq = ali.get(seq_name)

                        l_aln = len(aln_rep)
                        s_aln = [0] * l_aln
                        b = 0
                        while aln_seq[b] == "-":
                            b +=1
                        e = l_aln -1
                        while aln_seq[e] == "-":
                            e -=1
                        b_r = 0
                        while aln_rep[b_r] == "-":
                            b_r +=1
                        e_r = l_aln -1
                        while aln_rep[e_r] == "-":
                            e_r -=1

                        # SCORING
                        for i in range(b, e):

                            if aln_seq[i] == "-":
                                s_aln[i] = 2
                            else:
                                if aln_rep[i] != aln_seq[i]:
                                    s_aln[i] = 2
                            #print aln_rep[i] + " " + aln_seq[i] + ":" + str(s_aln[i])
                        # SCORE PARSING  and length measurement
                        ## make region:
                        lim_region = [0] + range(b_r, e_r+1, 25) + [e_r+1,l_aln]
                        len_lim_region = len(lim_region)
                        sp_l = [0] * len_lim_region
                        len_seq_l = [0] * len_lim_region

                        i_x = 0

                        while (i_x +1) < len_lim_region:
                            i=lim_region[i_x]
                            j=lim_region[i_x+1]
                            i_x+=1
                            len_seq_l[i_x] = len(aln_seq[i:j].replace("-",""))
                            if sum(s_aln[i:j]) > 2:
                                sp_l[i_x] = 1


                        # extend region with score > 10
                        e_sp_l = [0] * len_lim_region
                        for i_x in range(len_lim_region):
                            if sp_l[i_x]:
                                len_seq_to_add = len_seq_l[i_x]
                                e_sp_l[i_x] = 1
                                j = 0
                                while len_seq_to_add < 250 and j < (len_lim_region-1):
                                    j += 1
                                    if (i_x + j) in range(len_lim_region):
                                        e_sp_l[(i_x + j)] = 1
                                        len_seq_to_add += len_seq_l[i_x + j]
                                    if (i_x - j) in range(len_lim_region):
                                        e_sp_l[(i_x - j)] = 1
                                        len_seq_to_add += len_seq_l[i_x - j]

                        # detect and write sub sequences

                        subseq_to_keep = []
                        i_x = 0
                        while i_x < len_lim_region:
                            j_x = i_x
                            if e_sp_l[i_x]:
                                while e_sp_l[j_x] and j_x < (len_lim_region-1):
                                    j_x += 1
                                # write sub sequence:
                                i=lim_region[i_x]
                                j=lim_region[j_x]
                                seq_to_add = aln_seq[i:j].replace("-","")
                                subseq_to_keep.append(seq_to_add)
                            i_x = j_x
                            i_x += 1

                        seq_name_i = 0
                        for seq_i in subseq_to_keep:
                            seq_name_i +=1
                            nb_subseq +=1
                            seq_ok.append((seq_name+"_"+str(seq_name_i), seq_i))
                            seq_to_keep = []
        # Write new baitsequences file
        if seq_ok:
            string_seq_ok = []
            for (n,s) in seq_ok:
                string_seq_ok.append(">%s\n%s\n" %(n,s))
            ApytramNeeds.write_in_file("".join(string_seq_ok), Simplified_BaitSequences)

        print "#####"
        print "Avant:" + str(ApytramNeeds.count_sequences (self.BaitSequences))
        print "Avant:" + str(os.stat(self.BaitSequences).st_size)
        print "Après:" + str(ApytramNeeds.count_sequences (Simplified_BaitSequences))
        print "Après:" + str(os.stat(Simplified_BaitSequences).st_size)

        if nb_subseq > 0:
            Rep2Cluster_BaitSequences = "%s/BaitSequences.%d.s2.fasta" %(self.TmpDirName, self.CumulIteration)
            # Run cd-hit-est to remove redundancy
            start_fc = time.time()
            CdHitEstProcess = Aligner.CdHitEst(Simplified_BaitSequences,Rep2Cluster_BaitSequences, c=1)
            CdHitEstProcess.run()
            print((time.time() - start_fc))
            self.logger.info("reduce_final_redundancy_baitSequences clustering --- %s seconds ---" %(time.time() - start_fc))
            print "Après cdhit:" + str(ApytramNeeds.count_sequences (Rep2Cluster_BaitSequences))
            print "Après cdhit:" + str(os.stat(Rep2Cluster_BaitSequences).st_size)
            self.BaitSequences = Rep2Cluster_BaitSequences
        else:
            self.BaitSequences = Simplified_BaitSequences
        print self.BaitSequences
        print "#####"
        print "Reduce complexity:" + str(time.time() - start)
        self.logger.info("reduce complexity all--- %s seconds ---" %(time.time() - start))

    def measure_final_coverage(self):
        if not self.Aligned:
            self.align_query()
        # Use Mafft
        start = time.time()
        MafftProcess = Aligner.Mafft(self.AlignedQuery)
        MafftProcess.QuietOption = True
        MafftProcess.AutoOption = True
        #MafftProcess.AdjustdirectionOption = True
        MafftProcess.AddOption = self.FinalFastaFileName
        (self.MafftResult,err) = MafftProcess.get_output()

        self.logger.debug("Mafft --- %s seconds ---" %(time.time() - start))

        (StrictCoverage, LargeCoverage, self.DicPlotCov) = ApytramNeeds.calculate_coverage(self.MafftResult,self.ReferenceNames)

        self.logger.info("Strict Coverage: %s\tLarge Coverage: %s" %(StrictCoverage, LargeCoverage))
