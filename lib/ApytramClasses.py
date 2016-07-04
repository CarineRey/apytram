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


from lib import BlastPlus
from lib import Trinity
from lib import ApytramNeeds


#### Execution stats class

import pandas
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

class Exec_stats:
    def __init__(self,start_time):
        New_stat_dic = {"IterationTime": 0,
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
        "BlastTime": 0,
        "BlastdbcmdTime": 0,
        "TrinityTime": 0,
        "Exonerate1Time":0,
        "Exonerate2Time":0,
        "MafftTime":0,
        "PythonTime":0
        }
        self.StatsDict =  {0:New_stat_dic}
        self.StatsDict[0]["CumulTime"] = time.time() - start_time
        self.StatsDict[0]["PythonTime"] = time.time() - start_time


#### RNA sample class

class RNA_species:
    def __init__(self,start_time,logger):
        self.logger = logger
        self.Species = ""
        self.OutputPrefix = ""
        self.TmpDirName = ""

        self.Fasta = []
        self.Fastq = []
        self.InputFasta = ""
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

        # Intermediary files

        self.ReadNamesFile = ""
        self.TrinityFasta = ""
        self.FileteredTrinityFasta = ""

    def has_a_formated_database(self):
        CheckDatabase_BlastdbcmdProcess = BlastPlus.Blastdbcmd(self.DatabaseName,"","")
        return CheckDatabase_BlastdbcmdProcess.is_database()

    def set_TmpDir(self,TmpDirName):
        self.TmpDirName = TmpDirName + "/" + self.Species
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

    def build_database(self,FreeSpaceTmpDir,TmpDirName):
        logger.info("Database %s does not exist for the species: %s" % (self.DatabaseName, self.Species))
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
                    end(1,TmpDirName)
                # Format the fastq file in fasta
            self.InputFasta = "%s/input_fastq.fasta" %(self.TmpDirName)
            self.logger.info("Convert the fastq file in fasta format")
            start_convert = time.time()
            input_size = ApytramNeeds.get_size(self.Fastq)
            self.logger.debug("size of input files: %s " %input_size)
            if input_size*1.1 > FreeSpaceTmpDir:
                self.logger.error("Not enough available free space in %s to convert input files" %TmpDirName)
                end(1,TmpDirName)
            if input_size*2.5 > FreeSpaceDBDir:
                self.logger.error("Not enough available free space in %s to build the database" %DatabaseDirName)
                end(1,TmpDirName)
            (out,err) = ApytramNeeds.fastq2fasta(" ".join(self.Fastq),self.InputFasta)
            if err:
                self.logger.error(err)
                end(1,TmpDirName)
            self.logger.info("Convertion takes %s seconds" %(time.time() - start_convert))
        elif self.Fasta:
            for fasta in self.Fasta:
                if not os.path.isfile(fasta):
                    self.logger.error("The fasta file (%s) does not exist." %fasta)
                    end(1,TmpDirName)
            # Check enough available free space
            input_size = ApytramNeeds.get_size(self.Fasta)
            self.logger.debug("size of input files: %s " %input_size)
            if input_size*2.5 > FreeSpaceDBDir:
                self.logger.error("Not enough available free space in %s to build the database" %DatabaseDirName)
                end(1,TmpDirName)
            # Concatenate fasta files
            if len(self.Fasta) > 1:
                if input_size*1.1 > FreeSpaceTmpDir:
                    self.logger.error("Not enough available free space in %s to concatenate input files" %TmpDirName)
                    end(1,TmpDirName)
                self.InputFasta = "%s/input_fasta.fasta" %(self.TmpDirName)
                self.logger.info("Concatenate fasta files")
                start_convert = time.time()
                (out,err) = ApytramNeeds.cat_fasta(" ".join(args.fasta),InputFasta)
                if err:
                    self.logger.error(err)
                    end(1,TmpDirName)
                self.logger.info("Concatenation takes %s seconds" %(time.time() - start_convert))
            else:
                self.InputFasta = self.Fasta[0]

        #Check if the end of sequence name of paired data are 1 or 2
        if self.PairedData:
            BadReadName = ApytramNeeds.check_paired_data(self.InputFasta)
            if BadReadName:
                self.logger.error("Paired read names must finished by 1 or 2. %s is uncorrect" %(BadReadName))
                end(1,TmpDirName)
        #Build blast formated database from a fasta file
        if not os.path.isfile(self.InputFasta):
            self.logger.error("Error during concatenation or conversion of input files. %s is not a file" %(self.InputFasta))
            end(1,TmpDirName)

        # Database building
        self.logger.info(self.DatabaseName + " database building")
        MakeblastdbProcess = BlastPlus.Makeblastdb(self.InputFasta,self.DatabaseName)
        (out,err) = MakeblastdbProcess.launch()

        self.FormatedDatabase = self.has_a_formated_database()

        if not self.FormatedDatabase:
            self.logger.error("Problem in the database building.\nAre you sure of your input format?\nAre all read names unique?")
            self.logger.info("Database %s does not exist" % self.DatabaseName)
            end(1,TmpDirName)
        else:
            self.logger.info("Database %s exists" % self.DatabaseName)

    def new_iteration(self,iter_time):
        # Cleaning
        
        # Previous Iteration
        self.PreviousReadNamesFile = self.ReadNamesFile
        self.PreviousFileteredTrinityFasta = self.FileteredTrinityFasta

        # New names files
        self.CurrentIteration +=1
        self.Improvment = True
        self.CompletedIteration = True

        self.ReadNamesFile = "%s/ReadNames.%d.txt" %(self.TmpDirName,self.CurrentIteration)
        self.ReadsNumber = 0
        
        self.ReadNamesFile_Right = "%s/ReadNames.%d.1.txt" % (self.TmpDirName,self.CurrentIteration)
        self.ReadNamesFile_Left  = "%s/ReadNames.%d.2.txt" % (self.TmpDirName,self.CurrentIteration)
        self.ReadFasta_Right     = "%s/Reads.%d.1.fasta"   % (self.TmpDirName,self.CurrentIteration)
        self.ReadFasta_Left      = "%s/Reads.%d.2.fasta"   % (self.TmpDirName,self.CurrentIteration)
        
        self.ReadFasta      = "%s/Reads.%d.fasta"   % (self.TmpDirName,self.CurrentIteration)
                    
        self.TrinityFasta = "%s/Trinity_iter_%d" %(self.TmpDirName,self.CurrentIteration)
        self.FileteredTrinityFasta = "%s/Trinity_iter_%d.filtered.fasta" %(self.TmpDirName,self.CurrentIteration)

        self.TrinityExonerate = "%s/Trinity_iter_%d.exonerate_cdna2g" %(self.TmpDirName,self.CurrentIteration)
    
    def launch_Blastn(self,BaitSequences,Threads):
            BlastnProcess = BlastPlus.Blast("blastn", self.DatabaseName, BaitSequences)
            BlastnProcess.Evalue = self.Evalue
            BlastnProcess.Task = "blastn"
            BlastnProcess.Threads = Threads
            BlastnProcess.OutFormat = "6 sacc"
            
            (out,err) = BlastnProcess.launch(self.ReadNamesFile)

    def get_read_sequences_by_blasdbcmd(self, Threads, Memory):
        if self.DatabaseType in ["RF","FR"]:
            self.logger.info("Split read names depending on 1/ or 2/")
            (out, err) = ApytramNeeds.split_readnames_in_right_left(self.ReadNamesFile,self.ReadNamesFile_Right,self.ReadNamesFile_Left)
            if err:
                self.logger.error(err)
            StrandList = [".1",".2"]
        else:
            StrandList = [""]
        
        self.logger.info("Retrieve reads sequences")
        start_blastdbcmd_time = time.time()
        for strand in StrandList:
            ReadFasta = "%s/Reads.%d%s.fasta" %(self.TmpDirName,self.CurrentIteration,strand)
            ReadNamesFile = "%s/ReadNames.%d%s.txt" % (self.TmpDirName,self.CurrentIteration,strand)
            BlastdbcmdProcess = BlastPlus.Blastdbcmd(self.DatabaseName, ReadNamesFile, ReadFasta)
            if not os.path.isfile(ReadFasta):
                (out,err) = BlastdbcmdProcess.launch()
            else:
                self.logger.warn("%s has already been created, it will be used" %(ReadFasta) )
                
        #StatsDict[i]["BlastdbcmdTime"] += time.time() - start_blastdbcmd_time
        #self.logger.debug("blastdbcmd --- %s seconds ---" %(time.time() - start_blastdbcmd_time))
        
    def launch_Trinity(self,Threads, Memory):
        start_trinity_time = time.time()
        self.logger.info("Launch Trinity")
        ExitCode = 0
        if self.StrandedData:
            if self.DatabaseType in ["RF","FR"]:
                TrinityProcess = Trinity.Trinity(self.TrinityFasta, right = self.ReadFasta_Right,
                                                left = self.ReadFasta_Left)
            else:
                TrinityProcess = Trinity.Trinity(self.TrinityFasta, single = self.ReadFasta)
            
            TrinityProcess.SS_lib_type = self.DatabaseType
        else:
            TrinityProcess = Trinity.Trinity(self.TrinityFasta, single = self.ReadFasta)
            if self.PairedData:
                TrinityProcess.RunAsPaired = True
        # If there is a huge number of reads, remove duplicated reads
        if self.ReadsNumber > 1000:
            TrinityProcess.NormalizeReads = True
        
        TrinityProcess.CPU = Threads
        TrinityProcess.max_memory = Memory
        # Keep only contig with a length superior to MinLength
        TrinityProcess.MinLength = self.MinLength
        
        # Use the  --full_cleanup Trinity option to keep only the contig file
        TrinityProcess.FullCleanup = True
        if not os.path.isfile(self.TrinityFasta+".Trinity.fasta"):
            (out,err,ExitCode) = TrinityProcess.launch()
        else:
            self.logger.warn("%s has already been created, it will be used" %(self.TrinityFasta+".Trinity.fasta") )
        
        self.TrinityFasta = self.TrinityFasta + ".Trinity.fasta"
        
        if not os.path.isfile(self.TrinityFasta):
            if ExitCode == 2 or ExitCode == 0 : # Trinity exit 0 if "No butterfly assemblies to report"
               self.logger.debug("Trinity found nothing...\n[...]\n"+"\n".join(out.strip().split("\n")[-15:]))
               self.logger.warning("Trinity has assembled no contigs at the end of the iteration %s (ExitCode: %d)" %(i,ExitCode) )
            elif ExitCode != 0:
               self.logger.debug("Trinity found nothing...\n[...]\n"+"\n".join(out.strip().split("\n")[-15:]))
               self.logger.error("Trinity has crashed (ExitCode: %d). Are all dependencies satisfied?" %(ExitCode))
        
        #StatsDict[i]["TrinityTime"] = time.time() - start_trinity_time
        #self.logger.debug("trinity --- %s seconds ---" %(StatsDict[i]["TrinityTime"]))
    




#### query class

class Query:
    def __init__(self,Name,QueryPath):
        self.Name = Name
        self.RawQuery = QueryPath
        self.SequenceNb = ApytramNeeds.count_sequences(QueryPath)
        

        self.AlignedQuery = ""
        
        
        self.CumulIteration = 0
        self.AbsIteration = 0
        
        self.Improvment = True
        self.SpeciesWithoutImprovment = {0:[]}
        self.BaitSequences = ""
        self.PreviousBaitSequences = ""
        
    def continue_iter(self):
        NbSpeciesWithoutImprovment = len(self.SpeciesWithoutImprovment[self.AbsIteration])
        if NbSpeciesWithoutImprovment == self.NbSpecies:
            return(False)
        else:
            return(True)
    
    def add_BaitSequences(self, new_bait_sequences):
        self.PreviousBaitSequences = self.BaitSequences
        self.BaitSequences = "%s/BaitSequences.%d.fasta" %(self.TmpDirName, self.CumulIteration)
        ApytramNeeds.cat_fasta("%s %s" %(new_bait_sequences, self.PreviousBaitSequences), self.BaitSequences)
        
        



