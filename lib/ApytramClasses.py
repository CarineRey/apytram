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
from lib import Aligner
from lib import ApytramNeeds


#### Execution stats class

import pandas
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages



#### Sequence class     
        
class Sequence
	def __init__(self,ListLine):
		#TrinityExonerateProcess.Ryo = "%ti\t%qi\t%ql\t%tal\t%tl\t%tab\t%tae\t%s\t%pi\t%qab\t%qae\n"
		
		#t = target
		#q = query
		#i = id
		#l = length
		#a = aligned part
		#b = begin
		#e = end
		self.TrinityName = ListLine[0]
		self.ti = ListLine[1]
		self.ql =    float(ListLine[2])
		self.tal =   float(ListLine[3])
		self.tl =    float(ListLine[4])
		self.tae =   float(ListLine[5])
		self.score = float(ListLine[6])
		self.pi =    float(ListLine[7])
		self.qab =   float(ListLine[8])
		self.qae =   float(ListLine[9])
		
		self.BestSequence = ""


#### Execution statistics class
class Exec_stats:
    def __init__(self,start_time):
        New_Time_stat_dic = {"DatabaseBuilding": 0,
                             "Blast": 0,
                             "Blastdbcmd": 0,
                             "Trinity": 0,
                             "Exonerate_1":0,
                             "Exonerate_2":0,
                             "Mafft":0,
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
        self.FilteredTrinityFasta = ""
        
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
        start = time.time()
        self.logger.info("Database %s does not exist for the species: %s" % (self.DatabaseName, self.Species))
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
            self.add_time_statistic("DatabaseBuilding", start = start)
            self.logger.info("Database %s build in %s" %(self.DatabaseName,self.get_time_statistic("DatabaseBuilding")))
            
    def new_iteration(self):
        
        self.ExecutionStats.StartIterTime = time.time()
        
        # Cleaning
        
        # Previous Iteration
        self.PreviousReadNamesFile = self.ReadNamesFile
        self.PreviousFilteredTrinityFasta = self.FilteredTrinityFasta

        # New names files
        
        self.CurrentIteration +=1
        self.Improvment = True
        self.CompletedIteration = True
        
        self.ExecutionStats.new_iteration(self.CurrentIteration)
        
        self.ReadNamesFile = "%s/ReadNames.%d.txt" %(self.TmpDirName,self.CurrentIteration)
        self.ReadsNumber = 0
        
        self.ReadNamesFile_Right = "%s/ReadNames.%d.1.txt" % (self.TmpDirName,self.CurrentIteration)
        self.ReadNamesFile_Left  = "%s/ReadNames.%d.2.txt" % (self.TmpDirName,self.CurrentIteration)
        self.ReadFasta_Right     = "%s/Reads.%d.1.fasta"   % (self.TmpDirName,self.CurrentIteration)
        self.ReadFasta_Left      = "%s/Reads.%d.2.fasta"   % (self.TmpDirName,self.CurrentIteration)
        
        self.ReadFasta      = "%s/Reads.%d.fasta"   % (self.TmpDirName,self.CurrentIteration)
                    
        self.TrinityFasta = "%s/Trinity_iter_%d" %(self.TmpDirName,self.CurrentIteration)

        self.TrinityExonerateFile = "%s/Trinity_iter_%d.exonerate_cdna2g" %(self.TmpDirName,self.CurrentIteration)
        self.TrinityExonerateResult = ""
        
        self.ExonerateBetweenIterFile = "%s/iter_%d_%d.exonerate" %(self.TmpDirName,self.CurrentIteration -1, self.CurrentIteration)
    
    def launch_Blastn(self,BaitSequences,Threads):
        start = time.time()
        self.logger.info("Blast bait sequences on reads database")
        BlastnProcess = BlastPlus.Blast("blastn", self.DatabaseName, BaitSequences)
        BlastnProcess.Evalue = self.Evalue
        BlastnProcess.Task = "blastn"
        BlastnProcess.Threads = Threads
        BlastnProcess.OutFormat = "6 sacc"
        
        (out,err) = BlastnProcess.launch(self.ReadNamesFile)
        self.add_time_statistic("Blast", start = start)
        self.logger.debug("Blast --- %s seconds ---" %(self.get_time_statistic("Blast")))
        
    def get_read_sequences_by_blasdbcmd(self, Threads, Memory):
        start = time.time()
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
        
        self.add_time_statistic("Blastdbcmd", start = start)
        self.logger.debug("Blastdbcmd --- %s seconds ---" %(self.get_time_statistic("Blastdbcmd")))
        
    def launch_Trinity(self, Threads, Memory):
        start = time.time()
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
        
        self.add_time_statistic("Trinity", start = start)
        self.logger.debug("Trinity --- %s seconds ---" %(self.get_time_statistic("Trinity")))
           
    def compare_trinity_results_with_references(self,Query):
        # Use Exonerate
        start = time.time()
        TrinityExonerateProcess = Aligner.Exonerate(Query.RawQuery,self.TrinityFasta)
        # Keep only the best hit for each contig from Trinity
        TrinityExonerateProcess.Bestn = 1
        TrinityExonerateProcess.Model = "cdna2genome"
        # Customize the output format
        TrinityExonerateProcess.Ryo = self.TrinityExonerateRyo
        (out,err,self.TrinityExonerateResult) = TrinityExonerateProcess.get_output()
        # Write the result in a file
        if self.TrinityExonerateResult:
            ApytramNeeds.write_in_file(self.TrinityExonerateResult,self.TrinityExonerateFile)       
        else:
            logger.info("Reconstructed sequences but no homologous with references")
            logger.info("Try to get homologies with a more sensible model")
            ### Try to get homologies with a more sensible model
            self.TrinityExonerateFile = "%s/Trinity_iter_%d.exonerate_coding2g" %(self.TmpDirName,self.CurrentIteration)
            TrinityExonerateProcess = Aligner.Exonerate(Query.RawQuery,self.TrinityFasta)
            # Keep only the best hit for each contig from Trinity
            TrinityExonerateProcess.Bestn = 1
            TrinityExonerateProcess.Model = "coding2genome"
            # Customize the output format
            TrinityExonerateProcess.Ryo = self.TrinityExonerateRyo
            (out,err,self.TrinityExonerateResult) = TrinityExonerateProcess.get_output()
            # Write the result in a file
            ApytramNeeds.write_in_file(self.TrinityExonerateResult,self.TrinityExonerateFile)
        
        self.add_time_statistic("Exonerate_1", start = start)
        self.logger.debug("Exonerate_1 --- %s seconds ---" %(self.get_time_statistic("Exonerate_1")))
	

    def parse_and_filter_exonerate_results(self, final_iteration = False):
        self.FilteredTrinityFasta = "%s/Trinity_iter_%d.filtered.fasta" %(self.TmpDirName,self.CurrentIteration)
        if not final_iteration:
            # Keep only sequence with a identity percentage > MinIdentitypercentage on the whole hit
            (self.BestScoreNames, ReverseNames, self.TrinityExonerateResultsDict, StatsIter) = ApytramNeeds.parse_exonerate_results(self.TrinityExonerateResult,
                                                                                                                                              self.MinIdentityPercentage,
                                                                                                                                              minalilength = self.MinAliLength)
            self.logger.info("Filter sequence with a identity percentage superior to %d and a alignment len %d" %(self.MinIdentityPercentage, self.MinAliLength))                                                                                                                                 
        else:
            (self.BestScoreNames, ReverseNames, self.TrinityExonerateResultsDict, StatsIter) = ApytramNeeds.parse_exonerate_results(self.TrinityExonerateResult,
                                                                                                                                              self.FinalMinIdentityPercentage,
                                                                                                                                              minalilengthpercentage = self.FinalMinAliLength,
                                                                                                                                              minlengthpercentage = self.FinalMinLength)
            self.logger.info("Filter sequence with a identity percentage superior to %d and a percentage alignment len %d" %(self.FinalMinIdentityPercentage, self.FinalMinAliLength))
                #StatsDict[Reali].update(StatsIter)
        
        for (att,value) in StatsIter.items():
            self.add_iter_statistic(att,value)
             
        self.FilteredSequenceNames = self.TrinityExonerateResultsDict.keys()
        if self.FilteredSequenceNames: # If sequences pass the filter
            # Write filtered sequences in a file
            if not final_iteration:
                ExitCode = ApytramNeeds.filter_fasta(self.TrinityFasta, self.FilteredSequenceNames, self.FilteredTrinityFasta, ReverseNames = ReverseNames)
            else:
                ExitCode = ApytramNeeds.filter_fasta(self.PreviousFilteredTrinityFasta, self.FilteredSequenceNames, FilteredTrinityFasta, ReverseNames = ReverseNames)
    
    def compare_current_and_previous_iterations():
        self.logger.info("Refind the \"parent\" contig from the previous contig for each contig and check they are different")
        start = time.time()
        ExonerateProcess = Aligner.Exonerate(Species.FilteredTrinityFasta, Species.PreviousFilteredTrinityFasta)
        # Keep only the best hit for each contigs
        ExonerateProcess.Bestn = 1
        ExonerateProcess.Model =  "est2genome"
        # Customize the output format
        ExonerateProcess.Ryo = "%ti\t%qi\t%ql\t%qal\t%tal\t%tl\t%pi\n"
        (out,err,ExonerateResult) = ExonerateProcess.get_output()
        
        ApytramNeeds.write_in_file(ExonerateResult,self.ExonerateBetweenIterFile)

        AlmostIdenticalResults = ApytramNeeds.check_almost_identical_exonerate_results(ExonerateResult)
        
        if AlmostIdenticalResults:
            self.logger.info("Contigs are almost identical than the previous iteration (Same size (~98%), > 99% identity)")
            self.Improvment = False
            
        self.add_time_statistic("Exonerate_2", start = start)
        self.logger.debug("Exonerate_2 --- %s seconds ---" %(self.get_time_statistic("Exonerate_2")))     

    def measure_coverage(self,Query):
        # Use Mafft
        start = time.time()
        MafftProcess = Aligner.Mafft(Query.AlignedQuery)
        MafftProcess.QuietOption = True
        MafftProcess.AutoOption = True
        #MafftProcess.AdjustdirectionOption = True
        MafftProcess.AddOption = self.FilteredTrinityFasta
        (self.MafftResult,err) = MafftProcess.get_output()
        self.add_time_statistic("Mafft", start = start)
        
        self.logger.debug("Mafft --- %s seconds ---" %(self.get_time_statistic("Mafft"))) 
        
        (StrictCoverage, LargeCoverage, self.DicPlotCov) = ApytramNeeds.calculate_coverage(self.MafftResult)

        self.add_iter_statistic("StrictCoverage", StrictCoverage)
        self.add_iter_statistic("LargeCoverage", LargeCoverage)
        self.logger.info("Strict Coverage: %s\tLarge Coverage: %s" %(StrictCoverage, LargeCoverage))
                      
    def end_iteration(self):
        iter_time = time.time() - self.ExecutionStats.StartIterTime
        self.ExecutionStats.IterStatsDict[self.CurrentIteration]["IterationTime"] = iter_time
        self.ExecutionStats.IterStatsDict[self.CurrentIteration]["CumulTime"] += iter_time
        
        NoPythonTime = self.get_time_statistic("Blast") + \
                       self.get_time_statistic("Blastdbcmd") + \
                       self.get_time_statistic("Trinity") + \
                       self.get_time_statistic("Exonerate_1") + \
                       self.get_time_statistic("Exonerate_2")
        
        self.add_time_statistic("Python", inter = NoPythonTime )
        self.logger.debug("Python --- %s seconds ---" %(self.get_time_statistic("Python")))
            
    def write_fasta_apytram_outputs(self, Query, no_best_file = False, only_best_file = False):
        "Return fasta files with new names depending on the self.ExonerateResultsDict and options"
        self.logger.info("Write outputfiles for %s" %(self.Species))
        OutputFastaExtensions = []
        if not no_best_file:
            # Best sequences
            OutputFastaExtensions.append(".best.fasta")
        if not only_best_file:
            # Last iteration sequences
            OutputFastaExtensions.append(".fasta")
        
        if not OutputFastaExtensions:
            self.logger.error("No output files to write !")
        else:
            Header = self.TrinityExonerateRyo.replace('%',"").replace("\n","").split()
            df = pandas.DataFrame(self.TrinityExonerateResultsDict, index = Header)
            
            # read self.FilteredTrinity
            File = open(self.FilteredTrinityFasta,"r")
            Fasta = File.read().strip().split("\n")
            File.close()

            for OutFastaExtension in OutputFastaExtensions:
                ValidatedNames = []
                if (OutFastaExtension == ".fasta") :
                    ValidatedNames = self.FilteredSequenceNames
                    Message = self.Species + "_"
                else:
                    ValidatedNames = self.BestScoreNames.values()
                    Message = self.Species + "_best_"
                
                name = ""
                sequence = ""
                string_list = []
                i = 1
                for line in Fasta:
                    if re.match(">",line):
                        name = line.split()[0].replace(">","")
                        if name in ValidatedNames:
                            string_list.append(">APYTRAM_%s%d.len=%s.[%s]\n" %(Message,i,df[name]["ql"],df[name]["ti"]))
                            i+=1
                        else:
                            name = ""                   
                    elif name != "":
                        string_list.append(line + "\n")
                    else:
                        pass

                # Write sequences
                ApytramNeeds.write_in_file("".join(string_list),Query.OutPrefixName + "_" + self.Species + OutFastaExtension)
       
    def get_output_fasta(self, Query, fasta = "all"):
        assert fasta in ["all","best"], "fasta must be all or best"
        
        "Return fasta files with new names depending on the self.ExonerateResultsDict and options"
        self.logger.info("Prepare outputfiles for %s" %(self.Species))
        OutputFastaExtensions = []
        if fasta == "best":
            # Best sequences
            OutFastaExtension = ".best.fasta"
            Message = self.Species + "_best_"
            ValidatedNames = self.BestScoreNames.values()
            
        else:
            # Last iteration sequences
            OutFastaExtension = ".fasta"
            Message = self.Species + "_"
            ValidatedNames = self.FilteredSequenceNames
        
        Header = self.TrinityExonerateRyo.replace('%',"").replace("\n","").split()
        df = pandas.DataFrame(self.TrinityExonerateResultsDict, index = Header)
        
        # read self.FilteredTrinity
        File = open(self.FilteredTrinityFasta,"r")
        Fasta = File.read().strip().split("\n")
        File.close()
        
        name = ""
        sequence = ""
        string_list = []
        i = 1
        for line in Fasta:
            if re.match(">",line):
                name = line.split()[0].replace(">","")
                if name in ValidatedNames:
                    string_list.append(">APYTRAM_%s%d.len=%s.[%s]\n" %(Message,i,df[name]["ql"],df[name]["ti"]))
                    i+=1
                else:
                    name = ""                   
            elif name != "":
                string_list.append(line + "\n")
            else:
                pass

        # Write sequences
        return(string_list)
       
    def get_stats(self):
        df_time = pandas.DataFrame(self.ExecutionStats.TimeStatsDict)
        df_iter = pandas.DataFrame(self.ExecutionStats.IterStatsDict)
        df = pandas.concat([df_time,df_iter]).T
        df.insert(0, "Iteration", df.index)
        df.insert(0, "Species", self.Species)
        return([df])

    def parse_exonerate_results(self, final_iteration = False):
		"Return a list of Sequences if the identity percentage is superior to MinIdentityPercentage and the alignment length is superior to MinAliLen "
		
		if final_iteration:
			minidentypercentage =  self.FinalMinIdentityPercentage
			minalilength = 0
			minlengthpercentage = self.FinalMinLength
			minalilengthpercentag = self.FinalMinAliLength
		else:
			minidentypercentage = self.MinIdentityPercentage
			minalilength = self.MinAliLength
			minlengthpercentage = 0
			minalilengthpercentag = 0

		self.SequencesList = []
		BestScoreNames = {}
		
		List = self.ExonerateResult.strip().split("\n")
		
		for line in List:
			ListLine = line.split("\t")		
			(ti,qi) = ListLine[1]
			(ql,tal,tl,tab,tae,score,pi,qab,qae) = [float(x) for x in ListLine[2:]]
			
			if (pi >=  minidentypercentage) and (tal >= minalilength) and (ql >= minlengthpercentage*tl/100) and (tal >= minalilengthpercentage*tl/100) :
				# We keep this sequence
				# A same sequence can be present 2 time if the hit scores are identical.
				# We keep only the first
				ever_seen = False
				for Sequence in self.SequencesList:
					if Sequence.TrinityName == qi:
						ever_seen = True
						break
				if not ever_seen:
					new_sequence = Sequence(ListLine)					
					SequencesList.append(new_sequence)

					self.add_iter_statistic("TotalIdentity", pi, mode_add = True)
					if self.get_iter_statistic("BestIdentity") <= pi:
						self.add_iter_statistic("BestIdentity",pi)
						
					self.add_iter_statistic("TotalLength", ql, mode_add = True)
					if self.get_iter_statistic("BestLength") <= ql:
						self.add_iter_statistic("BestLength", ql)
						
					self.add_iter_statistic("TotalScore", score, mode_add = True)
					if self.get_iter_statistic("BestScore") <= score :
						self.add_iter_statistic("BestScore", score)
	
					IterStats["TotalLength"] += ql
					if IterStats["BestLength"] <= ql:
						IterStats["BestLength"] = ql
	
					IterStats["TotalScore"] += score
					if IterStats["BestScore"] <= score:
						IterStats["BestScore"] = score
						BestScoreNames[ti] = qi
						
					# Check if the seqeunce is reverse
					if ((qae-qab)*(tae-tab) < 0):
						new_sequence.reverse = True
		
		NbContigs = len(self.SequencesList)
		self.add_iter_statistic("NbContigs",NbContigs)
		
		if NbContigs:
			AverageIdentity = self.get_iter_statistic("TotalIdentity") / NbContigs 
			AverageLength = self.get_iter_statistic("TotalLength") / NbContigs
			AverageScore = self.get_iter_statistic("TotalScore") / NbContigs
			
			self.add_iter_statistic("AverageIdentity", pi)
			self.add_iter_statistic("AverageLength", pi)
			self.add_iter_statistic("AverageScore", pi)

		
		for (Target,Query) in BestScoreNames.items():
			for Sequence in SequencesList:
				if Query == Sequence.TrinityName:
					Sequence.BestSequence = Target
#   def new_query(self):
 
      


#### query class

class Query:
    def __init__(self,Name,QueryPath):
        self.Name = Name
        self.RawQuery = QueryPath
        self.SequenceNb = ApytramNeeds.count_sequences(QueryPath)
        

        self.AlignedQuery = ""
        
        
        self.CumulIteration = 0
        self.AbsIteration = 0
        
        self.Stop = False
        self.SpeciesWithoutImprovment = {0:[]}
        self.BaitSequences = ""
        self.PreviousBaitSequences = ""
        
        self.BestOutFileContent = []
        self.OutFileContent = []
        self.StatsFileContent = []
        
        
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
        SpeciesCurrentReconstructedSequencesFileList = [Species.FilteredTrinityFasta for Species in SpeciesList if Species.FilteredTrinityFasta]
        
        if self.AbsIteration == 1 :
            SpeciesCurrentReconstructedSequencesFileList.append(self.RawQuery)
            
        if SpeciesCurrentReconstructedSequencesFileList:
            ApytramNeeds.cat_fasta(" ".join(SpeciesCurrentReconstructedSequencesFileList), self.BaitSequences)

        
        



