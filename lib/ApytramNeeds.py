#!/usr/bin/python
# coding: utf-8

# File: ApytramNeeds.py
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


import os
import re
import shutil
import string
import numpy as np
import sys
import subprocess
import pandas
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
#matplotlib.style.use('ggplot')



# Class

#### Sequence class     
        
class Sequence:
    def __init__(self,):
        self.Name = ""
        self.OldName = ""
        
        self.Sequence = ""
        
        ### homology score with the best reference
        #t = target
        #q = query
        #i = id
        #l = length
        #a = aligned part
        #b = begin
        #e = end
        self.ti =    ""
        self.ql =    ""
        self.tal =   ""
        self.tl =    ""
        self.tae =   ""
        self.score = ""
        self.pi =    ""
        self.qab =   ""
        self.qae =   ""
        
        ### Caracteristics
        self.BestSequence = ""
        self.Reverse = False
    
    def __str__(self):
        return(">" + self.Name + "\n" + '\n'.join(self.Sequence[i:i+60] for i in range(0, len(self.Sequence), 60)))
           
    def copy(self):
        ns = Sequence()
        
        for (att,value) in self.__dict__.items():
            ns.__dict__[att] = value
        
        return(ns)
    
    def add_attribute_from_exonerate(self, ExonerateListLine):
        #TrinityExonerateProcess.Ryo = "%ti\t%qi\t%ql\t%tal\t%tl\t%tab\t%tae\t%s\t%pi\t%qab\t%qae\n"
        #"%ti\t%qi\t%ql\t%tal\t%tl\t%tab\t%tae\t%s\t%pi\t%qab\t%qae\n"
        self.ti     = ExonerateListLine[0]
        self.Name   = ExonerateListLine[1]
        self.ql     = float(ExonerateListLine[2])
        self.tal    = float(ExonerateListLine[3])
        self.tl     = float(ExonerateListLine[4])
        self.tab    = float(ExonerateListLine[5])
        self.tae    = float(ExonerateListLine[6])
        self.score  = float(ExonerateListLine[7])
        self.pi     = float(ExonerateListLine[8])
        self.qab    = float(ExonerateListLine[9])
        self.qae    = float(ExonerateListLine[10])
    
    def complete(self,newinfo_sequence):
        for (att,value) in self.__dict__.items():
            if not value and newinfo_sequence.__dict__[att]:
                self.__dict__[att] = newinfo_sequence.__dict__[att]
        


#### Fasta class

# import os,re
# f = Fasta()
# f.read_fasta(FastaFilename = "example/multi_ref.fasta")
# fi = f.filter_fasta("ENSMUSP00000002708_G01455")
# print fi
# d = {"ENSMUSP00000002708_G01455" : "totot"}
# fd = fi.rename_fasta(d)
# print fd
# fd.write_fasta("test")

class Fasta:
    def __init__(self):
        self.Sequences = []
        self.Names = []
        self.ReferenceSequences = []
    
    def __str__(self):
        string = []
        for s in self.Sequences:
            string.extend([str(s)])   
        return("".join(string))
        
    def append(self,new_sequence):
        assert isinstance(new_sequence, Sequence), "Sequence must belong to the Sequence class"
        if new_sequence.Reverse and new_sequence.Sequence:
            new_sequence.Sequence = reverse_complement(new_sequence.Sequence)
            new_sequence.Reverse = False
        
        self.Sequences.append(new_sequence)
        self.Names.append(new_sequence.Name)
        
    def read_fasta(self, FastaFilename = "" , String = ""):
        if String:
            Fasta = String.strip().split("\n")
        elif os.path.isfile(FastaFilename):
            File = open(FastaFilename,"r")
            Fasta = File.read().strip().split("\n")
            File.close()
        else:
			Fasta = []
        
        name = ""
        sequence = ""
        sequence_list = []
        
        for line in Fasta + [">"]:
            if re.match(">",line):
                # This is a new sequence write the previous sequence if it exists
                if sequence_list:
                    new_sequence = Sequence()
                    new_sequence.Name = name
                    new_sequence.Sequence = "".join(sequence_list)
                    self.append(new_sequence)
                    sequence_list = []
                    
                name = line.split()[0][1:] # remove the >
            
            elif name != "":
                sequence_list.append(line)
            else:
                pass
    
    def filter_fasta(self, SelectedNames):
        FilteredFasta = Fasta()
        for s in self.Sequences:
            if s.Name in SelectedNames:
                FilteredFasta.append(s)
        return FilteredFasta
    
    def complete_fasta(self,NewInfoFasta):
        assert isinstance(NewInfoFasta, Fasta), "NewInfoFasta must belong to the Fasta class"
        for i in range(len(self.Sequences)):
            if self.Sequences[i].Name in NewInfoFasta.Names:
                index = NewInfoFasta.Names.index(self.Sequences[i].Name)
                self.Sequences[i].complete(NewInfoFasta.Sequences[index])
                
        
    def rename_fasta(self, OldNamesNewNamesDict):
        assert isinstance(OldNamesNewNamesDict, dict), "OldNamesNewNamesDict must be a dict"
        RenamedFasta = Fasta()
        for s in self.Sequences:
            if s.Name in OldNamesNewNamesDict.keys():
                ns = s.copy()
                ns.Name = OldNamesNewNamesDict[s.Name]
                ns.OldName = s.Name
                RenamedFasta.append(ns)
        return RenamedFasta
    
    
    def write_fasta(self, OutFastaFile):
        # Write all sequences in the file
        write_in_file(str(self), OutFastaFile)


def end(exit_code, TmpDirName, keep_tmp = False, logger=""):
    "Functions to end apytram in removing temporary directory"
    ### Remove tempdir if the option --tmp have not been use
    if logger:
        logger.debug("Remove the temporary directory")
    #Remove the temporary directory :
    if not keep_tmp and "apytram" in TmpDirName:
        shutil.rmtree(TmpDirName)
    sys.exit(exit_code)
      
def fastq2fasta(FastqFile,FastaFile):
    ExitCode = 1
    command = """cat %s | awk 'NR%%4==1||NR%%4==2'  | tr "@" ">" > %s """ %(FastqFile, FastaFile)
    p = subprocess.Popen(command,
                          stdout = subprocess.PIPE,
                          stderr = subprocess.PIPE,
                          shell = True)
    out, err = p.communicate()
    return (out, err)

def cat_fasta(FastaFiles,CatFastaFile):
    ExitCode = 1
    command = """cat %s > %s""" %(FastaFiles, CatFastaFile)
    p = subprocess.Popen(command,
                          stdout = subprocess.PIPE,
                          stderr = subprocess.PIPE,
                          shell = True)
    (out, err) = p.communicate()
    return (out, err)
    
def reverse_complement(Sequence_str):
    intab = "ABCDGHMNRSTUVWXYabcdghmnrstuvwxy"
    outtab = "TVGHCDKNYSAABWXRtvghcdknysaabwxr"
    trantab = string.maketrans(intab, outtab)
    # Reverse
    Reverse = Sequence_str.replace("\n","")[::-1]
    # Complement 
    Complement = Reverse.translate(trantab)   
    Complement = '\n'.join(Complement[i:i+60] for i in range(0, len(Complement), 60))
    if not re.search("\n$",Complement):
        Complement += "\n"
    return Complement

def write_in_file(String,Filename,mode = "w"):
    if mode in ["w","a"]:
        File = open(Filename,mode)
        File.write(String)
        File.close()   
    
def add_paired_read_names(File, logger = ""):
    "Add paired read name to a read name list"
    if os.path.isfile(File):
        command1 = """awk '{ print $0; if (match($0,"1$")) sub("1$",2,$0); else if (match($0,"2$")) sub("2$",1,$0); print $0}' %s  | sort -u > %s """ %(File, File+".paired")
        command2 = "mv %s %s" %(File+".paired", File)
        
        p1 = subprocess.Popen(command1,
                          stdout = subprocess.PIPE,
                          stderr = subprocess.PIPE,
                          shell = True)
        (out1, err1) = p1.communicate()
        if not err1:
            p2 = subprocess.Popen(command2,
                          stdout = subprocess.PIPE,
                          stderr = subprocess.PIPE,
                          shell = True)
            (out2, err2) = p2.communicate()
            if err2 and logger:
                logger.error(err2)
        elif logger:
            logger.error(err1)
    elif logger:
        logger.error("%s is not a file" %(File))

def remove_duplicated_read_names(File,logger = ""):
    "remove duplicated read names"
    if os.path.isfile(File):
        command1 = """sort -u %s > %s""" %(File, File+".uniq")
        command2 = "mv %s %s" %(File+".uniq", File)
        p1 = subprocess.Popen(command1,
                              stdout = subprocess.PIPE,
                              stderr = subprocess.PIPE,
                              shell = True)
        (out1, err1) = p1.communicate()
        if not err1:
            p2 = subprocess.Popen(command2,
                              stdout = subprocess.PIPE,
                              stderr = subprocess.PIPE,
                              shell = True)
            (out2, err2) = p2.communicate()
            if err2 and logger:
                logger.error(err2)
        elif logger:
            logger.error(err1)

    elif logger:
        logger.error("File %s is not a valid path" %File)

def tmp_dir_clean_up(TmpDirName,i):
    if i == 1 :
        FilesToRemoves = ["%s/input_fastq.fasta" %TmpDirName,
                          "%s/input_fasta.fasta" %TmpDirName]
    else:
        FilesToRemoves = ["%s/ReadNames.%d.txt" %(TmpDirName,i),
                          "%s/ReadNames.%d.1.txt" %(TmpDirName,i),
                          "%s/ReadNames.%d.2.txt" %(TmpDirName,i),
                          "%s/Reads.%d.fasta" %(TmpDirName,i),
                          "%s/Reads.%d.1.fasta" %(TmpDirName,i),
                          "%s/Reads.%d.2.fasta" %(TmpDirName,i),
                          "%s/Trinity_iter_%d.exonerate_cdna2g" %(TmpDirName,i),
                          "%s/Trinity_iter_%d.exonerate_coding2g" % (TmpDirName, i),
                          "%s/Trinity_iter_%d.filtered.fasta" %(TmpDirName,i),
                          "%s/Trinity_iter_%d.Trinity.fasta" %(TmpDirName,i)]
    for f in FilesToRemoves:
        if os.path.isfile(f):
            os.remove(f)

def get_free_space(DirPath):
    if not DirPath:
        DirPath = "."
    free = 0
    if os.path.isdir(DirPath):
        st = os.statvfs(DirPath)
        free = (st.f_bavail * st.f_frsize)
    return free    

def get_size(Files):
    size = 0
    if not type(Files) == type([]):
        Files = [Files]
    for f in Files:
        if os.path.isfile(f):
            size += os.stat(f).st_size
    return size

def count_lines(Files):
    Number = 0
    if type(Files) != type([]):
        Files = [Files]
        
    for File in Files:
        if os.path.isfile(File):
            Exit = subprocess.check_output(["wc", "-l", File])  
            Number += int(Exit.split(" ")[0])
    return Number

def count_sequences(File):
    Number = 0
    if os.path.isfile(File):
        Exit = subprocess.check_output(["grep \"^>\" %s | wc -l" %File], shell=True) 
        Number = Exit.replace("\n","")
    return int(Number)

def split_readnames_in_right_left(File,Right,Left):
    if os.path.isfile(File):     
        command = """awk '{if (match($0,"1$")) print > "%s"; else print > "%s"}' %s""" %(Right,Left,File)
        p = subprocess.Popen(command,
                              stdout = subprocess.PIPE,
                              stderr = subprocess.PIPE,
                              shell = True)
        out, err = p.communicate()
        return (out, err)
    else:
        return("","File %s is not a valid path" %File)

def check_paired_data(FastaFile):
    BadReadName = ""
    if os.path.isfile(FastaFile):
        #Check First sequence name end with 1/ or 2/
        Exit = subprocess.check_output(["grep", "-e", "^>",FastaFile, "-m", "100"])
        Exit2 = subprocess.check_output(["tail %s -n 200 | grep -e \"^>\"" %FastaFile], shell=True)
        ReadNames = (Exit+Exit2).strip().split("\n")
        while (not BadReadName) and ReadNames:
            ReadName = ReadNames.pop()
            if not re.search("[12]$",ReadName):
                BadReadName = ReadName[1:]
    return BadReadName

def are_identical(File1,File2):
    Identical = False
    if os.path.isfile(File1) and  os.path.isfile(File2):
        diff = subprocess.call(["diff",File1,File2,"--brief"],
                                       stdout=open("/dev/null", "w"),
                                       stderr=open("/dev/null", "w"))
        if not diff:
            Identical = True
    return Identical


def check_almost_identical_exonerate_results(ExonerateResult):
    "Return True if all hit have a hit with 99% > id and a len = 98% len query "
    Result = False
    i = 0
    List = ExonerateResult.strip().split("\n")
    for line in List:
        ListLine = line.split("\t")
        #ExonerateProcess.Ryo = "%ti\t%qi\t%ql\t%qal\t%tal\t%tl\t%pi\t\n"
        pi = float(ListLine[6])
        ql = float(ListLine[2])
        tl = float(ListLine[5])
        pl = min( ql, tl) / max( ql, tl)
        if (pl >= 0.99) and (pi > 98) :
            i += 1

    if len(List) == i:
        Result = True
    return Result 

def write_stats(df_list, Output):
    df = pandas.concat(df_list)
    df.to_csv(Output, na_rep = "0.0",index=False)

def create_plot(TimeStatsDict, IterStatsDict, OutPrefixName):
    df_time = pandas.DataFrame(TimeStatsDict).T
    df_iter = pandas.DataFrame(IterStatsDict).T
  
    with PdfPages("%s.stats.pdf" %(OutPrefixName)) as pdf:
        ### Plot 1 ###
        df_iter.plot(subplots=True,
                     grid = True,
                     legend = "best",
                     figsize=(10, 20),
                     style = "-o")
        plt.xlabel('Iteration')
        pdf.savefig()
        plt.close()
        
        ### Plot 2 ###
        df_time.plot(kind ='bar',
                     grid = True,
                     figsize=(10, 10),
                     stacked=True,
                     legend = "best")
        plt.ylabel('seconds')
        plt.xlabel('Iteration')
        plt.title('Distribution of the execution time')
        pdf.savefig()
        plt.close()


def calculate_coverage(AlignmentString, RefSeq):
    assert AlignmentString, "AlignmentString must not be empty"
    assert isinstance(AlignmentString,str),  "AlignmentString must be str"
    assert isinstance(RefSeq,list),  "RefSeq must be list"
    # The first sequence must be the reference
    #Alignment reading
    Alignment = Fasta()
    Alignment.read_fasta(String = AlignmentString)
    
    #Build a reference consensus sequences
    
    Ref_Ali = Alignment.filter_fasta(RefSeq)
    ConsensusList = []
    AliLength = len(Ref_Ali.Sequences[0].Sequence)
    for i in range(AliLength):
        # get all nucl for reference for each position
        pos = {seq.Sequence[i] for seq in Ref_Ali.Sequences}
        if "-" in pos:
            pos.remove("-")
        pos = list(pos)
        if len(pos) == 1:
            ConsensusList.append(pos[0])
        elif len(pos) == 0:
            ConsensusList.append("-")
        else:
            ConsensusList.append("0")
    
    # Attribute a code for each position
    ## 0 Gap
    ## 1 All references identical
    ## 2 ambigous in reference
    ## 3 contig identical references
    ## 4 contig ok but ambigous reference
    ## 5 present in contig but not in reference
    ## 6 different in contig than in reference
    
    DicPlotCov = {}
    StrictRefCov = [0]*AliLength
    LargeRefCov = [0]*AliLength
    
    for Sequence in Alignment.Sequences:
        Ref = False
        Name = Sequence.Name
        seq = Sequence.Sequence
        DicPlotCov[Name] = {}

        if Name in RefSeq:
            Ref = True
            
        for i in range(AliLength):

            if Ref:
                if ConsensusList[i] == "0":
                    DicPlotCov[Name][i] = 2
                elif seq[i] == "-":
                    DicPlotCov[Name][i] = 0             
                else:
                    DicPlotCov[Name][i] = 1
            else:
                if seq[i] == "-":
                    DicPlotCov[Name][i] = 0
                elif ConsensusList[i] == "-":
                    DicPlotCov[Name][i] = 5
                    LargeRefCov[i] = 1
                elif ConsensusList[i] == seq[i]:
                    DicPlotCov[Name][i] = 3
                    StrictRefCov[i] = 1
                    LargeRefCov[i] = 1
                elif ConsensusList[i] == "0":
                    DicPlotCov[Name][i] = 4
                    StrictRefCov[i] = 1
                    LargeRefCov[i] = 1
                else:
                    DicPlotCov[Name][i] = 6
                    StrictRefCov[i] = 1
                    LargeRefCov[i] = 1
    
    # Get length of the consensus
    ConsensusLength = len([ p for p in ConsensusList if p != "-" ])

    # We count the number of base position in the alignment
    SumStrictRefCov = sum(StrictRefCov)
    SumLargeRefCov = sum(LargeRefCov)

    #Percentage of positon of the reference which is reprensented in contigs
    StrictCov=float(SumStrictRefCov)/float(ConsensusLength) *100
    #Percentage of positon in contigs on the number positon of the reference
    #It can be superior to 100 if the contigs are longer than the reference
    LargeCov=float(SumLargeRefCov)/float(ConsensusLength) *100
    return (StrictCov, LargeCov, DicPlotCov)

def create_plot_ali(DictPlotCov, OutPrefixName):
    ### Plot Ali ###
    "Create a png file containing a representation of an alignement. The first seqeunce must be the reference."
    df = pandas.DataFrame(DictPlotCov).T
    fig, ax = plt.subplots()
    ax.grid(False)
    cmap, norm = matplotlib.colors.from_levels_and_colors([0,0.5, 1.5, 2.5,3.5,4.5,5.5,6.5],
                                                          ["White","Blue","DarkViolet",'Darkgreen',"DarkViolet","Orange",'Darkred'])
    heatmap = ax.pcolor(df ,
                        edgecolors="white",  # put black lines between squares in heatmap
                        cmap=cmap,
                        norm=norm)
    # Format
    fig = plt.gcf()
    Width = 3 + float(df.shape[1])/12
    Heigth =  1+float(df.shape[0])/2

    # Maximum figure size (pixels)
    if Width > 32000:
        Width = 32000
    if Heigth > 32000:
        Heigth = 32000

    fig.set_size_inches(Width, Heigth)
    # put the major ticks at the middle of each cell
    plt.yticks(np.arange(len(df.index)) + 0.5, df.index, size = 12)
    ax.xaxis.tick_top()
    plt.xticks(np.arange(len(df.columns)) + 0.5, df.columns, rotation=90, size= 6)
    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.tick_params(bottom='off', top='off', left='off', right='off') 
    plt.xlabel('Position')
    fig.tight_layout()
    fig.savefig("%s.ali.png" % OutPrefixName)
    plt.close()
