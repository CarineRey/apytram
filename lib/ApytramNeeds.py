#!/usr/bin/python
# coding: utf-8

# File: ApytramNeeds.py
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


import os
import re
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

def fastq2fasta(FastqFile,FastaFile):
    ExitCode = 1
    command = """cat %s | awk 'NR%%4==1||NR%%4==2'  | tr "@" ">" > %s """ %(FastqFile, FastaFile)
    ExitCode = os.system(command)
    return ExitCode

def cat_fasta(FastaFiles,CatFastaFile):
    ExitCode = 1
    command = """cat %s > %s""" %(FastaFiles, CatFastaFile)
    ExitCode = os.system(command)
    return ExitCode
    
def reverse_complement(Sequence):
    intab = "ABCDGHMNRSTUVWXYabcdghmnrstuvwxy"
    outtab = "TVGHCDKNYSAABWXRtvghcdknysaabwxr"
    trantab = string.maketrans(intab, outtab)
    # Reverse
    Reverse = Sequence.replace("\n","")[::-1]
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
    
def add_paired_read_names(File):
    "Add paired read name to a read name list"
    if os.path.isfile(File):
        command1 = """awk '{ print $0; if (match($0,"1$")) sub("1$",2,$0); else if (match($0,"2$")) sub("2$",1,$0); print $0}' %s  | sort -u > %s """ %(File, File+".paired")
        command2 = "mv %s %s" %(File+".paired", File)
        os.system(command1)
        os.system(command2)
    return 0

def remove_duplicated_read_names(File):
    "remove duplicated read names"
    if os.path.isfile(File):
        command1 = """sort -u %s > %s""" %(File, File+".uniq")
        command2 = "mv %s %s" %(File+".uniq", File)
        os.system(command1)
        os.system(command2)
    return 0

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
        command1 = """awk '{if (match($0,"1$")) print > "%s"; else print > "%s"}' %s""" %(Right,Left,File)
        os.system(command1)

def check_paired_data(FastaFile):
    BadReadName = ""
    if os.path.isfile(FastaFile):
        #Check First sequence name end with 1/ or 2/
        Exit = subprocess.check_output(["grep", "-e", "^>",FastaFile, "-m", "100"])
        Exit2 = subprocess.check_output(["tail %s -n 200 | grep -e \"^>\"" %FastaFile], shell=True)
        ReadNames = (Exit+Exit2).strip().split("\n")
        while Res and ReadNames:
            ReadName = ReadNames.pop()
            if not re.search("[12]$",ReadName):
                BadReadName = ReadName
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

def parse_exonerate_results(ExonerateResult, MinIdentityPercentage,
                            minalilength = 0,
                            minlengthpercentage = 0,
                            minalilengthpercentage = 0):
    "Return Query names if the identity percentage is superior to MinIdentityPercentage and the alignment length is superior to MinAliLen "
        
    IterStats = {"AverageLength": 0,
                 "TotalLength": 0,
                 "BestLength":0,
                 "AverageScore": 0,
                 "TotalScore": 0,
                 "BestScore":0,
                 "AverageIdentity": 0,
                 "TotalIdentity": 0,
                 "BestIdentity":0,
                }
        
    BestScoreNames = {}
    ReverseNames = []
    ExonerateResultsDict = {}
    List = ExonerateResult.strip().split("\n")
    
    for line in List:
        ListLine = line.split("\t")
        #TrinityExonerateProcess.Ryo = "%ti\t%qi\t%ql\t%tal\t%tl\t%tab\t%tae\t%s\t%pi\t%qab\t%qae\n"
        ti = ListLine[0]
        qi = ListLine[1]
        ql = float(ListLine[2])
        tal = float(ListLine[3])
        tl = float(ListLine[4])
        tab = float(ListLine[5])
        tae = float(ListLine[6])
        qab = float(ListLine[9])
        qae = float(ListLine[10])

        score = float(ListLine[7])
        pi = float(ListLine[8])


        if (pi >=  MinIdentityPercentage) and (tal >= minalilength) and (ql >= minlengthpercentage*tl/100) and (tal >= minalilengthpercentage*tl/100) :
            # We keep this sequence
            if qi not in ExonerateResultsDict.keys():
            # A same sequence can be present 2 time if the hit scores are identical.
            # We keep only the first
                ExonerateResultsDict[qi] = ListLine
                IterStats["TotalIdentity"] += pi
                if IterStats["BestIdentity"] <= pi:
                    IterStats["BestIdentity"] = pi

                IterStats["TotalLength"] += ql
                if IterStats["BestLength"] <= ql:
                    IterStats["BestLength"] = ql

                IterStats["TotalScore"] += score
                if IterStats["BestScore"] <= score:
                    IterStats["BestScore"] = score
                    BestScoreNames[ti] = qi
                    
                # Check if the seqeunce is reverse
                if ((qae-qab)*(tae-tab) < 0):
                    ReverseNames.append(qi)
    
    NbContigs = len(ExonerateResultsDict.keys())
    if NbContigs:
        IterStats["AverageIdentity"] = IterStats["TotalIdentity"] / NbContigs 
        IterStats["AverageLength"] = IterStats["TotalLength"] / NbContigs
        IterStats["AverageScore"] = IterStats["TotalScore"] / NbContigs
    
    return BestScoreNames, ReverseNames, ExonerateResultsDict, IterStats

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
        
def filter_fasta(FastaFile, Names, OutFastaFile, ReverseNames = []):
    "Return a fasta file with only sequence in Names"
    File = open(FastaFile,"r")
    Fasta = File.read().strip().split("\n")
    File.close()
    name = ""
    sequence = ""
    sequence_list = []
    string_list = []
    for line in Fasta:
        if re.match(">",line):
            # This is a new sequence write the previous sequence if it exists
            if sequence_list:
                sequence = "".join(sequence_list)
                if name in ReverseNames:
                    sequence = reverse_complement(sequence)
                string_list.append(sequence)
                sequence_list = []
            name = line.split()[0].replace(">","")
            if name in Names:
                string_list.append( ">%s\n" % name)
            else:
                name = ""
        elif name != "":
            sequence_list.append(line + "\n")
        else:
            pass
    # Write the last sequence    
    if sequence_list:
        sequence = "".join(sequence_list)
        if name in ReverseNames:
            sequence = reverse_complement(sequence)
        string_list.append(sequence)
        sequence_list = []
                
    # Write all sequences in the file
    OutFile = open(OutFastaFile,"w")
    OutFile.write("".join(string_list))
    OutFile.close()        
    return 0

def write_apytram_output(FastaFile, ExonerateResultsDict, OutFastaFile, Header = None, Names = None, Message = ""):
    "Return a fasta file with new names depending on the ExonerateResultsDict"
    df = pandas.DataFrame(ExonerateResultsDict, index = Header)
    File = open(FastaFile,"r")
    Fasta = File.read().strip().split("\n")
    File.close()
    name = ""
    sequence = ""
    string_list = []
    i = 1

    for line in Fasta:
        if re.match(">",line):
            name = line.split()[0].replace(">","")
            if (Names):
                if name in Names:
                    string_list.append(">APYTRAM_%s%d.len=%s.[%s]\n" %(Message,i,df[name]["ql"],df[name]["ti"]))
                    i+=1
                else:
                    name = ""
            else:
                string_list.append(">APYTRAM_%s%d.len=%s.[%s]\n" %(Message,i,df[name]["ql"],df[name]["ti"]))
                i+=1
            
        elif name != "":
            string_list.append(line + "\n")
        else:
            pass
    # Write sequences
    OutFile = open(OutFastaFile,"w")
    OutFile.write("".join(string_list))
    OutFile.close()  
    return 0

def write_stats(StatsDict,OutPrefixName):
    df = pandas.DataFrame(StatsDict).T
    df.to_csv("%s.stats.csv" % OutPrefixName)

def create_plot(StatsDict, OutPrefixName):
    df = pandas.DataFrame(StatsDict).T
    # Categorize columns
    TimeColumns = []
    OtherColumns = []
    for col in df.columns:
        if "Time" in col and col not in ["CumulTime","IterationTime"]:
            TimeColumns.append(col)
        else:
            OtherColumns.append(col)
        
    with PdfPages("%s.stats.pdf" % OutPrefixName) as pdf:
        ### Plot 1 ###
        df[OtherColumns].plot(subplots=True,
                              grid = True,
                              legend = "best",
                              figsize=(10, 20),
                              style = "-o")
        plt.xlabel('Iteration')
        pdf.savefig()
        plt.close()
        
        ### Plot 2 ###
        df[TimeColumns].plot(kind ='bar',
                             grid = True,
                             figsize=(10, 10),
                             stacked=True,
                             legend = "best")
        plt.ylabel('seconds')
        plt.xlabel('Iteration')
        plt.title('Distribution of the execution time')
        pdf.savefig()
        plt.close()

def create_plot_ali(DictPlotCov, OutPrefixName):
    ### Plot Ali ###
    "Create a png file containing a representation of an alignement. The first seqeunce must be the reference."
    df = pandas.DataFrame(DictPlotCov).T
    fig, ax = plt.subplots()
    ax.grid(False)
    cmap, norm = matplotlib.colors.from_levels_and_colors([0,0.5, 1.5, 2.5,3.5,4.5],
                                                          ["White","Orange",'Darkgreen','Darkred',"Blue"])
    heatmap = ax.pcolor(df ,
                        edgecolors="white",  # put black lines between squares in heatmap
                        cmap=cmap,
                        norm=norm)
    # Format
    fig = plt.gcf()
    Width = 3+float(df.shape[1])/12
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

def calculate_coverage(Alignment, NbRefSeq = 1):
    # The first sequence must be the reference
    #Alignment reading
    Alignment= Alignment.split("\n")
    Dic = {}
    DicPlotCov = {}
    RefSequence = ""
    RefName = ""
    Sequence = ""
    Name = ""
    # Reference sequence reading
    while (len(Dic.keys()) != NbRefSeq):
        line = Alignment.pop(0).replace("\n","")
        if line == "":
            pass
        elif re.match("^>",line):
            if RefName != "":
                Dic[RefName] = RefSequence
            else:
                RefName = line.replace(">","")
        else:
            RefSequence+=line

    # We reput the last line which is the next sequence name in the list
    Alignment.insert(0,line)
    # Other contig reading
    for line in Alignment:
        line=line.replace("\n","")
        if line=="":
            pass
        elif re.match("^>",line):
            if Sequence != "":
                Dic[Name] = Sequence
                Sequence = ""
            Name=line.replace(">","")
        else:
            Sequence+=line

    if Sequence != "":
        Dic[Name]=Sequence

    # Coverage count
    RefLength = len(Dic[RefName])
    # Counters initilisation
    StrictRefCov = [0]*RefLength
    LargeRefCov = [0]*RefLength
    Ref = [0]*RefLength

    for Contig in Dic.keys():
        DicPlotCov[Contig] = {}
        if Contig == RefName:
            for i in range(0,RefLength):
                if RefSequence[i] != "-":
                    Ref[i]=1
                    DicPlotCov[Contig][i] = 4
        else:
            Sequence = Dic[Contig]
            for i in range(0,RefLength):
                if Sequence[i]!='-':
                    LargeRefCov[i] = 1
                    if RefSequence[i] != "-":
                        StrictRefCov[i] = 1
                        if RefSequence[i] == Sequence[i]:
                            DicPlotCov[Contig][i] = 2
                        else:
                            DicPlotCov[Contig][i] = 3
                    else:
                        DicPlotCov[Contig][i] = 1
                else:
                    DicPlotCov[Contig][i] = 0

    # We count the number of base position in the alignment
    SumStrictRefCov = sum(StrictRefCov)
    SumLargeRefCov = sum(LargeRefCov)
    SumRef = sum(Ref)

    #Percentage of positon of the reference which is reprensented in contigs
    StrictCov=float(SumStrictRefCov)/float(SumRef) *100
    #Percentage of positon in contigs on the number positon of the reference
    #It can be superior to 100 if the contigs are longer than the reference
    LargeCov=float(SumLargeRefCov)/float(SumRef) *100
    return StrictCov, LargeCov, DicPlotCov
