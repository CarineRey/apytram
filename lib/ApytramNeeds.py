import os
import re
import numpy as np
import sys
import subprocess
import pandas
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
#matplotlib.style.use('ggplot')

def fastq2fasta(FastqFile,FastaFile):
    ExitCode = 1
    command = """awk 'NR%%4==1||NR%%4==2' < %s | tr "@" ">" > %s """ %(FastqFile, FastaFile)
    os.system(command)
    return ExitCode

def add_paired_read_names(File):
    "Add paired read name to a read name list"
    if os.path.isfile(File):
        command1 = """awk '{ print $0; if (match($0,"1$")) sub("1$",2,$0); else if (match($0,"2$")) sub("2$",1,$0); print $0}' %s  | sort -u > %s """ %(File, File+".paired")
        command2 = "mv %s %s" %(File+".paired", File)
        os.system(command1)
        os.system(command2)
    return 0

def are_identical(File1,File2):
    Identical = False
    if os.path.isfile(File1) and  os.path.isfile(File2):
        diff = subprocess.call(["diff",File1,File2,"--brief"],
                                       stdout=open("/dev/null", "w"),
                                       stderr=open("/dev/null", "w"))
        if not diff:
            Identical = True
    return Identical

def parse_exonerate_results(ExonerateResult, MinIdentityPercentage):
    "Return Query names if the identity percentage is superior to MinIdentityPercentage"
        
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
        
    BestScoreNames = [""]
    ExonerateResultsDict = {}
    List = ExonerateResult.strip().split("\n")
    
    for line in List:
        ListLine = line.split("\t")
        #ExonerateProcess.Ryo = "%ti\t%qi\t%ql\t%tal\t%tl\t%tab\t%tae\t%s\t%pi\n"
        ti = ListLine[1]
        pi = float(ListLine[8])
        score = float(ListLine[7])
        tl = float(ListLine[4])
        if pi >=  MinIdentityPercentage:
            # We keep this sequence
            ExonerateResultsDict[ti] = ListLine

            IterStats["TotalIdentity"] += pi
            if IterStats["BestIdentity"] <= pi:
                IterStats["BestIdentity"] = pi
            
            IterStats["TotalLength"] += tl
            if IterStats["BestLength"] <= tl:
                IterStats["BestLength"] = tl
            
            IterStats["TotalScore"] += score
            if IterStats["BestScore"] <= score:
                IterStats["BestScore"] = score
                BestScoreNames[0] = ti
    
    NbContigs = len(ExonerateResultsDict.keys())
    IterStats["AverageIdentity"] = IterStats["TotalIdentity"] / NbContigs 
    IterStats["AverageLength"] = IterStats["TotalLength"] / NbContigs
    IterStats["AverageScore"] = IterStats["TotalScore"] / NbContigs
    
    return BestScoreNames, ExonerateResultsDict, IterStats

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
        
def filter_fasta(FastaFile, Names, OutFastaFile):
    "Return a fasta file with only sequence in Names"
    File = open(FastaFile,"r")
    Fasta = File.read().strip().split("\n")
    File.close()
    name = ""
    sequence = ""
    string = ""
    for line in Fasta:
        if re.match(">",line):
            name = line.split()[0].replace(">","")
            if name in Names:
                string += ">%s\n" % name
            else:
                name = ""
        elif name != "":
            string += line + "\n"
        else:
            pass
    # Write sequences
    OutFile = open(OutFastaFile,"w")
    OutFile.write(string)
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
    string = ""
    i = 1

    for line in Fasta:
        if re.match(">",line):
            name = line.split()[0].replace(">","")
            if (Names):
                if name in Names:
                    string += ">APYTRAM_%s%d[%s]\n" %(Message,i,df[name]["ti"])
                    i+=1
                else:
                    name = ""
            else:
                string += ">APYTRAM_%s%d[%s]\n" %(Message,i,df[name]["ti"])
                i+=1
            
        elif name != "":
            string += line + "\n"
        else:
            pass
    # Write sequences
    OutFile = open(OutFastaFile,"w")
    OutFile.write(string)
    OutFile.close()  
    return 0
def write_stats(StatsDict,OutPreffixName):
    df = pandas.DataFrame(StatsDict).T
    df.to_csv("%s.stats.csv" % OutPreffixName)

def create_plot(StatsDict,OutPreffixName):
    df = pandas.DataFrame(StatsDict).T
    with PdfPages("%s.stats.pdf" % OutPreffixName) as pdf:
        df.plot(subplots=True, grid = True,
                legend = "best",
                figsize=(10, 30), style = "-o")
        pdf.savefig()
        plt.close()
    
def calculate_coverage(Alignment):
    # The first sequence must be the reference
    #Alignment reading
    Alignment= Alignment.split("\n")
    dic = {}
    ref_sequence = ""
    ref_name = ""
    sequence = ""
    name = ""
    # Reference sequence reading
    while (len(dic.keys()) == 0):
        line = Alignment.pop(0).replace("\n","")
        if line == "":
            pass
        elif re.match("^>",line):
            if ref_name != "":
                dic[ref_name]=ref_sequence
            else:
                ref_name = line.replace(">","")
        else:
            ref_sequence+=line

    # We reput the last line which is the next sequence name in the list
    Alignment.insert(0,line)

    # Other contig reading
    for line in Alignment:
        line=line.replace("\n","")
        if line=="":
            pass
        elif re.match("^>",line):
            if sequence != "":
                dic[name]=sequence
                sequence = ""
            name=line.replace(">","")
        else:
            sequence+=line

    if sequence != "":
        dic[name]=sequence

    # Coverage count
    length_reference = len(dic[ref_name])
    # Counters initilisationwrite_stats
    cov_ref = [0]*length_reference
    cov_ref_ext = [0]*length_reference
    ref = [0]*length_reference

    for contig in dic.keys():
        if contig == ref_name:
            for i in range(0,length_reference):
                if ref_sequence[i] != "-":
                    ref[i]=1
        else:
            sequence = dic[contig]
            for i in range(0,length_reference):
                if sequence[i]!='-':
                    cov_ref_ext[i] =1
                    if ref_sequence[i] != "-":
                        cov_ref[i] = 1

    # We count the number of base position in the alignment
    sum_cov_ref_ext = 0
    sum_cov_ref = 0
    sum_ref = 0
    for i in range(0,length_reference):
        sum_cov_ref_ext += cov_ref_ext[i]
        sum_cov_ref += cov_ref[i]
        sum_ref += ref[i]

    #Percentage of positon of the reference which is reprensented in contigs
    p_cov=float(sum_cov_ref)/float(sum_ref) *100
    #Percentage of positon in contigs on the number positon of the reference
    #It can be superior to 100 if the contigs are longer than the reference
    p_cov_ext=float(sum_cov_ref_ext)/float(sum_ref) *100

    return p_cov, p_cov_ext