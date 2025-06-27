import os
import sys
import glob
import subprocess
import re
import gzip
import sortedcontainers
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from plotnine import *
import multiprocessing

def loadFile(filePath):
        openFile=open(filePath,"r")
        fileContents=[]
        for line in openFile:
                line=line.strip("\n").split("\t")
                fileContents.append(line)
        openFile.close()
        return fileContents


def simplifyDataVis(dataVisPath,simpleOutPath,distance):
    fileContents=loadFile(dataVisPath)
    fileDict={}

    for line in fileContents:
        if line[0] not in fileDict:
            fileDict[line[0]]=[line]
        else:
            fileDict[line[0]].append(line)

    newLines=[]

    for haplotig, SNPLines in fileDict.items():
        startPos=int(SNPLines[0][1])
        endPos=int(SNPLines[0][1])+1
        for line in SNPLines:
            if int(line[1])>endPos+distance:
                newLines.append([haplotig,str(startPos),str(endPos),line[2],line[3]])
                startPos=int(line[1])
                endPos=int(line[1])+1
            else:
                endPos=int(line[1])
        newLines.append([haplotig, str(startPos), str(endPos), line[2], line[3]])

    outFileText=""
    for line in newLines:
        outFileText+="\t".join(line)+"\n"

    outFile=open(simpleOutPath,"w")
    outFile.write(outFileText)
    outFile.close()


def generatePhasingVis(dataVisPath,outPath):
    #Figure out how supported (as a %) each base is for each position in the cluster

    outputPNG=outPath+"phasedVis.png"
    outputPDF=outPath+"phasedVis.pdf"
    outputSVG=outPath+"phasedVis.svg"

    tbl=pd.read_csv(dataVisPath,sep="\t",header=None)
    tbl.columns=["contigName","startPos","endPos","chr","yValue"]

    g=(ggplot(tbl,aes(y="contigName",yend="contigName",x="startPos",xend="endPos",color='contigName'))
       +geom_segment(size=1.5)
       +theme(legend_position="none")
       +theme(panel_grid_minor=element_blank(), strip_text=element_text(size=14))
       +facet_wrap("~chr",scales="free")
       +theme(axis_title_y=element_blank(),axis_text_y=element_blank(),axis_ticks_major_y=element_blank(),axis_title_x=element_text(size=15),axis_text_x=element_text(size=11))
       +xlab("Position (bp)"))

    ggsave(g,filename=outputSVG,width=18,height=10)
    ggsave(g,filename=outputPNG,width=18,height=10)
    ggsave(g,filename=outputPDF,width=18,height=10)

    return "Generated phased plots"


def giveMeFullData(clusters):
    clusterText=""
    sortedClusterLines=sortedcontainers.SortedList()
    clusterLines=[]
    for clusterName, cluster in clusters.items():
        for SNP in cluster:
            contig=SNP.split(":")[0]
            position=int(SNP.split(":")[1].split("=")[0])
            sortedClusterLines.add([position,clusterName,contig])
    i=1
    previouslySeen={sortedClusterLines[0][1]:i}
    for line in sortedClusterLines:
        if line[1] not in previouslySeen:
            i+=1
            previouslySeen[line[1]]=i
        clusterLines.append([line[1],line[0],line[2],previouslySeen[line[1]]])
    for line in clusterLines:
        clusterText+="\t".join([str(x) for x in line])+"\n"
    return clusterText

def loadConsensusCluster(file):
    clusters_consensus = dict()
    with open(file) as f:
        for line in f:
            line = line.strip("\n").split("\t")
            clusters_consensus.setdefault(line[0], set())
            snp = line[1]+":"+line[2]+"="+line[3]
            clusters_consensus[line[0]].add(snp)
    
    return clusters_consensus


def generateLongReadFastQFiles(haplotigReadNameFilePath,longReadFastQFilePath,outputPath): #Parallelize in a way that can save memory
    ##clusterStatsText=""
    haplotigReadDict={}
    allAllowedReads=set()
    haplotigReadNameFile=open(haplotigReadNameFilePath,"r")
    for line in haplotigReadNameFile:
        line=line.strip("\n").split("\t")
        haplotigName=line[0]
        readName="@"+line[1]
        haplotigReadDict.setdefault(readName,set())
        haplotigReadDict[readName].add(haplotigName)
        allAllowedReads.add(readName)

    haplotigReadNameFile.close()


    ##longReadMetaData={}

    fastQReadQueue = multiprocessing.Queue()

    fastQWriter_p = multiprocessing.Process(target=fastQWriter, args=((outputPath),(fastQReadQueue),))
    fastQWriter_p.daemon = True
    fastQWriter_p.start()

# 需要修改成读入多个fq，写入一个clusters 的fq
    readData=[]
    readDataBatch={}
    batchSize=2500
    try:
        longReadFastQFile=gzip.open(longReadFastQFilePath, "rt")
    except gzip.BadGzipFile:
        longReadFastQFile=open(longReadFastQFilePath,"r")

    i=0
    j=0
    for line in longReadFastQFile:
        line=line.strip("\n")
        if i%4==0:
            if readData!=[]:
                for haplotigName in haplotigNames:
                    readDataBatch.setdefault(haplotigName,[])
                    readDataBatch[haplotigName].append(readData)
                    j+=1
                if j>batchSize:
                    j=0
                    fastQReadQueue.put(readDataBatch)
                    readDataBatch={}
            line=line.split()
            readName=line[0]
            haplotigNames=[]
            readData=[]
            if readName in allAllowedReads:
                for haplotigName in haplotigReadDict[readName]:
                    haplotigNames.append(haplotigName)
                readData=[readName]
        else:
            if readName in allAllowedReads:
                ##if i%4==3:
                    ##longReadMetaData[readName]=getAvgFastQScore(line)
                readData.append(line)
        i+=1
    longReadFastQFile.close()

    if readDataBatch!={}:
        fastQReadQueue.put(readDataBatch)

    fastQReadQueue.put("end")

    fastQWriter_p.join()

    return "Successfully generated phased FastQ files"

def fastQWriter(outputPath,queue):
    while True:
        readDataBatch=queue.get()
        if readDataBatch == "end":
            return
        for haplotigName, haplotigReads in readDataBatch.items():
            haplotigFastQText=""
            for readData in haplotigReads:
                readName=readData[0]
                readString=readName+"\n"+"\n".join(readData[1:])+"\n"
                haplotigFastQText+=readString
            haplotigFile=gzip.open(os.path.join(outputPath,haplotigName+".fastq.gz"),'at')
            haplotigFile.write(haplotigFastQText)
            haplotigFile.close()


def Plot(consensusClusterPath, readClusterPath,outPath,longReadFastQFilesPath,phasedName):
    plotOutdir = os.path.join(outPath,"Plots")
    fastQOut = os.path.join(outPath,"Fastq")
    allPath = [plotOutdir, fastQOut]
    for path in allPath:
        os.makedirs(path, exist_ok=True)

    # generateLongReadFastQFiles(readClusterPath, longReadFastQFilePath, fastQOut)

    # phasedName = "phased"+"_"+str(min_ovl)+"_"+str(min_sim)+"_"+str(max_ID)+"_"+str(min_len)+"_"+str(min_overlap)+"_"+str(min_ovl_len)
    clustersConsensus = loadConsensusCluster(consensusClusterPath)
    visDataTextFull=giveMeFullData(clustersConsensus)
    visDataFilePath=os.path.join(plotOutdir,phasedName+"_phasedDataFull.tsv")
    visDataFile=open(visDataFilePath,"w")
    visDataFile.write(visDataTextFull)
    visDataFile.close()

    #Simplify datavis
    dataVisPath=os.path.join(plotOutdir,phasedName+"_phasedDataFull.tsv")
    simpleOutPath=os.path.join(plotOutdir,phasedName+"_phasedDataSimple.tsv")
    simplifyDataVis(dataVisPath,simpleOutPath,1000)

    ################
    #Generate plots#
    ################

    #Phased
    datavisPath=os.path.join(plotOutdir,phasedName+"_")
    generatePhasingVis(simpleOutPath,datavisPath)

    return None


def WHPPPlot(WHPPPath,vcfName,outPath):
    consensusClusterPath = os.path.join(WHPPPath,vcfName+".vcf")
    plotOutdir = os.path.join(outPath,"Plots")
    # fastQOut = os.path.join(PhaseResultFolder, "Phased","Fastq")
    allPath = [plotOutdir]
    for path in allPath:
        os.makedirs(path, exist_ok=True)

    # generateLongReadFastQFiles(readClusterPath, longReadFastQFilePath, fastQOut)

    clustersConsensus = loadConsensusCluster(consensusClusterPath)
    visDataTextFull=giveMeFullData(clustersConsensus)
    visDataFilePath=os.path.join(plotOutdir,vcfName+"_phasedDataFull.tsv")
    visDataFile=open(visDataFilePath,"w")
    visDataFile.write(visDataTextFull)
    visDataFile.close()

    #Simplify datavis
    dataVisPath=os.path.join(plotOutdir,vcfName+"_phasedDataFull.tsv")
    simpleOutPath=os.path.join(plotOutdir,vcfName+"_phasedDataSimple.tsv")
    simplifyDataVis(dataVisPath,simpleOutPath,1000)

    ################
    #Generate plots#
    ################

    #Phased
    datavisPath=os.path.join(plotOutdir,vcfName+"_")
    generatePhasingVis(simpleOutPath,datavisPath)

    return None

if __name__ == "__main__":
    consensusClusterPath = sys.argv[1]
    readClusterPath = sys.argv[2]
    outdir = sys.argv[3]
    longReadFastQFilesPath = sys.argv[4]
    out_name = sys.argv[5]
    Plot(consensusClusterPath, readClusterPath,outdir,longReadFastQFilesPath,out_name)
