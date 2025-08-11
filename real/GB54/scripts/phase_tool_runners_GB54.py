import sys
import os
import subprocess
import timeit
import psutil
import time
import glob

def updateLog(logFilePath,logText):
    logFile=open(logFilePath,"a")
    logFile.write(logText)
    logFile.close()

def loadBenchmarkParameterFile(parameterFilePath):
    benchmarkParameters = {"mainPath": "",
                           "longReads": [],
                           "shortReads": [],
                           "reference": {},
                           "genomeSize": 0,
                           "threads": 0,
                           "max_ID": 0,
                           "min_len": 0,
                           "min_ovl": 0,
                           "min_sim": 0,
                           "min_overlap": 0,
                           "min_ovl_len": 0}
    benchmarkParameterFile=open(parameterFilePath,"r")
    for line in benchmarkParameterFile:
        line=line.strip("\n").split("\t")
        if line[0] in {"mainPath","genomeSize","threads","max_ID","min_len","min_ovl","min_sim","min_overlap","min_ovl_len"}:
            benchmarkParameters[line[0]]=line[1]
        elif line[0] in {"longReads","shortReads"}:
            for dataPoint in line[1:]:
                benchmarkParameters[line[0]].append(dataPoint)
        elif line[0] in {"reference"}:
            for dataPoint in line[1:]:
                dataPoint=dataPoint.split(":")
                benchmarkParameters[line[0]][dataPoint[0]]=dataPoint[1]
        else:
            print("Ignoring '","\t".join(line),"' because it is not recognized",sep="")
    benchmarkParameterFile.close()
    return benchmarkParameters

def launchFunction(functionName,parameters):
    if functionName=="ntchap":
        mappedLR=parameters[0]
        vcfFile=parameters[1]
        PhasingPath=parameters[2]
        referenceFilePath=parameters[3]
        threads=parameters[4]
        testName=parameters[5]
        outputLog, systemMessage=launchnTChap(mappedLR, vcfFile, PhasingPath, referenceFilePath, threads, testName)
        print(systemMessage)
        updateLog(fullLogPath, outputLog)
    elif functionName=="plot":
        consensusClusterPath = parameters[0]
        readClusterPath = parameters[1]
        outdir = parameters[2]
        longReadFastQFilesPath = parameters[3]
        outName = parameters[4]
        outputLog, systemMessage=launchPlot(consensusClusterPath, readClusterPath,outdir, longReadFastQFilesPath,outName)
        print(systemMessage)
        updateLog(fullLogPath, outputLog)
    elif functionName=="Phasing_clean_chr":
        PhaseResultFolder=parameters[0]
        readSNPsPath=parameters[1]
        max_ID = parameters[2]
        min_len = parameters[3] 
        min_ovl = parameters[4] 
        min_sim = parameters[5]
        min_overlap = parameters[6]
        min_ovl_len = parameters[7]
        percentKept = parameters[8]
        deduplicate = parameters[9]
        outputLog, systemMessage=launchPhasingCleaningChr(PhaseResultFolder, readSNPsPath, max_ID, min_len, min_ovl, min_sim, min_overlap, min_ovl_len,percentKept,deduplicate)
        print(systemMessage)
        updateLog(fullLogPath, outputLog)

def launchPlot(consensusClusterPath, readClusterPath,outdir, longReadFastQFilesPath,outName):
    start = timeit.default_timer()
    outputLog = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()) + "\t"

    command_py = "/home/gaoyun/poly/code/phasing/final_nTChap/process_data/real/GB54/scripts/plot.py"
    p = subprocess.Popen(["python3", command_py, consensusClusterPath, readClusterPath,outdir, longReadFastQFilesPath,outName], stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    outputLog += "COMMAND: " + " ".join(["python3", command_py, consensusClusterPath, readClusterPath,outdir, longReadFastQFilesPath,outName]) + "\n\n"

    maxRSS_mem = 0
    maxVMS_mem = 0

    try:
        poll = p.poll()
    except:
        poll = "Done"
    while poll is None:
        monitorP = psutil.Process(p.pid)
        children = list(monitorP.children(recursive=True))
        children = children + [monitorP]
        RSS_mem = 0
        VMS_mem = 0
        for child in children:
            mem_info = child.memory_info()
            RSS_mem += mem_info[0]
            VMS_mem += mem_info[1]
            if RSS_mem > maxRSS_mem:
                maxRSS_mem = RSS_mem
            if VMS_mem > maxVMS_mem:
                maxVMS_mem = VMS_mem
        time.sleep(1)
        try:
            poll = p.poll()
        except:
            poll = "Done"

    if p.stderr != "" or p.stdout != "":
        outputLog += "STDERR:\n\n" + p.stderr.read() + "\n\nSTDOUT:\n\n" + p.stdout.read() + "\n\n"

    print("Max memory usage, RSS:", maxRSS_mem / 1024 / 1024 / 1024, "VMS:", maxVMS_mem / 1024 / 1024 / 1024)

    stop = timeit.default_timer()
    totalRunTime = stop - start
    print("This took", totalRunTime, "seconds to run.")
    command = " ".join(["python3", command_py, consensusClusterPath, readClusterPath,outdir, longReadFastQFilesPath,outName])
    performanceLine = "\t".join([str(x) for x in [totalRunTime, maxRSS_mem, maxVMS_mem, testName,command]]) + "\n"
    pFile = open(performanceMetricFile, "a")
    pFile.write(performanceLine)
    pFile.close()
    return outputLog, "Plot ran successfully"


def launchPhasingCleaningChr(PhaseResultFolder, readSNPsPath, max_ID, min_len, min_ovl, min_sim, min_overlap, min_ovl_len,percentKept,deduplicate):
    start = timeit.default_timer()
    outputLog = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()) + "\t"

    command_py = "/home/gaoyun/poly/code/phasing/final_process_data/real/GB54/scripts/cleaning_phase_chr_nx.py"
    p = subprocess.Popen(["python3",command_py, PhaseResultFolder, readSNPsPath, max_ID, min_len, min_ovl, min_sim, min_overlap, min_ovl_len,percentKept,deduplicate], stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    outputLog += "COMMAND: " + " ".join(["python3", command_py, PhaseResultFolder, readSNPsPath, max_ID, min_len, min_ovl, min_sim, min_overlap, min_ovl_len]) + "\n\n"

    maxRSS_mem = 0
    maxVMS_mem = 0

    try:
        poll = p.poll()
    except:
        poll = "Done"
    while poll is None:
        monitorP = psutil.Process(p.pid)
        children = list(monitorP.children(recursive=True))
        children = children + [monitorP]
        RSS_mem = 0
        VMS_mem = 0
        for child in children:
            mem_info = child.memory_info()
            RSS_mem += mem_info[0]
            VMS_mem += mem_info[1]
            if RSS_mem > maxRSS_mem:
                maxRSS_mem = RSS_mem
            if VMS_mem > maxVMS_mem:
                maxVMS_mem = VMS_mem
        time.sleep(1)
        try:
            poll = p.poll()
        except:
            poll = "Done"

    if p.stderr != "" or p.stdout != "":
        outputLog += "STDERR:\n\n" + p.stderr.read() + "\n\nSTDOUT:\n\n" + p.stdout.read() + "\n\n"

    print("Max memory usage, RSS:", maxRSS_mem / 1024 / 1024 / 1024, "VMS:", maxVMS_mem / 1024 / 1024 / 1024)

    stop = timeit.default_timer()
    totalRunTime = stop - start
    print("This took", totalRunTime, "seconds to run.")
    command = " ".join(["python3", command_py, PhaseResultFolder, readSNPsPath, max_ID, min_len, min_ovl, min_sim, min_overlap, min_ovl_len])
    performanceLine = "\t".join([str(x) for x in [totalRunTime, maxRSS_mem, maxVMS_mem, testName,command]]) + "\n"
    pFile = open(performanceMetricFile, "a")
    pFile.write(performanceLine)
    pFile.close()
    return outputLog, "Phasing Cleaning ran successfully"



def launchnTChap(mappedLR, vcfFile, PhasingPath, referenceFilePath, threads, testName):
    start = timeit.default_timer()
    outputLog = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()) + "\t"

    p = subprocess.Popen(["ntchap", "--bam", mappedLR, "--vcf", vcfFile, "--output", PhasingPath, "--reference", referenceFilePath, "-t", threads, "--sample", testName], stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    outputLog += "COMMAND: " + " ".join(["ntchap", "--bam", mappedLR, "--vcf", vcfFile, "--output", PhasingPath, "--reference", referenceFilePath, "-t", threads, "--sample", testName]) + "\n\n"

    maxRSS_mem = 0
    maxVMS_mem = 0

    try:
        poll = p.poll()
    except:
        poll = "Done"
    while poll is None:
        monitorP = psutil.Process(p.pid)
        children = list(monitorP.children(recursive=True))
        children = children + [monitorP]
        RSS_mem = 0
        VMS_mem = 0
        for child in children:
            mem_info = child.memory_info()
            RSS_mem += mem_info[0]
            VMS_mem += mem_info[1]
            if RSS_mem > maxRSS_mem:
                maxRSS_mem = RSS_mem
            if VMS_mem > maxVMS_mem:
                maxVMS_mem = VMS_mem
        time.sleep(1)
        try:
            poll = p.poll()
        except:
            poll = "Done"

    if p.stderr != "" or p.stdout != "":
        outputLog += "STDERR:\n\n" + p.stderr.read() + "\n\nSTDOUT:\n\n" + p.stdout.read() + "\n\n"

    print("Max memory usage, RSS:", maxRSS_mem / 1024 / 1024 / 1024, "VMS:", maxVMS_mem / 1024 / 1024 / 1024)

    stop = timeit.default_timer()
    totalRunTime = stop - start
    print("This took", totalRunTime, "seconds to run.")
    command = " ".join(["ntchap", "--bam", mappedLR, "--vcf", vcfFile, "--output", PhasingPath, "--reference", referenceFilePath, "-t", threads, "--sample", testName])
    performanceLine = "\t".join([str(x) for x in [totalRunTime, maxRSS_mem, maxVMS_mem, testName,command]]) + "\n"
    pFile = open(performanceMetricFile, "a")
    pFile.write(performanceLine)
    pFile.close()
    return outputLog, "ntchap ran successfully"


if __name__ == "__main__":
    # Load parameters
    benchmarkParameterFilePath = sys.argv[1]
    benchmarkParameters = loadBenchmarkParameterFile(benchmarkParameterFilePath)
    phasingMethods=sys.argv[2:]

    # Handle paths
    mainPath = benchmarkParameters["mainPath"]

    speciesName=benchmarkParameters["reference"]["speciesName"]
    refName=benchmarkParameters["reference"]["refName"]
    shortReadsName = benchmarkParameters["shortReads"]
    longReadsName = benchmarkParameters["longReads"]

    predictionPath=os.path.join(mainPath, "Phasing/")
    hybridRawReads=os.path.join(mainPath, "rawData/")
    hybridLongReads=os.path.join(hybridRawReads, "longReads/")
    referencePath=os.path.join(mainPath, "reference/", speciesName+"/")
    referenceFilePath=os.path.join(referencePath, refName)
    mappedLongReads=os.path.join(mainPath, "Mapped/longReads/")
    variantCalledShortReads=os.path.join(mainPath, "Variants/shortReads/")
    performanceMetricPath=os.path.join(mainPath, "performanceMetrics/")

    allPaths=[predictionPath,performanceMetricPath]

    logPath=os.path.join(mainPath,"log/")
    allPaths.append(logPath)

    for path in allPaths:
        os.makedirs(path, exist_ok=True)

    #Prepare log

    fullLogPath=os.path.join(logPath, "phasingLog.txt")
    logText=""

    #Prepare performance metrics file
    performanceMetricFile=os.path.join(performanceMetricPath, "timeAndMemoryMetrics.txt")
    pFile = open(performanceMetricFile, "a")
    pFile.write("#time (seconds)\tRSS (bytes)\tVMS (bytes)\ttoolName\ttestName\n")
    pFile.close()

    threads=benchmarkParameters["threads"]

    # run our method
    if "ntchap" in phasingMethods:
        max_ID = benchmarkParameters["max_ID"]
        min_len = benchmarkParameters["min_len"] 
        min_ovl = benchmarkParameters["min_ovl"]  
        min_sim = benchmarkParameters["min_sim"] 
        min_overlap = benchmarkParameters["min_overlap"] 
        min_ovl_len = benchmarkParameters["min_ovl_len"] 
        for short_read in shortReadsName:
            for long_read in longReadsName:
                testName=speciesName+long_read
                vcfFile=os.path.join(variantCalledShortReads, short_read+".vcf")
                mappedLR=os.path.join(mappedLongReads, long_read+".sorted.sam")
                launchFunction("ntchap", [mappedLR, vcfFile, predictionPath, referenceFilePath, threads, testName])
        print("Done with all ntchap")


    if "plot" in phasingMethods:
        max_ID = benchmarkParameters["max_ID"]
        min_len = benchmarkParameters["min_len"] 
        min_ovl = benchmarkParameters["min_ovl"]  
        min_sim = benchmarkParameters["min_sim"] 
        min_overlap = benchmarkParameters["min_overlap"] 
        min_ovl_len = benchmarkParameters["min_ovl_len"]
        phasedName = "phased"+"_"+str(min_ovl)+"_"+str(min_sim)+"_"+str(max_ID)+"_"+str(min_len)+"_"+str(min_overlap)+"_"+str(min_ovl_len)
        for short_read in shortReadsName:
            for long_read in longReadsName:
                testName=speciesName+long_read
                longReadFastQFileDir = os.path.join(mainPath,"rawData/longReads",long_read+".fastq.gz")
                outdir = os.path.join(predictionPath, testName,"Phased")
                consensusClusterPath = os.path.join(outdir, phasedName+"_variants.tsv")
                readClusterPath = os.path.join(outdir, phasedName+"_clustered_read_name.tsv")
                launchFunction("plot", [consensusClusterPath, readClusterPath,outdir, longReadFastQFileDir,phasedName])


    if "Phasing_clean_chr" in phasingMethods:
        max_ID = benchmarkParameters["max_ID"]
        min_len = benchmarkParameters["min_len"] 
        min_ovl = benchmarkParameters["min_ovl"]  
        min_sim = benchmarkParameters["min_sim"] 
        min_overlap = benchmarkParameters["min_overlap"] 
        min_ovl_len = benchmarkParameters["min_ovl_len"]
        for short_read in shortReadsName:
            for long_read in longReadsName:
                testName=speciesName+long_read
                PhaseResultFolder = os.path.join(predictionPath,testName,"Phased")
                readSNPsPath = os.path.join(predictionPath,testName, "Variants", "Reads", testName+".hetPositions.SNPxLongReads.vcf.chr.tsv")
                percentKept="99"
                deduplicate="1"
                launchFunction("Phasing_clean_chr", [PhaseResultFolder, readSNPsPath, max_ID, min_len, min_ovl, min_sim, min_overlap, min_ovl_len,percentKept,deduplicate])
        print("Done with all Phasing cleaning_chr")


    if "plot_clean_chr" in phasingMethods:
        clean_ways = ["filter","dup"]
        max_ID = benchmarkParameters["max_ID"]
        min_len = benchmarkParameters["min_len"] 
        min_ovl = benchmarkParameters["min_ovl"]  
        min_sim = benchmarkParameters["min_sim"] 
        min_overlap = benchmarkParameters["min_overlap"] 
        min_ovl_len = benchmarkParameters["min_ovl_len"]
        phasedName = "phased"+"_"+str(min_ovl)+"_"+str(min_sim)+"_"+str(max_ID)+"_"+str(min_len)+"_"+str(min_overlap)+"_"+str(min_ovl_len)
        for short_read in shortReadsName:
            for long_read in longReadsName:
                testName=speciesName+long_read
                longReadFastQFileDir = os.path.join(mainPath,"rawData/longReads",long_read+".fastq.gz")
                outdir = os.path.join(predictionPath, testName,"Phased", "Cleaned_chr")
                for clean_way in clean_ways:
                    consensusClusterPath = os.path.join(outdir, phasedName+"_cleaned_"+clean_way+"_variants.tsv")
                    readClusterPath = os.path.join(outdir, phasedName+"_cleaned_"+clean_way+"_clustered_read_name.tsv")
                    launchFunction("plot", [consensusClusterPath, readClusterPath,outdir, longReadFastQFileDir,phasedName+"_cleaned_"+clean_way])

