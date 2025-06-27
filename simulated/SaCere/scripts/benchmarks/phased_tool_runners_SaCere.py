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
                           "longReads": {},
                           "longReadMethod": "",
                           "shortReads": {},
                           "coverages": [],
                           "heterozygosityRates": [],
                           "reference": {},
                           "genomeSize": 0,
                           "strainLists": [],
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
        if line[0] in {"mainPath","longReadMethod","genomeSize","threads","max_ID","min_len","min_ovl","min_sim","min_overlap","min_ovl_len"}:
            benchmarkParameters[line[0]]=line[1]
        elif line[0] in {"coverages","heterozygosityRates"}:
            for dataPoint in line[1:]:
                benchmarkParameters[line[0]].append(dataPoint)
        elif line[0] in {"longReads","shortReads","reference"}:
            for dataPoint in line[1:]:
                dataPoint=dataPoint.split(":")
                benchmarkParameters[line[0]][dataPoint[0]]=dataPoint[1]
        elif line[0] in {"strainLists"}:
            for dataPoint in line[1:]:
                dataPoint = dataPoint.split(",")
                benchmarkParameters[line[0]].append(dataPoint)
        else:
            print("Ignoring '","\t".join(line),"' because it is not recognized",sep="")
    benchmarkParameterFile.close()
    return benchmarkParameters

def launchFunction(functionName,parameters):
    if functionName=="nphase":
        strainName=parameters[0]
        referenceFilePath=parameters[1]
        outputFolder=parameters[2]
        longReadFile=parameters[3]
        vcfFile=parameters[4]
        mappedLR=parameters[5]
        threadNumber=parameters[6]
        outputLog, systemMessage=launchnPhase(strainName, referenceFilePath, outputFolder, longReadFile, vcfFile, mappedLR, threadNumber)
        print(systemMessage)
        updateLog(fullLogPath, outputLog)
    elif functionName=="whatshap polyphase":
        ploidy=parameters[0]
        whatsHapPolyphaseParameter=parameters[1]
        outPath=parameters[2]
        referenceFilePath=parameters[3]
        threads=parameters[4]
        vcfFile=parameters[5]
        mappedLR=parameters[6]
        indelBool=parameters[7]
        localStrainName=parameters[8]
        outputLog, systemMessage=launchWhatshapPolyphase(ploidy,whatsHapPolyphaseParameter,outPath,referenceFilePath,threads,
                                                         vcfFile,mappedLR,indelBool,localStrainName)
        print(systemMessage)
        updateLog(fullLogPath, outputLog)
    elif functionName=="flopp":
        mappedLR=parameters[0]
        vcfFile=parameters[1]
        ploidy=parameters[2]
        outPath=parameters[3]
        readPartionPath = parameters[4]
        threads=parameters[5]
        localStrainName=parameters[6]
        outputLog, systemMessage=launchFlopp(mappedLR,vcfFile,ploidy,outPath,readPartionPath,threads,localStrainName)
        print(systemMessage)
        updateLog(fullLogPath, outputLog)
    elif functionName=="ntchap":
        mappedLR=parameters[0]
        vcfFile=parameters[1]
        PhasingPath=parameters[2]
        referenceFilePath=parameters[3]
        threads=parameters[4]
        testName=parameters[5]
        outputLog, systemMessage=launchnTChap(mappedLR, vcfFile, PhasingPath, referenceFilePath, threads, testName)
        print(systemMessage)
        updateLog(fullLogPath, outputLog)


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


def launchFlopp(mappedLR,vcfFile,ploidy,outPath,threads,localStrainName):
    start = timeit.default_timer()
    outputLog = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()) + "\t"

    p = subprocess.Popen(["flopp","-b",mappedLR, "-c", vcfFile, "-p", ploidy, "-o", outPath, "-t", threads], stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    outputLog += "COMMAND: " + " ".join(["flopp","-b",mappedLR, "-c", vcfFile, "-p", ploidy, "-o", outPath, "-t", threads]) + "\n\n"

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
    command = " ".join(["flopp","-b",mappedLR, "-c", vcfFile, "-p", ploidy, "-o", outPath, "-t", threads])
    performanceLine = "\t".join([str(x) for x in [totalRunTime, maxRSS_mem, maxVMS_mem, localStrainName, command]]) + "\n"
    pFile = open(performanceMetricFile, "a")
    pFile.write(performanceLine)
    pFile.close()
    return outputLog, "flopp ran successfully"

def launchnPhase(strainName,referenceFilePath,outputFolder,longReadFile,vcfFile,mappedLR,threadNumber):
    start = timeit.default_timer()
    outputLog = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()) + "\t"

    p = subprocess.Popen(["nphase", "partial", "--sampleName", strainName, "--reference", referenceFilePath, "--output",
                        outputFolder, "--longReads", longReadFile, "--vcf", vcfFile, "--mappedLongReads",
                        mappedLR, "--threads", threadNumber], stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["nphase partial --sampleName",strainName,"--reference",referenceFilePath,"--output",
                                     outputFolder,"--longReads",longReadFile,"--vcf",vcfFile,"--mappedLongReads",
                                     mappedLR,"--threads", threadNumber]) + "\n\n"

    maxRSS_mem=0
    maxVMS_mem=0

    try:
        poll = p.poll()
    except:
        poll="Done"
    while poll is None:
        monitorP=psutil.Process(p.pid)
        children=list(monitorP.children(recursive=True))
        children=children+[monitorP]
        RSS_mem=0
        VMS_mem=0
        for child in children:
            mem_info=child.memory_info()
            RSS_mem+=mem_info[0]
            VMS_mem+=mem_info[1]
            if RSS_mem>maxRSS_mem:
                maxRSS_mem=RSS_mem
            if VMS_mem>maxVMS_mem:
                maxVMS_mem=VMS_mem
        time.sleep(1)
        try:
            poll=p.poll()
        except:
            poll="Done"


    if p.stderr!="" or p.stdout!="":
        outputLog+="STDERR:\n\n"+p.stderr.read()+"\n\nSTDOUT:\n\n"+p.stdout.read()+"\n\n"

    print("Max memory usage, RSS:",maxRSS_mem/1024/1024/1024,"VMS:",maxVMS_mem/1024/1024/1024)

    stop = timeit.default_timer()
    totalRunTime = stop - start
    print("This took", totalRunTime, "seconds to run.")
    command = " ".join(["nphase partial --sampleName",strainName,"--reference",referenceFilePath,"--output",
                                     outputFolder,"--longReads",longReadFile,"--vcf",vcfFile,"--mappedLongReads",
                                     mappedLR,"--threads", threadNumber])
    performanceLine="\t".join([str(x) for x in [totalRunTime,maxRSS_mem,maxVMS_mem,strainName,command]])+"\n"
    pFile=open(performanceMetricFile,"a")
    pFile.write(performanceLine)
    pFile.close()
    return outputLog, "nPhase ran successfully"

def launchWhatshapPolyphase(ploidy,whatsHapPolyphaseParameter,outPath,referenceFilePath,threads,vcfFile,mappedLR,indelBool,localStrainName):
    start = timeit.default_timer()

    outputLog = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()) + "\t"

    if indelBool==False:
        p = subprocess.Popen(["whatshap","polyphase","--ploidy",ploidy,"--ignore-read-groups","--block-cut-sensitivity",
                              whatsHapPolyphaseParameter ,"-o",outPath,"--reference",referenceFilePath,"--include-haploid-sets","--threads",threads,
                              vcfFile,mappedLR], stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
        outputLog += "COMMAND: " + " ".join(["whatshap","polyphase","--ploidy",ploidy,"--ignore-read-groups","--block-cut-sensitivity",
                                             whatsHapPolyphaseParameter ,"-o",outPath,"--reference",referenceFilePath,"--include-haploid-sets","--threads",
                                             threads,vcfFile,mappedLR]) + "\n\n"
    else:
        p = subprocess.Popen(["whatshap", "polyphase", "--ploidy", ploidy, "--indels","--ignore-read-groups",
                              "--block-cut-sensitivity",whatsHapPolyphaseParameter, "-o", outPath, "--reference",
                              referenceFilePath,"--include-haploid-sets", "--threads", threads,vcfFile, mappedLR], stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
        outputLog += "COMMAND: " + " ".join(["whatshap", "polyphase", "--ploidy", ploidy,"--indels", "--ignore-read-groups",
                                             "--block-cut-sensitivity",whatsHapPolyphaseParameter, "-o", outPath, "--reference",
                                             referenceFilePath,"--include-haploid-sets", "--threads",threads, vcfFile, mappedLR]) + "\n\n"

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
    print("This took", totalRunTime / 60 / 60, "hours to run.")
    command = " ".join(["whatshap","polyphase","--ploidy",ploidy,"--ignore-read-groups","--block-cut-sensitivity",
                                             whatsHapPolyphaseParameter ,"-o",outPath,"--reference",referenceFilePath,"--include-haploid-sets","--threads",
                                             threads,vcfFile,mappedLR])
    performanceLine="\t".join([str(x) for x in [totalRunTime,maxRSS_mem,maxVMS_mem,localStrainName,command]])+"\n"
    pFile=open(performanceMetricFile,"a")
    pFile.write(performanceLine)
    pFile.close()
    return outputLog, "whatshap polyphase ran successfully"


if __name__ == "__main__":
    # Load parameters
    benchmarkParameterFilePath = sys.argv[1]
    benchmarkParameters = loadBenchmarkParameterFile(benchmarkParameterFilePath)
    phasingMethods=sys.argv[2:]

    # Handle paths
    mainPath = benchmarkParameters["mainPath"]

    speciesName=benchmarkParameters["reference"]["speciesName"]

    predictionPath=os.path.join(mainPath, "phasingPredictions/")
    hybridPath=os.path.join(mainPath, "virtualPolyploids/")
    hybridRawReads=os.path.join(hybridPath, "rawReads/")
    hybridLongReads=os.path.join(hybridRawReads, "longReads/")
    referencePath=os.path.join(mainPath, "reference/", speciesName+"/")
    referenceFilePath=os.path.join(referencePath, speciesName+".fasta")
    hybridPath=os.path.join(mainPath, "virtualPolyploids/")
    mappedLongReads=os.path.join(hybridPath, "Mapped/longReads/")
    variantCalledShortReads=os.path.join(hybridPath, "Variants/shortReads/")
    performanceMetricPath=os.path.join(mainPath, "performanceMetrics/")

    allPaths=[predictionPath,performanceMetricPath]

    if "flopp" in phasingMethods:
        floppPath = os.path.join(predictionPath, "flopp/")
        allPaths.append(floppPath)
    if "nphase" in phasingMethods:
        nPhasePath=os.path.join(predictionPath, "nPhase/")
        allPaths.append(nPhasePath)
    if "whatshap-polyphase" in phasingMethods:
        whatsHapPolyphasePath=os.path.join(predictionPath, "whatsHapPolyphase/")
        allPaths.append(whatsHapPolyphasePath)
    if "ntchap" in phasingMethods:
        PhasingPath = os.path.join(predictionPath, "nTChap/")
        allPaths.append(PhasingPath)

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

    #Run all flopp
    if "flopp" in phasingMethods:
        for strainList in benchmarkParameters["strainLists"]:
            strainName="_".join(strainList)
            ploidy = str(len(strainList))
            for coverage in benchmarkParameters["coverages"]:
                for heterozygosityRate in benchmarkParameters["heterozygosityRates"]:
                    testName = strainName + "_" + str(coverage) + "X_"+heterozygosityRate
                    vcfFile=os.path.join(variantCalledShortReads, testName+".SNPs.vcf")
                    mappedLR = os.path.join(mappedLongReads, strainName + "_" + str(coverage) + "X.sorted.bam")
                    outPath = os.path.join(floppPath,"flopp_"+testName+".SNPs.vcf")
                    launchFunction("flopp",[mappedLR, vcfFile, ploidy, outPath, threads,testName])

        print("Done with all flopp")

    #Run all whatshap polyphase
    if "whatshap-polyphase" in phasingMethods:
        for strainList in benchmarkParameters["strainLists"]:
            strainName="_".join(strainList)
            ploidy = str(len(strainList))
            for coverage in benchmarkParameters["coverages"]:
                for heterozygosityRate in benchmarkParameters["heterozygosityRates"]:
                    whatsHapPolyphaseParameter = "4"
                    testName = strainName + "_" + str(coverage) + "X_" + heterozygosityRate
                    vcfFile=os.path.join(variantCalledShortReads, testName+".SNPs.vcf")
                    mappedLR = os.path.join(mappedLongReads, strainName + "_" + str(coverage) + "X.sorted.bam")
                    outPath = os.path.join(whatsHapPolyphasePath, "WHP_"+testName+"_B"+whatsHapPolyphaseParameter+".SNPs.vcf")

                    indelBool = False
                    launchFunction("whatshap Polyphase",[ploidy,whatsHapPolyphaseParameter,outPath,referenceFilePath,threads,vcfFile,
                                                         mappedLR,indelBool,testName])

        print("Done with all whatshap polyphase")

    #Run all nPhase
    if "nphase" in phasingMethods:
        for strainList in benchmarkParameters["strainLists"]:
            strainName="_".join(strainList)
            ploidy = str(len(strainList))
            for coverage in benchmarkParameters["coverages"]:
                for heterozygosityRate in benchmarkParameters["heterozygosityRates"]:
                    testName=strainName+"_"+str(coverage)+"X_"+heterozygosityRate
                    longReadFile=os.path.join(hybridLongReads, strainName+"_"+str(coverage)+"X.fastq.gz")
                    vcfFile=os.path.join(variantCalledShortReads, testName+".SNPs.vcf")
                    mappedLR=os.path.join(mappedLongReads, strainName+"_"+str(coverage)+"X.sorted.sam")
                    launchFunction("nphase", [testName,referenceFilePath,nPhasePath,longReadFile,vcfFile,mappedLR,threads])

        print("Done with all nPhase")

    # run our method
    if "ntchap" in phasingMethods:
        max_ID = benchmarkParameters["max_ID"]
        min_len = benchmarkParameters["min_len"] 
        min_ovl = benchmarkParameters["min_ovl"]  
        min_sim = benchmarkParameters["min_sim"] 
        min_overlap = benchmarkParameters["min_overlap"] 
        min_ovl_len = benchmarkParameters["min_ovl_len"] 
        for strainList in benchmarkParameters["strainLists"]:
            strainName="_".join(strainList)
            ploidy = str(len(strainList))
            for coverage in benchmarkParameters["coverages"]:
                for heterozygosityRate in benchmarkParameters["heterozygosityRates"]:
                    testName=strainName+"_"+str(coverage)+"X_"+heterozygosityRate
                    vcfFile=os.path.join(variantCalledShortReads, testName+".SNPs.vcf")
                    mappedLR=os.path.join(mappedLongReads, strainName+"_"+str(coverage)+"X.sorted.sam")
                    launchFunction("ntchap", [mappedLR, vcfFile, PhasingPath, referenceFilePath, threads, testName])

        print("Done with all ntchap")

    print("All datasets phased with all tools")
