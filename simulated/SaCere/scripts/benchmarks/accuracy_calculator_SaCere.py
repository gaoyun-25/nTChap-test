import sys
import os
import random
import subprocess
import pysam
import glob

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

def loadVCFInformation(vcfFilePath):
    VCFInformation={}
    vcfFile=open(vcfFilePath,"r")
    for line in vcfFile:
        line=line.strip("\n")
        if "#" not in line:
            line=line.split("\t")
            chromosome=line[0]
            position=line[1]
            SNPs=[line[3]]+line[4].split(",")
            VCFInformation.setdefault(chromosome,{})
            VCFInformation[chromosome][position]=SNPs
    vcfFile.close()
    return VCFInformation

def loadClusteredReadName(clusterReadFile):
    # cid:rid
    ClusteredReads = dict()
    with open(clusterReadFile) as f:
        for line in f:
            line = line.strip().split("\t")
            ClusteredReads.setdefault(line[0], set())
            ClusteredReads[line[0]].add(line[1])

    return ClusteredReads

def loadnPhaseClusteredReadName(clusterReadFile):
    # cid:rid
    ClusteredReads = dict()
    with open(clusterReadFile) as f:
        for line in f:
            line = line.strip().split("\t")
            ClusteredReads.setdefault(line[0], set())
            ClusteredReads[line[0]].add(line[1].split("_VCSTTP")[0])

    return ClusteredReads

def floppResultTranslator(floppResultFilePath,outputFilePath,VCFInformation):
    translatedFloppResults=""
    floppFile=open(floppResultFilePath,"r")
    chromosomeName=None
    haplotigIndex=0
    haplotigNames=set()
    for line in floppFile:
        line=line.strip("\n")
        if line=="*****":
            haplotigIndex=max(haplotigNames)
        elif line[0:2]=="**" and line[-2:]=="**":
            chromosomeName=line[2:-2]
        else:
            line=line.split("\t")
            SNPPosition=line[0].split(":")[1]
            SNPCodes=line[1:1+(len(line)-1)//2] #There will always be position+nSNPs+nSNPs in the results
            localIndex=0
            for SNPCode in SNPCodes:
                localIndex += 1
                currentHaplotigName = "flopp_" + str(haplotigIndex + localIndex)
                haplotigNames.add(haplotigIndex + localIndex)
                if int(SNPCode)==-1:
                    SNP="*" #-1 means no read coverage, which means indel.
                else:
                    SNP=VCFInformation[chromosomeName][SNPPosition][int(SNPCode)]
                translatedFloppResults+="\t".join([currentHaplotigName,chromosomeName,SNPPosition,SNP])+"\n"
    floppFile.close()
    outputFile=open(outputFilePath,"w")
    outputFile.write(translatedFloppResults)
    outputFile.close()

def whatsHapPolyphaseResultTranslator(whatsHapPolyphaseFilePath,outputFilePath):
    translatedWhatsHapPolyphaseResults=""
    whatsHapPolyphaseFile=open(whatsHapPolyphaseFilePath,"r")
    for line in whatsHapPolyphaseFile:
        line=line.strip("\n")
        if line[0]=="#":
            pass
        else:
            line=line.split("\t")
            if "PS" in line[8] and "HS" in line[8]:
                chromosomeName=line[0]
                SNPPosition=line[1]
                SNPs=[line[3]]+line[4].split(",")
                infoBlock=line[9]
                GTBlock=infoBlock.split(":")[0]
                # phaseBlocks=infoBlock.split(":")[-1].split(",")
                phaseBlocks = infoBlock.split(":")[-2]
                localPhase=0
                for SNPCode in GTBlock.split("|"):
                    localPhase+=1
                    # currentHaplotigName=chromosomeName+"_"+phaseBlocks[localPhase-1]+"_"+str(localPhase)
                    currentHaplotigName=chromosomeName+"_"+phaseBlocks+"_"+str(localPhase)
                    SNP=SNPs[int(SNPCode)]
                    translatedWhatsHapPolyphaseResults+="\t".join([currentHaplotigName,chromosomeName,SNPPosition,SNP])+"\n"
    whatsHapPolyphaseFile.close()
    outputFile=open(outputFilePath,"w")
    outputFile.write(translatedWhatsHapPolyphaseResults)
    outputFile.close()


def whatsHapPolyphaseResultVariantReads(bam_file, phased_bam, reference_file, ploidy, phased_vcf, variant_reads_file):
    phased_vcf_gz = phased_vcf + ".gz"
    # # gzip vcf
    p = subprocess.run(["bgzip","-o",phased_vcf_gz,phased_vcf],stderr=subprocess.PIPE,stdout=subprocess.PIPE,universal_newlines=True)
    # index vcf
    p = subprocess.run(["tabix","-p","vcf",phased_vcf_gz],stderr=subprocess.PIPE,stdout=subprocess.PIPE,universal_newlines=True)
    # tag reads by haplotype
    p = subprocess.run(["whatshap","haplotag","--ignore-read-groups","--ploidy",ploidy,"-o",phased_bam,"--reference",reference_file,phased_vcf_gz,bam_file],stderr=subprocess.PIPE,stdout=subprocess.PIPE,universal_newlines=True)

    # output clusterReadNames.tsv
    bf = pysam.AlignmentFile(phased_bam, "rb")
    with open(variant_reads_file, "w") as f:
        for r in bf:
            flag = r.flag
            if flag == 4:
                continue
            tags =r.tags
            reference = r.reference_name
            read_name = r.query_name
            rf_pos = r.pos
            HP = None
            PS = None
            for tag in tags:
                if tag[0] == 'HP':
                    HP = tag[1]
                if tag[0] == 'PS':
                    PS = tag[1]
            if HP != None and PS != None:
                cid = reference + "_" + str(PS)+"_"+str(HP)
                f.write(cid+"\t"+read_name+"\n")


#Now a function to generate the "true" phases using the original VCFs
def getTruePhases(originalVCFs):
    referenceDict={}
    truePhaseDict={}
    heterozygosityChecker={}
    sampleNames=[]
    for originalVCFFilePath, sampleName in originalVCFs:
        sampleNames.append(sampleName)
        originalVCFFile=open(originalVCFFilePath,"r")
        for line in originalVCFFile:
            line=line.strip("\n")
            if line[0]=="#":
                pass
            else:
                line=line.split("\t")
                chromosome=line[0]
                SNPPosition=line[1]
                refSNP=line[3]
                SNPs=line[4].split(",")
                for SNP in SNPs:
                    referenceDict.setdefault(chromosome,{})
                    referenceDict[chromosome][SNPPosition]=chromosome+":"+SNPPosition+"="+refSNP
                    truePhaseDict.setdefault(chromosome,{})
                    truePhaseDict[chromosome].setdefault(SNPPosition, {})
                    truePhaseDict[chromosome][SNPPosition].setdefault(sampleName,set())
                    fullSNP=chromosome+":"+SNPPosition+"="+SNP
                    truePhaseDict[chromosome][SNPPosition][sampleName].add(fullSNP)
                    heterozygosityChecker.setdefault(chromosome, {})
                    heterozygosityChecker[chromosome].setdefault(SNPPosition, set())
                    heterozygosityChecker[chromosome][SNPPosition].add(SNP)
        originalVCFFile.close()
    #Must include the ref allele for positions present in one of the VCFs but not others (Keep track of ref)
    for chromosome, SNPPositions in referenceDict.items():
        for SNPPosition, refHaplotype in SNPPositions.items():
            for sampleName in sampleNames:
                if sampleName not in truePhaseDict[chromosome][SNPPosition].keys():  # 有菌株未覆盖此pos，则用ref base作为其pos base 补上，这样可以保证ref上的每一个snp pos都有4种菌株cov
                    truePhaseDict[chromosome][SNPPosition][sampleName]={refHaplotype}
                    heterozygosityChecker[chromosome][SNPPosition].add(refHaplotype)
    #Only want to include heterozygous SNPs
    heterozygousPhaseDict={}
    for chromosome, chromosomeData in truePhaseDict.items():
        heterozygousPhaseDict.setdefault(chromosome,{})
        for SNPPosition, haplotypeData in chromosomeData.items():
            for sampleName, sampleSNPs in haplotypeData.items():
                heterozygousPhaseDict[chromosome].setdefault(sampleName,set())
                for fullSNP in sampleSNPs:
                    if len(heterozygosityChecker[chromosome][SNPPosition])>1:
                        heterozygousPhaseDict[chromosome][sampleName].add(fullSNP)
    return heterozygousPhaseDict

#Now a function to load translated results
def loadPhasingResults(phasingResultFilePath):
    phasingResults={}
    phasingResultFile=open(phasingResultFilePath,"r")
    for line in phasingResultFile:
        line=line.strip("\n").split("\t")
        haplotigName=line[0]
        chromosome=line[1]
        position=line[2]
        SNP=line[3]
        fullSNP=chromosome+":"+position+"="+SNP
        phasingResults.setdefault(chromosome, {})
        phasingResults[chromosome].setdefault(haplotigName,set())
        phasingResults[chromosome][haplotigName].add(fullSNP)
    phasingResultFile.close()
    return phasingResults

#Now a function to compare & calculate, when given a standardized file
def calculateAccuracyMetrics(groundTruthInfo,predictedInfo,PhasingReadsResults,longReadsTag,testInfo,accuracyMetricFile):
    accuracyData={"bestScore":{"Total":0},"totalTruePositives":{},"totalFalsePositives":{},"totalFalseNegatives":{},"results":{}, "readsClusterAccuracy":{}}
    #Initializing
    haplotypeNames=set()
    nContigs=0
    for chromosome, haplotypes in groundTruthInfo.items():
        for haplotypeName, SNPs in haplotypes.items():
            nContigs+=1
            haplotypeNames.add(haplotypeName)
            accuracyData["bestScore"].setdefault(haplotypeName,0)
            accuracyData["bestScore"][haplotypeName]+=len(SNPs)
            accuracyData["bestScore"]["Total"]+=len(SNPs)
            accuracyData["totalTruePositives"].setdefault(haplotypeName, 0)
            accuracyData["totalFalsePositives"].setdefault(haplotypeName,0)
            accuracyData["totalFalseNegatives"].setdefault(haplotypeName,set())
            SNPPos=set([fullSNP.split("=")[0] for fullSNP in SNPs])
            accuracyData["totalFalseNegatives"][haplotypeName]=accuracyData["totalFalseNegatives"][haplotypeName].union(SNPPos)
            accuracyData["readsClusterAccuracy"].setdefault(haplotypeName,dict())
            accuracyData["readsClusterAccuracy"][haplotypeName]["Total"] = 0
            accuracyData["readsClusterAccuracy"][haplotypeName]["nTrueReads"] = 0

    nHaplotigs=0
    for chromosome, haplotypes in predictedInfo.items():
        for predictionName, SNPs in haplotypes.items():
            nHaplotigs+=1  # 针对predict中每条染色体每个块中的snps pos都找他们跟ground truth中哪个菌株更相似，用TP/块中snps数
            closestHaplotypeName,TP,FP,FNPos=getBestMatch(SNPs,groundTruthInfo)
            if closestHaplotypeName=="None":
                closestHaplotypeName=random.choice(list(haplotypeNames))
            accuracyData["totalTruePositives"][closestHaplotypeName]+=TP
            accuracyData["totalFalsePositives"][closestHaplotypeName]+=FP
            accuracyData["totalFalseNegatives"][closestHaplotypeName].difference_update(FNPos)
            closestHaplotypeReadsPrefix = longReadsTag[closestHaplotypeName]
            if predictionName in PhasingReadsResults:
                readNameSet = PhasingReadsResults[predictionName]
                for readName in readNameSet:
                    accuracyData["readsClusterAccuracy"][closestHaplotypeName]["Total"] += 1
                    readNameHap = readName.split(".")[0]
                    if readNameHap == closestHaplotypeReadsPrefix:
                        accuracyData["readsClusterAccuracy"][closestHaplotypeName]["nTrueReads"] += 1

    accuracyData["results"]["nHaplotigs"]=nHaplotigs
    accuracyData["results"]["nContigs"]=nContigs
    fullTP=0
    fullFP=0
    fullFN=0
    readsAccurTotal = 0
    readsAccurTrue = 0
    for haplotypeName in haplotypeNames:
        accuracyData["results"].setdefault(haplotypeName,{"TP":0,"FP":0,"FN":0})
        TP=accuracyData["totalTruePositives"][haplotypeName]
        FP=accuracyData["totalFalsePositives"][haplotypeName]
        FN=len(accuracyData["totalFalseNegatives"][haplotypeName])
        total=TP+FP+FN
        fullTP+=TP
        fullFP+=FP
        fullFN+=FN
        accuracyData["totalFalseNegatives"][haplotypeName]=len(accuracyData["totalFalseNegatives"][haplotypeName])
        accuracyData["results"][haplotypeName]["TP"]=round(TP/total*100,2)
        accuracyData["results"][haplotypeName]["FP"]=round(FP/total*100,2)
        accuracyData["results"][haplotypeName]["FN"]=round(FN/total*100,2)
        readsAccurTotal += accuracyData["readsClusterAccuracy"][haplotypeName]["Total"]
        readsAccurTrue += accuracyData["readsClusterAccuracy"][haplotypeName]["nTrueReads"]
    fullTotal=fullTP+fullFP+fullFN
    accuracyData["results"]["all"]={}
    accuracyData["results"]["all"]["TP"]=round(fullTP/fullTotal*100,2)
    accuracyData["results"]["all"]["FP"]=round(fullFP/fullTotal*100,2)
    accuracyData["results"]["all"]["FN"]=round(fullFN/fullTotal*100,2)

    readsClusterAccuracy = round(readsAccurTrue/readsAccurTotal*100,2)
    resultTable=testInfo+[accuracyData["results"]["all"]["TP"],accuracyData["results"]["all"]["FP"],accuracyData["results"]["all"]["FN"],accuracyData["results"]["nHaplotigs"],accuracyData["results"]["nContigs"],readsClusterAccuracy]
    accuracyLine="\t".join([str(element) for element in resultTable])
    pFile = open(accuracyMetricFile, "a")
    pFile.write(accuracyLine+"\n")
    pFile.close()
    print(accuracyLine)

def getBestMatch(predictedSNPs,groundTruthInfo):
    bestScore=0
    bestScoreDetail=("None",0,0,set())
    for SNP in predictedSNPs:
        chromosome=SNP.split(":")[0]
        break
    for haplotypeName, trueSNPs in groundTruthInfo[chromosome].items():
        TPScore=len(predictedSNPs&trueSNPs)
        posSNPs=set([fullSNP.split("=")[0] for fullSNP in predictedSNPs])
        posTrueSNPs=set([fullSNP.split("=")[0] for fullSNP in trueSNPs])
        commonPos=posSNPs&posTrueSNPs
        FPScore=0
        FPCandidates=predictedSNPs.difference(trueSNPs)
        for FPCandidate in FPCandidates:
            if FPCandidate.split("=")[0] in commonPos:
                FPScore+=1
        commonPosNumber=len(commonPos)
        noiseScore = len(posSNPs.difference(posTrueSNPs))  # Sometimes there are predictions that aren't in the ground truth, not sure how to handle this
        missingScore=commonPosNumber+noiseScore-(TPScore+FPScore) #Adding noiseScore since it's not a TP, or a FP in the sense we're using it, or missing.
        localScore=TPScore/(TPScore+FPScore+missingScore)  # 分母是针对predict的base数
        if localScore>bestScore:
            bestScore=localScore
            bestScoreDetail=(haplotypeName,TPScore,FPScore,commonPos)
    return bestScoreDetail

def subsetTruePhases(groundTruthDict,subsetFile,hybridVCFPrefix):
    VCFInfo=loadVCFInformation(hybridVCFPrefix+subsetFile)
    subsetPositions=set()  # 'NC_001136.10:1318054'
    for chromosome, chromosomeData in VCFInfo.items():  # 模拟数据的vcf snp pos
        for SNPPosition in chromosomeData.keys():
            subsetPositions.add(chromosome+":"+SNPPosition)
    truePhaseSubset={}
    for chromosome, haplotypeData in groundTruthDict.items():  # 以模拟vcf pos为准，筛选ground truth中与其共有的pos，记录下这些pos的ground truth情况
        truePhaseSubset.setdefault(chromosome,{})
        for haplotypeName, SNPSet in haplotypeData.items():
            truePhaseSubset[chromosome].setdefault(haplotypeName,set())
            for fullSNP in SNPSet:
                if fullSNP.split("=")[0] in subsetPositions:  #fullSNP='NC_001133.9:33512=G'
                    truePhaseSubset[chromosome][haplotypeName].add(fullSNP)
    return truePhaseSubset


def loadFloppReadsResults(FloppClusteredReadsFiles):
    ClusteredReads = dict()
    for file in FloppClusteredReadsFiles:
        with open(file) as f:
            for line in f:
                line = line.strip().split("\t")
                if len(line) == 1:
                    temp = int(line[0][1:])
                    cid = "flopp_"+str(temp+1)
                    ClusteredReads.setdefault(cid, set())
                elif len(line) > 1:
                    ClusteredReads[cid].add(line[0])
                else:
                    continue
    return ClusteredReads

if __name__ == "__main__":
    # Load parameters
    benchmarkParameterFilePath = sys.argv[1]
    benchmarkParameters = loadBenchmarkParameterFile(benchmarkParameterFilePath)
    phasingMethods=sys.argv[2:]

    # Handle paths
    mainPath=benchmarkParameters["mainPath"]
    speciesName=benchmarkParameters["reference"]["speciesName"]
    longReadsTag = benchmarkParameters["longReads"]
    referencePath=os.path.join(mainPath, "reference/", speciesName+"/")
    referenceFilePath=os.path.join(referencePath, speciesName+".fasta")
    accuracyMetricPath=os.path.join(mainPath, "accuracyMetrics/")
    accuracyMetricFile=os.path.join(accuracyMetricPath, "accuracyMetrics.tsv")
    VCFPrefix=os.path.join(mainPath, "groundTruth/Variants/shortReads/")
    accuracyCalculationPath=os.path.join(mainPath, "accuracyCalculations/")
    WHPPrefix=os.path.join(mainPath, "phasingPredictions/whatsHapPolyphase/")
    floppPrefix=os.path.join(mainPath, "phasingPredictions/flopp/")
    nPhasePrefix=os.path.join(mainPath, "phasingPredictions/nPhase/")
    PhasingPrefix=os.path.join(mainPath, "phasingPredictions/nTChap")
    hybridVCFPrefix=os.path.join(mainPath, "virtualPolyploids/Variants/shortReads/")
    hybridLongMapPrefix = os.path.join(mainPath, "virtualPolyploids/Mapped/longReads")
    # print(PhasingPrefix)
    allPaths=[accuracyMetricPath,accuracyCalculationPath]

    logPath=os.path.join(mainPath, "log/")
    allPaths.append(logPath)

    for path in allPaths:
        os.makedirs(path, exist_ok=True)

    fullLogPath=os.path.join(logPath, "phasingLog.txt")
    logText=""

    pFile = open(accuracyMetricFile, "a")
    pFile.write("#tool\tploidy\theterozygosity level (%)\tcoverage\tTrue Positives (%)\tFalse Positives (%)\tmissing (%)\tnHaplotigs\tnContigs\treadsClusteredAccuracyRate\n")
    pFile.close()

    max_ID = benchmarkParameters["max_ID"]
    min_len = benchmarkParameters["min_len"] 
    min_ovl = benchmarkParameters["min_ovl"]  
    min_sim = benchmarkParameters["min_sim"] 
    min_overlap = benchmarkParameters["min_overlap"] 
    min_ovl_len = benchmarkParameters["min_ovl_len"] 

    truePhaseDict={} #Will hold the ground truth of individual virtual strains
    allSubsets={} #Will hold the subset ground truth of individual virtual strains
    for strainList in benchmarkParameters["strainLists"]:
        virtualStrainName="_".join(strainList)
        truePhaseGenotypes=[]
        for strainName in strainList:
            strainVCF=os.path.join(VCFPrefix, strainName+".vcf")
            truePhaseGenotypes.append((strainVCF,strainName))
        truePhaseDict[virtualStrainName]=getTruePhases(truePhaseGenotypes) #The phases we get here are full phases but we generate subsets
        #So we want it to only consider the subsets
        for heterozygosityRate in benchmarkParameters["heterozygosityRates"]:
            for coverageLevel in benchmarkParameters["coverages"]:
                subsetFileName=virtualStrainName+"_"+coverageLevel+"X_"+heterozygosityRate+".SNPs.vcf"
                # subsetFileNameIndels=virtualStrainName+"_"+coverageLevel+"X_"+heterozygosityRate+".vcf"
                allSubsets[subsetFileName]=subsetTruePhases(truePhaseDict[virtualStrainName],subsetFileName,hybridVCFPrefix)
                # allSubsets[subsetFileNameIndels]=subsetTruePhases(truePhaseDict[virtualStrainName],subsetFileNameIndels,hybridVCFPrefix)


    if "whatshap-polyphase" in phasingMethods:
        #Translate WH-PP results
        for strainList in benchmarkParameters["strainLists"]:
            virtualStrainName = "_".join(strainList)
            for heterozygosityRate in benchmarkParameters["heterozygosityRates"]:
                for coverage in benchmarkParameters["coverages"]:
                    VCFFileName = virtualStrainName + "_" + coverage + "X_" + heterozygosityRate + ".SNPs.vcf"
                    mappedLongFile = os.path.join(hybridLongMapPrefix, virtualStrainName + "_" + str(coverage) + "X.sorted.bam")
                    phasingMethod="WHP-PP"
                    ploidy=str(len(strainList))
                    if ".SNPs." in VCFFileName:
                        indelStatus="noIndels"
                    else:
                        indelStatus="Indels"
                    testInfo=[phasingMethod,ploidy,heterozygosityRate,coverage]
                    whatsHapPolyphaseResultTranslator(os.path.join(WHPPrefix, "WHP_"+VCFFileName), os.path.join(accuracyCalculationPath, "Translated_WHP_"+VCFFileName))
                    WHPPhasingResults=loadPhasingResults(os.path.join(accuracyCalculationPath, "Translated_WHP_"+VCFFileName))
                    phasedBam = os.path.join(accuracyCalculationPath, "Translated_WHP_"+virtualStrainName + "_" + coverage + "X_" + heterozygosityRate+".sorted.bam")
                    TranslatedClusteredReadsFile = os.path.join(accuracyCalculationPath, "Translated_WHP_"+virtualStrainName + "_" + coverage + "X_" + heterozygosityRate+"_clustered_read_name.tsv")
                    whatsHapPolyphaseResultVariantReads(mappedLongFile, phasedBam, referenceFilePath, ploidy, os.path.join(WHPPrefix, "WHP_"+VCFFileName), TranslatedClusteredReadsFile)
                    WHPPhasingReadsResults = loadClusteredReadName(TranslatedClusteredReadsFile)
                    calculateAccuracyMetrics(allSubsets[VCFFileName],WHPPhasingResults,WHPPhasingReadsResults,longReadsTag,testInfo,accuracyMetricFile)

        print("Done calculating WhatsHap polyphase accuracy")

    if "flopp" in phasingMethods:
        #Translate all flopp results
        for strainList in benchmarkParameters["strainLists"]:
            virtualStrainName = "_".join(strainList)
            for heterozygosityRate in benchmarkParameters["heterozygosityRates"]:
                for coverage in benchmarkParameters["coverages"]:
                    VCFFileName = virtualStrainName + "_" + coverage + "X_" + heterozygosityRate + ".SNPs.vcf"
                    phasingMethod="flopp"
                    ploidy=str(len(strainList))
                    testInfo=[phasingMethod,ploidy,heterozygosityRate,coverage]
                    VCFInformation=loadVCFInformation(os.path.join(hybridVCFPrefix,VCFFileName))
                    floppResultTranslator(os.path.join(floppPrefix, "flopp_"+VCFFileName),os.path.join(accuracyCalculationPath,"Translated_flopp_"+VCFFileName),VCFInformation)
                    floppPredictions=loadPhasingResults(os.path.join(accuracyCalculationPath, "Translated_flopp_"+VCFFileName))
                    floppClusteredReadsPath = os.path.join(floppPrefix, virtualStrainName + "_" + coverage + "X_" + heterozygosityRate)
                    floppClusteredReadsFiles = glob.glob(os.path.join(floppClusteredReadsPath, "*_part.txt"))
                    floppPhasingReadsResults = loadFloppReadsResults(floppClusteredReadsFiles)
                    calculateAccuracyMetrics(allSubsets[VCFFileName], floppPredictions,floppPhasingReadsResults,longReadsTag,testInfo,accuracyMetricFile)

        print("Done calculating flopp accuracy")

    if "nphase" in phasingMethods:
        nPhaseDefaultSuffix="_0.1_0.01_0.05_0_variants.tsv"
        nPhaseDefaultReadsSuffix="_0.1_0.01_0.05_0_clusterReadNames.tsv"
        for strainList in benchmarkParameters["strainLists"]:
            virtualStrainName = "_".join(strainList)
            for heterozygosityRate in benchmarkParameters["heterozygosityRates"]:
                for coverage in benchmarkParameters["coverages"]:
                    exactTestName=virtualStrainName + "_" + coverage + "X_" + heterozygosityRate
                    VCFFileName = exactTestName+nPhaseDefaultSuffix
                    phasingMethod="nPhase"
                    ploidy=str(len(strainList))
                    testInfo=[phasingMethod,ploidy,heterozygosityRate,coverage]
                    nPhaseResultPath = os.path.join(nPhasePrefix, exactTestName, "Phased/", VCFFileName)
                    nPhaseReadsResultPath = os.path.join(nPhasePrefix, exactTestName, "Phased/", exactTestName+nPhaseDefaultReadsSuffix)
                    subsetFileName=exactTestName+".SNPs.vcf"
                    nPhasePhasingResults = loadPhasingResults(nPhaseResultPath)
                    nPhasePhasingReadsResults = loadnPhaseClusteredReadName(nPhaseReadsResultPath)
                    calculateAccuracyMetrics(allSubsets[subsetFileName],nPhasePhasingResults,nPhasePhasingReadsResults,longReadsTag,testInfo,accuracyMetricFile)

        print("Done calculating nPhase accuracy")

    if "ntchap" in phasingMethods:
        PhasingDefaultSuffix="_"+str(min_ovl)+"_"+str(min_sim)+"_"+str(max_ID)+"_"+str(min_len)+"_"+str(min_overlap)+"_"+str(min_ovl_len)+"_variants.tsv"
        PhasingDefaultReadsSuffix="_"+str(min_ovl)+"_"+str(min_sim)+"_"+str(max_ID)+"_"+str(min_len)+"_"+str(min_overlap)+"_"+str(min_ovl_len)+"_clustered_read_name.tsv"
        for strainList in benchmarkParameters["strainLists"]:
            virtualStrainName = "_".join(strainList)
            for heterozygosityRate in benchmarkParameters["heterozygosityRates"]:
                for coverage in benchmarkParameters["coverages"]:
                    exactTestName=virtualStrainName + "_" + coverage + "X_" + heterozygosityRate
                    """ if run ./phased_tool_runners_SaCere.py, 
                    PhasingDefaultSuffix="_"+testName+"_variants.tsv"
                    PhasingDefaultReadsSuffix="_"+testName+"_clustered_read_name.tsv" """
                    VCFFileName = "phased"+PhasingDefaultSuffix
                    phasingMethod="Phasing_"+str(min_ovl)+"_"+str(min_sim)+"_"+str(max_ID)+"_"+str(min_len)+"_"+str(min_overlap)+"_"+str(min_ovl_len)
                    ploidy=str(len(strainList))
                    testInfo=[phasingMethod,ploidy,heterozygosityRate,coverage]
                    PhasingResultPath = os.path.join(PhasingPrefix, exactTestName, "Phased/", VCFFileName)
                    PhasingReadsResultPath = os.path.join(PhasingPrefix, exactTestName, "Phased/", "phased"+PhasingDefaultReadsSuffix)
                    subsetFileName=exactTestName+".SNPs.vcf"
                    PhasingResults = loadPhasingResults(PhasingResultPath)
                    PhasingReadsResults = loadClusteredReadName(PhasingReadsResultPath)
                    calculateAccuracyMetrics(allSubsets[subsetFileName],PhasingResults,PhasingReadsResults,longReadsTag,testInfo,accuracyMetricFile)

        print("Done calculating Phasing accuracy")

    if "Phasing_clean" in phasingMethods:
        PhasingDefaultSuffix="_"+str(min_ovl)+"_"+str(min_sim)+"_"+str(max_ID)+"_"+str(min_len)+"_"+str(min_overlap)+"_"+str(min_ovl_len)+"_variants.tsv"
        PhasingDefaultReadsSuffix="_"+str(min_ovl)+"_"+str(min_sim)+"_"+str(max_ID)+"_"+str(min_len)+"_"+str(min_overlap)+"_"+str(min_ovl_len)+"_clustered_read_name.tsv"
        # PhasingDefaultSuffixIndel="_Indels_0.1_0.01_0.05_0_variants.tsv"
        # clean_ways = ["short", "short_filter", "short_filter_stitch"]
        clean_ways = ["consensus_short"]
        for strainList in benchmarkParameters["strainLists"]:
            virtualStrainName = "_".join(strainList)
            for heterozygosityRate in benchmarkParameters["heterozygosityRates"]:
                for coverage in benchmarkParameters["coverages"]:
                    for clean_way in clean_ways:
                        exactTestName=virtualStrainName + "_" + coverage + "X_" + heterozygosityRate
                        print(exactTestName)
                        ploidy=str(len(strainList))
                        subsetFileName=exactTestName+".SNPs.vcf"
                        Clean_VCFFileName = "phased_cleaned_"+clean_way+PhasingDefaultSuffix
                        Clean_phasingMethod="Phasing_"+clean_way+str(min_ovl)+"_"+str(min_sim)+"_"+str(max_ID)+"_"+str(min_len)+"_"+str(min_overlap)+"_"+str(min_ovl_len)
                        Clean_testInfo=[Clean_phasingMethod,ploidy,heterozygosityRate,coverage,"noIndels"]
                        Clean_PhasingResultPath = os.path.join(PhasingPrefix, exactTestName, "Phased/Cleaned", Clean_VCFFileName)
                        Clean_PhasingReadsResultPath = os.path.join(PhasingPrefix, exactTestName, "Phased/Cleaned", "phased_cleaned_"+clean_way+PhasingDefaultReadsSuffix)
                        Clean_PhasingResults = loadPhasingResults(Clean_PhasingResultPath)
                        Clean_PhasingReadsResults = loadClusteredReadName(Clean_PhasingReadsResultPath)
                        calculateAccuracyMetrics(allSubsets[subsetFileName],Clean_PhasingResults,Clean_PhasingReadsResults,longReadsTag,Clean_testInfo,accuracyMetricFile)

        print("Done calculating Phasing cleaning accuracy")

    if "PhasingAndClean" in phasingMethods:
        PhasingDefaultSuffix="_"+str(min_ovl)+"_"+str(min_sim)+"_"+str(max_ID)+"_"+str(min_len)+"_"+str(min_overlap)+"_"+str(min_ovl_len)+"_variants.tsv"
        PhasingDefaultReadsSuffix="_"+str(min_ovl)+"_"+str(min_sim)+"_"+str(max_ID)+"_"+str(min_len)+"_"+str(min_overlap)+"_"+str(min_ovl_len)+"_clustered_read_name.tsv"
        # clean_ways = ["short", "short_filter", "short_filter_stitch"]
        clean_ways = ["consensus_short"]
        # PhasingDefaultSuffixIndel="_Indels_0.1_0.01_0.05_0_variants.tsv"
        for strainList in benchmarkParameters["strainLists"]:
            virtualStrainName = "_".join(strainList)
            for heterozygosityRate in benchmarkParameters["heterozygosityRates"]:
                for coverage in benchmarkParameters["coverages"]:
                    exactTestName=virtualStrainName + "_" + coverage + "X_" + heterozygosityRate
                    Phasing_VCFFileName = "phased"+PhasingDefaultSuffix
                    # VCFFileNameIndels = exactTestName +PhasingDefaultSuffixIndel
                    Phasing_phasingMethod="Phasing_"+str(min_ovl)+"_"+str(min_sim)+"_"+str(max_ID)+"_"+str(min_len)+"_"+str(min_overlap)+"_"+str(min_ovl_len)
                    ploidy=str(len(strainList))
                    #Without indels
                    Phasing_testInfo=[Phasing_phasingMethod,ploidy,heterozygosityRate,coverage,"noIndels"]
                    PhasingResultPath = os.path.join(PhasingPrefix, exactTestName, "Phased/", Phasing_VCFFileName)
                    PhasingReadsResultPath = os.path.join(PhasingPrefix, exactTestName, "Phased/", "phased"+PhasingDefaultReadsSuffix)
                    subsetFileName=exactTestName+".SNPs.vcf"
                    PhasingResults = loadPhasingResults(PhasingResultPath)
                    PhasingReadsResults = loadClusteredReadName(PhasingReadsResultPath)
                    calculateAccuracyMetrics(allSubsets[subsetFileName],PhasingResults,PhasingReadsResults,longReadsTag,Phasing_testInfo,accuracyMetricFile)

                    for clean_way in clean_ways:
                        Clean_VCFFileName = "phased_cleaned_"+clean_way+PhasingDefaultSuffix
                        Clean_phasingMethod="Phasing_"+clean_way+str(min_ovl)+"_"+str(min_sim)+"_"+str(max_ID)+"_"+str(min_len)+"_"+str(min_overlap)+"_"+str(min_ovl_len)
                        Clean_testInfo=[Clean_phasingMethod,ploidy,heterozygosityRate,coverage,"noIndels"]
                        Clean_PhasingResultPath = os.path.join(PhasingPrefix, exactTestName, "Phased/Cleaned", Clean_VCFFileName)
                        Clean_PhasingReadsResultPath = os.path.join(PhasingPrefix, exactTestName, "Phased/Cleaned", "phased_cleaned_"+clean_way+PhasingDefaultReadsSuffix)
                        Clean_PhasingResults = loadPhasingResults(Clean_PhasingResultPath)
                        Clean_PhasingReadsResults = loadClusteredReadName(Clean_PhasingReadsResultPath)
                        calculateAccuracyMetrics(allSubsets[subsetFileName],Clean_PhasingResults,Clean_PhasingReadsResults,longReadsTag,Clean_testInfo,accuracyMetricFile)

                    # # 计算每次迭代的准确度，ncontigs等指标
                    # itr_dir_name = itr_dir_name = "phased"+"_"+str(min_ovl)+"_"+str(min_sim)+"_"+str(max_ID)+"_"+str(min_len)+"_"+str(min_overlap)+"_"+str(min_ovl_len)
                    # itr_files = glob.glob(os.path.join(PhasingPrefix, exactTestName, "Phased", itr_dir_name, itr_dir_name+"_*"))
                    # itr_num = int(len(itr_files)/2)
                    # for i in range(itr_num):
                    #     PhasingItrSuffix = "_"+str(min_ovl)+"_"+str(min_sim)+"_"+str(max_ID)+"_"+str(min_len)+"_"+str(min_overlap)+"_"+str(min_ovl_len)+"_itr"+str(i+1)+"_variants.tsv"
                    #     PhasingItrReadsSuffix="_"+str(min_ovl)+"_"+str(min_sim)+"_"+str(max_ID)+"_"+str(min_len)+"_"+str(min_overlap)+"_"+str(min_ovl_len)+"_itr"+str(i+1)+"_clustered_read_name.tsv"
                    #     Phasing_VCFFileName = "phased"+PhasingItrSuffix
                    #     # VCFFileNameIndels = exactTestName +PhasingDefaultSuffixIndel
                    #     Phasing_phasingMethod="Phasing_"+str(min_ovl)+"_"+str(min_sim)+"_"+str(max_ID)+"_"+str(min_len)+"_"+str(min_overlap)+"_"+str(min_ovl_len)+"_itr"+str(i+1)
                    #     #Without indels
                    #     Phasing_testInfo=[Phasing_phasingMethod,ploidy,heterozygosityRate,coverage,"noIndels"]
                    #     PhasingResultPath = os.path.join(PhasingPrefix, exactTestName, "Phased/", itr_dir_name, Phasing_VCFFileName)
                    #     PhasingReadsResultPath = os.path.join(PhasingPrefix, exactTestName, "Phased/", itr_dir_name, "phased"+PhasingItrReadsSuffix)
                    #     PhasingResults = loadPhasingResults(PhasingResultPath)
                    #     PhasingReadsResults = loadClusteredReadName(PhasingReadsResultPath)
                    #     calculateAccuracyMetrics(allSubsets[subsetFileName],PhasingResults,PhasingReadsResults,longReadsTag,Phasing_testInfo,accuracyMetricFile)

                    #     for clean_way in clean_ways:
                    #         Clean_VCFFileName = "phased_cleaned_"+clean_way+PhasingItrSuffix
                    #         Clean_phasingMethod="Phasing_"+clean_way+" "*9
                    #         Clean_testInfo=[Clean_phasingMethod,ploidy,heterozygosityRate,coverage,"noIndels"]
                    #         Clean_PhasingResultPath = os.path.join(PhasingPrefix, exactTestName, "Phased/Cleaned", itr_dir_name, Clean_VCFFileName)
                    #         Clean_PhasingReadsResultPath = os.path.join(PhasingPrefix, exactTestName, "Phased/Cleaned", itr_dir_name,"phased_cleaned_"+clean_way+PhasingItrReadsSuffix)
                    #         Clean_PhasingResults = loadPhasingResults(Clean_PhasingResultPath)
                    #         Clean_PhasingReadsResults = loadClusteredReadName(Clean_PhasingReadsResultPath)
                    #         calculateAccuracyMetrics(allSubsets[subsetFileName],Clean_PhasingResults,Clean_PhasingReadsResults,longReadsTag,Clean_testInfo,accuracyMetricFile)


        print("Done calculating Phasing and Cleaning accuracy")
    
    print("Done calculating accuracy metrics")