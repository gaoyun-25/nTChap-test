import os
import sys
import subprocess
import json

import addGTtoVCF

def loadBenchmarkParameterFile(parameterFilePath):
    benchmarkParameters = {"mainPath": "",
                           "longReadMethod": [],
                           "shortReadMethod": "",
                           "coverages": [],
                           "reference": {},
                           "haplogeneratorModel": [],
                           "Ploidy": [],
                           "mutation":[],
                           "dosage":[],
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
        if line[0] in {"mainPath","shortReadMethod", "genomeSize","threads","max_ID","min_len","min_ovl","min_sim","min_overlap","min_ovl_len"}:
            benchmarkParameters[line[0]]=line[1]
        elif line[0] in {"longReadMethod","coverages","Ploidy","mutation","dosage","haplogeneratorModel"}:
            for dataPoint in line[1:]:
                benchmarkParameters[line[0]].append(dataPoint)
        elif line[0] in {"reference"}:
            for dataPoint in line[1:]:
                dataPoint=dataPoint.split(":")
                benchmarkParameters[line[0]][dataPoint[0]]=dataPoint[1]
        else:
            print("Ignoring '","\t".join(line),"' because it is not recognized")
    benchmarkParameterFile.close()
    return benchmarkParameters


execute = True
def call(s,execute = execute, check_code = True):
    if execute:
        code = subprocess.call(s,shell=True)
        if code > 0 and check_code:
            error_string = "Return code from command %s is not zero. Exiting!" %(s)
            print(error_string)
            exit()
    else:
        print(s)


if __name__ == "__main__":
    benchmarkParameterFilePath = sys.argv[1]
    benchmarkParameters = loadBenchmarkParameterFile(benchmarkParameterFilePath)
    my_code_path = sys.argv[2]

    main_path=benchmarkParameters["mainPath"]
    ref_name=benchmarkParameters["reference"]["speciesName"]
    threads=benchmarkParameters["threads"]
    Ploidy=benchmarkParameters["Ploidy"]
    haplogeneratorModels = benchmarkParameters["haplogeneratorModel"]
    mutations = benchmarkParameters["mutation"]
    dosages = benchmarkParameters["dosage"]
    coverages = benchmarkParameters["coverages"]
    longReadMethod = benchmarkParameters["longReadMethod"]

    hybridPath=os.path.join(main_path, "virtualPolyploids/")
    GTPrefix = os.path.join(main_path, "groundTruth/")
    referencePath = os.path.join(main_path,"reference/",ref_name)

    variantCalledShortReads=os.path.join(hybridPath, "Variants/")
    variantCalledShortReadsFilter=os.path.join(hybridPath, "Variants_filter/")
    hybridRawReads=os.path.join(hybridPath, "rawReads/")
    hybridLongReads=os.path.join(hybridRawReads, "longReads/")
    hybridShortReads=os.path.join(hybridRawReads, "shortReads/")
    mappedLongReads=os.path.join(hybridPath, "Mapped/longReads/")
    mappedShortReads=os.path.join(hybridPath, "Mapped/shortReads/")

    all_paths = [hybridPath, GTPrefix,referencePath, variantCalledShortReads,variantCalledShortReadsFilter,hybridLongReads,hybridShortReads,mappedLongReads,mappedShortReads]
    for path in all_paths:
        os.makedirs(path, exist_ok=True)

    ref_path=os.path.join(referencePath,ref_name+".fa")

    ## Random standard mutation dict
    mut_dict = "\"{'A':'C','G':'A','C':'T','T':'G'}\""
    ## Polyalleic case
    mut_dict_poly = "\"{'A':'C', 'A':'G' ,'G':'A','C':'T','T':'G'}\""

    s = 'samtools faidx %s' %(ref_path)
    call(s)

    s = 'bwa index %s' %(ref_path)
    call(s)

    get_hap_files = True
    simulate_reads = True
    gen_nanopore_reads = True
    get_vcf = True

    for haplogeneratorModel in haplogeneratorModels:
        if haplogeneratorModel == "Poisson":
            for mut in mutations:
                for index in range(len(Ploidy)):
                    ploidy = int(Ploidy[index])
                    virtualName = ref_name+"_p"+str(ploidy)+"_mt"+str(mut)
                    GTPath = os.path.join(GTPrefix, virtualName)

                    all_paths = [GTPath]
                    for path in all_paths:
                        os.makedirs(path, exist_ok=True)

                    if get_hap_files:
                        dosage = dosages[index]
                        s_mut = "\"["+str(mut)+",0,0"+"]\""
                        mut_dict_str = "\""+json.dumps(mut_dict)+"\""
                        # 用shell脚本跑主要是conda切换 haplogenerator.py(py2.7--Haplosim_py2)
                        # -o out_name -s 此处是lognormal的mean，--sdlog lognormal的STD
                        # hap_out_name = os.path.join(GTPath,virtualName)
                        simulate_haplotypes_code = os.path.join(my_code_path, "run_haplogenerator.sh")
                        s = "%s %s %s/%s \"[%s,0,0]\" %s %s %s" % (simulate_haplotypes_code,ref_path,GTPath,virtualName,mut,ploidy,mut_dict,dosage)
                        call(s)
                        os.system("sh "+simulate_haplotypes_code+" "+ref_path+" "+os.path.join(GTPath,virtualName)+" "+s_mut+" "+str(ploidy)+" "+mut_dict_str+" "+dosage)

                        # 获得由ground truth得到的vcf文件，而不是由short reads call所得
                        get_vcf_from_haplotypes = os.path.join(my_code_path, "get_vcf_from_haplotypes.py")
                        vcf_outfile = os.path.join(variantCalledShortReads, virtualName+".vcf")
                        s = 'python %s %s %s %s/*.fa' % (get_vcf_from_haplotypes, ref_path, vcf_outfile, GTPath)
                        call(s)

                    if simulate_reads:
                        # 需要改表头的名称，方便知道生成的reads是来自哪个chr和haplotype的
                        for p in range(ploidy):
                        # 改变表头，影响short reads name,long reads name后的备注会有hap_i#chr,
                        # 改一下模拟的haplotypes表头名字，方便后续模拟reads名都按照各自haplotype名字，便于计算clustering accuracy
                            haplo_ref = os.path.join(GTPath, virtualName+"_hap"+str(p+1)+".fa")
                            haplo_ref_temp = os.path.join(GTPath, virtualName+"_hap"+str(p+1)+".fa.Temp")
                            with open(haplo_ref, "r") as fr, open(haplo_ref_temp, "w") as fw:
                                for line in fr:
                                    if line.startswith(">"):
                                        header = line[1:]
                                        fw.write(">hap_"+str(p+1)+"#"+header)
                                    else:
                                        fw.write(line)
                            os.remove(haplo_ref)
                            os.rename(haplo_ref_temp, haplo_ref)
                            s = "samtools faidx %s" %(haplo_ref)
                            call(s)

                        # 模拟长reads
                    if gen_nanopore_reads:
                        for coverage in coverages:
                            longfqgz=""
                            for p in range(ploidy):
                                haplo_ref = os.path.join(GTPath, virtualName+"_hap"+str(p+1)+".fa")
                                hap_long_reads = os.path.join(hybridLongReads, virtualName+"_hap"+str(p+1)+"_"+str(coverage)+"X_ont.fastq.gz")
                                s = 'badread-runner.py simulate --reference %s --quantity %sx  | gzip > %s' % (haplo_ref,coverage,hap_long_reads)
                                call(s)

                                longfqgz = longfqgz + " " + hap_long_reads

                            longReadsFq = os.path.join(hybridLongReads,virtualName+"_"+str(coverage)+"X_ont.fastq.gz")
                            s = "zcat %s | gzip - > %s" % (longfqgz,longReadsFq)
                            call(s)
                            for p in range(ploidy):
                                hap_long_reads = os.path.join(hybridLongReads, virtualName+"_hap"+str(p+1)+"_"+str(coverage)+"X_ont.fastq.gz")
                                os.remove(hap_long_reads)

                            longReadsMapPrefix = os.path.join(mappedLongReads,virtualName+"_"+str(coverage)+"X_ont")
                            s = "minimap2 -ax map-ont %s %s -o %s.sam " % (ref_path, longReadsFq,longReadsMapPrefix)
                            call(s)
                            s = "samtools view -@ %s -bS %s.sam -o %s.bam " % (threads, longReadsMapPrefix,longReadsMapPrefix)
                            call(s)
                            s = "samtools sort -@ %s %s.bam -o %s.sorted.bam" % (threads,longReadsMapPrefix,longReadsMapPrefix)
                            call(s)
                            s = "samtools index %s.sorted.bam" %(longReadsMapPrefix)
                            call(s)
                            s = "samtools view -@ %s %s.sorted.bam -o %s.sorted.sam" % (threads, longReadsMapPrefix,longReadsMapPrefix)
                            call(s)

                    if get_vcf:
                        caller = "lofreq"
                        variantsFaVCF = os.path.join(variantCalledShortReads, virtualName+".vcf")
                        for coverage in coverages:
                            longReadsMapPrefix = os.path.join(mappedLongReads,virtualName+"_"+str(coverage)+"X_ont")
                            variantsVCF = os.path.join(variantCalledShortReadsFilter,virtualName+"_"+str(coverage)+"X_ont_"+caller+".vcf")
                            os.system("lofreq call -B -f "+ref_path+" "+longReadsMapPrefix+".sorted.bam"+" -o "+variantsVCF+" --force-overwrite")
                            if caller == "lofreq":
                                variantsFilterSNPVCF=os.path.join(variantCalledShortReadsFilter,virtualName+"_"+str(coverage)+"X_ont_"+caller+"_filter.vcf")
                                addGTtoVCF.addGTtoVCF(variantsFaVCF, variantsVCF, variantsFilterSNPVCF)

