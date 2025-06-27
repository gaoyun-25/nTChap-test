
def loadVCF(vcf_file):
    variants_pos = set()
    # vcf_dict["headers"] = list()
    with open(vcf_file) as f:
        for line in f:
            if line[0] != "#":
                # vcf_dict["headers"].append(line)

                line = line.strip("\n").split("\t")
                chromosome=line[0]
                position=line[1]
                SNPs=[line[3]]+line[4].split(",")
                # if len(line) > 7:
                #     GT_tag = line[9]
                #     vcf_dict.setdefault(chromosome,{})
                #     vcf_dict[chromosome][position] = {"SNPs":SNPs,"GT":GT_tag}
                # else:
                variants_pos.add(chromosome+":"+position)

    return variants_pos

def addGTtoVCF(vcf1,lofreq_vcf,out_vcf):
    # 根据lofreq提取的位点筛选haplotypes比较得到的vcf，相同chr，pos，ref，alt保留
    lofreq_pos = loadVCF(lofreq_vcf)
    with open(vcf1) as f,open(out_vcf,"w") as g:
        for line in f:
            if line[0] == "#":
                g.write(line)
            else:
                item = line.strip("\n").split("\t")
                chromosome=item[0]
                position=item[1]
                if chromosome+":"+position in lofreq_pos:
                    g.write(line)

