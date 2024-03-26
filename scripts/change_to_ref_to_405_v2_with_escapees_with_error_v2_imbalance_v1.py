import pysam
import matplotlib.pyplot as plt
import random

samfile = pysam.AlignmentFile("possorted_baq_ne_filtered_genome_bam.bam", "rb")
vcffile = pysam.VariantFile("extracted_403_to_405.bcf", "rb")
escapeefile = open("escapees_sorted.txt","r")
outfile = open("modified_reads.fastq", "w")
outfile_barcode = open("barcodes.fastq", "w")
logfile = open("change_log.txt", "w")
logfile_2 = open("cell_stat.log", "w")

SEQUENCING_ERROR_RATE = 0.0005
FATHER_ODDS = 0.5

escape_pos_list =[] 
for line in escapeefile.readlines():
    sp = line.strip().split("\t")
    tag = random.randrange(3)
    if len(sp)>=4 and tag == 0 and sp[3] == "confirmed":
        escape_pos_list.append([sp[0], int(sp[1]), int(sp[2])])

num_change_to_father = 0
num_change_to_mother = 0
base_dict = {}
cell_log_dict = {}
i = 0
for read in samfile.fetch():
    #if i>10:
    #    break
    ref_name = read.reference_name
    #print(ref_name)
    ref_read = read.get_reference_sequence()

    # Add random sequencing error on reads
    rand_mut_pos_i = 0
    rand_mut_pos = 0
    for rand_mut in ref_read:
        result = random.choices([0,1], weights=[SEQUENCING_ERROR_RATE, 1-SEQUENCING_ERROR_RATE], k=1)
        # if it is sequencing error
        if result[0] == 0:            
            base_choosed = random.choice(["A","C","G","T"])
            ref_read = ref_read[:rand_mut_pos_i]+base_choosed+ref_read[rand_mut_pos_i+1:]
        rand_mut_pos = rand_mut_pos_i
        rand_mut_pos_i = rand_mut_pos_i + 1

    #print(read)
    #print(ref_read)
    pos_list = read.get_reference_positions()
    #print(pos_list)
    if read.has_tag("CB"):
        cell_barcode = read.get_tag("CB")
    else:
        continue
    #print("Cell barcode: "+cell_barcode)

    if cell_barcode in cell_log_dict.keys():
        cell_log_dict[cell_barcode][0] = cell_log_dict[cell_barcode][0] + 1
    else:
        cell_log_dict[cell_barcode] = [1, 0]

    for snp in vcffile.fetch(ref_name, pos_list[0], pos_list[-1]):
        # make sure it is all single base mutation
        single_base_flag = 1
        alleles = snp.alleles
        for allele in alleles:
            if len(allele) != 1:
                single_base_flag = 0
                break
        if single_base_flag == 0:
            continue

        # begin to check and change
        zero_based_pos = snp.pos - 1
        if zero_based_pos in pos_list:
            mut_pos = zero_based_pos - pos_list[0]
            if mut_pos == rand_mut_pos:
                continue
            father_tuple = snp.samples["HG00403"]["GT"]
            mother_tuple = snp.samples["HG00404"]["GT"]
            child_tuple = snp.samples["HG00405"]["GT"]

            homozygous_condition = (((len(father_tuple)==1) or (len(father_tuple) == 2 and father_tuple[0] == father_tuple[1])) and (len(mother_tuple) == 2 and mother_tuple[0] == mother_tuple[1]) and len(child_tuple) == 2) and father_tuple[0] != mother_tuple[0] and ((child_tuple[0] == father_tuple[0] and child_tuple[1] == mother_tuple[0]) or (child_tuple[0] == mother_tuple[0] and child_tuple[1] == father_tuple[0]))

            father_set = set()
            if len(father_tuple) == 1:
                father_set.add(father_tuple[0])
            elif len(father_tuple) == 2:
                father_set.add(father_tuple[0])
                father_set.add(father_tuple[1])
            mother_set = set()
            mother_set.add(mother_tuple[0])
            mother_set.add(mother_tuple[1])
            symmetric_diffrence_set = father_set.symmetric_difference(mother_set)
            heterozygous_condition_1 = (child_tuple[0] in father_set.intersection(symmetric_diffrence_set) and child_tuple[1] in mother_set.intersection(symmetric_diffrence_set))
            heterozygous_condition_2 = (child_tuple[0] in mother_set.intersection(symmetric_diffrence_set) and child_tuple[1] in father_set.intersection(symmetric_diffrence_set))

            if heterozygous_condition_1 or heterozygous_condition_2:
            # check whether the snp is on an escapee gene
                escapee_flag = 0
                for pos_range in escape_pos_list:
                    if zero_based_pos >= pos_range[1]-1 and zero_based_pos <= pos_range[2]-1:
                        escapee_flag = 1
                        break
                # if is escapee gene, randomly inherited from father or mother
                if escapee_flag:
                    inherited_from = random.choice(["father", "mother"])
                    print("Pos in escapee gene: "+str(zero_based_pos))
                else:
                    if cell_barcode in base_dict.keys():
                        inherited_from = base_dict[cell_barcode]
                    else:
                        random_result = random.choices(["father", "mother"], weights=[FATHER_ODDS, 1-FATHER_ODDS], k=1)
                        inherited_from = random_result[0]
                        base_dict[cell_barcode] = inherited_from

                cell_log_dict[cell_barcode][1] = cell_log_dict[cell_barcode][1] + 1

                
                if inherited_from == "father" and snp.alleles[father_tuple[0]] != "*" and snp.alleles[mother_tuple[0]] != "*":
                    if heterozygous_condition_1:
                        ref_read = ref_read[:mut_pos]+snp.alleles[child_tuple[0]]+ref_read[mut_pos+1:]
                        print("Pos "+str(zero_based_pos)+" - Father alleles: "+str(father_tuple)+", mother alleles: "+str(mother_tuple) + ", child alleles: "+str(child_tuple)+". Changed to father.")
                        if len(father_tuple) == 1:
                            logfile.write(cell_barcode+"\t"+str(zero_based_pos)+"\t"+snp.alleles[father_tuple[0]]+"\t"+snp.alleles[mother_tuple[0]]+snp.alleles[mother_tuple[1]]+"\t"+snp.alleles[child_tuple[0]]+"\n")
                        elif len(father_tuple) == 2:
                            logfile.write(cell_barcode+"\t"+str(zero_based_pos)+"\t"+snp.alleles[father_tuple[0]]+snp.alleles[father_tuple[1]]+"\t"+snp.alleles[mother_tuple[0]]+snp.alleles[mother_tuple[1]]+"\t"+snp.alleles[child_tuple[0]]+"\n")
                    elif heterozygous_condition_2:
                        ref_read = ref_read[:mut_pos]+snp.alleles[child_tuple[1]]+ref_read[mut_pos+1:]
                        print("Pos "+str(zero_based_pos)+" - Father alleles: "+str(father_tuple)+", mother alleles: "+str(mother_tuple) + ", child alleles: "+str(child_tuple)+". Changed to father.")
                        if len(father_tuple) == 1:
                            logfile.write(cell_barcode+"\t"+str(zero_based_pos)+"\t"+snp.alleles[father_tuple[0]]+"\t"+snp.alleles[mother_tuple[0]]+snp.alleles[mother_tuple[1]]+"\t"+snp.alleles[child_tuple[1]]+"\n")
                        elif len(father_tuple) == 2:
                            logfile.write(cell_barcode+"\t"+str(zero_based_pos)+"\t"+snp.alleles[father_tuple[0]]+snp.alleles[father_tuple[1]]+"\t"+snp.alleles[mother_tuple[0]]+snp.alleles[mother_tuple[1]]+"\t"+snp.alleles[child_tuple[1]]+"\n")
                    num_change_to_father = num_change_to_father + 1
                elif inherited_from == "mother" and snp.alleles[father_tuple[0]] != "*" and snp.alleles[mother_tuple[0]] != "*":
                    if heterozygous_condition_1:
                        ref_read = ref_read[:mut_pos]+snp.alleles[child_tuple[1]]+ref_read[mut_pos+1:]
                        print("Pos "+str(zero_based_pos)+" - Father alleles: "+str(father_tuple)+", mother alleles: "+str(mother_tuple) + ", child alleles: "+str(child_tuple)+". Changed to mother.")
                        if len(father_tuple) == 1:
                            logfile.write(cell_barcode+"\t"+str(zero_based_pos)+"\t"+snp.alleles[father_tuple[0]]+"\t"+snp.alleles[mother_tuple[0]]+snp.alleles[mother_tuple[1]]+"\t"+snp.alleles[child_tuple[1]]+"\n")
                        elif len(father_tuple) == 2:
                            logfile.write(cell_barcode+"\t"+str(zero_based_pos)+"\t"+snp.alleles[father_tuple[0]]+snp.alleles[father_tuple[1]]+"\t"+snp.alleles[mother_tuple[0]]+snp.alleles[mother_tuple[1]]+"\t"+snp.alleles[child_tuple[1]]+"\n")
                    elif heterozygous_condition_2:
                        ref_read = ref_read[:mut_pos]+snp.alleles[child_tuple[0]]+ref_read[mut_pos+1:]
                        print("Pos "+str(zero_based_pos)+" - Father alleles: "+str(father_tuple)+", mother alleles: "+str(mother_tuple) + ", child alleles: "+str(child_tuple)+". Changed to mother.")
                        if len(father_tuple) == 1:
                            logfile.write(cell_barcode+"\t"+str(zero_based_pos)+"\t"+snp.alleles[father_tuple[0]]+"\t"+snp.alleles[mother_tuple[0]]+snp.alleles[mother_tuple[1]]+"\t"+snp.alleles[child_tuple[0]]+"\n")
                        elif len(father_tuple) == 2:
                            logfile.write(cell_barcode+"\t"+str(zero_based_pos)+"\t"+snp.alleles[father_tuple[0]]+snp.alleles[father_tuple[1]]+"\t"+snp.alleles[mother_tuple[0]]+snp.alleles[mother_tuple[1]]+"\t"+snp.alleles[child_tuple[0]]+"\n")
                    num_change_to_mother = num_change_to_mother + 1
    outfile.write("@"+read.query_name+"\n")
    outfile.write(ref_read.upper()+"\n")
    outfile.write("+\n")
    outfile.write("F"*len(ref_read)+"\n")

    umi = read.get_tag("UR")
    outfile_barcode.write("@"+read.query_name+"\n")
    outfile_barcode.write(cell_barcode.upper()[:-2]+umi.upper()+"\n")
    outfile_barcode.write("+\n")
    outfile_barcode.write("F"*(len(umi)+len(cell_barcode)-2)+"\n")

    i=i+1
print("Base changed to father allele: " + str(num_change_to_father))
print("Base changed to mother allele: " + str(num_change_to_mother))

def plotxy(x, y):
    plt.plot(x, y, 'og')
    plt.xlabel("Reads number in cell")
    plt.ylabel("Number of SNPs")
    plt.title("SNPs vs reads number")
    plt.savefig("SNPs_vs_reads_number.pdf")
    plt.savefig("SNPs_vs_reads_number.png")

x = []
y = []

for key in cell_log_dict.keys():
    x.append(cell_log_dict[key][0])
    y.append(cell_log_dict[key][1])
    logfile_2.write(key+"\t"+str(cell_log_dict[key][0])+"\t"+str(cell_log_dict[key][1])+"\n")

plotxy(x,y)

samfile.close()
vcffile.close()
escapeefile.close()
outfile.close()
outfile_barcode.close()
logfile.close()
logfile_2.close()
