#encoding:utf-8

# 从每个样本检测到的TR结果里面取出拷贝数，support_reads, 与target motif的编辑距离
# 输出clean data

# 2022-0727 vcf文件的坐标已经改成1-based

import sys
import numpy as np
import argparse
import logging
import random
import numpy as np
from collections import OrderedDict

def get_base_score(base1, base2, match=0, mismatch=1):
    '''
        若两个碱基相同，则返回match得分，若不同则返回mismtch得分
    '''
    return match if base1 == base2 else mismatch

def cal_edit_distance(seq1, seq2, match=0, mismatch=1, indel=1):
    '''
        仅在score_matrix[j,i]这一行与序列比对不同，在这里是取最小值
        可参考该链接：https://users.monash.edu/~lloyd/tildeAlgDS/Dynamic/Edit/#:~:text=The%20simple%20edit%20distance%20algorithm%20would%20normally%20be,a%20distance%20between%20DNA%20sequences%20%28strings%20over%20%7BA%2CC%5D
    '''
    seq1_len = len(seq1)
    seq2_len = len(seq2)
    score_matrix = np.zeros((seq2_len+1, seq1_len+1))

    # 将第一行与第一列置为indel的得分
    for i in range(1, seq1_len+1):
        score_matrix[0, i] = score_matrix[0, i-1] + indel

    for i in range(1, seq2_len+1):
        score_matrix[i, 0] = score_matrix[i-1, 0] + indel

    # 计算score_matrix的其他分值
    for j in range(1, seq2_len+1):
        for i in range(1, seq1_len+1):
            # print(i,j)
            base1 = seq1[i-1]
            base2 = seq2[j-1]
            base_score = get_base_score(base1, base2, match, mismatch)
            mis_up = score_matrix[j-1,i] + indel
            mis_left = score_matrix[j,i-1] + indel
            score_matrix[j,i] = min(mis_up, mis_left, base_score+score_matrix[j-1,i-1])

    final_score = score_matrix[-1,-1] # 两条序列的最佳匹配得分

    return int(final_score)

logging.basicConfig(level = logging.INFO, 
                    format = '%(asctime)s %(levelname)s %(filename)s %(funcName)s[line:%(lineno)d] : %(message)s', 
                    datefmt = '%Y-%m-%d %H:%M:%S'
                    )

contigs = ["chr1", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
        "chr2", "chr20", "chr21", "chr22", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chrX", "chrY"]

genome_fai = "/public/home/fan_lab/wangjie/genome/hg38/hg38_22_XYM.fa.fai"
contig_length = {}
with open(genome_fai, 'r') as ifl:
    for l in ifl.readlines():
        line = l.strip().split()
        contig_length[line[0]] = line[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sample_id', help='sample id')
    parser.add_argument('-i', '--input_file', help='bed file of Straglr')
    parser.add_argument('-tsv', '--tsv_file', help='tsv file of Straglr')
    parser.add_argument('-o', '--output', help='output vcf file, uncompressed')
    parser.add_argument('-tr', '--tr_region', help='reference tr region file. \
                                It must have 7 columns which are: chrom, start, end, ref_motif, motif_length, copy_number, tr_pos. \
                                tr_pos format: chrom:start_end. \
                                This file should be sorted by chrom and position')
    parser.add_argument('--min_support', default=5, type=int, help="Minimum number of reads supporting an allele")

    if len(sys.argv)==1:   
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    clean_bed_file(args)

def cal_std_cv(read_copy_numbers):
    '''
        read_copy_numbers: a list of float numbers
    '''
    read_copy_numbers = np.array(read_copy_numbers)
    std = np.std(read_copy_numbers, ddof = 1)
    ave = np.mean(read_copy_numbers)
    cv = std / ave

    return round(std, 3), round(cv, 3)

def list_to_str(alist, sep=","):
    '''
        alist: a list of numbers
        return: a string seperated by sep    
    '''
    return sep.join([str(x) for x in alist])

def clean_tsv_file(raw_tsv_file):
    '''
        对tsv文件进行clean
        返回clean后的字典，其中key为tr_pos, value为所有reads拷贝数的标准差和变异系数，以及每条reads的拷贝数(字符串形式)
    '''
    per_read_result = {}

    counter = 0
    p_itr = 100000
    logging.info("Cleaning raw tsv file: %s"%(raw_tsv_file))
    with open(raw_tsv_file, 'r') as ifl:
        for l in ifl.readlines():
            if l.startswith("#"):
                continue
            
            line = l.strip().split("\t")
            chrom, start, end = line[0:3]
            genotype, read_copy_number = line[4], line[6]
            tr_pos = "%s:%s_%s"%(chrom, start, end)

            if ";" in genotype: # 找到了两个allele
                # 需要分别存储两个allele对应的reads的拷贝数
                locus_type = "het" # 杂合
                alleles = [float(g.split('(')[0]) for g in genotype.split(';')]
                read_allele = float(line[10])

                if read_allele not in alleles:
                    raise Exception("Error! Copy number of this read is not in genotype. Locus is %s"%(tr_pos))

                if tr_pos not in per_read_result:
                    per_read_result[tr_pos] = [locus_type, [], []]

                per_read_result[tr_pos][alleles.index(read_allele) + 1].append(float(read_copy_number))

            else: # 仅找到了一个allele
                locus_type = "homo" # 杂合

                if tr_pos not in per_read_result:
                    per_read_result[tr_pos] = [locus_type, []]

                per_read_result[tr_pos][1].append(float(read_copy_number))
            
            counter += 1
            if counter % p_itr == 0:
                logging.info("%s done"%(counter))

    copy_number_std_cvs = {} # 存储每个TR所有reads拷贝数的标准差和变异系数
    logging.info("Calculate std and cv of read copy numbers")
    for tr_pos, results in per_read_result.items():
        locus_type = results[0]

        if locus_type == "het":
            allele1_cns, allele2_cns = results[1:3] # cns represents copy_numbers
            std_cvs = [cal_std_cv(x) for x in [allele1_cns, allele2_cns]]
            std = ','.join([str(x[0]) for x in std_cvs])
            cv = ','.join([str(x[1]) for x in std_cvs])
            allele_cns_str = ';'.join([list_to_str(x) for x in [allele1_cns, allele2_cns]])

        elif locus_type == "homo":
            allele_cns = results[1]
            std, cv = [str(x) for x in cal_std_cv(allele_cns)]
            std = "%s,%s"%(std, std)
            cv = "%s,%s"%(cv, cv)
            allele_cns_str = list_to_str(allele_cns)

        else:
            raise Exception("Wrong locus type %s. Must be one of het and homo"%(locus_type))

        copy_number_std_cvs[tr_pos] = [std, cv, allele_cns_str]
    logging.info("Done. Total %s TRs are calculated"%(len(copy_number_std_cvs)))

    return copy_number_std_cvs

def clean_bed_file(args):
    input_file = args.input_file
    output_file = args.output
    tr_region = args.tr_region
    sample_id = args.sample_id
    raw_tsv_file = args.tsv_file
    min_support_allepe = args.min_support
    
    # read trs
    tr_region_dict = OrderedDict()
    with open(tr_region, 'r') as ifl:
        for l in ifl.readlines():
            line = l.strip().split('\t')
            if len(line) != 7:
                raise Exception("Error! Wrong format of tr_region file %s. It must have 7 columns which are: chrom, start, end, ref_motif, motif_length, copy_number, tr_pos "%(tr_region))

            chrom, start, end, ref_motif, motif_length, copy_number, tr_pos = line
            # todo 可以检查一下 tr_pos 是不是chrom:start_end_motif格式
            tr_region_dict[tr_pos] = [chrom, start, end, ref_motif, motif_length, copy_number]
    logging.info("Total %s loci of %s are loaded"%(len(tr_region_dict), tr_region))

    copy_number_std_cvs = clean_tsv_file(raw_tsv_file)

    # 读入bed文件里vntr的结果
    vntr_results = {}
    with open(input_file, 'r') as ifl:
        ifl.readline()
        for l in ifl.readlines():
            line = l.strip().split('\t')
            chrom, start, end, qry_motif = line[0:4]
            tr_pos = "%s:%s_%s"%(chrom, start, end)
            if len(line) == 7:
                logging.warning("This TR(%s) is not been genotyped. Will skip it"%(tr_pos))
                continue

            # sanity check
            if tr_pos not in tr_region_dict:
                raise Exception("Error! Your tr_region file: %s does not contain this TR: %s. Maybe you supplied a wrong tr_region file?"%(tr_region, tr_pos))

            copy_number_A1, num_support_reads_A1 = line[5], line[6]
            copy_number_A2, num_support_reads_A2 = line[8], line[9]
            vntr_results[tr_pos] = [qry_motif, copy_number_A1, copy_number_A2, num_support_reads_A1, num_support_reads_A2]

    logging.info("Total %s vntrs of %s are loaded"%(len(vntr_results), input_file))

    logging.info("Begin to write clean results to %s"%(output_file))
    # header = 'chrom start end ref_motif ref_motif_len ref_motif_copy_number tr_pos qry_motif dis_from_ref_motif copy_number_A1 copy_number_A2 num_support_reads_A1 num_support_reads_A2\n'.replace(' ', '\t')
    with open(output_file, 'w') as ofl:
        # write header
        ofl.write("##fileformat=VCFv4.2\n")

        # write contigs
        for contig in contigs:
            contig_len = contig_length[contig]
            ofl.write('##contig=<ID=%s,length=%s>\n' % (contig, contig_len))
        
        # write info and format
        ofl.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate of the VNTR">\n')
        ofl.write('##INFO=<ID=RefMotif,Number=1,Type=String,Description="Reference motif of the VNTR">\n')
        ofl.write('##INFO=<ID=MotifLen,Number=1,Type=Float,Description="Motif length of the VNTR">\n')
        ofl.write('##INFO=<ID=CopyNumber,Number=1,Type=Float,Description="Copy number of the VNTR">\n')
        ofl.write('##FILTER=<ID=PASS,Description="All filters passed">\n')

        ofl.write('##FORMAT=<ID=Motif,Number=1,Type=String,Description="Detected motif(query motif)">\n')
        ofl.write('##FORMAT=<ID=CNB,Number=2,Type=Float,Description="Copy numbers">\n')
        ofl.write('##FORMAT=<ID=ED,Number=1,Type=Float,Description="Edit distance between detected motif and reference motif">\n')
        ofl.write('##FORMAT=<ID=NSR,Number=2,Type=Integer,Description="Number support reads of two alleles">\n')
        ofl.write('##FORMAT=<ID=RSTD,Number=2,Type=Float,Description="Standard deviation(STD) of copy numbers of each reads(R) supporting this TR">\n')
        ofl.write('##FORMAT=<ID=RCV,Number=2,Type=Float,Description="Coefficient of variation(CV) of copy numbers of each reads(R) supporting this TR">\n')
        ofl.write('##FORMAT=<ID=RCNS,Number=1,Type=String,Description="Copy numbers(CNS) of each reads(R) supporting this TR">\n')

        # write header of data lines
        ofl.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n"%(sample_id))
        # ofl.write(header)
        skipped_vntr = []
        for tr_pos in vntr_results.keys():
            chrom, start, end, ref_motif, motif_length, copy_number = tr_region_dict[tr_pos]
            info = "END=%s;RefMotif=%s;MotifLen=%s;CopyNumber=%s" % (end, ref_motif, motif_length, copy_number)
            tr_format = "Motif:ED:CNB:NSR:RSTD:RCV:RCNS" # 

            qry_motif, copy_number_A1, copy_number_A2, num_support_reads_A1, num_support_reads_A2 = vntr_results[tr_pos]
            
            if copy_number_A1 == '-' and copy_number_A2 == '-':
                # skipped_vntr.append([tr_pos, 'empty_copy_number'])
                logging.warning('Copy numbers of tr %s are empty! It will be skipped.'%(tr_pos))
                continue
            else:
                # 处理检测到的vntr
                dis_from_ref_motif = cal_edit_distance(qry_motif, ref_motif)
                
                if copy_number_A2 != '-': # 找到了两个allele
                    # 判断num_support是否小于某个阈值
                    if int(num_support_reads_A1) < min_support_allepe and int(num_support_reads_A2) < min_support_allepe: # 两个allele均不可信
                        logging.info("Number of support reads of both alleles(%s) are below %s. It will be discarded"%(tr_pos, min_support_allepe))
                        # skipped_vntr.append([tr_pos, 'Number_of_support_reads_both_alleles_below_%s'%(min_support_allepe)])
                        continue
                    elif int(num_support_reads_A1) < min_support_allepe: # allele1 不可信，认为这里是纯合的
                        copy_number_A1 = copy_number_A2
                        num_support_reads_A1 = num_support_reads_A2
                        logging.info("Number of support reads of allele1(%s) is below %s. It will be discarded"%(tr_pos, min_support_allepe))
                        # skipped_vntr.append([tr_pos, 'Number_of_support_reads_allele1_below_%s'%(min_support_allepe)])
                    elif int(num_support_reads_A2) < min_support_allepe: # allele2 不可信，认为这里是纯合的
                        copy_number_A2 = copy_number_A1
                        num_support_reads_A2 = num_support_reads_A1
                        # skipped_vntr.append([tr_pos, 'Number_of_support_reads_allele2_below_%s'%(min_support_allepe)])
                        logging.info("Number of support reads of allele2(%s) is below %s. It will be discarded"%(tr_pos, min_support_allepe))
                    elif float(copy_number_A1) < float(copy_number_A2): # 两个allele均可信，需要互换位置，将大的allele放到前面
                        copy_number_A1 = copy_number_A2
                        copy_number_A2 = copy_number_A1
                        num_support_reads_A1, num_support_reads_A2 = num_support_reads_A2, num_support_reads_A1
                    else: # 无需改动allele的拷贝数与num_support
                        pass
                else:
                    if int(num_support_reads_A1) < min_support_allepe:
                        # skipped_vntr.append([tr_pos, 'Number_of_support_reads_allele_below_%s'%(min_support_allepe)])
                        logging.info("Number of support reads of this allele(%s) is below %s. It will be discarded"%(tr_pos, min_support_allepe))
                        continue

                    copy_number_A2 = copy_number_A1
                    num_support_reads_A2 = num_support_reads_A1

                rstd, rcv, allele_cns_str = copy_number_std_cvs[tr_pos]
                tr_format_value = "%s:%s:%s,%s:%s,%s:%s:%s:%s" % (qry_motif, dis_from_ref_motif, 
                                                                copy_number_A1, copy_number_A2, 
                                                                num_support_reads_A1, num_support_reads_A2, 
                                                                rstd, rcv, allele_cns_str)

            # 将最终结果写出
            tr_pos = "%s:%s_%s"%(chrom, int(start)+1, end)
            ofl.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(
                    chrom, int(start)+1, tr_pos, "N", ".", ".", "PASS", info, tr_format, tr_format_value)
                    )

    # logging.info('Total %s VNTRs are skipped'%(len(skipped_vntr)))
    # logging.info('Skipped VNTRs and reason are listed below')
    # for vntr in skipped_vntr:
    #     print('%s\t%s'%(vntr[0], vntr[1]))
    logging.info('Done. Output file is %s'%(output_file))


if __name__ == '__main__':
    main()
