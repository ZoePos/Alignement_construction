#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import collections
import copy
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def main():

    # main algo
    # 1 collect all data from file to dictionnaries
    # 2 create new seq for each sample by applying variants to ref with some criteria
    # 3 export new seq and trace per sample and per gene 

    # vocabulary :
    # gene = sequence name in fasta file (= seqid)
    # site = pair (gene, pos) : this is used for variants

    ################
    # 1 collect all data from file to dictionnaries
    ################
    
    # collect pops
    pops = ('E', 'W1', 'W2', 'W3')
    samples, sample2pop, pop2samples = collect_pop(pops)
    # data structures
    # sample => list of all samples
    # sample2pop => dictionnary to get get pop corresponding to a sample : sample2pop[sample] = pop
    # pop2samples => dictionnary to get sample list for one pop : pop2samples[pop] = set of samples
    # restrict samples for test purpose
    # samples = samples[:1]
    print("{} samples".format(len(samples)))
        
    # read ref
    ref = read_ref(args.ref)
    # data structure
    # dictionary of sequence string : ref[gene] = seq_string
    
    # reads vcf
    # the vcf file should be a clean vcf (only filtered snps) since these snps will be used to create sample sequce by applying snps to ref sequence
    gene2pos, sites = read_vcf(args.clean_vcf_file, samples)
    # data structure 
    # sites[site] = {'ref_all' : nuc, 'alt_all' : nuc_list } ==> data for each site
    # gene2pos = {'gene1' : (pos1, pos2 ...), 'gene2' : (pos1, pos2 ...), } ==> list of pos for each gene (or gene)

    # list of genes found in vcf file
    genes = sorted(gene2pos.keys())
    # restrict to some genes
    # genes = genes[:2]

    # collect genotypes and allele depth from GT and AD FORMAT fields
    # FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    # FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
    # data structure 
    #  vcf_format[key][(sample, site)] = value 
    vcf_format = dict()
    vcf_format ['GT'] = read_vcf_FORMAT('GT', "./.")
    vcf_format ['AD'] = read_vcf_FORMAT('AD', ".")

    ################
    # 2 create new seq for each sample by applying variants to ref with some criteria
    ################

    # for each sample, create new seqs and traces from ref + vcf
    # choose to collect new sequences and traces in data structures instead of writing files
    # so that we can later create various sequence and trace files, per sample and per gene 
    # new_seq[sample][gene] = sequence
    # traces[sample][gene][pos] = message text
    new_seqs = dict()
    traces = dict()
    for sample in samples:
        create_new_seq(ref, sample, sites, vcf_format, new_seqs, traces)

    ################
    # 3 export new seq and trace per sample and per gene 
    ################

    # write new seq and trace per sample
    for sample in samples:
        pop = sample2pop[sample]
        per_sample_seq_file = "{}/{}_{}.fasta".format(args.seq_per_sample_dir, pop, sample)
        per_sample_trace_file = "{}/{}_{}.trace".format(args.seq_per_sample_dir, pop, sample)
        print("write per sample seq file {} and trace file {}".format(per_sample_seq_file, per_sample_trace_file))
        with open(per_sample_seq_file, 'w') as f_out_seq, open(per_sample_trace_file, 'w') as f_out_trace:
            for gene in genes:
                my_seq = Seq("".join(new_seqs[sample][gene]))
                my_id = "{}_{}_{}".format(pop, sample, gene) # new seqence ID = POP_SAMPLE_GENE
                record = SeqRecord(my_seq, id=my_id, description='')
                if args.verbose:
                    print("write seq {} for sample {} gene {}".format(my_id, sample, gene))
                SeqIO.write(record, f_out_seq, "fasta")
                for gene in genes:
                    if gene in traces[sample]:
                        for pos in gene2pos[gene]:
                            if pos in traces[sample][gene]:
                                f_out_trace.write("{}_{}\t{}\t{}\t{}\n".format(pop, sample, gene, pos, traces[sample][gene][pos]))

    # write new seq and trace per gene
    for gene in genes:
        per_gene_seq_file = "{}/{}.fasta".format(args.seq_per_gene_dir, gene)
        print("write per gene seq file {}".format(per_gene_seq_file))
        with open(per_gene_seq_file, 'w') as f_out_seq:
            for sample in samples:
                pop = sample2pop[sample]
                my_seq = Seq("".join(new_seqs[sample][gene]))
                my_id = "{}_{}_{}".format(pop, sample, gene)  # new seqence ID = POP_SAMPLE_GENE
                record = SeqRecord(my_seq, id=my_id, description='')
                if args.verbose:
                    print("write seq {} for sample {} gene {}".format(my_id, sample, gene))
                SeqIO.write(record, f_out_seq, "fasta")

        per_gene_trace_file = "{}/{}.trace".format(args.seq_per_gene_dir, gene)
        print("write per gene trace file {}".format(per_gene_trace_file))
        with open(per_gene_trace_file, 'w') as f_out_trace:
            for pos in gene2pos[gene]:
                for sample in samples:
                    pop = sample2pop[sample]
                    if gene in traces[sample]:
                        if pos in traces[sample][gene]:
                            f_out_trace.write("{}_{}\t{}\t{}\t{}\n".format(pop, sample, gene, pos, traces[sample][gene][pos]))


def create_new_seq(ref, sample, sites, vcf_format, new_seqs, traces):
    """
    create new seqs and traces from ref + vcf
    data structure in :
    - ref : dictionary of sequence string : ref[gene] = seq_string
    - sample : one sample
    - sites : data for each site : sites[site] = {'ref_all' : nuc, 'alt_all' : nuc_list }
    - vcf_format : genotype ('GT' key) and allele depth ('AD' key) for each site and each sample : vcf_format[key][(sample, site)] = value
    data structure out : 
    - new_seq[sample][gene] = new sequence string
    - traces[sample][gene][pos] = message text

    choose to collect new sequences and traces in data structures instead of writing files
    so that we can later create various sequence and trace files, per sample and per gene 
    """
    
    new_seqs[sample] = copy.deepcopy(ref) # copy sequence from ref
    traces[sample] = dict()
    
    for site in sites:
        gene, pos = site
        ref_all = sites[site]['ref_allele']
        alt_all = sites[site]['alt_alleles']
        alleles = [ ref_all ] + alt_all
        #print("{} {} : ref_all={} alt_all={} all_all={}".format(sample, site, ref_all, alt_all, alleles))

        # set !!! to trace variable so taht we can check at the end that trace is OK
        ref_allele = "!!! REF_ALLELE"
        alt_allele = "!!! ALT_ALLELE"
        GT = "!!! GENOTYPE"
        A1 = "!!! ALLELE1"
        A2 = "!!! ALLELE2"
        AD1 = "!!! AD1"
        AD2 = "!!! AD2"
        AD_tot = "!!! AD_tot"
        site_type = "!!! SITE_TYPE"
        allele_type = "!!! ALLELE_TYPE"
        nuc = "!!! NUCLEOTIDE"
        msg = "!!! MESSAGE"
        
        if (sample, site) not in vcf_format['GT']:  # no genotype for this site and this sample
            nuc = "N"
            site_type = "no_reads"
        else:
            GT = vcf_format['GT'][(sample, site)]
            # alleles
            A1, A2 = [int(i) for i in GT.split('/')]
            # alleles depth
            AD = [int(i) for i in vcf_format['AD'][(sample, site)].split(',')]
            AD1 = AD[A1]
            AD2 = AD[A2]
            AD_tot = AD1 + AD2
            # alleles nucleotides
            N1 = alleles[A1]
            N2 = alleles[A2]
            if A1 == A2:
                site_type = "Hom"
            else:
                site_type = "Het"

            if (AD_tot) < args.AD_thr:
                # alleles depth too small, => N
                nuc = 'N'
                msg = 'AD_tot {} < {}'.format(AD_tot, args.AD_thr)
            elif A1 == A2:
                # Hom, use any allele
                nuc = alleles[A1]
                msg = ""
            else:
                # Het, use major allele
                if AD1 == AD2:
                    # same count => mask with N
                    nuc = 'N'
                    msg = 'AD1 == AD2'
                elif AD1 > AD2:
                    nuc = alleles[A1]
                    msg = 'AD1 >= AD2'
                else:
                    nuc = alleles[A2]
                    msg = 'AD1 < AD2'

        if nuc == '*': # stand for indel
            nuc = 'N'
            msg = "Alert : nuc=* => N"

        # modif sequence
        # pos is 1-based, seq has been converted to string list which is 0-based
        new_seqs[sample][gene][pos - 1] = nuc

        # format trace line
        if nuc == "N":
            allele_type = "UNKNOWN"
        elif nuc == ref_all:
            allele_type = "REF"
        else :
            allele_type = "ALT"

        if site_type == 'no_reads':
            trace = "ref={} alt={} ==> nuc={} [no reads]".format(ref_all, alt_all, nuc)
        else:
            if site_type == 'Hom':
                count = "{}x{}:{}".format(AD1, A1, N1)
            else:
                count = "{}x{}:{}, {}x{}:{}".format(AD1, A1, N1, AD2, A2, N2)
            trace = "ref={} alt={} GT={} {} ==> nuc={} [{} {}] {}".format(ref_all, alt_all, GT, count, nuc, site_type, allele_type, msg)
        if "!!!" in trace:  # check trace is OK
            print("bad trace {}".format(trace))
            exit(0)
                
        if sample not in traces:
            traces[sample] = dict()
        if gene not in traces[sample]:
            traces[sample][gene] = dict()
        traces[sample][gene][pos] = trace

def read_vcf(in_file, samples):
    # read vcf file to collect
    # - ref and alt alleles
    # - variants as a list of sites (ie pair gene, pos) 

    # the vcf file should be a clean vcf (only filtered snps) since these snps will be used to create sample sequce by applying snps to ref sequence
    # build data structures 
    # sites[site] = {'ref_all' : nuc, 'alt_all' : nuc_list } ==> data for each site
    # gene2pos = {'gene1' : (pos1, pos2 ...), 'gene2' : (pos1, pos2 ...), } ==> list of pos for each gene (or gene)

    # also check that no sample is missing in vcf file
    
    sites = dict()
    gene2pos = dict()
    with open(in_file) as f_in:
        for l in f_in:
            tab = l.strip().split()
            if l.startswith("##"): # skip header lines
                continue
            if l.startswith("#CHROM"): # check sample list
                vcf_samples = tab[9:]
                missing_samples = set(samples) - set(vcf_samples)
                if len(missing_samples) > 0:
                    print("{} samples from pops files, {} vcf_samples, {} missing_sample".format(len(samples), len(vcf_samples), len(missing_samples)))
                    exit(0)
                continue
            gene= tab[0]
            pos = int(tab[1])
            site = (gene, pos)
            
            ref_allele = tab[3] # allele in reference
            alt_alleles = tab[4].split(',') # list of alternative alleles (* means indels)
            sites[site] = {'ref_allele' : ref_allele, 'alt_alleles' : alt_alleles} 

            # add pos for this gene
            if gene not in gene2pos:
                gene2pos[gene] = list()
            gene2pos[gene].append(pos)

    print("{} genes, {} sites from clean_vcf file {}".format(len(gene2pos), len(sites), in_file))
    return gene2pos, sites

def read_vcf_FORMAT(key, missing_value):
    """ 
    extract one FORMAT field from vcf to tsv file using vcftools
    then parse this tsv file
    """
    cmd = "vcftools --vcf {} --extract-FORMAT-info {} --out {} 2> /dev/null".format(args.clean_vcf_file, key, args.out_basename)
    os.system(cmd)
    res_file = "{}.{}.FORMAT".format(args.out_basename, key)
    return read_vcf_FORMAT_file(res_file, missing_value)
    
def read_vcf_FORMAT_file(in_file, missing_value):
    """
    reads tsv file with one FORMAT field extracted using vcftools
    build data structre
    - data[(sample, site)] = value
    """
    
    data = dict()
    with open(in_file) as f_in:
        for l in f_in:
            tab = l.strip().split()
            if l.startswith("CHROM"): # collect sample names for each col
                samples = tab[2:]
                continue
            gene = tab[0]
            pos = int(tab[1])
            site = (gene, pos)
            values = tab[2:]
            for sample, value in zip(samples, values):
                #print("site {} sample {} value {}".format(site, sample, value))
                key = (sample, site)
                if value == missing_value:
                    continue
                data[key] = value
    print("{:,} non missing values in {}".format(len(data), in_file))
    return(data)

def read_ref(in_file):
    """ 
    read reference fasta file and build a dictionnary ref[gene] = sequence_string 
    use sequence string instead of sequence object because to facilitate sequences modification
    """
    
    ref = dict()
    n = 0 # number of sequence
    n_bp = 0 # numbep of base pair
    with open(in_file) as f_in:
        for record in SeqIO.parse(f_in, "fasta") :
            ref[record.id] = list(str(record.seq))
            n += 1
            n_bp += len(record)
    print("{:,} seq {:,} bp from {}".format(n, n_bp, in_file))
    return ref

def collect_pop(pops):
    """
    read each pop file pop_E, pop_W1, pop_W2 ... containing one sample per line
    build 3 data structures
    - sample => list of all samples
    - sample2pop => dictionnary to get get pop corresponding to a sample : sample2pop[sample] = pop
    - pop2samples => dictionnary to get sample list for one pop : pop2samples[pop] = set of samples
    """
    pop2samples = dict()
    sample2pop = dict()
    samples = list()
    for pop in pops:
        pop2samples[pop] = set()
        with open("pop_{}".format(pop)) as f_in:
            for l in f_in:
                sample = l.strip()
                samples.append(sample)
                pop2samples[pop].add(sample)
                sample2pop[sample] = pop
        print("pop {} : {} samples".format(pop, len(pop2samples[pop])))
    return samples, sample2pop, pop2samples
                            
if __name__ == '__main__':
    description = """ 
    create seq per ind or per gene from a reference in fasta format and variants in vcf format
    per ind and site :
        set N if depth < threshold
        set N if both alleles have same depth
        modif ref seq with major allele

    uses v4 file tree for collectig ref.fasta, clean_vcf and for writing results
    output : fasta file + trace file in both folders v4/res/new_seq_per_sample and v4/res/new_seq_per_gene
    """

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-v", "--verbose", action="store_const", const="-v", default="")
    parser.add_argument("--name", help="name to be used for collecting data from files and creating results", required=True)
    parser.add_argument("--AD_thr", type=int, help="threshold for Allele depth, default=5", default=5)
    
    args = parser.parse_args()
    args.ref = "v4/ref/{}.fasta".format(args.name)
    args.clean_vcf_file = "v4/vcf/{}/clean5q_all_snp.recode.vcf".format(args.name)
    args.out_basename = args.clean_vcf_file.replace('.vcf', '')

    args.seq_per_sample_dir = "v4/new_seq_per_sample/{}".format(args.name)
    args.seq_per_gene_dir = "v4/new_seq_per_gene/{}".format(args.name)
    
    if not os.path.isdir(args.seq_per_sample_dir):
        print("create dir {}".format(args.seq_per_sample_dir))
        os.makedirs(args.seq_per_sample_dir)
    if not os.path.isdir(args.seq_per_gene_dir):
        print("create dir {}".format(args.seq_per_gene_dir))
        os.makedirs(args.seq_per_gene_dir)

    print(args)
    main() 
