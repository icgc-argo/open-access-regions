#!/usr/bin/env python3
"""
Authors:
  Linda Xiang <linda.xiang@oicr.on.ca>

Copyright (c) 2021, Ontario Institute for Cancer Research (OICR).
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.
You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.

Acknowledgements:
Some functions and code are adapted and modified from GTFtools: 
http://www.genemine.org/gtftools.php 

"""

from argparse import ArgumentParser
import gzip
from copy import deepcopy
import os
import shutil
import urllib.request as request
from contextlib import closing
from urllib.error import URLError
from datetime import date

def get_promoter(tcx, glen, w=200):
    promoter = []
    for each in tcx:
        if each[3] == "+":
            iregion = (each[0],max(int(each[1])-w, 0),int(each[1])-1,each[3])
        elif each[3] == '-':
            iregion = (each[0],int(each[2])+1,min(int(each[2])+w, int(glen['chr'+each[0]])),each[3])
        else:
            pass
        promoter.append(iregion)
    return promoter

def get_splice_site(intron):
    splice_site = []
    for each in intron:
        strand=each[3]
        if strand=='+':
            donor_left=each[1]-3
            donor_right=each[1]+6

            acceptor_left=each[2]-20
            acceptor_right=each[2]+3
        elif strand=='-':
            donor_left=each[2]-6
            donor_right=each[2]+3

            acceptor_left=each[1]-3
            acceptor_right=each[1]+20
      
        joined_donor=(each[0],donor_left,donor_right,each[3])
        joined_acceptor=(each[0],acceptor_left,acceptor_right,each[3])
        
        splice_site.append(joined_donor)
        splice_site.append(joined_acceptor)
    return splice_site

def get_UTR(utr, cds):
    # determine 5UTR and 3UTR
    UTR5, UTR3 = [], []
    listCDS = deepcopy(cds)
    listCDS.sort(key = lambda x: (x[0],x[1]))
    if not listCDS or not utr: 
        return UTR5, UTR3
    first_cds = listCDS[0]
    last_cds = listCDS[-1]
    for origUTR in utr:
        thisUTR = deepcopy(origUTR)
        if origUTR[1] < first_cds[1]:
            if origUTR[2] >= first_cds[1]:
                thisUTR = (origUTR[0], origUTR[1], first_cds[1]-1, origUTR[3])
            if thisUTR[3] == '+':
                UTR5.append(thisUTR)
            else:
                UTR3.append(thisUTR)
        elif origUTR[2] > last_cds[2]:
            if origUTR[1] <= last_cds[2]:
                thisUTR = (origUTR[0], last_cds[2]+1, origUTR[2], origUTR[3])
            if thisUTR[3] == '+':
                UTR3.append(thisUTR)
            else:
                UTR5.append(thisUTR)
        else:
            pass
    return UTR5, UTR3

def exon2intron(featureRange):
    intron=[]
    nRange=len(featureRange)
    featureRange.sort(key = lambda x: (x[0],x[1]))	
    if nRange > 1:
        for i in range(1,nRange):
            exon0=featureRange[i-1]
            exon1=featureRange[i]
            thisintron=(exon0[0],int(exon0[2])+1,int(exon1[1])-1,exon1[3])
            intron.append(thisintron)
    return intron

def get_conjugate(featureRange, glen):
    conjugate=[]
    nRange=len(featureRange)
    featureRange.sort(key = lambda x: (x[0],x[1]))	
    if nRange > 1:
        first = (featureRange[0][0],1,int(featureRange[0][1])-1,featureRange[0][3])
        conjugate.append(first)

        for i in range(1,nRange):
            feature0=featureRange[i-1]
            feature1=featureRange[i]
            if feature0[0] == feature1[0]:
                inbetween=(feature0[0],int(feature0[2])+1,int(feature1[1])-1,feature1[3])
                conjugate.append(inbetween)
            else:
                chr_last = (feature0[0],int(feature0[2])+1,int(glen['chr'+feature0[0]]),feature0[3])
                conjugate.append(chr_last)
                nextchr_first = (feature1[0],1,int(feature1[1])-1,feature1[3])
                conjugate.append(nextchr_first)

        last = (featureRange[-1][0],int(featureRange[-1][2])+1,int(glen['chr'+featureRange[-1][0]]),featureRange[-1][3])
        conjugate.append(last)

    return conjugate

def neighbor_merge(range1,range2,gene_id=None):
    if range1[0] == range2[0] and range2[1]<=range1[2]:
        if gene_id:
            merged = [(range1[0],range1[1],max(range1[2],range2[2]),range1[3],gene_id)]
        else:
            merged = [(range1[0],range1[1],max(range1[2],range2[2]),range1[3])]
    else:
        if gene_id:
            merged = [(range1[0],range1[1],range1[2],range1[3],gene_id), (range2[0],range2[1],range2[2],range2[3],gene_id)]
        else: 
            merged=[(range1[0], range1[1], range1[2], range1[3]),(range2[0], range2[1], range2[2], range2[3])]
    return merged

def bedmerge(featureRange, gene_id=None):
	# featureRange: a list of ranges in bed format, such as [(1,1000,2000,+),(1,2200,3000,-)]
    featureRange.sort(key = lambda x: (x[0],x[1]))	
    merged=[]
    nRange=len(featureRange)
    if nRange == 0:
        return merged
    elif nRange == 1:
        if gene_id:
            merged=[featureRange[0]+(gene_id,)]
        else:
            merged=featureRange
    elif nRange == 2:
        imerge=neighbor_merge(featureRange[0],featureRange[1],gene_id)
        for each in imerge:
            merged.append(each)
    else:
        i = 2
        imerge=neighbor_merge(featureRange[0],featureRange[1],gene_id)
        n_imerge=len(imerge)
        while n_imerge > 0:
            if n_imerge == 2:
                merged.append(imerge[0])
			
            imerge=neighbor_merge(imerge[n_imerge-1],featureRange[i],gene_id)
            n_imerge=len(imerge)
            if i == nRange-1:
                for each in imerge:
                    merged.append(each)
                n_imerge = -1
            i+=1
    return merged


# convert GENCODE to ENSEMBL format
def gencode2ensembl(gtf1,gtf2):
    #gtf1: ensembl GTF
    #gtf2: gencode GTF
    fn1_open = gzip.open if gtf1.endswith(".gz") else open
    with fn1_open(gtf1) as fn1:
        with open(gtf2, 'w') as fn2:
            for line in fn1:
                line_str = line.decode('utf8')
                if line_str[0:3] == 'chr':
                    t = line_str.strip().split('\t')
                    # update
                    cname = t[0].split('chr')[1]
                    if cname == 'M': cname = 'MT'
                    # replace
                    t[0] = cname
                    fn2.write('\t'.join(t)+'\n')
                else:
                    fn2.write(line_str)
	

# check format of GTF: ensembl or gencode
def gtf_format_check(gtf):
    fn_open = gzip.open if gtf.endswith(".gz") else open
    with fn_open(gtf) as f:
        for line in f:
            line_str = line.decode('utf8')
            if line_str.startswith("#"): continue
            chrom=line_str.split('\t')[0]
            if len(chrom) >= 4:
                return 'GENCODE'
            else:
                return 'ENSEMBL'

def get_gene_dict(GTFfile,feature_type=['CDS','exon','UTR','transcript'],chroms=[str(r) for r in range(1,23)]+['X','Y']):	
    gene_dict={}   
    with open(GTFfile, 'r') as f:
        for line in f:
            if line[0] == '#': continue
            table = line.split('\t')
            if not table[0] in chroms: continue
            if table[2] in feature_type:
                ensid  = line.split('gene_id')[1].split('"')[1]
                transcript_type = line.split('transcript_type')[1].split('"')[1]			
                tcx = line.split('transcript_id')[1].split('"')[1]			
                record = (table[0],int(table[3]),int(table[4]),table[6])
                if not ensid in gene_dict: gene_dict[ensid] = {}
                if not tcx in gene_dict[ensid]: gene_dict[ensid][tcx] = {}
                if not table[2] in gene_dict[ensid][tcx]: gene_dict[ensid][tcx][table[2]] = []
                gene_dict[ensid][tcx][table[2]].append(record)
                if not 'transcript_type' in gene_dict[ensid][tcx]: gene_dict[ensid][tcx]['transcript_type'] = transcript_type
                            
    return gene_dict

def create_region_bed(bed_region_list, bed_dir, bed_fname, bedtype='0') :
    if not os.path.exists(bed_dir):
        os.makedirs(bed_dir)
    bed_filename = os.path.join(bed_dir, bed_fname+'.bed.gz')
	
    bed_region_list.sort(key = lambda x: (0,int(x[0]),x[1]) if x[0].isdigit() else (1,x[0],x[1]))
    with gzip.open(bed_filename,'wb') as f:
        for item in bed_region_list:
            item_new = [e for e in item] 
            item_new[1] = int(item[1])-1 if bedtype == '0' else item[1] 
            item_string = 'chr'+'\t'.join(str(e) for e in item_new)+'\n'
            f.write(item_string.encode('utf-8'))
    return(bed_region_list)

def download_file(url, local_filename):
    try:
        with closing(request.urlopen(url)) as r:
            with open(local_filename, 'wb') as f:
                shutil.copyfileobj(r, f)
    except URLError as e:
        if e.reason.find('No such file or directory') >= 0:
            raise Exception('FileNotFound')
        else:
            raise Exception(f'Something else happened. "{e.reason}"')

    return local_filename

element_type = {
	'protein_coding': ['protein_coding'],
    'lncRNA': ['lncRNA'],
	'smallRNA': ['Mt_rRNA', 'Mt_tRNA', 'miRNA', 'misc_RNA', 'rRNA', 'scRNA', 'snRNA', 'snoRNA', 'ribozyme', 'sRNA', 'scaRNA']
}

def main():
    parser = ArgumentParser()
    parser.add_argument('-g','--gtf_dir', dest="gtf_dir", type=str, default="../data/hg38/annotation", help="local directory for GTF file: only ENSEMBL or GENCODE GTF file accepted")
    parser.add_argument('-u','--gtf_url', dest="gtf_url", type=str, default="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human", help="specify URL to download GTF file. Use GENCODE by default")
    parser.add_argument('-v','--gtf_version', dest="gtf_version", type=str, default="37", help="specify GTF release version to download")
    parser.add_argument('-l','--glen', dest='glen',type=str, default="../data/hg38/hg38.chrom.sizes",help="file for genome size")
    parser.add_argument('-w','--window',dest='window_size',type=int,default=200,help="the value of w in calculating TSS regions. Default: 200")
    parser.add_argument('-b','--bed_dir',dest='bed_dir',default="../data/hg38/bed",help="directory for output bed files")
    parser.add_argument('-r', '--regions',dest="regions",type=str, nargs='+',default=['utr5', 'utr3', 'cds', 'exon', 'intron', 'protein_coding','protein_coding_promoter', \
        'protein_coding_splice_site', 'lncRNA', 'lncRNA_promoter', 'lncRNA_splice_site', 'smallRNA', 'smallRNA_promoter', 'smallRNA_splice_site'],\
            help="specify regions to generate bed files" )
    parser.add_argument('-o', '--open_access_regions',dest="open_access_regions",type=str, nargs='+',default=['utr5', 'utr3', 'cds', 'exon', 'protein_coding_promoter', \
        'protein_coding_splice_site', 'lncRNA_promoter', 'lncRNA_splice_site', 'smallRNA', 'smallRNA_promoter'],\
            help="specify regions to be included into open access tier" )

    args = parser.parse_args()
    regions = args.regions
    open_access_regions = args.open_access_regions

    # download the GENCODE annotation file if it does not exist locally
    url = args.gtf_url+ '/release_'+args.gtf_version+'/gencode.v'+args.gtf_version+'.annotation.gtf.gz'
    if not os.path.exists(args.gtf_dir):
        os.makedirs(args.gtf_dir)
    GTFfile = os.path.join(args.gtf_dir, url.split('/')[-1])
    if not os.path.exists(GTFfile):
        download_file(url, GTFfile)

    ftype = gtf_format_check(GTFfile)
    print('GTF format is '+ftype)
    if ftype == 'GENCODE':
        print('Converting GTF to ENSEMBL format')
        GTFfile_ensembl=GTFfile+'.ensembl'
        if not os.path.exists(GTFfile_ensembl): gencode2ensembl(GTFfile,GTFfile_ensembl)
        GTFfile = GTFfile_ensembl
    print('------> Analyzing '+GTFfile)

    # create output bed dir if not existing
    bed_dir = os.path.join(args.bed_dir, ftype.lower()+'.v'+args.gtf_version)
    if not os.path.exists(bed_dir):
        os.makedirs(bed_dir)

    gene_dict = get_gene_dict(GTFfile)

    if args.glen:
        genomeLen = {}
        with open(args.glen, 'r') as f:
            for line in f:
                table = line.strip().split('\t')
                genomeLen[table[0]] = table[1]
    
    gencode_regions = {}
    for region in regions:
        gencode_regions[region] = []

    for gene_id in gene_dict:
        gene_region = {}
        for region in regions:
            gene_region[region] = []
            gene_region['merged_'+region] = []
        
        for transcript_id in gene_dict[gene_id]:
            tcx_type = gene_dict[gene_id][transcript_id]['transcript_type']
            icds = gene_dict[gene_id][transcript_id].get('CDS', [])
            iexon = gene_dict[gene_id][transcript_id].get('exon', [])
            iutr = gene_dict[gene_id][transcript_id].get('UTR', [])
            itcx = gene_dict[gene_id][transcript_id].get('transcript', [])
            
            iintron = exon2intron(iexon)
            iutr5, iutr3 = get_UTR(iutr, icds)
            ipromoter = get_promoter(itcx, genomeLen, args.window_size)
            isplice_site = get_splice_site(iintron)

            gene_region['cds'] += icds
            gene_region['exon'] += iexon
            gene_region['intron'] += iintron
            gene_region['utr5'] += iutr5
            gene_region['utr3'] += iutr3
            for iregion in ['protein_coding', 'lncRNA', 'smallRNA']:
                if tcx_type in element_type[iregion]:
                    gene_region[iregion] += itcx
                    gene_region[iregion+'_promoter'] += ipromoter
                    gene_region[iregion+'_splice_site'] += isplice_site

        for region in regions:
            if not gene_region[region]: continue
            gene_region['merged_'+region] = bedmerge(gene_region[region], gene_id)
            gencode_regions[region] += gene_region['merged_'+region]

    date_str = date.today().strftime("%Y%m%d")
    open_access_region_list = []
    for region in gencode_regions:
        create_region_bed(gencode_regions[region], bed_dir, region+'.'+date_str)
        if not region in open_access_regions: continue
        open_access_region_list += gencode_regions[region]
    merged_open_access_region_list = bedmerge(open_access_region_list)

    control_access_region_list = get_conjugate(merged_open_access_region_list, genomeLen)
    create_region_bed(merged_open_access_region_list, bed_dir, 'open_access.'+date_str)
    create_region_bed(control_access_region_list, bed_dir, 'control_access.'+date_str)        
        
if __name__ == "__main__":
    main()