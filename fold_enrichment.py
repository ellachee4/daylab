import pandas as pd
import gffpandas.gffpandas as gffpd
from pybedtools import BedTool
import pybedtools

# bedtools intersect given a single genomic element bed file
def intersect_func(genomic_element, rq):
    intersection = BedTool(genomic_element).sort().intersect(BedTool(rq).sort())
    intersection_df = pd.read_table(intersection.fn, names=['seq_id', 'start', 'end'])
    num_overlaps = intersection_df.shape[0]
    return num_overlaps


# calculate actual overlaps for each genomic element
genomic_elements = ['c-elegans-cds.bed', 'c-elegans-promoter.bed', 'c-elegans-exon.bed', 'c-elegans-five-prime-utr.bed', 'c-elegans-three-prime-utr.bed', 'c-elegans-intron.bed', 'c-elegans-genes.bed', 'c-elegans-mrna.bed', 'c-elegans-ncrna.bed', 'c-elegans-tss.bed',]

genomic_df = pd.DataFrame(genomic_elements, columns=['genomic element'])
genomic_df['overlaps'] = genomic_df['genomic element'].map(lambda x: intersect_func(x, 'BG4_L1_HQ.bed'))

# save dataframe to csv file
genomic_df.to_csv('actual_overlaps.csv', index=False)

# shuffles
def shuffle_intersect_func(genomic_element):
    rq = BedTool('BG4_L1_HQ.bed').sort()
    shuffle1 = gen.shuffle(genome='ce11')
    shuffle2 = shuffle1.shuffle(genome='ce11').sort()
    shuffle3 = shuffle2.shuffle(genome='ce11').sort()
    return intersect_func(shuffle3)

genomic_df['shuffle_overlaps'] = genomic_df['genomic element'].map(shuffle_intersect_func)

# save dataframe to csv file
genomic_df.to_csv('all_overlaps.csv', index=False)