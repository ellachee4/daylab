import pandas as pd
import csv

# File paths and names
chambers_path = '/work/daylab/AV_for_coops/Chambers_data/'
chambers_peak_path = 'work/daylab/fall-2024-coops/Chambers_data/'

chambers_files = ['GSM3003540_Homo_all_w15_th-1_minus.hits.max.PDS.w50.35.bed',
                  'GSM3003540_Homo_all_w15_th-1_minus.PDS.bw',
                  'GSM3003540_Homo_all_w15_th-1_plus.hits.max.PDS.w50.35.bed',
                  'GSM3003540_Homo_all_w15_th-1_plus.PDS.bw']

# Calculate peak length for each peak in single bed file
def calc_peak_lengths(input_file_name):
    '''Calculate peak length and add column for each peak in bed file'''

    output_file_name = input_file_name.replace('.bed', '_lengths.bed')
    output_file_name.replace('AV_for_coops', 'fall-2024-coops')
    with open(input_file_name, 'r') as infile, open(output_file_name, 'w') as outfile:
        writer = csv.writer(outfile)
        for row in infile:
            columns = row.strip().split('\t')
            start = int(columns[1])
            end = int(columns[2])
            difference = end - start
            # add difference column
            writer.writerow(row+[difference])
    
    print(f'Calculated peak length for {input_file_name}')
    return output_file_name

# Calculate average peak length for single bed file
def calc_average_peak_length(input_file_name):
    '''Removes lines from a G4catchall BED file where the start site is greater than the end site'''

    with open(input_file_name, 'r') as infile:
        length = len(infile.readlines())
        sum_peaks = 0
        for row in infile:
            columns = row.strip().split('\t')
            difference = int(columns[4])
            sum_peaks += difference
    
    average_peak_length = sum_peaks / length

    print(f'Average peak length for {input_file_name}: {average_peak_length} bp')
    return average_peak_length

# Calculate average peak length for each bed file in list
def calc_all_avg_peaks(i):
    '''Calculate average peak length for given bed file'''
    bed_file = chambers_files[i]
    bed_file_name = chambers_path + bed_file

    # Calculate peak lengths for each peak in file
    peak_lengths_file = calc_peak_lengths(bed_file_name)
    
    # Calculate average peak length for file
    avg_peak_length = calc_average_peak_length(peak_lengths_file)
    return avg_peak_length

# For all files in list
avg_peak_lengths = {}
for i in range(len(chambers_files)):
    avg_peak_length = calc_all_avg_peaks(i)
    bed_file = chambers_files[i]
    avg_peak_lengths.update({bed_file: avg_peak_length})
avg_peak_lengths.to_csv('work/daylab/fall-2024-coops/Chambers_data/average_peak_lengths_chambers_data.csv')
print('Average peak length of each file saved to average_peak_lengths_chambers_data.csv')