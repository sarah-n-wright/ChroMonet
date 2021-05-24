import os
import pandas as pd
import numpy as np
import argparse

import plot_karyogram_LB.PlotKaryogram as PlotKaryogram 

def convert_tsv_to_beds(filename, genome_tsv_file, output_dir):
    """
    Convert genome tsv file (simulated or otherwise) to bed files
    Original: 
    POS     A1      A2      POP1    POP2
    0       C       T       0       0
    
    Output:
    Copy 1: CHR    START_POS  END_POS  POP1
    Copy 2: CHR    START_POS  END_POS  POP2
    """
    df = pd.read_csv(genome_tsv_file, sep=('\t'))
    
    # Insert chromosome column
    df.insert(0,'CHR','21')
    
    # Insert stop position for bed file
    df.insert(2,'STOP_POS',df['POS'] + 1)
    df.rename(columns={'POS':'START_POS'}, inplace=True)
    
    # Split into chromosome copy
    copy1 = df[['CHR', 'START_POS', 'STOP_POS', 'POP1']]
    copy2 = df[['CHR', 'START_POS', 'STOP_POS', 'POP2']]
    
    # Write to new files
    copy1_filename = filename.split('.')[0] + '_copy1.bed'
    copy1.to_csv(output_dir + copy1_filename, sep='\t', index=False)
    
    copy2_filename = filename.split('.')[0] + '_copy2.bed'
    copy2.to_csv(output_dir + copy2_filename, sep='\t', index=False) 
    
    return(copy1_filename, copy2_filename)

def merge_regions(bed_file, chrom, output_prefix, pop_header):
    """
    Merge positions where the adjacent population labels are the same.
    Output bed file formatted for chromosome painting processing.
    
    Input:
    chr21 1  2 AFR
    chr21 5  6 AFR
    chr21 8  9 EUR 
    Output:
    chr21 1 6 AFR
    chr21 8 9 EUR
    """
    df = pd.read_csv(bed_file, sep='\t')
    
    output_file = output_prefix + ".merged.bed"
    out = open(output_prefix + ".merged.bed", 'w')
    
    pop_map = {0: "AFR", 1:"EUR"}
    pop_range = []
    pop=0
    for index,row in df.iterrows():
        #print("Index", index, df.shape[0])
        if index < df.shape[0] - 1: 
            # If current position has the same population, add to range
            if row[pop_header] == pop:
                pop_range.append(row['START_POS'])
            else:
                # Write previous pop range to merged bed file
                if pop_range != []:
                    output_line = [str(chrom), str(pop_range[0]), str(pop_range[-1] + 1), str(pop_map[pop])]
                    out.write('\t'.join(output_line) + '\n')

                # Start new range and switch pop
                if pop == 0:
                    pop = 1
                elif pop == 1:
                    pop = 0
                pop_range = []

                # Add current position
                pop_range.append(row['START_POS'])

        else:
            if row[pop_header] == pop:
                pop_range.append(row['START_POS'])
                output_line = [str(chrom), str(pop_range[0]), str(pop_range[-1] + 1), str(pop_map[pop])]
                out.write('\t'.join(output_line) + '\n')
            else:
                # Write previous pop range to merged bed file
                if pop_range != []:
                    output_line = [str(chrom), str(pop_range[0]), str(pop_range[-1] + 1), str(pop_map[pop])]
                    out.write('\t'.join(output_line) + '\n')
                
                # Start new range and switch pop
                if pop == 0:
                    pop = 1
                elif pop == 1:
                    pop = 0
                pop_range = []

                # Add current position
                pop_range.append(row['START_POS'])
                output_line = [str(chrom), str(pop_range[0]), str(pop_range[-1] + 1), str(pop_map[pop])]
                out.write('\t'.join(output_line) + '\n')
    return output_file 

def plot_karyogram(bed_a_formatted, bed_b_formatted, plot_filename_title, outdir, pop_list):
    """
    Plot Chromosome painting/kayrogram
    Input: 1 Bedfile for each copy of the chromosome: 
    CHR StartPos StopPos Population StartCM StopCM
    """
    outprefix = outdir + plot_filename_title + '.png'
    PlotKaryogram(bed_a_formatted, bed_b_formatted, plot_filename_title, pop_list, outprefix)

def format_bed_files_for_karyogram(bed_a, bed_b, bed_a_formatted, bed_b_formatted, chromosome_length):
    """
    Format bed file for Karyogram script by adding 'cM' start and stop positions 
    by dividing the start and stop position by the total length of chromosome 21
    then multiplying by 64.6, the length of the chromosome in cM
    """
    # Chromosome21: 48129895 bps
    #chromosome_length = 48129895
    
    # Generate fake cM values knowing the 
    # chromosome 21 centromere range from the CENTROMERE file: 
    # 0.075990715724  64.6258250955
    # First Copy of Chromosome
    bed_a_df = pd.read_csv(bed_a, sep='\t', header=None)
    bed_a_df['start_cM'] = (64.6 * (bed_a_df[1] + 1) / chromosome_length)
    bed_a_df['stop_cM'] = (64.6 * (bed_a_df[2] + 1) / chromosome_length)
    bed_a_df=bed_a_df.astype(str)
    display(bed_a_df.head())
    bed_a_df.to_csv(bed_a_formatted, header=False, sep='\t', index=False)

    # Second Copy of Chromosome
    bed_b_df = pd.read_csv(bed_b, sep='\t', header=None)
    bed_b_df['start_cM'] = (64.6 * (bed_b_df[1] + 1) / chromosome_length)
    bed_b_df['stop_cM'] = (64.6 * (bed_b_df[2] + 1)/ chromosome_length)
    bed_b_df=bed_b_df.astype(str)
    display(bed_b_df.head())
    bed_b_df.to_csv(bed_b_formatted, header=False, sep='\t', index=False)


def main():

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--input_dir', help='Input directory of tsv file from HMM', required=True)
    parser.add_argument('--output_dir', help='Output dir for intermediate files and plots', required=True)
    parser.add_argument('--tsv_file', help='TSV file name (not including directory)', required=True)
    parser.add_argument('--plot_tite', help='Title string for chromosome painting', required=True)
    parser.add_argument('--chromosome_length', help='Length of the chromosome e.g. synthetic = 100, actual Chr21 = 48129895', required=True)
    parser.add_argument('--pop_list', nargs='+', help='Population list, default = ["AFR","EUR"]', default =["AFR","EUR"], required=True)
    args = parser.parse_args()

    # 1. Convert TSV files to BED files
    genome_tsv_file = args.input_dir + args.tsv_file
    filename = os.path.basename(genome_tsv_file)

    (copy1_filename, copy2_filename) = convert_tsv_to_beds(filename, genome_tsv_file, args.output_dir)

    # 2. Merge the file regions for each separate chromosome copy 
    copy1_bed = merge_regions(args.output_dir + copy1_filename, '21', args.output_dir + copy1_filename.split('.')[0], 'POP1')
    copy2_bed = merge_regions(args.output_dir + copy2_filename, '21', args.output_dir + copy2_filename.split('.')[0], 'POP2') 

    # 3. Format bed file for karyogram
    # Example bed line: 21      9411245 48119700        EUR
    # File names with fake centimorgans to the input bed files
    copy_1_formatted = args.output_dir + os.path.basename(copy1_bed)
    copy_2_formatted = args.output_dir + os.path.basename(copy2_bed) 

    # 4. Format bed files for chromosome painting input
    format_bed_files_for_karyogram(copy1_bed, copy2_bed, copy_1_formatted, copy_2_formatted, args.chromosome_length)

    # 5. Plot Chromosome 21 Chromosome Painting
    #plot_filename_title = "Admixed Chromosome 21 Painting"
    #pop_list = ["AFR","EUR"]
    plot_filename_title = args.plot_title
    pop_list =  args.pop_list

    plot_karyogram(copy_1_formatted, copy_2_formatted, plot_filename_title, args.output_dir, pop_list) 

    print(args.pop_list)


if __name__ == "__main__":
    main()
    
