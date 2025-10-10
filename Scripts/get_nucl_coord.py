import pysam
import pandas as pd
import sys
import argparse

def extract_coordinates(input_file, out_file):
    try:
        mt_data = []
        with pysam.AlignmentFile(input_file, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                cigar = read.cigartuples
                if cigar:
                    match, del_count, ins_count = 0, 0, 0
                    for i in range(len(cigar)):
                        if cigar[i][0] == 0:
                            match += cigar[i][1]
                start = read.reference_start + 1
                end = start + match + 1
                if read.has_tag('SA'):
                    mt_data.append({
                        "MT_HEADER": "MT",
                        "MT_START": start,
                        "MT_END": end,
                        "READ_ID": read.query_name,
                        "READ_SA": read.get_tag('SA')
                    })
        sa_df = pd.DataFrame(mt_data)
        sa_temp_data = sa_df['READ_SA'].str.split(';', expand=True)
        sa_temp_data.columns = [f'READ_SA_{i+1}' for i in range(sa_temp_data.shape[1])]
        result = pd.concat([sa_df, sa_temp_data], axis = 1)
        result.drop(columns=["READ_SA"], inplace=True)
        result = result.melt(id_vars = ['MT_HEADER', 'MT_START', 'MT_END', 'READ_ID'], var_name = 'SA_TAG', value_name = 'SA_VALUE')
        result = result[['MT_HEADER', 'MT_START', 'MT_END', 'READ_ID', 'SA_VALUE']]
        result = result.dropna(subset=['SA_VALUE'])
        result = result[result['SA_VALUE'].astype(str).str.strip() != '']
        result = result[~result['SA_VALUE'].str.contains('MT')]
        chr_data = result['SA_VALUE'].str.split(',', n=2, expand=True)
        result["CHR"] = chr_data[0]
        result["START"] = chr_data[1].astype(int)-1
        result["END"] = chr_data[1].astype(int)
        result.drop(columns=["SA_VALUE"], inplace=True)
        result = result.sort_values(by='CHR')
        result = result[['CHR', 'START', 'END', 'READ_ID', 'MT_HEADER', 'MT_START', 'MT_END']]
        result.to_csv(out_file, sep="\t", index=False, header=False)

    except Exception as e:
        print(f"Error occured: {e}")
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract Supplementary Alignment coordinates from BAM/SAM")
    parser.add_argument("-i", "--input", type=str, required=True, help="Input BAM/SAM file")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output file")
    args = parser.parse_args()

    extract_coordinates(args.input, args.output)