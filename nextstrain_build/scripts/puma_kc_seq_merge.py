import pandas as pd
import numpy as np
import argparse


def load_metadata(file):
    '''
    Loads metadata tsv as df.
    '''
    with open(file) as tfile:
        metadata = pd.read_csv(tfile, sep = '\t')
    return metadata


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description= adds PUMA to King County sequences ,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--metadata', type=str, required=True, help= downloaded and unzipped metadata file )
    parser.add_argument('--puma_raw', type=str, required=True, help= tsv file containing additional puma metadata for king county seqs )
    parser.add_argument('--output', type=str, required=True, help= new metadata file )

    args = parser.parse_args()


metadata = load_metadata(args.metadata)
additional = load_metadata(args.puma_raw)


new_merged = pd.merge(metadata, additional, how= left , on=  strain )

new_merged.location[new_merged.COUNTY == King ] =  King County 

north_kc = [11601.0, 11602.0, 11603.0, 11604.0, 11605.0, 11606.0, 11607.0 , 11608.0 , 11609.0 , 11616.0]
south_kc = [11610.0, 11611.0, 11612.0, 11613.0, 11614.0, 11615.0]

new_merged.Puma[(new_merged.division == 'Washington') & (new_merged.location == 'King County') & (new_merged.Puma.isna())] = 'King County'
new_merged.Puma[(new_merged.division == 'Washington') & (new_merged.Puma.isna())] = 'Washington'
new_merged.Puma[(new_merged.Puma.isna())] = new_merged.region

new_merged['ns_kc'] = np.nan

new_merged.ns_kc[(new_merged.location == 'King County') & (new_merged.Puma.isin(north_kc))] = 'North_King_County'
new_merged.ns_kc[(new_merged.location == 'King County') & (new_merged.Puma.isin(south_kc))] = 'South_King_County'
new_merged.ns_kc[(new_merged.division == 'Washington') & (new_merged.location == 'King County') & ((new_merged.ns_kc.isna()))] = 'Other_King_County'
new_merged.ns_kc[(new_merged.division == 'Washington') & ((new_merged.ns_kc.isna()))] = 'Washington'
new_merged.ns_kc[(new_merged.ns_kc.isna())] = new_merged.region

new_merged['ns_kc'] = new_merged['ns_kc'].astype(str)

with open(args.output, 'w') as f:
    new_merged.to_csv(f, sep = '\t')
