import yaml
import os
import pandas as pd
import subprocess

def remove_primers(sample_name, primers_df, seq_file, config):
    sample_data = primers_df.loc[primers_df['sample_name'] == sample_name].reset_index(drop=True)
    forward = sample_data['forward'][0]
    reverse = sample_data['reverse'][0]
    min = sample_data['min'][0]
    max = sample_data['max'][0]

    subprocess.run(
        f'bash virtual_pcr.sh -n {sample_name} -o {seq_file} -n {sample_name} -f {forward} -r {reverse} -i {seq_file} -m {min} -M {max} -t {config["THREADS"]}',
        shell=True,
        executable='/bin/bash'
    )