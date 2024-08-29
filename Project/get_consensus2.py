import subprocess
import re
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import pymysql.cursors
import concurrent.futures
import yaml
from copy import deepcopy

from remove_primers import remove_primers


def read_yaml(file_path):
	with open(file_path, "r") as f:
		return yaml.safe_load(f)

def add_clusters_keys(alignment_info, clusters2process):
	for i in range(clusters2process):
		alignment_info[f'Sequências Cluster {i+1}'] = []
		alignment_info[f'% Cluster {i+1}'] = []	

	return alignment_info


def read_primers_file(config):
	for file in os.listdir(config['INPUT_PATH']):
		if '.tsv' in file:
			primers = pd.read_csv(f"{config['INPUT_PATH']}{file}", sep='\t', names=['sample_name', 'forward', 'reverse', 'min', 'max'])

	return primers

def processs_files(config, files_report, alignment_info):
	for file in os.listdir(config['INPUT_PATH']):
		if '.fasta' in file or '.fa' in file or '.fas' in file:
			file_path = f"{config['INPUT_PATH']}{file}"
			files_report, alignment_info = process_file(file, file_path, alignment_info, config)

	return files_report, alignment_info


def process_file(file, file_path, alignment_info, config):
	file_name = file.replace('.fa', '')
	file_name = file_name.replace('.fasta', '')
	file_name = file_name.replace('.fas', '')

	if not os.path.exists(f"{config['OUTPUT_PATH']}{file_name}"):
		os.makedirs(f"{config['OUTPUT_PATH']}{file_name}")
	
	with open(file_path, 'r') as file:
		lines = file.readlines()

	sequences = []
	seq = ''
	count = 1

	for line in lines:
		if line.startswith('>'):
			sequences.append(f'>{count}\n')
			count += 1
		else:
			sequences.append(line)

	with open(file_path, 'w') as file:
		file.write(''.join(sequences))

	sequences = SeqIO.to_dict(SeqIO.parse(file_path, "fasta"))
	total_sequences = len(sequences.keys())
	
	subprocess.run(
		f"source {config['CONDA_SOURCE']} && \
		conda activate cdhit-env && \
		cd-hit -T 10 -i {file_path} -o {config['OUTPUT_PATH']}{file_name}/clustering.fasta -n 3 -aL {config['CLUSTERING_PARAMS']['aL']} -aS {config['CLUSTERING_PARAMS']['aS']} -c {config['CLUSTERING_PARAMS']['c']} && \
		conda deactivate",
		shell=True,
		executable='/bin/bash'
	)


	with open(f"{config['OUTPUT_PATH']}{file_name}/clustering.fasta.clstr", 'r') as f:
		clustering = f.read()

		clusters = re.findall(r'>Cluster \d*', clustering)
		clusters_content = re.split(r'>Cluster \d*\n', clustering)
		clusters_content.remove('')
		clusters_sizes = {}
		for pos, cluster in enumerate(clusters):
			cluster_number = int(cluster[-1])
			cluster_content = clusters_content[cluster_number]
			cluster_seq_count = len(cluster_content.split('\n')) - 1
			if cluster_seq_count >= int(config['MIN_SEQS_FOR_CONSENSUS']):
				clusters_sizes[cluster_number] = cluster_seq_count


		clusters_sizes = {k: v for k, v in sorted(clusters_sizes.items(), reverse=True, key=lambda item: item[1])}

	files_report[file_name] = {'clusters_content': clusters_content, 'clusters_sizes': clusters_sizes, 'sequences': sequences}

	
	alignment_info['Amostra'].append(file_name)
	alignment_info['Total de Sequências'].append(total_sequences)
	alignment_info['Clusters'].append(len(clusters_sizes.keys()))

	for i in range(config['CLUSTERS2PROCESS']):
		if len(list(clusters_sizes.keys())) > i:
			nseqs_cluster = clusters_sizes[list(clusters_sizes.keys())[i]]
			alignment_info[f'% Cluster {i + 1}'].append(nseqs_cluster / total_sequences * 100)
			alignment_info[f'Sequências Cluster {i + 1}'].append(nseqs_cluster)
		else:
			alignment_info[f'% Cluster {i + 1}'].append(0)
			alignment_info[f'Sequências Cluster {i + 1}'].append(0)

	return files_report, alignment_info


def get_consensus(config, files_report, primers_df):
	all_sequences = []
	clusters_reports = [{'Amostra': [], 'Consenso': [], 'Comprimento': []} for i in range(config['CLUSTERS2PROCESS'])]
	for cluster_report_pos in range(len(clusters_reports)):
		for file, report in files_report.items():
			clusters_content = report['clusters_content']
			clusters_sizes = report['clusters_sizes']
			sequences = report['sequences']
			if cluster_report_pos < len(clusters_sizes):
				cluster = list(clusters_sizes.keys())[cluster_report_pos]
				
				cluster_content = clusters_content[cluster]
				scores = re.findall(r'>(.*)\.\.\. (.*)', cluster_content)
				sequences2align = []
				for score in scores:
					sequences2align.append(score[0])

				cluster_sequences = {item[0]: item[1] for item in sequences.items() if item[0] in sequences2align}
				with open(f"{config['OUTPUT_PATH']}{file}/cluster{cluster_report_pos + 1}_sequences2align.fasta", 'w') as handle:
					SeqIO.write(cluster_sequences.values(), handle, 'fasta')

				if len(cluster_sequences) > 2:
					subprocess.run(
						f"mafft --quiet --thread 10 {config['OUTPUT_PATH']}{file}/cluster{cluster_report_pos + 1}_sequences2align.fasta > {config['OUTPUT_PATH']}{file}/cluster{cluster_report_pos + 1}_alignment.fasta",
						shell=True,
						executable='/bin/bash'
					)

					subprocess.run(
						f"source {config['CONDA_SOURCE']} && \
						conda activate emboss-env && \
						cons {config['OUTPUT_PATH']}{file}/cluster{cluster_report_pos + 1}_alignment.fasta {config['OUTPUT_PATH']}{file}/cluster{cluster_report_pos + 1}_consensus.fasta && \
						conda deactivate",
						shell=True,
						executable='/bin/bash'
					)


					consensus_sequences = []
					with open(f"{config['OUTPUT_PATH']}{file}/cluster{cluster_report_pos + 1}_consensus.fasta", 'r') as consensus_file:
						for record in SeqIO.parse(consensus_file, "fasta"):
							record.seq = Seq(str(record.seq).replace('\n', ''))
							record.seq = Seq(str(record.seq).replace('n', ''))
							record.seq = Seq(str(record.seq).replace('N', ''))
							record.id = f"cluster{cluster_report_pos + 1}@{file}"
							consensus_sequences.append(record)

					with open(f"{config['OUTPUT_PATH']}{file}/cluster{cluster_report_pos + 1}_consensus.fasta", 'w') as consensus_file:
						SeqIO.write(consensus_sequences, consensus_file, "fasta")

					if config['REMOVE_PRIMERS']:
						remove_primers(file, primers_df, f"{config['OUTPUT_PATH']}{file}/cluster{cluster_report_pos + 1}_consensus.fasta", config)
					
					with open(f"{config['OUTPUT_PATH']}{file}/cluster{cluster_report_pos + 1}_consensus.fasta", 'r') as consensus_file:
						for record in SeqIO.parse(consensus_file, "fasta"):
							all_sequences.append(record)

							clusters_reports[cluster_report_pos]['Amostra'].append(file)
							clusters_reports[cluster_report_pos]['Consenso'].append(record.seq)
							clusters_reports[cluster_report_pos]['Comprimento'].append(len(record.seq))
							


	return clusters_reports, all_sequences


def run_blast(config, all_sequences, clusters_reports):
	clusters_reports_copy = deepcopy(clusters_reports)

	with open(f"{config['OUTPUT_PATH']}all_sequences.fasta", 'w') as file:
		SeqIO.write(all_sequences, file, "fasta")


	subprocess.run(
		f'export BLASTDB={config["TAX_DATABASE_PATH"]} && \
		blastn \
		-db {config["BLAST_DATABASE_PATH"]} \
		-query {config["OUTPUT_PATH"]}all_sequences.fasta \
		-outfmt "6 qseqid length score bitscore pident nident evalue gapopen gaps qcovs qcovhsp stitle sscinames mismatch qstart qend sstart send" \
		-out {config["OUTPUT_PATH"]}blast.tsv \
		-num_threads 6 \
		-max_target_seqs 10',
		shell=True,
		executable='/bin/bash'
	)


	alignment_info_df = pd.DataFrame(alignment_info)
	blast_df = pd.read_csv(f"{config['OUTPUT_PATH']}blast.tsv", sep='\t', names=['Amostra', 'length', 'score', 'bitscore', 'pident', 'nident', 'evalue', 'gapopen', 'gaps', 'qcovs', 'qcovhsp', 'stitle', 'sscinames', 'mismatch', 'qstart', 'qend', 'sstart', 'send'])
	blast_df[['Cluster', 'Amostra']] = blast_df['Amostra'].str.split('@', n=1, expand=True)

	return blast_df, alignment_info_df


def get_organism_taxid(organism, cursor):
	cursor.execute(f"SELECT * FROM organisms WHERE tax_name = '{organism}' LIMIT 1")
	result = cursor.fetchone()
	if result is None:
		cursor.execute(f"SELECT * FROM names WHERE name_txt = '{organism}' LIMIT 1")
		result = cursor.fetchone()
		if result is None:
			cursor.execute(f"SELECT * FROM names WHERE MATCH(name_txt) AGAINST('{organism}' IN NATURAL LANGUAGE MODE)")
			result = cursor.fetchone()

		taxid = result['tax_id']
		name_txt = result['name_txt']
	else:
		taxid = result['tax_id']
		name_txt = result['tax_name']
 
	return taxid, name_txt


def get_taxonomy(config, especie):
	# Connect to the database
	connection = pymysql.connect(
		host=config['HIERARCHY_DATABASE']['HOST'],
		user=config['HIERARCHY_DATABASE']['USER'],
		password=config['HIERARCHY_DATABASE']['PASSWORD'],
		database=config['HIERARCHY_DATABASE']['DATABASE'],
		cursorclass=pymysql.cursors.DictCursor

	)
	taxonomies = {'sscinames': None, 'Gênero': None, 'Família': None, 'Ordem': None, 'Classe': None, 'Filo': None, 'Super-Reino': None}
	with connection:
		with connection.cursor() as cursor:
			try:
				taxid, name_txt = get_organism_taxid(especie, cursor)
				
				cursor.execute(f"SELECT * FROM organisms WHERE tax_id = {taxid} LIMIT 1")
				# cursor.execute(f"SELECT * FROM organisms WHERE MATCH(tax_name) AGAINST('{especie}' IN NATURAL LANGUAGE MODE)")
				result = cursor.fetchone()

				taxonomies['sscinames'] = especie
				taxonomies['Gênero'] = result['genus']
				taxonomies['Família'] = result['family']
				taxonomies['Ordem'] = result['_order']
				taxonomies['Classe'] = result['class']
				taxonomies['Filo'] = result['phylum']
				taxonomies['Super-Reino'] = result['superkingdom']
				
			except Exception as e:				
				taxonomies['sscinames'] = especie
				taxonomies['Gênero'] = 'unclassified'
				taxonomies['Família'] = 'unclassified'
				taxonomies['Ordem'] = 'unclassified'
				taxonomies['Classe'] = 'unclassified'
				taxonomies['Filo'] = 'unclassified'
				taxonomies['Super-Reino'] = 'unclassified'
				
		
	return taxonomies


def build_excel_sheet(writer, sheet_name, df):
	df.to_excel(writer, sheet_name=sheet_name, startrow=1, header=False, index=False)

	worksheet = writer.sheets[sheet_name]

	(max_row, max_col) = df.shape
	column_settings = [{'header': column} for column in df.columns]
	worksheet.add_table(0, 0, max_row, max_col - 1, {'columns': column_settings})
	worksheet.set_column(0, max_col - 1, 15)

	return (writer, worksheet)


def build_reports_sheets(config, clusters_reports, blast_df, alignment_info_df):
	clusters_reports = [pd.DataFrame(cluster_report) for cluster_report in clusters_reports]

	clusters_blast = blast_df.groupby('Cluster')
	writer = pd.ExcelWriter(f"{config['OUTPUT_PATH']}{config['OUTPUT_FILENAME']}.xlsx", engine='xlsxwriter')
	writer, worksheet = build_excel_sheet(writer, f'INFO', alignment_info_df)

	for group in clusters_blast.groups:
		cluster_blast_df = clusters_blast.get_group(group)
		cluster_pos = int(group.replace('cluster', '')) - 1
		samples_statistics = cluster_blast_df.groupby('Amostra')[['pident', 'qcovs', 'gaps']].mean().reset_index()
		samples_statistics.rename(columns={'pident': 'pident mean', 'qcovs': 'qcovs mean', 'gaps': 'gaps mean'}, inplace=True)

		most_frequent_organisms = cluster_blast_df.groupby('Amostra')['sscinames'].agg(pd.Series.mode).reset_index()
		most_frequent_organisms.rename(columns={'sscinames': 'sscinames mode'}, inplace=True)

		samples_statistics = samples_statistics.merge(most_frequent_organisms, on='Amostra')

		first_organisms = cluster_blast_df.groupby('Amostra').first().reset_index()
		first_organisms.rename(columns={'sscinames': 'Primeiro Resultado', 'pident': 'Primeiro Resultado pident', 'qcovs': 'Primeiro Resultado qcovs', 'gaps': 'Primeiro Resultado gaps'}, inplace=True)
		first_organisms = first_organisms[['Amostra', 'Primeiro Resultado', 'Primeiro Resultado pident', 'Primeiro Resultado qcovs', 'Primeiro Resultado gaps']]

		clusters_reports[cluster_pos] = clusters_reports[cluster_pos].merge(alignment_info_df[['Amostra', f'% Cluster {cluster_pos + 1}', f'Sequências Cluster {cluster_pos + 1}']], on='Amostra', how='left')
		seq_count = clusters_reports[cluster_pos].pop(f'Sequências Cluster {cluster_pos + 1}')
		clusters_reports[cluster_pos].insert(1, f'Sequências Cluster {cluster_pos + 1}', seq_count)
		seq_perc = clusters_reports[cluster_pos].pop(f'% Cluster {cluster_pos + 1}')
		clusters_reports[cluster_pos].insert(2, f'% Cluster {cluster_pos + 1}', seq_perc)

		clusters_reports[cluster_pos] = clusters_reports[cluster_pos].merge(first_organisms, on='Amostra')
		clusters_reports[cluster_pos] = clusters_reports[cluster_pos].merge(samples_statistics, on='Amostra')
	

		unique_sscinames = cluster_blast_df['sscinames'].unique()

		with concurrent.futures.ThreadPoolExecutor(max_workers=24) as executor:
			futures = [executor.submit(get_taxonomy, config, ssciname) for ssciname in unique_sscinames]
			results = [future.result() for future in futures]
			taxonomies_df = pd.json_normalize(results)

		cluster_blast_df = cluster_blast_df.merge(taxonomies_df, on='sscinames', how='left').fillna('unclassified')


		cluster_blast_df = cluster_blast_df[['Amostra', 'pident', 'qcovs', 'gaps', 'stitle', 'sscinames', 'Gênero', 'Família', 'Ordem', 'Classe', 'Filo', 'Super-Reino', 'length', 'score', 'bitscore', 'nident', 'evalue', 'gapopen', 'qcovhsp', 'mismatch', 'qstart', 'qend', 'sstart', 'send']]

		writer, worksheet = build_excel_sheet(writer, f'BLAST Cluster {cluster_pos + 1}', cluster_blast_df)
		writer, worksheet = build_excel_sheet(writer, f'Consenso Cluster {cluster_pos + 1}', clusters_reports[cluster_pos])

	writer.close()


if __name__ == "__main__":
	config = read_yaml("config.yaml")
	all_sequences = []
	alignment_info = {'Amostra': [], 'Total de Sequências': [], 'Clusters': []}
	files_report = {}

	primers_df = None
	if config['REMOVE_PRIMERS']:
		primers_df = read_primers_file(config)

	alignment_info = add_clusters_keys(alignment_info, config['CLUSTERS2PROCESS'])
	files_report, alignment_info = processs_files(config, files_report, alignment_info)
	clusters_reports, all_sequences = get_consensus(config, files_report, primers_df)
	blast_df, alignment_info_df = run_blast(config, all_sequences, clusters_reports)
	build_reports_sheets(config, clusters_reports, blast_df, alignment_info_df)
