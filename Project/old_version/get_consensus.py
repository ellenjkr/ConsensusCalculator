import subprocess
import re
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import pymysql.cursors
import concurrent.futures


all_sequences = []
alignment_info = {'Amostra': [], 'Consenso': [], 'Comprimento': [], 'Clusters': [], 'Total de Sequências': [], 'Sequências Cluster 1': [], '% Cluster 1': [], 'Sequências Cluster 2': [], '% Cluster 2': [], 'Sequências Cluster 3': [], '% Cluster 3': []}

with open('output/log.txt', 'w') as log:
	for file in os.listdir('input/lib7/bac/'):
		file_path = f'./input/lib7/bac/{file}'
		# original_file_path = f'./fasta_files/{file}'
		file_name = file.replace('.fa', '')
		file_name = file_name.replace('.fasta', '')
		file_name = file_name.replace('.fas', '')

		if not os.path.exists(f"output/{file_name}"):
			os.makedirs(f"output/{file_name}")

		
		
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

		subprocess.run(
			f'source /home/bioinfo/miniconda3/etc/profile.d/conda.sh && \
			conda activate cdhit-env && \
			cd-hit -T 10 -i {file_path} -o ./output/{file_name}/clustering.fasta -n 3 -aL 0.8 -aS 0.8 -c 0.7 && \
			conda deactivate',
			shell=True,
			executable='/bin/bash'
		)


		with open(f'output/{file_name}/clustering.fasta.clstr', 'r') as f:
			clustering = f.read()

			clusters = re.findall(r'>Cluster \d*', clustering)
			clusters_content = re.split(r'>Cluster \d*\n', clustering)
			clusters_content.remove('')
			n_seq = 0
			clusters_sizes = {}
			for pos, cluster in enumerate(clusters):
				cluster_number = int(cluster[-1])
				cluster_content = clusters_content[cluster_number]
				cluster_seq_count = len(cluster_content.split('\n')) - 1
				if cluster_seq_count >= 5:
					clusters_sizes[cluster_number] = cluster_seq_count

			if len(clusters_sizes.keys()) > 0:
				clusters_sizes = {k: v for k, v in sorted(clusters_sizes.items(), reverse=True, key=lambda item: item[1])}
				n_seq = clusters_sizes[list(clusters_sizes.keys())[0]]
				main_cluster = clusters_content[list(clusters_sizes.keys())[0]]

				scores = re.findall(r'>(.*)\.\.\. (.*)', main_cluster)
				sequences2align = []
				for score in scores:
					sequences2align.append(score[0])

		if len(clusters_sizes.keys()) > 0:
			sequences = SeqIO.to_dict(SeqIO.parse(file_path, "fasta"))
			total_sequences = len(sequences.keys())

			sequences = {item[0]: item[1] for item in sequences.items() if item[0] in sequences2align}
			with open(f'./output/{file_name}/sequences2align.fasta', 'w') as handle:
				SeqIO.write(sequences.values(), handle, 'fasta')

			aligned_sequences = len(sequences) + 1

			if aligned_sequences > 1:
				subprocess.run(
					f'mafft --quiet --thread 10 ./output/{file_name}/sequences2align.fasta > ./output/{file_name}/alignment.fasta',
					shell=True,
					executable='/bin/bash'
				)


				subprocess.run(
					f'source /home/bioinfo/miniconda3/etc/profile.d/conda.sh && \
					conda activate emboss-env && \
					cons ./output/{file_name}/alignment.fasta ./output/{file_name}/consensus.fasta && \
					conda deactivate',
					shell=True,
					executable='/bin/bash'
				)

				sequences = []
				with open(f'output/{file_name}/consensus.fasta', 'r') as file:
					for record in SeqIO.parse(file, "fasta"):
						record.seq = Seq(str(record.seq).replace('\n', ''))
						record.seq = Seq(str(record.seq).replace('n', ''))
						record.seq = Seq(str(record.seq).replace('N', ''))
						record.id = file_name
						sequences.append(record)
						all_sequences.append(record)

						nseqs_cluster1 = clusters_sizes[list(clusters_sizes.keys())[0]]
						if len(list(clusters_sizes.keys())) > 1:
							nseqs_cluster2 = clusters_sizes[list(clusters_sizes.keys())[1]]
							alignment_info['Sequências Cluster 2'].append(nseqs_cluster2)
							alignment_info['% Cluster 2'].append(nseqs_cluster2 / total_sequences * 100)
						else:
							alignment_info['Sequências Cluster 2'].append(0)
							alignment_info['% Cluster 2'].append(0)

						if len(list(clusters_sizes.keys())) > 2:
							nseqs_cluster3 = clusters_sizes[list(clusters_sizes.keys())[2]]
							alignment_info['Sequências Cluster 3'].append(nseqs_cluster3)
							alignment_info['% Cluster 3'].append(nseqs_cluster3 / total_sequences * 100)
						else:
							alignment_info['Sequências Cluster 3'].append(0)
							alignment_info['% Cluster 3'].append(0)

						alignment_info['Amostra'].append(file_name)
						alignment_info['Consenso'].append(record.seq)
						alignment_info['Comprimento'].append(len(record.seq))
						alignment_info['Clusters'].append(len(clusters_sizes.keys()))
						alignment_info['Total de Sequências'].append(total_sequences)
						alignment_info['Sequências Cluster 1'].append(nseqs_cluster1)
						alignment_info['% Cluster 1'].append(nseqs_cluster1 / total_sequences * 100)
						

				with open(f'output/{file_name}/consensus.fasta', 'w') as file:
					SeqIO.write(sequences, file, "fasta")


		else:
			log.write(f'Error in file {file}\tNot enough sequences to align\n')



with open(f'output/all_sequences.fasta', 'w') as file:
	SeqIO.write(all_sequences, file, "fasta")

subprocess.run(
	f'export BLASTDB=/media/bioinfo/6tb_hdd/04_Blast_Databases/taxdb && \
	blastn \
	-db /media/bioinfo/6tb_hdd/04_Blast_Databases/16sDB/16S_ribosomal_RNA \
	-query ./output/all_sequences.fasta \
	-outfmt "6 qseqid length score bitscore pident nident evalue gapopen gaps qcovs qcovhsp stitle sscinames mismatch qstart qend sstart send" \
	-out ./output/blast.tsv \
	-num_threads 6 \
	-max_target_seqs 10',
	shell=True,
	executable='/bin/bash'
)


alignment_info_df = pd.DataFrame(alignment_info)
blast_df = pd.read_csv('output/blast.tsv', sep='\t', names=['Amostra', 'length', 'score', 'bitscore', 'pident', 'nident', 'evalue', 'gapopen', 'gaps', 'qcovs', 'qcovhsp', 'stitle', 'sscinames', 'mismatch', 'qstart', 'qend', 'sstart', 'send'])


samples_statistics = blast_df.groupby('Amostra')[['pident', 'qcovs', 'gaps']].mean().reset_index()
samples_statistics.rename(columns={'pident': 'pident mean', 'qcovs': 'qcovs mean', 'gaps': 'gaps mean'}, inplace=True)

most_frequent_organisms = blast_df.groupby('Amostra')['sscinames'].agg(pd.Series.mode).reset_index()
most_frequent_organisms.rename(columns={'sscinames': 'sscinames mode'}, inplace=True)

samples_statistics = samples_statistics.merge(most_frequent_organisms, on='Amostra')

first_organisms = blast_df.groupby('Amostra').first().reset_index()
first_organisms.rename(columns={'sscinames': 'Primeiro Resultado', 'pident': 'Primeiro Resultado pident', 'qcovs': 'Primeiro Resultado qcovs', 'gaps': 'Primeiro Resultado gaps'}, inplace=True)
first_organisms = first_organisms[['Amostra', 'Primeiro Resultado', 'Primeiro Resultado pident', 'Primeiro Resultado qcovs', 'Primeiro Resultado gaps']]

alignment_info_df = alignment_info_df.merge(samples_statistics, on='Amostra')

alignment_info_df = alignment_info_df.merge(first_organisms, on='Amostra')


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


def get_taxonomy(especie):
	# Connect to the database
	connection = pymysql.connect(
		host='localhost',
		user='root',
		password='Amora#1000',
		database='ncbi_data',
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

unique_sscinames = blast_df['sscinames'].unique()

with concurrent.futures.ThreadPoolExecutor(max_workers=24) as executor:
	futures = [executor.submit(get_taxonomy, ssciname) for ssciname in unique_sscinames]
	results = [future.result() for future in futures]
	taxonomies_df = pd.json_normalize(results)

blast_df = blast_df.merge(taxonomies_df, on='sscinames', how='left').fillna('unclassified')


blast_df = blast_df[['Amostra', 'pident', 'qcovs', 'gaps', 'stitle', 'sscinames', 'Gênero', 'Família', 'Ordem', 'Classe', 'Filo', 'Super-Reino', 'length', 'score', 'bitscore', 'nident', 'evalue', 'gapopen', 'qcovhsp', 'mismatch', 'qstart', 'qend', 'sstart', 'send']]


def build_excel_sheet(writer, sheet_name, df):
	df.to_excel(writer, sheet_name=sheet_name, startrow=1, header=False, index=False)

	worksheet = writer.sheets[sheet_name]

	(max_row, max_col) = df.shape
	column_settings = [{'header': column} for column in df.columns]
	worksheet.add_table(0, 0, max_row, max_col - 1, {'columns': column_settings})
	worksheet.set_column(0, max_col - 1, 15)

	return (writer, worksheet)


writer = pd.ExcelWriter(f'output/lib7_ani_85perc.xlsx', engine='xlsxwriter')
writer, worksheet = build_excel_sheet(writer, 'BLAST', blast_df)
writer, worksheet = build_excel_sheet(writer, 'Consenso', alignment_info_df)

writer.close()

