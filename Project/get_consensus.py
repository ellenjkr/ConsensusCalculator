import subprocess
import re
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


all_sequences = []
alignment_info = {'Amostra': [], 'Consenso': [], 'Sequências Utilizadas': [], '%': []}

with open('output/log.txt', 'w') as log:
	for file in os.listdir('corrected_orientations'):
		file_path = f'./corrected_orientations/{file}'
		original_file_path = f'./fasta_files/{file}'
		file_name = file.replace('.fa', '')
		file_name = file_name.replace('.fasta', '')
		file_name = file_name.replace('.fas', '')

		if not os.path.exists(f"output/{file_name}"):
			os.makedirs(f"output/{file_name}")

		subprocess.run(
			f'source /home/bioinfo/miniconda3/etc/profile.d/conda.sh && \
			conda activate cdhit-env && \
			cd-hit -i {file_path} -o ./output/{file_name}/clustering.fasta -n 3 -aL 0.8 -aS 0.8 -c 0.7 && \
			conda deactivate',
			shell=True,
			executable='/bin/bash'
		)


		with open(f'output/{file_name}/clustering.fasta.clstr', 'r') as f:
			clustering = f.read()

			clusters = re.findall(r'>Cluster \d*', clustering)

			n_seq = 0
			main_cluster = None
			for cluster in clusters:
				cluster_content = re.findall(rf'({cluster}[\s\S]*)(?=>Cluster)', clustering)
				if len(cluster_content) == 0:
					cluster_content = re.findall(rf'({cluster}[\s\S]*)', clustering)
				
				if len(cluster_content[0]) > n_seq:
					n_seq = len(cluster_content[0])
					main_cluster = cluster_content[0]

			scores = re.findall(r'>(.*)\.\.\. (.*)', main_cluster)
			reference = None
			sequences2align = []
			for score in scores:
				if score[1] == '*':
					reference = score[0]
				else:
					sequences2align.append(score[0])
			# clusters = re.findall(r'(>Cluster .*[\s\S]*)(?=>Cluster)', clustering)
			# if len(main_cluster) == 0:
			# 	main_cluster = re.findall(r'>Cluster 0[\s\S]*', clustering)[0]
			# else:
			# 	main_cluster = main_cluster[0]
			
			# scores = re.findall(r'>(.*)\.\.\. (.*)', main_cluster)

			# reference = None
			# sequences2align = []
			# for score in scores:
			# 	if score[1] == '*':
			# 		reference = score[0]
			# 	else:
			# 		sequences2align.append(score[0])


		sequences = SeqIO.to_dict(SeqIO.parse(file_path, "fasta"))
		reference_seq = sequences[reference]
		SeqIO.write([reference_seq], f'./output/{file_name}/reference.fasta', "fasta")

		sequences = {item[0]: item[1] for item in sequences.items() if item[0] in sequences2align}

		with open(f'./output/{file_name}/sequences2align.fasta', 'w') as handle:
			SeqIO.write(sequences.values(), handle, 'fasta')

		aligned_sequences = len(sequences) + 1
		total_sequences = len(SeqIO.to_dict(SeqIO.parse(original_file_path, "fasta")).keys())

		if aligned_sequences > 1:
			subprocess.run(
				f'mafft --6merpair --addfragments ./output/{file_name}/sequences2align.fasta ./output/{file_name}/reference.fasta > ./output/{file_name}/alignment.fasta',
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
					alignment_info['Amostra'].append(file_name)
					alignment_info['Consenso'].append(record.seq)
					alignment_info['Sequências Utilizadas'].append(aligned_sequences)
					alignment_info['%'].append(aligned_sequences / total_sequences * 100)

			with open(f'output/{file_name}/consensus.fasta', 'w') as file:
				SeqIO.write(sequences, file, "fasta")


		else:
			log.write(f'Error in file {file}\tNo sequences to align\n')



with open(f'output/all_sequences.fasta', 'w') as file:
	SeqIO.write(all_sequences, file, "fasta")

subprocess.run(
	f'export BLASTDB=/media/bioinfo/6tb_hdd/04_Blast_Databases/taxdb && \
	blastn \
	-db /media/bioinfo/6tb_hdd/04_Blast_Databases/BLAST_DB_nt/nt \
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

alignment_info_df = alignment_info_df.merge(samples_statistics, on='Amostra')

def build_excel_sheet(writer, sheet_name, df):
	df.to_excel(writer, sheet_name=sheet_name, startrow=1, header=False, index=False)

	worksheet = writer.sheets[sheet_name]

	(max_row, max_col) = df.shape
	column_settings = [{'header': column} for column in df.columns]
	worksheet.add_table(0, 0, max_row, max_col - 1, {'columns': column_settings})
	worksheet.set_column(0, max_col - 1, 15)

	return (writer, worksheet)


writer = pd.ExcelWriter(f'output/lib3_consensus_statistics.xlsx', engine='xlsxwriter')
writer, worksheet = build_excel_sheet(writer, 'BLAST', blast_df)
writer, worksheet = build_excel_sheet(writer, 'Consenso', alignment_info_df)

writer.save()


# verificar se o alinhamento sem o cdhit nao é melhor. verificar tambem alinhamento pelo muscle


# melhores:
# normal com reverse complement
# reverse com complement

# é uma possibilidade mas não muito boa:
# normal do 1 com normal do 5
# reverse com reverse
# complement com complement
# reverse complement com reverse complement