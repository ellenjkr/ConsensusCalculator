# import pandas as pd
# import matplotlib.pyplot as plt

# df = pd.read_excel('output/lib4/lib4_ani1.xlsx', sheet_name='Consenso')
# # Criando boxplots
# columns_to_plot = ['Clusters', '% Cluster 1', 'Comprimento', 'Primeiro Resultado pident', 'Primeiro Resultado qcovs', 'Primeiro Resultado gaps']

# plt.figure(figsize=(12, 8))

# for i, column in enumerate(columns_to_plot, 1):
#     plt.subplot(2, 3, i)
#     plt.boxplot(df[column])
#     plt.title(f'{column}')
#     plt.grid()
#     #  plt.xlabel(column)
#     plt.xticks([1], ['ANI1'])  # Define o rótulo no eixo x

# plt.tight_layout()
# plt.show()


import pandas as pd
import matplotlib.pyplot as plt

# Lista de nomes dos arquivos Excel
file_paths = [
    'output/LIB6/lib6_ani_consensus.xlsx'
]

xticks = ['ANI']

# Lista para armazenar os DataFrames de cada arquivo Excel
dfs = []

# Lê cada arquivo Excel e armazena os DataFrames na lista dfs
for file_path in file_paths:
    df = pd.read_excel(file_path, sheet_name='Consenso')
    dfs.append(df)

# Criando boxplots
columns_to_plot = ['Clusters', '% Cluster 1', 'Comprimento', 'Primeiro Resultado pident', 'Primeiro Resultado qcovs', 'Primeiro Resultado gaps']
# columns_to_plot = ['Primeiro Resultado pident', 'Primeiro Resultado qcovs', 'Primeiro Resultado gaps']

plt.figure(figsize=(18, 12))

for i, column in enumerate(columns_to_plot, 1):
    plt.subplot(2, 3, i)
    data_to_plot = [df[column] for df in dfs]
    plt.boxplot(data_to_plot)
    column = column.replace('pident', 'Identidade')
    column = column.replace('qcovs', 'Query Coverage')
    column = column.replace('gaps', 'Gaps')
    plt.title(f'{column}')
    plt.grid()
    plt.xticks(range(1, len(file_paths) + 1), xticks)  # Define os rótulos no eixo x

plt.tight_layout(pad=5.0)
plt.show()
