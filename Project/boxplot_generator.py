import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_excel('output/lib3_consensus_statistics.xlsx', sheet_name='Consenso')
print(df.columns)
# Criando boxplots
columns_to_plot = ['%', 'Comprimento', 'Primeiro Resultado pident', 'Primeiro Resultado qcovs', 'Primeiro Resultado gaps', 'Clusters']

plt.figure(figsize=(12, 8))

for i, column in enumerate(columns_to_plot, 1):
    plt.subplot(2, 3, i)
    plt.boxplot(df[column])
    plt.title(f'{column}')
    #  plt.xlabel(column)
    plt.xticks([1], ['BAC'])  # Define o rótulo no eixo x

plt.tight_layout()
plt.show()
