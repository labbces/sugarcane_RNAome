import matplotlib.pyplot as plt

# Valores originais e transformados
original = [0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00]
transformado = [x - 0.7 for x in original]

plt.figure(figsize=(10, 6))
plt.plot(original, transformado, marker='o', linestyle='-', color='b')
plt.title('Transformação das Correlações de Pearson')
plt.xlabel('Correlações Originais')
plt.ylabel('Correlações Transformadas')
plt.grid(True)
plt.show()