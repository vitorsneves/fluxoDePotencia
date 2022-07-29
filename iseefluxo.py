from functools import reduce

def obter_matriz_nula(linhas, colunas):
  return [[0 for i in range(linhas)] for j in range(colunas)]

# Dados do problema

numero_barras = 9

## Dados de linhas e transformadores
impedancia_serie = obter_matriz_nula(numero_barras + 1, numero_barras + 1)
admitancia_shunt = obter_matriz_nula(numero_barras + 1, numero_barras + 1)

### Entre 1 e 4
impedancia_serie[1][4] = 0 + 0.0576j
impedancia_serie[4][1] = 0 + 0.0576j

### Entre 2 e 7
impedancia_serie[2][7] = 0 + 0.0625j
impedancia_serie[7][2] = 0 + 0.0625j

### Entre 3 e 9
impedancia_serie[3][9] = 0 + 0.0586j
impedancia_serie[9][3] = 0 + 0.0586j

### Entre 4 e 5
impedancia_serie[4][5] = 0.0100 + 0.0850j
impedancia_serie[5][4] = 0.0100 + 0.0850j
admitancia_shunt[4][5] = 0.0176j
admitancia_shunt[5][4] = 0.0176j

### Entre 4 e 6
impedancia_serie[4][6] = 0.0170 + 0.0920j
impedancia_serie[6][4] = 0.0170 + 0.0920j
admitancia_shunt[4][6] = 0.1580j
admitancia_shunt[6][4] = 0.1580j

### Entre 5 e 7
impedancia_serie[5][7] = 0.0320 + 0.1610j
impedancia_serie[7][5] = 0.0320 + 0.1610j
admitancia_shunt[5][7] = 0.3060j
admitancia_shunt[7][5] = 0.3060j

### Entre 6 e 9
impedancia_serie[6][9] = 0.0390 + 0.1700j
impedancia_serie[9][6] = 0.0390 + 0.1700j
admitancia_shunt[6][9] = 0.3580j
admitancia_shunt[9][6] = 0.3580j

### Entre 7 e 8
impedancia_serie[7][8] = 0.0085 + 0.0720j
impedancia_serie[8][7] = 0.0085 + 0.0720j
admitancia_shunt[7][8] = 0.1490j
admitancia_shunt[8][7] = 0.1490j

### Entre 8 e 9
impedancia_serie[8][9] = 0.0119 + 0.1010j
impedancia_serie[9][8] = 0.0119 + 0.1010j
admitancia_shunt[8][9] = 0.2090j
admitancia_shunt[9][8] = 0.2090j

def calcular_matriz_admitancia(impedancia_serie, admitancia_shunt):
  matriz_admitancia = obter_matriz_nula(numero_barras + 1, numero_barras + 1)

  def inverte(a):
    return 0 if a == 0 else 1 / a

  for i in range(numero_barras + 1):
    for j in range(numero_barras + 1):

      if(i != j):
        matriz_admitancia[i][j] = -1 * inverte(impedancia_serie[i][j])
        continue

      if(i == j):     
        def soma(a, b):
          return a + b

        admitancias_serie_barra = complex(reduce(soma, map(inverte, impedancia_serie[i])))
        admitancias_shunt_barra = complex(reduce(soma, admitancia_shunt[i]))

        matriz_admitancia[i][j] = admitancias_serie_barra + admitancias_shunt_barra
  
  return matriz_admitancia

for linha in calcular_matriz_admitancia(impedancia_serie, admitancia_shunt):
  print(linha)