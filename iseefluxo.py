from functools import reduce
import numpy as np
import math

# Funções auxiliares


def obter_matriz_nula(linhas, colunas):
    return [[0 for i in range(colunas)] for j in range(linhas)]


def soma(a, b):
    return a + b


# Dados do problema

total_barras = 9
barras_PV = 2
barras_PQ = 6
tolerancia = 0.001

## Dados de linhas e transformadores
impedancia_serie = obter_matriz_nula(total_barras, total_barras)
admitancia_shunt = obter_matriz_nula(total_barras, total_barras)

### Entre 0 e 3
impedancia_serie[0][3] = 0 + 0.0576j
impedancia_serie[3][0] = 0 + 0.0576j

### Entre 1 e 6
impedancia_serie[1][6] = 0 + 0.0625j
impedancia_serie[6][1] = 0 + 0.0625j

### Entre 2 e 8
impedancia_serie[2][8] = 0 + 0.0586j
impedancia_serie[8][2] = 0 + 0.0586j

### Entre 3 e 4
impedancia_serie[3][4] = 0.0100 + 0.0850j
impedancia_serie[4][3] = 0.0100 + 0.0850j
admitancia_shunt[3][4] = 0.0176j
admitancia_shunt[4][3] = 0.0176j

### Entre 3 e 5
impedancia_serie[3][5] = 0.0170 + 0.0920j
impedancia_serie[5][3] = 0.0170 + 0.0920j
admitancia_shunt[3][5] = 0.1580j
admitancia_shunt[5][3] = 0.1580j

### Entre 4 e 6
impedancia_serie[4][6] = 0.0320 + 0.1610j
impedancia_serie[6][4] = 0.0320 + 0.1610j
admitancia_shunt[4][6] = 0.3060j
admitancia_shunt[6][4] = 0.3060j

### Entre 5 e 8
impedancia_serie[5][8] = 0.0390 + 0.1700j
impedancia_serie[8][5] = 0.0390 + 0.1700j
admitancia_shunt[5][8] = 0.3580j
admitancia_shunt[8][5] = 0.3580j

### Entre 6 e 7
impedancia_serie[6][7] = 0.0085 + 0.0720j
impedancia_serie[7][6] = 0.0085 + 0.0720j
admitancia_shunt[6][7] = 0.1490j
admitancia_shunt[7][6] = 0.1490j

### Entre 7 e 8
impedancia_serie[7][8] = 0.0119 + 0.1010j
impedancia_serie[8][7] = 0.0119 + 0.1010j
admitancia_shunt[7][8] = 0.2090j
admitancia_shunt[8][7] = 0.2090j


def calcular_matriz_admitancia(impedancia_serie, admitancia_shunt):
    matriz_admitancia = obter_matriz_nula(total_barras, total_barras)

    def inverte(a):
        return 0 if a == 0 else 1 / a

    for i in range(total_barras):
        for j in range(total_barras):

            if i != j:
                matriz_admitancia[i][j] = -1 * inverte(impedancia_serie[i][j])
                continue

            if i == j:
                admitancias_serie_barra = complex(
                    reduce(soma, map(inverte, impedancia_serie[i]))
                )
                admitancias_shunt_barra = complex(reduce(soma, admitancia_shunt[i]))

                matriz_admitancia[i][j] = (
                    admitancias_serie_barra + admitancias_shunt_barra
                )

    return matriz_admitancia


matriz_admitancia = calcular_matriz_admitancia(impedancia_serie, admitancia_shunt)
tensao = [1.04, 1.025, 1.025, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]
fase = [0] * (total_barras)


def calcular_potencia_ativa(num_barra):
    potencia_ativa = 0

    for j in range(total_barras):

        condutancia = matriz_admitancia[num_barra][j].real
        susceptancia = matriz_admitancia[num_barra][j].imag
        teta = fase[j] - fase[num_barra]

        potencia_ativa += (
            tensao[num_barra]
            * tensao[j]
            * (condutancia * math.cos(teta) - susceptancia * math.sin(teta))
        )

    return potencia_ativa


def calcular_potencia_reativa(num_barra):
    potencia_reativa = 0

    for j in range(total_barras):

        condutancia = matriz_admitancia[num_barra][j].real
        susceptancia = matriz_admitancia[num_barra][j].imag
        teta = fase[j] - fase[num_barra]

        potencia_reativa -= (
            tensao[num_barra]
            * tensao[j]
            * (susceptancia * math.cos(teta) + condutancia * math.sin(teta))
        )

    return potencia_reativa


def calcular_submatriz_h():
    submatriz_h = obter_matriz_nula(total_barras - 1, total_barras - 1)

    for i in range(len(submatriz_h)):
        for j in range(len(submatriz_h[0])):

            condutancia = matriz_admitancia[i][j].real
            susceptancia = matriz_admitancia[i][j].imag
            teta = fase[j] - fase[i]

            if i != j:

                submatriz_h[i][j] = (
                    -1
                    * tensao[i]
                    * tensao[j]
                    * (susceptancia * math.cos(teta) + condutancia * math.sin(teta))
                )
                continue

            if i == j:
                submatriz_h[i][j] = (
                    -1 * calcular_potencia_reativa(i) - (tensao[i] ** 2) * susceptancia
                )

    return submatriz_h


def calcular_submatriz_n():
    submatriz_n = obter_matriz_nula(total_barras - 1, barras_PQ)

    for i in range(len(submatriz_n)):
        for j in range(len(submatriz_n[0])):

            condutancia = matriz_admitancia[i][j].real
            susceptancia = matriz_admitancia[i][j].imag
            teta = fase[j] - fase[i]

            if i != j:
                submatriz_n[i][j] = (
                    tensao[i]
                    * tensao[j]
                    * (condutancia * math.cos(teta) - susceptancia * math.sin(teta))
                )
                continue

            if i == j:
                submatriz_n[i][j] = (
                    calcular_potencia_ativa(i) + (tensao[i] ** 2) * condutancia
                )

    return submatriz_n


def calcular_submatriz_j():
    submatriz_j = obter_matriz_nula(barras_PQ, total_barras - 1)

    for i in range(len(submatriz_j)):
        for j in range(len(submatriz_j[0])):

            condutancia = matriz_admitancia[i][j].real
            susceptancia = matriz_admitancia[i][j].imag
            teta = fase[j] - fase[i]

        if i != j:
            submatriz_j[i][j] = (
                tensao[i]
                * tensao[j]
                * (susceptancia * math.sin(teta) - condutancia * math.cos(teta))
            )
            continue

        if i == j:
            submatriz_j[i][j] = (
                calcular_potencia_ativa(i) - (tensao[i] ** 2) * condutancia
            )

    return submatriz_j


def calcular_submatriz_l():
    submatriz_l = obter_matriz_nula(barras_PQ, barras_PQ)

    for i in range(len(submatriz_l)):
        for j in range(len(submatriz_l[0])):

            condutancia = matriz_admitancia[i][j].real
            susceptancia = matriz_admitancia[i][j].imag
            teta = fase[j] - fase[i]

            if i != j:
                submatriz_l[i][j] = (
                    -1
                    * tensao[i]
                    * tensao[j]
                    * (susceptancia * math.cos(teta) + condutancia * math.sin(teta))
                )
                continue

            if i == j:
                submatriz_l[i][j] = (
                    calcular_potencia_reativa(i) - (tensao[i] ** 2) * susceptancia
                )

    return submatriz_l


def calcular_jacobiano():
    jacobiano = obter_matriz_nula(barras_PQ * 2 + barras_PV, barras_PQ * 2 + barras_PV)

    submatriz_h = calcular_submatriz_h()
    submatriz_n = calcular_submatriz_n()
    submatriz_j = calcular_submatriz_j()
    submatriz_l = calcular_submatriz_l()

    for i in range(len(jacobiano)):
        for j in range(len(jacobiano[0])):

            if i < len(submatriz_h) and j < len(submatriz_h[0]):
                jacobiano[i][j] = submatriz_h[i][j]
                continue

            if i < len(submatriz_h) and j >= len(submatriz_h[0]):
                jacobiano[i][j] = submatriz_n[i][j - len(submatriz_h[0])]
                continue

            if i >= len(submatriz_h) and j < len(submatriz_h[0]):
                jacobiano[i][j] = submatriz_j[i - len(submatriz_h)][j]
                continue

            if i >= len(submatriz_h) and j >= len(submatriz_h[0]):
                jacobiano[i][j] = submatriz_l[i - len(submatriz_h)][
                    j - len(submatriz_h[0])
                ]

    return jacobiano


def obter_potencias_calculadas():
    potencias_calculadas = [0] * (barras_PQ * 2 + barras_PV)
    posicao_atual = 0

    # Cálculo de potências ativas
    # A barra slack foi desconsiderada
    for i in range(1, total_barras):
        potencias_calculadas[posicao_atual] = calcular_potencia_ativa(i)
        posicao_atual += 1

    # Cálculo de potências reativas
    # A barra slack e as barras pv foram desconsideradas
    for i in range(3, total_barras):
        potencias_calculadas[posicao_atual] = calcular_potencia_reativa(i)
        posicao_atual += 1

    return potencias_calculadas


def subtrair_vetores(vetor1, vetor2):
    resultado = [0] * len(vetor1)

    for i in range(len(resultado)):
        resultado[i] = vetor1[i] - vetor2[i]

    return resultado


def resultado_esta_bom(delta_pot):
    resultado_bom = True

    for i in range(len(delta_pot)):
        if delta_pot[i] > tolerancia:
            resultado_bom = False
            break

    return resultado_bom


potencias_esperadas = [
    1.63,
    0.85,
    0.00,
    -1.25,
    -0.90,
    0.00,
    -1.00,
    0.00,
    0.00,
    -0.50,
    -0.30,
    0.00,
    -0.35,
    0.00,
]

while True:
    potencias_calculadas = obter_potencias_calculadas()
    delta_pot = subtrair_vetores(potencias_esperadas, potencias_calculadas)

    if resultado_esta_bom(delta_pot):
        break

    jacobiano = np.array(calcular_jacobiano())

    jacobiano_inverso = np.linalg.inv(jacobiano)

    delta_x = jacobiano_inverso * delta_pot

    print(delta_x)
    break
