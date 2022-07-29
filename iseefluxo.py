from functools import reduce
import numpy as np
import math


def obter_matriz_nula(linhas, colunas):
    return [[0 for i in range(colunas)] for j in range(linhas)]


# Dados do problema

## Parametros gerais

tolerancia = 0.00000001

## Dados das barras

### Quantidade
total_barras = 4
barras_PV = 1
barras_PQ = 2

### Tensões e fases iniciais supostas
tensao = [1, 1.05, 1, 1]
fase = [0] * (total_barras)

potencias_esperadas = [0.6, -0.5, -0.6, -0.2, -0.1]


## Dados de linhas e transformadores
impedancia_serie = obter_matriz_nula(total_barras, total_barras)
admitancia_shunt = obter_matriz_nula(total_barras, total_barras)

impedancia_serie[1][0] = 1 / (0.5 - 5j)
impedancia_serie[0][1] = 1 / (0.5 - 5j)

impedancia_serie[1][2] = 1 / (0.5 - 5j)
impedancia_serie[2][1] = 1 / (0.5 - 5j)

impedancia_serie[1][3] = 1 / (1 - 3j)
impedancia_serie[3][1] = 1 / (1 - 3j)

impedancia_serie[2][3] = 1 / (1 - 3j)
impedancia_serie[3][2] = 1 / (1 - 3j)

impedancia_serie[3][0] = 1 / (0.2 - 3j)
impedancia_serie[0][3] = 1 / (0.2 - 3j)

# Fim dos dados do problema


def soma(a, b):
    return a + b


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


# Calculo da matriz admitancia
matriz_admitancia = calcular_matriz_admitancia(impedancia_serie, admitancia_shunt)


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

    i_inicial = 1
    j_inicial = 1

    for i in range(1, total_barras):
        for j in range(1, total_barras):

            condutancia = matriz_admitancia[i][j].real
            susceptancia = matriz_admitancia[i][j].imag
            teta = fase[j] - fase[i]

            if i != j:

                submatriz_h[i - i_inicial][j - j_inicial] = (
                    -1
                    * tensao[i]
                    * tensao[j]
                    * (susceptancia * math.cos(teta) + condutancia * math.sin(teta))
                )
                continue

            if i == j:
                submatriz_h[i - i_inicial][j - j_inicial] = (
                    -1 * calcular_potencia_reativa(i) - (tensao[i] ** 2) * susceptancia
                )

    return submatriz_h


def calcular_submatriz_n():
    submatriz_n = obter_matriz_nula(total_barras - 1, barras_PQ)

    i_inicial = 1
    j_inicial = 1 + barras_PV

    for i in range(i_inicial, total_barras):
        for j in range(j_inicial, total_barras):

            condutancia = matriz_admitancia[i][j].real
            susceptancia = matriz_admitancia[i][j].imag
            teta = fase[j] - fase[i]

            if i != j:
                submatriz_n[i - i_inicial][j - j_inicial] = (
                    tensao[i]
                    * tensao[j]
                    * (condutancia * math.cos(teta) - susceptancia * math.sin(teta))
                )
                continue

            if i == j:
                submatriz_n[i - i_inicial][j - j_inicial] = (
                    calcular_potencia_ativa(i) + (tensao[i] ** 2) * condutancia
                )

    return submatriz_n


def calcular_submatriz_j():
    submatriz_j = obter_matriz_nula(barras_PQ, total_barras - 1)

    i_inicial = 1 + barras_PV
    j_inicial = 1

    for i in range(i_inicial, total_barras):
        for j in range(1, total_barras):

            condutancia = matriz_admitancia[i][j].real
            susceptancia = matriz_admitancia[i][j].imag
            teta = fase[j] - fase[i]

            if i != j:
                submatriz_j[i - i_inicial][j - j_inicial] = (
                    tensao[i]
                    * tensao[j]
                    * (susceptancia * math.sin(teta) - condutancia * math.cos(teta))
                )
                continue

            if i == j:
                submatriz_j[i - i_inicial][j - j_inicial] = (
                    calcular_potencia_ativa(i) - (tensao[i] ** 2) * condutancia
                )

    return submatriz_j


def calcular_submatriz_l():
    submatriz_l = obter_matriz_nula(barras_PQ, barras_PQ)

    i_inicial = 1 + barras_PV
    j_inicial = 1 + barras_PV

    for i in range(i_inicial, total_barras):
        for j in range(j_inicial, total_barras):

            condutancia = matriz_admitancia[i][j].real
            susceptancia = matriz_admitancia[i][j].imag
            teta = fase[j] - fase[i]

            if i != j:
                submatriz_l[i - i_inicial][j - j_inicial] = (
                    -1
                    * tensao[i]
                    * tensao[j]
                    * (susceptancia * math.cos(teta) + condutancia * math.sin(teta))
                )
                continue

            if i == j:
                submatriz_l[i - i_inicial][j - j_inicial] = (
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
    for i in range(barras_PV + 1, total_barras):
        potencias_calculadas[posicao_atual] = calcular_potencia_reativa(i)
        posicao_atual += 1

    return potencias_calculadas


def subtrair_vetores(vetor1, vetor2):
    resultado = [0] * len(vetor1)

    for i in range(len(resultado)):
        resultado[i] = vetor1[i] - vetor2[i]

    return resultado


def resultado_esta_bom(delta_pot):

    for i in range(len(delta_pot)):
        if abs(delta_pot[i]) > tolerancia:
            return False

    return True


def atualizar_fase_e_tensao(variacao_fase_e_tensao):

    # Atualização da fase
    # Pula a barra slack
    for i in range(1, total_barras):
        fase[i] += variacao_fase_e_tensao[i - 1]

    # Atualização tensão
    # Começa na primeira tensão de barra PQ
    # Pula a slack e as barras PV
    posicao_tensao = barras_PV + 1

    for i in range(total_barras, len(variacao_fase_e_tensao) + 1):

        tensao[posicao_tensao] += variacao_fase_e_tensao[i - 1]
        posicao_tensao += 1


quantidade_iteracoes = 0

# Laço de repetição responsável por calcular os valores de tensão nodal
# As fases também são calculadas
while True:
    potencias_calculadas = obter_potencias_calculadas()
    delta_pot = subtrair_vetores(potencias_esperadas, potencias_calculadas)

    if resultado_esta_bom(delta_pot):
        break

    quantidade_iteracoes += 1

    jacobiano = np.array(calcular_jacobiano())

    jacobiano_inverso = np.linalg.inv(jacobiano)

    variacao_fase_e_tensao = np.dot(jacobiano_inverso, delta_pot)

    atualizar_fase_e_tensao(variacao_fase_e_tensao)

print(tensao)
print(fase)
