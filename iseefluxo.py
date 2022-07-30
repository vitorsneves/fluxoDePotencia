from functools import reduce
import numpy as np
import math
import cmath


def obter_matriz_nula(linhas, colunas):
    return [[0 for i in range(colunas)] for j in range(linhas)]


# Dados do problema

## Parametros gerais

tolerancia = 0.001
potencia_complexa_base_MVA = 100
tensoes_base_kV = [16.5, 18, 13.8, *([230] * 6)]

## Dados das barras

### Quantidade
total_barras = 9
barras_PV = 2
barras_PQ = 6

### Tensões e fases iniciais supostas
tensao = [1.04, 1.025, 1.025, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]
fase = [0] * (total_barras)

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


# Fim dos dados do problema


def soma(a, b):
    return a + b


def inverte(a):
    return 0 if a == 0 else 1 / a


def calcular_matriz_admitancia(impedancia_serie, admitancia_shunt):
    matriz_admitancia = obter_matriz_nula(total_barras, total_barras)

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


def norma_vetor(vetor):
    norma = 0

    for i in range(len(vetor)):
        norma += vetor[i] ** 2

    return math.sqrt(norma)


def resultado_esta_bom(delta_pot):
    erro = norma_vetor(delta_pot)

    esta_bom = erro < tolerancia

    return esta_bom, erro


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

        tensao[posicao_tensao] = (-1 * tensao[posicao_tensao]) / (
            variacao_fase_e_tensao[i - 1] - 1
        )
        posicao_tensao += 1


quantidade_iteracoes = 0
historico_erros_maximos = []

# Laço de repetição responsável por calcular os valores de tensão nodal
# As fases também são calculadas
while True:
    potencias_calculadas = obter_potencias_calculadas()
    delta_pot = subtrair_vetores(potencias_esperadas, potencias_calculadas)

    esta_bom, erro_maximo = resultado_esta_bom(delta_pot)

    historico_erros_maximos.append(erro_maximo)
    if esta_bom:
        break

    quantidade_iteracoes += 1

    jacobiano = np.array(calcular_jacobiano())
    jacobiano_inverso = np.linalg.inv(jacobiano)

    variacao_fase_e_tensao = np.dot(jacobiano_inverso, delta_pot)
    atualizar_fase_e_tensao(variacao_fase_e_tensao)


def calcular_fluxo_potencia(barra_1, barra_2):
    tensao_1 = cmath.rect(tensao[barra_1], fase[barra_1])
    tensao_2 = cmath.rect(tensao[barra_2], fase[barra_2])

    corrente_12 = (tensao_1 - tensao_2) * inverte(impedancia_serie[barra_1][barra_2])
    corrente_12 += tensao_1 * admitancia_shunt[barra_1][barra_2]

    fluxo_potencia_12 = tensao_1 * np.conjugate(corrente_12)

    return fluxo_potencia_12


def calcular_todos_fluxos_potencia_MVA():
    todos_fluxos_potencia = obter_matriz_nula(total_barras, total_barras)

    for i in range(total_barras):
        for j in range(total_barras):
            todos_fluxos_potencia[i][j] = (
                calcular_fluxo_potencia(i, j) * potencia_complexa_base_MVA
            )

    return todos_fluxos_potencia


def calcular_perdas_MVA(todos_fluxos_potencia):
    perdas_MVA = 0 + 0j

    for i in range(total_barras):
        for j in range(total_barras):
            perdas_MVA += todos_fluxos_potencia[i][j]

    return perdas_MVA


todos_fluxos_potencia = calcular_todos_fluxos_potencia_MVA()
perdas_MVA = calcular_perdas_MVA(todos_fluxos_potencia)


# Daqui em diante está as funções responsáveis por
# exibir os dados na tela


def exibe_quantidade_iteracoes():
    print(f"\ntolerância = {tolerancia}")

    for i in range(len(historico_erros_maximos)):
        print(f"Iteração {i} => erro máximo = {'%.4f' % historico_erros_maximos[i]}")

    print(f"\nQuantidade total de iterações: {quantidade_iteracoes}\n")


def padding_numero(numero):
    if numero >= 0:
        return " %.4f" % numero
    else:
        return "%.4f" % numero


def exibe_tensoes_nodais():
    def tensao_pu(numero_barra):
        return padding_numero(tensao[numero_barra])

    def angulo_rad(numero_barra):
        return padding_numero(fase[numero_barra])

    def tensao_kV(numero_barra):
        return padding_numero(tensao[numero_barra] * tensoes_base_kV[numero_barra])

    def angulo_graus(numero_barra):
        angulo = fase[numero_barra] * (180 / math.pi)
        return padding_numero(angulo)

    print("\n#######################  Tensões Nodais  #######################")
    print("################################################################")
    print("    Barra     V[pu]      Fase[rad]     V[kV]     Fase[grau]    ")
    for i in range(len(tensao)):
        print(
            f"      {i + 1}     {tensao_pu(i)}      {angulo_rad(i)}    { tensao_kV(i)}      {angulo_graus(i)}"
        )
    print("################################################################\n")


def exibe_fluxo_potencia():
    def imprime_linha_tabela(barra_1, barra_2):
        pot_complexa = todos_fluxos_potencia[barra_1][barra_2]

        pot_ativa = padding_numero(pot_complexa.real)
        pot_reativa = padding_numero(pot_complexa.imag)

        print(
            f"      {barra_1 + 1}  ---->  {barra_2 + 1}          {pot_ativa}           {pot_reativa}"
        )

    print("\n#####################  Fluxos de Potência  #####################")
    print("################################################################")
    print("    From       To            P[MW]             Q[MVAr]       ")
    for i in range(total_barras):
        for j in range(i + 1, total_barras):
            if todos_fluxos_potencia[i][j] != 0:
                imprime_linha_tabela(i, j)
                imprime_linha_tabela(j, i)
    print("################################################################\n")


def exibe_perdas_totais():
    print(f"\nPerdas de potência ativa totais [MW]: {padding_numero(perdas_MVA.real)}")
    print(
        f"\nPerdas de potência reativa totais [MVAr]: {padding_numero(perdas_MVA.imag)}\n"
    )


print("----------------------------------------------------------------")
exibe_quantidade_iteracoes()
print("----------------------------------------------------------------")
exibe_tensoes_nodais()
print("----------------------------------------------------------------")
exibe_fluxo_potencia()
print("----------------------------------------------------------------")
exibe_perdas_totais()
print("----------------------------------------------------------------")
