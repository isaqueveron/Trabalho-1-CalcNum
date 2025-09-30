#CODIGO DA ELIOMINACAO DE GAUSS COM PIVOTEAMENTO IMPLEMENTADOS DENTRO DA FUNCAO DE ELIMINACAO

import numpy as np
import time

def EliGaussPivotParc(Matriznaosingular,vetorcolunab): #n^2
    print("Matriz Inicial:\n",Matriznaosingular)
    n = len(Matriznaosingular) #numero de linhas
    j = 0 

    while j < (n-1):
        ####PIVOTEAMENTO
        valor_max = 0.0 #variavel para armazenar o maior valor
        p = j
        for k in range(j,n): #para cada linha, na mesma coluna
            valor_abs_atual = abs(Matriznaosingular[k][j]) #ve o valor alocado nessa posicao
            if valor_abs_atual>valor_max: #se for maior (em modulo) que o atual adiciona na variavel de valor max
                valor_max = valor_abs_atual
                p = k
        if p != j: #se o valor de maior modulo nao for o valor do pivo, troca as linhas
            Matriznaosingular[[k, p]] = Matriznaosingular[[p, k]] #troca as linhas de lugar (operacao do numpy)
        ####FIMDOPIVOTEAMENTO

        i = j + 1
        while i < n:
            ####ELIMINACAODEGAUSS
            mij = Matriznaosingular[i][j]/Matriznaosingular[j][j] #(valor da linha i e coluna do pivo)/pivo
            Matriznaosingular[i] = Matriznaosingular[i] - mij*Matriznaosingular[j] #Linha [i] <- Linha [i] - mij*linha [i-1]
            vetorcolunab[i] = vetorcolunab[i] - mij*vetorcolunab[j] #mesma coisa para o vetor b
            i+=1
        j+=1
    return  Matriznaosingular, vetorcolunab
    
def ResolverMatriz_U(Matriz_U_b_vetorcolunab:tuple):
    """
    Args: 
        Matriz_U_b: Matriz quadratida triangular superior + ultima coluna para b, em uma unica matriz
    #Formato das Matrizes: Matriz[linhas][colunas]
    """

    Matriz_U_b = Matriz_U_b_vetorcolunab[0]
    vetorcolunab = Matriz_U_b_vetorcolunab[1]
    linhas = len(Matriz_U_b)
    vetor_x = np.zeros(linhas, dtype=float)
    for n in range(linhas - 1, -1, -1): #faz a operacao de baixo para cima
        #print('n:',n)
        b = vetorcolunab[n][0]
        for m in range(n + 1, linhas):
            #print('m',m)
            b = b - Matriz_U_b[n][m] * vetor_x[m] #a partir da segunda iteracao, ja temos valores para usar em vetor_x
        pivo = Matriz_U_b[n][n]
        vetor_x[n] = b/pivo
    return vetor_x

####Questao2

matriz7_7 = np.array([
                        [-1, 0, 0, 0, 0, 0, 2],
                        [1, -1, 0, 0, 0, 0, 0],
                        [0, 1, -1, 0, 0, 0, 0],
                        [0, 0, 1, -1, 0, 0, 0],
                        [0, 0, 0, 1, -1, 0, 0],
                        [0, 0, 0, 0, 0, 1, 0],
                        [0, 0, 0, 0, 1, 1, -1],
                    ])

matrizB = np.array([[0],
                    [-130],
                    [102],
                    [22],
                    [-1],
                    [29],
                    [38]])

tempo = time.time()
print(ResolverMatriz_U(EliGaussPivotParc(matriz7_7,matrizB)))
print("tempo para resolver: ",time.time()-tempo)
