#CODIGO IMPLEMENTADO COM PIVOTEAMENTO OCORRENDO ANTES DA ELIMINACAO DE GAUSS

import numpy as np
import time

def Pivotemanto(MatrizA_b: np.ndarray) -> np.ndarray: #n^2
    """
    Funcao para fazer o pivotemanto antes de iniciar a eliminacao de gauss
    """

    n = len(MatrizA_b)

    for k in range(n-1):
        max_valor = 0.0
        p = k
        for i in range(k, n):
            valor_absoluto_atual = abs(MatrizA_b[i][k])
            if valor_absoluto_atual > max_valor:
                max_valor = valor_absoluto_atual
                p = i

        if p != k:
            MatrizA_b[[k, p]] = MatrizA_b[[p, k]] #troca as linhas de lugar (operacao do numpy)
            
    return MatrizA_b

def EliminacaoGauss(MatrizA_b: np.ndarray) -> np.ndarray: #n^3
    """
    Args:
        MatrizA_b: Matriz quadratida + ultima coluna para b, em uma unica matriz
    #Formato das Matrizes: Matriz[linhas][colunas]
    """

    MatrizA_b=Pivotemanto(MatrizA_b) #Faz o pivoteamento

    linhas = len(MatrizA_b)
    colunas = len(MatrizA_b[0])
    for n in range(linhas-1):
        for m in range(colunas-1):
            if (n==m and MatrizA_b[n][m-1]==0) or n==m==0:
                pivo = MatrizA_b[n][m]
                for i in range(linhas-1):
                    try: 
                        MatrizA_b[n+i+1] = MatrizA_b[n+i+1] - MatrizA_b[n]*(MatrizA_b[n+i+1][m]/pivo)
                        print("Matriz apos troca:\n", MatrizA_b) #debug
                    except: pass
    return MatrizA_b

def ResolverMatriz_U(Matriz_U_b: np.ndarray) -> np.ndarray: #n^2
    """
    Args:
        Matriz_U_b: Matriz quadratida triangular superior + ultima coluna para b, em uma unica matriz
    #Formato das Matrizes: Matriz[linhas][colunas]
    """

    linhas = len(Matriz_U_b)
    vetor_x = np.zeros(linhas, dtype=float)
    for n in range(linhas - 1, -1, -1): #faz a operacao de baixo para cima
        #print('n:',n) #debug
        b = Matriz_U_b[n][-1]
        for m in range(n + 1, linhas):
            #print('m',m) #debug
            b = b - Matriz_U_b[n][m] * vetor_x[m] #a partir da segunda iteracao, ja temos valores para usar em vetor_x
        pivo = Matriz_U_b[n][n]
        vetor_x[n] = b/pivo
    return vetor_x

####Resolvendo a questao Questao 2

matriz7_7_B = np.array([
                        [-1, 0, 0, 0, 0, 0, 2, 0],
                        [1, -1, 0, 0, 0, 0, 0, -130],
                        [0, 1, -1, 0, 0, 0, 0, 102],
                        [0, 0, 1, -1, 0, 0, 0, 22],
                        [0, 0, 0, 1, -1, 0, 0, -1],
                        [0, 0, 0, 0, 0, 1, 0, 29],
                        [0, 0, 0, 0, 1, 1, -1, 38],
                    ])

Matriz_resolvida_U_b =  EliminacaoGauss(matriz7_7_B)

tempo = time.time()
print("Matriz 7X7|B resolvida por algoritimo desenvolvido: \n",Matriz_resolvida_U_b)
print("Vetor x = ",ResolverMatriz_U(Matriz_resolvida_U_b))
print("tempo para resolver: ",time.time()-tempo)
