import numpy as np
import matplotlib.pyplot as plt

# Arrays extraidos do trabalho usando Gemini IA
# Array para a coluna "Tempo (ms)"
tempo_ms = [
    0.00, 0.50, 1.00, 1.50, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50, 
    5.00, 5.50, 6.00, 6.50, 7.00, 7.50, 8.00, 8.50, 9.00, 9.50, 
    10.00
]

# Array para a coluna "x_1(t)"
x1_t = [
    2.0000, 1.8377, 1.5080, 1.0408, 0.3228, -0.8418, -1.8678, -2.0077, -1.7759, -1.4131, 
    -0.8997, -0.0883, 1.1594, 1.9594, 1.9382, 1.6603, 1.2561, 0.6639, -0.3052, -1.5577, 
    -2.0068
]

# Array para a coluna "x_2(t)"
x2_t = [
    0.0000, -0.5346, -0.7804, -1.1247, -1.8337, -2.6765, -1.0843, 0.2562, 0.6106, 0.8494, 
    1.2514, 2.0850, 2.5563, 0.5355, -0.3999, -0.6822, -0.9552, -1.4788, -2.4492, -1.9799, 
    -0.0507
]

def AproPoli(x,y,grau):
    """
    retorna o vetor x com os coeficientes do polinomio para aproximar
    uma funcao amostrada pelo conjunto de dados x,y
    """
    matriz_A = [] #array para armazenar as linhas da matriz A
    matriz_Y = []   #array para armazenar as linhas do vetor coluna Y
    for linha in range(len(x)):
        matriz_A.append([]) # cria uma nova linha na matriz
        matriz_Y.append([y[linha]]) # armazena o valor de y[n] na linha do vetor coluna Y
        for i in range(grau+1):
            matriz_A[linha].append(x[linha]**i) # adiciona na linha da matriz a os valores de x^0 ate x^(grau)
    A = np.array(matriz_A) #transformamos as matrizes para o formato np.array
    Y = np.array(matriz_Y)
    # o formato np.array permite realizar operacoes com matrizes
    # a linha abaixo faz a conta dos minimos quadrados x = (At.A)'.At.y
    return np.dot(np.dot(np.linalg.inv(np.dot(np.transpose(A),A)),np.transpose(A)),Y)


def aprox(coeficientes,x):
    """
    python func que vai aplicar o polinomio no tempo x
    essa funcao automatiza a aplicacao no ponto para n coeficientes
    """
    sum = 0
    for i in range(len(coeficientes)):
        sum+=coeficientes[i][0]*x**i
    return sum

def EQM(coef_poli, x_real, tempo):
    sum = 0
    for i in range(len(x_real)):
        sum+= (x_real[i]-(aprox(coef_poli,tempo[i])))**2
    return sum/len(x_real)

for i in [3,5,15]:
    grau = i
    coeficientesx1 = AproPoli(tempo_ms,x1_t,grau)
    print("coef grau "+str(i)+" p/ x1_t [c(0),..., c(n)]' =\n"+str(coeficientesx1))
    eqm1 = EQM(coeficientesx1, x1_t, tempo_ms)
    coeficientesx2 = AproPoli(tempo_ms,x2_t,grau)
    print("coef grau "+str(i)+" p/ x2_t [c(0),..., c(n)]' =\n"+str(coeficientesx2))
    eqm2 = EQM(coeficientesx2, x2_t, tempo_ms)

    vectempo = []
    vecx1 = []
    vecx2 = []
    i = 0
    while i <= 10:
        vectempo.append(i)
        vecx1.append(aprox(coeficientesx1,i))
        vecx2.append(aprox(coeficientesx2,i))
        i+=0.01

    plt.subplot(1,2,1)
    plt.title("x1(t) -- EQM = "+str(eqm1))
    plt.xlabel('tempo [segundos]')
    plt.ylabel('tensao [volts]')
    plt.plot(tempo_ms,x1_t,"*",label = 'dados coletados')
    plt.plot(vectempo,vecx1,label = 'Polinomio ajustado de grau '+str(grau))
    plt.legend()

    plt.subplot(1,2,2)
    plt.title("x2(t) -- EQM = "+str(eqm2))
    plt.xlabel('tempo [segundos]')
    plt.ylabel('variacao de tensao [dV/dt]')
    plt.plot(tempo_ms,x2_t,"*",label = 'dados coletados')
    plt.plot(vectempo,vecx2,label = 'Polinomio ajustado de grau '+str(grau))
    plt.legend()
    plt.show()