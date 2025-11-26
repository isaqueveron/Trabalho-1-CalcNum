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

def Subst_regressiva(Lt,Y):

    Matriz_L_transposta = Lt
    vetorcolunaY = Y
    linhas = len(Matriz_L_transposta)
    vetor_x = np.zeros(linhas, dtype=float)
    for n in range(linhas - 1, -1, -1): #faz a operacao de baixo para cima
        #print('n:',n)
        b = vetorcolunaY[n][0]
        for m in range(n + 1, linhas):
            #print('m',m)
            b = b - Matriz_L_transposta[n][m] * vetor_x[m] #a partir da segunda iteracao, ja temos valores para usar em vetor_x
        pivo = Matriz_L_transposta[n][n]
        vetor_x[n] = b/pivo
    return vetor_x

def Subst_progressiva(L,Y):

    Matriz_L = L
    vetorcolunaY = Y
    linhas = len(Matriz_L)
    vetor_x = np.zeros(linhas, dtype=float)
    for n in range(linhas): #faz a operacao de cima para baixo
        "print('n:',n)"
        b = vetorcolunaY[n][0]
        for m in range(n):
            #print('m',m)
            b = b - Matriz_L[n][m] * vetor_x[m] #a partir da segunda iteracao, ja temos valores para usar em vetor_x
        pivo = Matriz_L[n][n]
        vetor_x[n] = b/pivo
    return vetor_x.reshape(-1,1)

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

    L = np.linalg.cholesky(np.dot(np.transpose(A),A))
    Atb = np.dot(np.transpose(A),Y)
    y = Subst_progressiva(L,Atb)
    x = Subst_regressiva(np.transpose(L),y)
    return x

def aprox(coeficientes,x):
    """
    python func que vai aplicar o polinomio no tempo x
    essa funcao automatiza a aplicacao no ponto para n coeficientes
    """
    sum = 0
    for i in range(len(coeficientes)):
        sum+=coeficientes[i]*x**i
    return sum

def EQM(coef_poli, x_real, tempo):
    sum = 0
    for i in range(len(x_real)):
        sum+= (x_real[i]-(aprox(coef_poli,tempo[i])))**2
    return sum/len(x_real)


grau = 3
coeficientesx1 = AproPoli(tempo_ms,x1_t,grau)
print("coef grau "+str(grau)+" p/ x1_t [c(0),..., c(n)]' =\n"+str(coeficientesx1))
print("EQM x1_t grau "+str(grau)+" = "+str(EQM(coeficientesx1, x1_t, tempo_ms)))
coeficientesx2 = AproPoli(tempo_ms,x2_t,grau)
print("coef grau "+str(grau)+" p/ x2_t [c(0),..., c(n)]' =\n"+str(coeficientesx2))
print("EQM x2_t grau "+str(grau)+" = "+str(EQM(coeficientesx2, x2_t, tempo_ms)))
vectempo = []
vecx1_3 = []
vecx2_3 = []
i = 0
while i <= 10:
    vectempo.append(i)
    vecx1_3.append(aprox(coeficientesx1,i))
    vecx2_3.append(aprox(coeficientesx2,i))
    i+=0.01

grau = 5
coeficientesx1 = AproPoli(tempo_ms,x1_t,grau)
print("coef grau "+str(grau)+" p/ x1_t [c(0),..., c(n)]' =\n"+str(coeficientesx1))
print("EQM x1_t grau "+str(grau)+" = "+str(EQM(coeficientesx1, x1_t, tempo_ms)))
coeficientesx2 = AproPoli(tempo_ms,x2_t,grau)
print("coef grau "+str(grau)+" p/ x2_t [c(0),..., c(n)]' =\n"+str(coeficientesx2))
print("EQM x2_t grau "+str(grau)+" = "+str(EQM(coeficientesx2, x2_t, tempo_ms)))
vectempo = []
vecx1_5 = []
vecx2_5 = []
i = 0
while i <= 10:
    vectempo.append(i)
    vecx1_5.append(aprox(coeficientesx1,i))
    vecx2_5.append(aprox(coeficientesx2,i))
    i+=0.01

grau = 13
coeficientesx1 = AproPoli(tempo_ms,x1_t,grau)
print("coef grau "+str(grau)+" p/ x1_t [c(0),..., c(n)]' =\n"+str(coeficientesx1))
print("EQM x1_t grau "+str(grau)+" = "+str(EQM(coeficientesx1, x1_t, tempo_ms)))
coeficientesx2 = AproPoli(tempo_ms,x2_t,grau)
print("coef grau "+str(grau)+" p/ x2_t [c(0),..., c(n)]' =\n"+str(coeficientesx2))
print("EQM x2_t grau "+str(grau)+" = "+str(EQM(coeficientesx2, x2_t, tempo_ms)))
vectempo = []
vecx1_13 = []
vecx2_13 = []
i = 0
while i <= 10:
    vectempo.append(i)
    vecx1_13.append(aprox(coeficientesx1,i))
    vecx2_13.append(aprox(coeficientesx2,i))
    i+=0.01

plt.subplot(1,2,1)
plt.title("x1(t)")
plt.xlabel('tempo [segundos]')
plt.ylabel('tensao [volts]')
plt.plot(tempo_ms,x1_t,"*",label = 'dados coletados')
plt.plot(vectempo,vecx1_3,label = 'Polinomio ajustado de grau '+str(3))
plt.plot(vectempo,vecx1_5,label = 'Polinomio ajustado de grau '+str(5))
plt.plot(vectempo,vecx1_13,label = 'Polinomio ajustado de grau '+str(13))
plt.legend()

plt.subplot(1,2,2)
plt.title("x2(t)")
plt.xlabel('tempo [segundos]')
plt.ylabel('variacao de tensao [dV/dt]')
plt.plot(tempo_ms,x2_t,"*",label = 'dados coletados')
plt.plot(vectempo,vecx2_3,label = 'Polinomio ajustado de grau '+str(3))
plt.plot(vectempo,vecx2_5,label = 'Polinomio ajustado de grau '+str(5))
plt.plot(vectempo,vecx2_13,label = 'Polinomio ajustado de grau '+str(13))
plt.legend()
plt.show()