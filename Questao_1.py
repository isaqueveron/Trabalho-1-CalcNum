from math import *
import matplotlib.pyplot as plt

def Funcao(angulo):
    """
    Args: algulo [radianos], velocidade inicial [m/s], gravidade [m/s^2]
        
    Returns: alcance horizontal [metros]
    """
    return ((100**2)/9.8)*(sin(angulo)+0.5*sin(2*angulo))-1000

def DerivadaReal(angulo):
    """
    Args: algulo [radianos], velocidade inicial [m/s], gravidade [m/s^2]
        
    Returns: derivada no ponto [metros]/[s]
    """
    return ((100**2/9.8)*(cos(angulo)+cos(2*angulo)))

#PLOT DO MODELO VARIANDO ANGULO DE 0 graus a 90 graus
vetor_angulo = []
vetor_alcance = []
i=0
while i <= pi/2 :
    vetor_angulo.append(i)
    vetor_alcance.append(Funcao(i)+1000)
    i+=0.01
plt.plot(vetor_angulo,vetor_alcance,vetor_angulo,[1000]*len(vetor_angulo))
plt.title("Alcance x Angulo")
plt.xlabel('angulo [radianos]')
plt.ylabel('Alcance [metros]')
#plt.show()
##################################################

#PARAMETROS DO CRITERIO DE PARADA
epislon = 1e-5
max_iter = 100
##################################

def Bissecao(x1 = float, x2 = float):
    iter = 0

    ########VETORES VAZIOS PARA PLOT
    xk_vector   = []
    yk_vector   = []
    iter_vector = []
    ################################

    while True:  
        iter+=1
        iter_vector.append(iter)

        ########################METODO
        xn = 0.5*(x1+x2) #calcula o valor medio do intervalo (x1 e x2) e armazena em xn
        xk_vector.append(xn)
        fxn = Funcao(xn) #calcula a funcao no ponto xn
        yk_vector.append(fxn)
        if (abs(fxn)<epislon )or(iter>=100): break #testa se o valor da funcao eh suficientemente perto de zero
        else:                                      #ou se o numero maximo de iteracoes foi atingido.
            if fxn*Funcao(x1)>0: x1 = xn #descobre o sinal da funcao no ponto para 
            else: x2 = xn                #substituir o limite do intervalo respectivo.
        ###############################

    #####PLOTS
    plt.subplot(1,2,1)
    plt.plot(vetor_angulo,vetor_alcance,xk_vector,[i+1000 for i in yk_vector], "*")
    plt.grid()
    plt.title("Alcance x Angulo (azul)\n Alcance em cada Iteracao (laranja '*')\n (BISSECAO)")
    plt.xlabel('angulo [radianos]')
    plt.ylabel('alcance [metros]')
    plt.subplot(1,2,2)
    plt.plot(iter_vector,yk_vector)
    plt.grid()
    plt.title("Iteracoes x erro \n(BISSECAO)")
    plt.xlabel('iteracoes [n]')
    plt.ylabel('erro [metros]')
    plt.show()
    ###########

    return xn, iter

result = Bissecao(0.5,1.5)
print("Metodo da Bissecao: raiz = ",result[0], "Com "+ str(result[1])+" iteracoes") 

############################################################

def FalsaPos(x1 = float, x2 = float):
    iter = 0
    
    ########VETORES VAZIOS PARA PLOT
    xk_vector   = []
    yk_vector   = []
    iter_vector = []
    ################################

    while True:  
        iter+=1
        iter_vector.append(iter)

        ########################METODO      #calcula uma media ponderada dos valores do intervalo (x1 e x2)
        xn = (x1*Funcao(x2) - x2*Funcao(x1))/ (Funcao(x2) - Funcao(x1)) 
        xk_vector.append(xn)
        fxn = Funcao(xn)
        yk_vector.append(fxn)
        if (abs(fxn)<epislon )or(iter>=100): break #mesma logica da anterior
        else: 
            if fxn*Funcao(x1)>0: x1 = xn #mesma logica da anterior
            else: x2 = xn
        ###############################

    #####PLOTS        
    plt.subplot(1,2,1)
    plt.plot(vetor_angulo,vetor_alcance,xk_vector,[i+1000 for i in yk_vector], "*")
    plt.grid()
    plt.title("Alcance x Angulo (azul)\n Alcance em cada Iteracao (laranja '*')\n (FALSA POSICAO)")
    plt.xlabel('angulo [radianos]')
    plt.ylabel('alcance [metros]')
    plt.subplot(1,2,2)
    plt.plot(iter_vector,yk_vector)
    plt.grid()
    plt.title("Iteracoes x erro \n(FALSAPOS)")
    plt.xlabel('iteracoes [n]')
    plt.ylabel('erro [metros]')
    plt.show()
    ##########################3

    return xn,iter

result = FalsaPos(0.5,1.5)
print("Metodo da Falsa Posicao: raiz = ", result[0], "Com "+ str(result[1])+" iteracoes")

############################################################

def Newton(xn = float):
    iter = 0
    xk_vector   = []
    yk_vector   = []
    iter_vector = []

    while True:  
        iter+=1
        iter_vector.append(iter)

        ########################METODO
        xn = xn - Funcao(xn)/DerivadaReal(xn) #calcula o proximo valor usando o metodo de newton 
        xk_vector.append(xn)                  #onde DerivadaReal() eh a derivada analitica de Funcao()
        fxn = Funcao(xn)
        yk_vector.append(fxn)
        if (abs(fxn)<epislon )or(iter>=100): break #mesma logica das anteriores
        ###############################

    #####PLOTS        
    plt.subplot(1,2,1)
    plt.plot(vetor_angulo,vetor_alcance,xk_vector,[i+1000 for i in yk_vector], "*")
    plt.grid()
    plt.title("Alcance x Angulo (azul)\n Alcance em cada Iteracao (laranja '*')\n (NEWTON)")
    plt.xlabel('angulo [radianos]')
    plt.ylabel('alcance [metros]')
    plt.subplot(1,2,2)
    plt.plot(iter_vector,yk_vector)
    plt.grid()
    plt.title("Iteracoes x erro \n(NEWTON)")
    plt.xlabel('iteracoes [n]')
    plt.ylabel('erro [metros]')
    plt.show()
    ##########

    return xn, iter

result = Newton(0.8)
print("Metodo de Newton: raiz = ",result[0], "Com "+ str(result[1])+" iteracoes")

############################################################


def Secante(x1 = float, x2 = float):
    iter = 0
    xk_vector   = []
    yk_vector   = []
    iter_vector = []

    while True:  
        iter+=1
        iter_vector.append(iter)

        ########################METODO
        Derivada_aprox = (Funcao(x2) - Funcao(x1))/(x2-x1) #calcula a variacao da funcao no intervalo dado
        xn = x2 - Funcao(x2)/Derivada_aprox #de maneira similar ao metodo de newton, calcula a proxima possivel raiz
        xk_vector.append(xn)
        fxn = Funcao(xn)
        yk_vector.append(fxn)
        if (abs(fxn)<epislon )or(iter>=100): break #mesma logica das anteriores
        else: 
            if fxn*Funcao(x1)>0: x1 = xn #mesma logica da bissecao e falsapos
            else: x2 = xn
        ###############################

    #####PLOTS        
    plt.subplot(1,2,1)
    plt.plot(vetor_angulo,vetor_alcance,xk_vector,[i+1000 for i in yk_vector], "*")
    plt.grid()
    plt.title("Alcance x Angulo (azul)\n Alcance em cada Iteracao (laranja '*')\n (SECANTE)")
    plt.xlabel('angulo [radianos]')
    plt.ylabel('alcance [metros]')
    plt.subplot(1,2,2)
    plt.plot(iter_vector,yk_vector)
    plt.grid()
    plt.title("Iteracoes x erro \n(SECANTE)")
    plt.xlabel('iteracoes [n]')
    plt.ylabel('erro [metros]')
    plt.show()
    ##########

    return xn, iter

result = Secante(0.7,0.8)
print("Metodo da Secante: raiz = ", result[0], "Com "+ str(result[1])+" iteracoes")


