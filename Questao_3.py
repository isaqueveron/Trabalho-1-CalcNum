import matplotlib.pyplot as plt
import numpy as np

#para essa implementacao foi utilizado matrizes em formato np.ndarray

##------PARAMETROS
L       = 32000 #[m] nao usado
D       = 0.6 #[m] nao usado
As      = 60319 #[m^2]
Ad      = 0.2827 #[m^2]
dmdt    = 44 #[kg/s]
Tsource = 10 # [graus]
Tamb    = 4 # [graus]
rho     = 950 # [kg/m^3]
v       = 0.1638 # [m/s]
Cp      = 2100 # [J/kg.K]
U       = 0.510 #[W/m^2.K] nao usado
nbomba  = 0.75 # [-]
Celet   = 0.17 # [USD/kW.h]
Ccalor  = 0.05 # [USD/kW.h]
euller  = np.exp(1)
################

#PARAMETROS DO CRITERIO DE PARADA
epislon = 1e-5
maxiter = 100
#################################

def Funcao(vetorx:np.ndarray) -> np.ndarray:
    """
    Args: vetorx: um vetor coluna 4x1 (matriz com 1 coluna de indice 0 e 4 linhas, indices 0 a 3)
    
    Returns: Vetor Funcao aplicada nos pontos dados pelo vetorx. Um vetor coluna 4x1 (matriz com 1 coluna de indice 0 e 4 linhas, indices 0 a 3)
    """

    deltaP  = vetorx[0][0] #como vetorx eh um vetor coluna, O valor associado a variavel esta na linha n coluna 0
    Tin     = vetorx[1][0]
    Tout    = vetorx[2][0]
    mimed   = vetorx[3][0]

    f1      = deltaP*1e6-(32*mimed*L*v)/D**2
    f2      = dmdt*Cp*(Tin-Tout)-U*As*((Tin-Tamb)/2+(Tout-Tamb)/2)
    f3      = 6100-(((deltaP*1e3*dmdt)/(rho*nbomba)))*Celet*24-((dmdt*Cp*(Tin-Tsource))/1000)*Ccalor*24
    f4      = mimed-1100*np.exp(-0.1*(Tin+Tout)/2)

    return np.array([[f1],
                     [f2],
                     [f3],
                     [f4]]) #retorna um vetor coluna

def Jacobiano(vetorx:np.ndarray) -> np.ndarray:
    """
    Args: vetorx: um vetor coluna 4x1 (matriz com 1 coluna de indice 0 e 4 linhas, indices 0 a 3)
    
    Returns: uma matriz 4x4 representando o jacobiano da funcao aplicado nos pontos do vetorx
    """
    #-----LINHA 1
    df1ddeltaP = 1e6
    df1dTin = 0
    df1dTout = 0
    df1dmimed = -(32*L*v)/D**2
    ##########

    #-----LINHA 2
    df2ddeltaP = 0
    df2dTin = dmdt*Cp - (U*As/2)
    df2dTout = -dmdt*Cp - (U*As/2)
    df2dmimed = 0
    ##########

    #-----LINHA 3 #VER
    df3ddeltaP = (((-1e3*dmdt)/(rho*nbomba)))*Celet*24
    df3dTin = ((-dmdt*Cp)/1000)*Ccalor*24
    df3dTout = 0
    df3dmimed = 0
    ##########

    #-----LINHA 4 VER
    df4ddeltaP = 0
    df4dTin = 55*np.exp(-0.05*(vetorx[1][0]+vetorx[2][0]))
    df4dTout = 55*np.exp(-0.05*(vetorx[1][0]+vetorx[2][0]))
    df4dmimed = 1
    ##########

    return np.array([
                    [df1ddeltaP,df1dTin,df1dTout,df1dmimed], 
                    [df2ddeltaP,df2dTin,df2dTout,df2dmimed],
                    [df3ddeltaP,df3dTin,df3dTout,df3dmimed], 
                    [df4ddeltaP,df4dTin,df4dTout,df4dmimed]
                    ])

def Newton(vetorx:np.ndarray) -> np.ndarray:
    """
    Args: vetorx: um vetor coluna 4x1 representando o chute inicial;
    
    Returns: um vetor coluna 4x1.
    """

    iter            = 0

    ###-----VETORES VAZIOS PARA PLOTS
    iter_vector     = []
    deltap_vector   = []
    tin_vector      = []
    tout_vector     = []
    mimed_vector    = []
    f1_vec          = []
    f2_vec          = []
    f3_vec          = []
    f4_vec          = []
    #################################

    while True:  
        iter+=1
        """if iter == 1:
            print("(Jacobiano(vetorx))\n",Jacobiano(vetorx)) #debug
            print("Funcao(vetorx)\n",Funcao(vetorx)) #debug"""

        ########################-----METODO NEWTON-RAPHSON
        vetorx = vetorx - np.dot(np.linalg.inv(Jacobiano(vetorx)),Funcao(vetorx))
        #print("vetorxT\n",vetorx) #debug
        ##############################

        ####################-----VETORES PARA PLOTS
        iter_vector.append(iter)

        deltap_vector.append(vetorx[0][0])
        tin_vector.append(vetorx[1][0])
        tout_vector.append(vetorx[2][0])
        mimed_vector.append(vetorx[3][0])
        
        fxn = Funcao((vetorx))
        #print("Funcao(vetorx):\n",fxn)
        f1_vec.append(fxn[0][0])
        f2_vec.append(fxn[1][0])
        f3_vec.append(fxn[2][0])
        f4_vec.append(fxn[3][0])
        #########################################

        ###############-----CRITERIO DE PARADA
        absfnx = (fxn[0]**2+fxn[1]**2+fxn[2]**2+fxn[3]**2)**(1/2) #calcula o modulo do vetor Funcao(vetorx)
        if (absfnx<epislon)or(iter>=maxiter): break #criterio de parada
        ###############################

    ###################-----PLOTS 
    plt.subplot(1,2,1) 
    plt.plot(iter_vector,deltap_vector,label = '∆p [Pa]')
    plt.plot(iter_vector,mimed_vector,label = 'µmed[Pa.s]')
    plt.grid()
    plt.legend()
    plt.title("Queda de pressao no oleoduto &\nViscosidade media do oleo extra pesado no trecho",loc = 'left')
    plt.xlabel('iteracoes [n]',loc='right')
    plt.subplot(1,2,2) 
    plt.plot(iter_vector,tin_vector,label = 'Entrada [°C]')
    plt.plot(iter_vector,tout_vector,label = 'Saida [°C]')
    plt.grid()
    plt.legend()
    plt.title("Temperatura do oleoduto [graus]",loc = 'left')
    plt.xlabel('iteracoes [n]',loc='right')
    plt.ylabel('Temperatura [graus]')
    plt.show()
    #############################

    return vetorx,iter

chute_inicial = np.array([[1],
                          [12],
                          [10],
                          [0.05]])

result = Newton(chute_inicial)

print("Metodo de Newton: raiz = \n",result[0], "Com "+ str(result[1])+" iteracoes")

print("Funcao aplicada em Vetor encontrado:\n", Funcao(result[0]))

chute_inicial_2 = np.array([[15],[25],[10],[1]])

result2 = Newton(chute_inicial_2)

print("Metodo de Newton: raiz = \n",result2[0], "Com "+ str(result2[1])+" iteracoes")

print("Funcao aplicada em Vetor encontrado:\n", Funcao(result[0]))
