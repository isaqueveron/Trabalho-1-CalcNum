import matplotlib.pyplot as plt

def Liq_Profit_negative(qinj):
    return -(-300*(qinj**3) + 8170.4*(qinj**2) - 53112*qinj + 111800)

def first_dFdt_Liq_Profit_negative(qinj):
    return -(-3*300*(qinj**2) + 2*8170.4*(qinj) - 53112)

def second_dFdt_Liq_Profit_negative(qinj):
    return -(-2*3*300*(qinj) + 2*8170.4)

#plotando a funcao originao (nao negativa)
vecqinj = []
vecliqprofit = []
vecliqprofit1d = []
vecliqprofit2d = []
for qinj in range(20):
    vecqinj.append(qinj)
    vecliqprofit.append(-Liq_Profit_negative(qinj))
    vecliqprofit1d.append(-first_dFdt_Liq_Profit_negative(qinj))
    vecliqprofit2d.append(-second_dFdt_Liq_Profit_negative(qinj))
plt.figure('L')
plt.plot(vecqinj,vecliqprofit,label = "L(qinj)")
plt.plot(vecqinj,vecliqprofit1d,label = "L'(qinj)")
plt.plot(vecqinj,vecliqprofit2d,label = "L''(qinj)")
plt.title("Func L(qinj) e suas derivadas")
plt.xlabel('qinj [milhoes de pes cubicos por dia]')
plt.ylabel('Lucro [USD/dia]')
plt.legend()
plt.grid()
plt.show()

def Newton(x0 = float, epislon = float, max_iter = float, passo = 1):
    iter = 0

    ########VETORES VAZIOS PARA PLOT
    xk_vector   = []
    yk_vector   = []
    iter_vector = []
    ################################

    while not(iter>=max_iter):
        iter+=1
        iter_vector.append(iter)
        xk_vector.append(x0)
        yk_vector.append(first_dFdt_Liq_Profit_negative(x0))
        
        x0 = x0 - first_dFdt_Liq_Profit_negative(x0)/second_dFdt_Liq_Profit_negative(x0)
        if (abs(first_dFdt_Liq_Profit_negative(x0))<=epislon): break

    #####PLOTS
    plt.figure(chute)
    plt.subplot(1,2,1)
    plt.plot(vecqinj,vecliqprofit,xk_vector,[-Liq_Profit_negative(i) for i in xk_vector], "*")
    plt.grid()
    plt.title("Avanco do metodo em L(qinj)\nchute inicial = "+str(chute))
    plt.xlabel('qinj [milhoes de pes cubicos por dia]')
    plt.ylabel('lucro [USD/dia]')
    
    plt.subplot(1,2,2)
    plt.plot(iter_vector,yk_vector)
    plt.grid()
    plt.title("Iteracoes x erro \n(Newton sem Backtracking)")
    plt.xlabel('iteracoes [n]')
    plt.ylabel('erro')
    plt.show()
    ###########


    return x0,iter

def NewtonB_not_conv(x0 = float, epislon = float, max_iter = float, passo = 1):
    iter = 0
    chute_inicial = x0

    ########VETORES VAZIOS PARA PLOT
    xk_vector   = []
    yk_vector   = []
    iter_vector = []
    ################################

    while not(iter>=max_iter):
        iter+=1
        iter_vector.append(iter)
        xk_vector.append(x0)
        yk_vector.append(first_dFdt_Liq_Profit_negative(x0))
        
        if second_dFdt_Liq_Profit_negative(x0)>0:
            v = -first_dFdt_Liq_Profit_negative(x0)/second_dFdt_Liq_Profit_negative(x0)
        else:
            v = -first_dFdt_Liq_Profit_negative(x0)
        
        while Liq_Profit_negative(x0+v)>first_dFdt_Liq_Profit_negative(x0)+passo*first_dFdt_Liq_Profit_negative(x0)*v:
            v = v/2
        
        x0 = x0 + v
        
        if (abs(first_dFdt_Liq_Profit_negative(x0))<=epislon): break

    #####PLOTS
    plt.subplot(1,2,1)
    plt.plot(vecqinj,vecliqprofit,xk_vector,[-Liq_Profit_negative(i) for i in xk_vector], "*")
    plt.grid()
    plt.title("Avanco do metodo\nchute inicial = "+str(chute))
    plt.xlabel('qinj [milhoes de pes cubicos por dia]')
    plt.ylabel('lucro [USD/dia]')
    
    plt.subplot(1,2,2)
    plt.plot(iter_vector,yk_vector)
    plt.grid()
    plt.title("Iteracoes x erro \n(Newton com Backtracking\n NonConvex)")
    plt.xlabel('iteracoes [n]')
    plt.ylabel('erro')
    plt.show()
    ###########


    return x0,iter

for chute in [6,11]:
    N = Newton(chute,10e-5,100)
    NB = NewtonB_not_conv(chute,10e-5,100)
    print("Sem backtracking e chute inicial "+str(chute)+" = "+str(N[0])+" com "+str(N[1])+" iteracoes")
    print("Com backtracking e chute inicial "+str(chute)+" = "+str(NB[0])+" com "+str(NB[1])+" iteracoes")
