import matplotlib.pyplot as plt

u = 1

def dx1dt(x1,x2):
    return x2

def dx2dt(x1,x2):
    return u*(1-x1**2)*x2-x1

def Euller(x_inicial,h,tempo):
    tms = 0
    x = x_inicial
    t_vec = []
    x1_vec = []
    x2_vec = []
    while tms<=tempo:
        t_vec.append(tms)
        x1_vec.append(x[0])
        x2_vec.append(x[1])
        x[0] = x[0] + h*dx1dt(x[0],x[1])
        x[1] = x[1] + h*dx2dt(x[0],x[1])
        tms += h
    print('')
    return t_vec,x1_vec,x2_vec

def RK4(x_inicial,h,tempo):
    tms = 0
    x = x_inicial
    t_vec = []
    x1_vec = []
    x2_vec = []
    while tms<=tempo:
        t_vec.append(tms)
        x1_vec.append(x[0])
        x2_vec.append(x[1])

        k1_1 = h*dx1dt(x[0],x[1])
        k1_2 = h*dx2dt(x[0],x[1])

        k2_1 = h*dx1dt(x[0]+k1_1/2,x[1]+k1_2/2)
        k2_2 = h*dx2dt(x[0]+k1_1/2,x[1]+k1_2/2)

        k3_1 = h*dx1dt(x[0]+k2_1/2,x[1]+k2_2/2)
        k3_2 = h*dx2dt(x[0]+k2_1/2,x[1]+k2_2/2)

        k4_1 = h*dx1dt(x[0]+k3_1,x[1]+k3_2)
        k4_2 = h*dx2dt(x[0]+k3_1,x[1]+k3_2)

        x[0] = x[0] + (k1_1 + 2*k2_1 + 2*k3_1 + k4_1)/6
        x[1] = x[1] + (k1_2 + 2*k2_2 + 2*k3_2 + k4_2)/6
        tms += h

    return t_vec,x1_vec,x2_vec

for h in [0.2 , 0.05]:

    plt.figure('h = '+str(h))

    t_vec,x1_vec,x2_vec = Euller([2,0],h,20)
    t_vec_rk4,x1_vec_rk4,x2_vec_rk4 = RK4([2,0],h,20)

    plt.subplot(2,1,1)
    plt.title('x1(t) = V(t) [Volts] vs tempo [ms], com h = '+str(h))
    plt.plot(t_vec,x1_vec,label='Euller')
    plt.plot(t_vec,x1_vec_rk4,label='Rk4')
    plt.legend()
    plt.grid()

    plt.subplot(2,1,2)
    plt.title('x2(t) = dV(t)/dt vs tempo [ms], com h = '+str(h))
    plt.plot(t_vec,x2_vec,label='Euller')
    plt.plot(t_vec,x2_vec_rk4,label='Rk4')
    plt.legend()
    plt.grid()

    plt.figure('retrato de fase com h = '+str(h))

    plt.title('retrato de fase com h = '+str(h))
    plt.plot(x1_vec,x2_vec,label = 'Euller')
    plt.plot(x1_vec_rk4,x2_vec_rk4,label = 'RK4')
    plt.xlabel('x1_vec')
    plt.ylabel('x2_vec')
    plt.legend()
    plt.grid()

plt.show()