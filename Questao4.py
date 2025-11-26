import math
import matplotlib.pyplot as plt

def funcao(s):
    return math.sin(s**2)

def R(func,a,b,n):
    h = (b-a)/n
    sum = 0
    x = a
    j = 1
    while j<=n:
        sum += func((2*x+h)/2)
        x += h
        j += 1
    return h*sum

def T(func,a,b,n):
    h = (b-a)/n
    sum = 0
    x = a
    j = 1
    while j<=n:
        sum += (func(x)+func(x+h))
        x += h
        j += 1
    return (h/2)*sum

def S(func,a,b,n):
    vec_Y = []
    vec_L = []
    h = (b-a)/n
    sum = 0
    x = a
    j = 1
    while j<=n+1:
        if j == 1 or j == n+1:
            coef = 1
        elif j%2 == 0:
            coef = 4
        else:
            coef = 2
        sum += (h/3)*coef*func(x)
        vec_Y.append(sum)
        vec_L.append(x)
        x += h
        j += 1
    return sum,vec_Y,vec_L


I = 0.7782378
for n in [10,100]:
    r = R(funcao,0,1.5,n)
    t = T(funcao,0,1.5,n)
    s,vec_Y,vec_L = S(funcao,0,1.5,n)
    print('regra dos retangulos com N = '+str(n)+' = '+str(r)+' com erro = '+str(I-r))
    print('regra dos trapezios com  N = '+str(n)+' = '+str(t)+' com erro = '+str(I-t))
    print('regra de simpson com     N = '+str(n)+' = '+str(s)+' com erro = '+str(I-s))

plt.plot(vec_L,vec_Y)
plt.title('Integral de Fresnel')
plt.xlabel('distancia percorrida L [km]')
plt.ylabel('Y(L) [km]')
plt.grid()
plt.show()


