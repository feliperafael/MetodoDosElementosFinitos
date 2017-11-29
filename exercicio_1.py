import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

inicio = 0
fim = 1

def du_dx(x):
    return np.cos(2*np.pi*x)
def exata(x):
    return np.sin(2*np.pi*x)

def plot(x):
    fig = plt.figure()
    ax = plt.subplot(111)
    plt.xlabel('x')
    plt.ylabel('u')
    ax.plot(x,exata(x),'r')
    plt.show()

def solucao_aproximada(x):
    print("a")

# t = ponto de gauss
def peso_gauss(n,t):
    if n == 2:
        return 1.0
    elif n == 3:
        if t == 0:
            return 8.0/9.0
        else:
            return 5.0/9.0
    print "ERROR: invalid value of n"
    return -1

def pontos_gauss(n):
    pontos = np.zeros(n)
    if n == 2:
        pontos[0] = -(np.sqrt(3.0)/3.0)
        pontos[1] = -pontos[0]
    if n == 3:
        pontos[0] = 0
        pontos[1] = -np.sqrt(3.0/5.0)
        pontos[2] = -pontos[1]    
    return pontos

#def converte_espaco(x1,x2):
#    t =  np.linspace(-1, 1, num=num_elementos)
#    
#    F = (t*((x2-x1)/2.0) + (x2+x1)/2.0)* (x2+x1/2.0)
#    return F
    
def quadratura_de_gauss(n):
    print("quadratura de gauss")

# onde grau eh o grau do polinomio
# i e qual funcao phi
# t e o valor para avaliar
def phi_function(grau,i,t):
    if grau == 1: #linear
        if i == 1 :
            return 0.5*(1-t)
        if i == 2 :
            return 0.5*(1+t)
        else:
            print("PHI FUNCTION RETURN ZERO")
            return 0
    if grau == 2: #quadratico
        if i == 1 :
            return 0.5*t*(t-1)
        if i == 2 :
            return -1*(t-1)*(t+1)
        if i == 3:
            return 0.5*t*(t+1)
        else:
            print("PHI FUNCTION RETURN ZERO")
            return 0
    print("ERROR PHI FUNCTION")
       
# k = ordem_polinomio
# n = numero de nos da integracao gauss
# num_elementos = num_elementos da malha de integracao
def calc_matriz_local(k,n,num_elementos):
    dimensao_matriz_global = k*num_elementos + 1
    ke = np.zeros(((k+1),(k+1)))
    K = np.zeros(((dimensao_matriz_global),(dimensao_matriz_global)))
    F = np.zeros(dimensao_matriz_global)
    pontos_de_gauss = pontos_gauss(n)
    h = 2.0/(num_elementos-1) # (1 - (-1))/num_elementos
    t =  np.linspace(-1, 1, num=num_elementos)
    numero_pontos_de_gauss = k+1

    for a in range(dimensao_matriz_global - k): #dimensao - grau
        #constroi matriz local
        ke = np.zeros(((k+1),(k+1)))
        fe = np.zeros(k+1)
        for i in range(k+1):
            t_pontos_de_gauss = pontos_gauss(numero_pontos_de_gauss)
            if i == 0:
                fe[i] += 1#0.5*(1-t_pontos_de_gauss[i]) # incremento vetor local
            if i == 1:
                fe[i] += 1#0.5*(1+t_pontos_de_gauss[i])
            #fe[i] += du_dx(fim)
            for j in range(k+1):
                for p in range(numero_pontos_de_gauss):
                    #print t_pontos_de_gauss
                    #print phi_function(k,phi_de_gauss+1,t_pontos_de_gauss)
                    if k == 1: #linear
                        if i == 0 :
                            if j == 0:
                                ke[i][j] += 0.5*(1-t_pontos_de_gauss[p])
                            if j == 1:
                                ke[i][j] -= 0.5*(1+t_pontos_de_gauss[p])
                        if i == 1 :
                            if j == 0:
                                ke[i][j] += -0.5
                            if j == 1:
                                ke[i][j] += 0.5
                        else:
                            #print("PHI FUNCTION RETURN ZERO")
                            ke[i][j] += 0
                    if k == 2: #quadratico
                        if i == 0 :
                            ke[i][j] += 0.5*t_pontos_de_gauss[p]*(t_pontos_de_gauss[p]-1)
                        if i == 1 :
                            ke[i][j] += -1*(t_pontos_de_gauss[p]-1)*(t_pontos_de_gauss[p]+1)
                        if i == 2:
                            ke[i][j] += 0.5*t_pontos_de_gauss[p]*(t_pontos_de_gauss[p]+1)
                    #print ke
                    print("vetor local")
                    print fe
        #constroi vetor local
        print h
       # for i in range(k+1):
       #     fe[i] += 1
       #     print fe
        
        
        for i in range(k+1):
            for j in range(k+1):
                    K[a+i][a+j] += ke[i][j]

        for i in range(k+1):
            F[i+a] += fe[i]
    #print ke
    print K
    print F
    K = K*(1.0/h)
    F = F*(h/2.0)

    #aplico condicoes de contorno
    K[0][1] = 0
    K[1][0] = 0
    
    K[dimensao_matriz_global-2][dimensao_matriz_global-1] = 0
    K[dimensao_matriz_global-1][dimensao_matriz_global-2] = 0
    
    F[0] = 0
    F[1] += 0/h
    F[dimensao_matriz_global-1] = 2*np.pi
    F[dimensao_matriz_global-2] += 2*np.pi/h

    print K
    print F
    print LA.solve(K,F)
    
if __name__ == "__main__":
    num_elementos = 4 #na real eh numero de pontos 
    k = 1 #ordem_polinomio = 1

    x = np.linspace(inicio, fim, num=num_elementos)
    i = 0
    calc_matriz_local(k,k+1,num_elementos)
