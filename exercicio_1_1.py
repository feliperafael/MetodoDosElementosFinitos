import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from numpy.polynomial.legendre import leggauss

inicio = 0
fim = 1

def du2_dx(x):
    return 4*np.pi*np.pi*np.sin(2*np.pi*x)
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
def peso_gauss(grau,t):
    if grau == 1:
        return 1.0
    elif grau == 3:
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
# t ponto de gauss
def phi_function(i,t):
	grau = 1
    if grau == 1: #linear
        if i == 0 :
            return 0.5*(1-t)
        if i == 1 :
            return 0.5*(1+t)
        else:
            print("PHI FUNCTION RETURN ZERO")
            return 0
    if grau == 2: #quadratico
        if i == 0 :
            return 0.5*t*(t-1)
        if i == 1 :
            return -1*(t-1)*(t+1)
        if i == 2:
            return 0.5*t*(t+1)
        else:
            print("PHI FUNCTION RETURN ZERO")
            return 0
    print("ERROR PHI FUNCTION")


def phi_function(i,t):
	grau = 1
    if grau == 1: #linear
        if i == 0 :
            return 0.5*(1-t)
        if i == 1 :
            return 0.5*(1+t)
        else:
            print("PHI FUNCTION RETURN ZERO")
            return 0
    if grau == 2: #quadratico
        if i == 0 :
            return 0.5*t*(t-1)
        if i == 1 :
            return -1*(t-1)*(t+1)
        if i == 2:
            return 0.5*t*(t+1)
        else:
            print("PHI FUNCTION RETURN ZERO")
            return 0
    print("ERROR PHI FUNCTION")

def dphi_function(i,t):
	grau = 1	
	if grau == 1: #linear
        if i == 0 :
            return 0.5*(1-t)
        if i == 1 :
            return 0.5*(1+t)
        else:
            print("PHI FUNCTION RETURN ZERO")
            return 0
    if grau == 2: #quadratico
        if i == 0 :
            return 0.5*t*(t-1)
        if i == 1 :
            return -1*(t-1)*(t+1)
        if i == 2:
            return 0.5*t*(t+1)
        else:
            print("PHI FUNCTION RETURN ZERO")
            return 0
    print("ERROR PHI FUNCTION")

def solucao_aproximada(nint,grau,nel):
	qtd_de_pontos_de_gauss = 2
	p,w = leggauss(qtd_de_pontos_de_gauss)
	dimensao_matriz_global = grau*nel + 1
	K = np.zeros(((dimensao_matriz_global),(dimensao_matriz_global)))
    F = np.zeros(dimensao_matriz_global)

	for elemento in range(nel): #loop em elementos
		ke = np.zeros(((k+1),(k+1)))
        fe = np.zeros(k+1)
		for l in range(nint): #loop em ponto de gaus
			xx = 0
			for i in range(grau+1):
				fe[i] += du2_dx(0.5*(1))
				for j in range(grau+1):
					ke[i][j] += w[l]* dfi[i,l] * dfi[j,l]
    
if __name__ == "__main__":
    num_elementos = 4 #na real eh numero de pontos 
    #k = 1 #ordem_polinomio = 1
	grau_polinomio = 1	
	nint = 2 # numero de pontos de integracao de gauss
    x = np.linspace(inicio, fim, num=num_elementos)
    i = 0
	solucao_aproximada(nint,grau_polinomio,num_elementos)	
	#calc_matriz_local(k,k+1,num_elementos)
