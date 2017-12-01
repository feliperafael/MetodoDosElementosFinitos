import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from numpy.polynomial.legendre import leggauss

inicio = 0
fim = 1

def du(x):
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
#	t =  np.linspace(-1, 1, num=num_elementos)
#	
#	F = (t*((x2-x1)/2.0) + (x2+x1)/2.0)* (x2+x1/2.0)
#	return F
	
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
		if i == 0:
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
			return -0.5
			return 0.5*(1-t)
		if i == 1 :
			return 0.5
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
	x = np.linspace(inicio, fim, num=nel+1)
	h = x[1]-x[0]
	det = h/2.0
	ii = 0
	for elemento in range(nel): #loop em elementos
		for l in range(nint): #loop em ponto de gaus
			
			fe = np.zeros(grau+1)
			ke = np.zeros(((grau+1),(grau+1)))
			if l == 0:				
				ponto = -1
			if l == 1:
				ponto = 1
			for i in range(grau+1):
				ponto = 0.5*(x[ii]+x[ii+grau] + p[l]*h)
				print("\nponto"+str(ponto))	
				print("phi(ponto) = "+str(phi_function(i,ponto)))			
				fe[i] += du(ponto) * phi_function(i,p[l]) * w[l]
				for j in range(grau+1):
					ke[i,j] += w[l]* dphi_function(i,p[l]) * dphi_function(i,p[l])
		print fe*(h/2.0)
		print ke/h
		ii = ii + grau
		#atualiza matriz global   
		for i in range(grau+1):
			for j in range(grau+1):
				if i!= j:
					K[elemento+i][elemento+j] -= ke[i][j]
				else:				
					K[elemento+i][elemento+j] += ke[i][j]
		#atualiza vetor local
		for i in range(grau+1):
			F[i+elemento] += fe[i]
	F *=h/2.0
	K *= 1.0/h
	#aplico condicoes de contorno
	K[0][1] = 0
	K[1][0] = 0
	
	#K[dimensao_matriz_global-2][dimensao_matriz_global-1] = 0
	#K[dimensao_matriz_global-1][dimensao_matriz_global-2] = 0
	
	F[0] = 0
	F[1] += 0/h
	F[dimensao_matriz_global-1] = 2*np.pi
	#F[dimensao_matriz_global-2] += 2*np.pi/h

	print K
	print F
	solucao = LA.solve(K,F)
   
	#solucao = np.insert(solucao, num_elementos-1, 2*np.pi)
	#solucao[4] = 2*np.pi
	print("Solucao :",solucao)

	print F
	print K
if __name__ == "__main__":
	num_elementos = 4 #na real eh numero de pontos 
	#k = 1 #ordem_polinomio = 1
	grau_polinomio = 1	
	nint = 2 # numero de pontos de integracao de gauss
	x = np.linspace(inicio, fim, num=num_elementos)
	i = 0
	solucao_aproximada(nint,grau_polinomio,num_elementos)	
	#calc_matriz_local(k,k+1,num_elementos)

#[[ 0.88729833  0.5	0.11270167]
#[ 0.11270167  0.5	0.88729833]]
	#p,w = leggauss(3)
	#for i in range(3):
	#	print phi_function(0,p[i])
	#	print phi_function(1,p[i])

