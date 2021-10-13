#%%
import matplotlib
import matplotlib.pyplot as plt 
import numpy
import numpy.polynomial.hermite as Herm
import math

#Basitleştirilmiş olarak kullanıldı.
m=1.  #kütle
w=1.  #Açısal frekan
hbar=1. #Planck sabiti

#Discretized space
dx = 0.05
x_lim = 12
x = numpy.arange(-x_lim,x_lim,dx)

def hermite(x, n):     #Harmonik Salınıcının normlanmış kararlı durumlarında Hermite Polinomları bulunur. Burada bu polinomları ifade eden kodlar var.
    xi = numpy.sqrt(m*w/hbar)*x
    herm_coeffs = numpy.zeros(n+1)
    herm_coeffs[n] = 1
    return Herm.hermval(xi, herm_coeffs)


def stationary_state(x,n): #Harmonik salınıcının normlanmış kararlı durumu.
    xi = numpy.sqrt(m*w/hbar)*x
    prefactor = 1./math.sqrt(2.**n * math.factorial(n)) * (m*w/(numpy.pi*hbar))**(0.25)
    psi = prefactor * numpy.exp(- xi**2 / 2) * hermite(x,n)
    return psi
  


plt.figure(figsize=(10, 8))
plt.subplot(3,2,1)
plt.plot(x, numpy.conjugate(stationary_state(x,0))*stationary_state(x,0), label="n=0")
plt.title(r"Harmonik Salınıcının ilk 5 Kararlı Durumu")

plt.legend()
plt.subplot(3,2,2)
plt.plot(x, numpy.conjugate(stationary_state(x,1))*stationary_state(x,1), label="n=1")


plt.legend()
plt.subplot(3,2,3)
plt.plot(x, numpy.conjugate(stationary_state(x,2))*stationary_state(x,2), label="n=2")


plt.legend()
plt.subplot(3,2,4)
plt.plot(x, numpy.conjugate(stationary_state(x,3))*stationary_state(x,3), label="n=3")


plt.legend()
plt.subplot(3,2,5)
plt.plot(x, numpy.conjugate(stationary_state(x,4))*stationary_state(x,4), label="n=4")


plt.legend()
plt.subplot(3,2,6)
plt.plot(x, numpy.conjugate(stationary_state(x,5))*stationary_state(x,5), label="n=5")


plt.legend()
plt.show()
# %%

# %%
