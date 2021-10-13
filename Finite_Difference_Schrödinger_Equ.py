#%%

import numpy as np
import matplotlib.pyplot as plt

def Vpot(x):      ##Potansiyel Fonksiyonu Tanımlayalım
    return x**2

a = float(input('Alt Sınır: '))
b = float(input('Üst sınır: '))
N = int(input('Grid noktalarının sayısını belirtin: '))


x = np.linspace(a,b,N) ## a ile b değerleri arasında N tane değer oluştur. (Array olarak)
h = x[1]-x[0] #arrayın 1. değerinden 0. değeri çıkart



T = np.zeros((N-2)**2).reshape(N-2,N-2) #Kinetik enerji matrisi yaratılıyor. N-2'ye N-2'lik bir matris
for i in range(N-2):
    for j in range(N-2):
        if i==j:
            T[i,j]= -2  ## Matrisin Diagonal elemanlarını -2 yap. 
        elif np.abs(i-j)==1:
            T[i,j]=1
        else:
            T[i,j]=0




V = np.zeros((N-2)**2).reshape(N-2,N-2)      # Potansiyel enerji matrisi. N-2'ye N-2'lik bir matris.
for i in range(N-2):
    for j in range(N-2):
        if i==j:
            V[i,j]= Vpot(x[i+1])    #Diagonal elemanları buraya göre yaz.
        else:
            V[i,j]=0

H = -T/(2*h**2) + V       # İki matrisin toplamından HAMİLTONİAN matrisini oluştur.


val,vec=np.linalg.eig(H)         # Hamiltonian'ın özdeğerlerini ve özvektörlerini hesaplatan komut
z = np.argsort(val)       #özdeğerleri sıralama komutu 

z = z[0:5] # ilk 5 değeri kaydet
energies=(val[z]/val[z][0]) #ilk 5 terimin her birini ilk terime böl. Enerji özdeğerlerini verecek
print(energies)


plt.figure(figsize=(12,10))
for i in range(len(z)):
    y = []
    y = np.append(y,vec[:,z[i]])
    y = np.append(y,0)
    y = np.insert(y,0,0)
    plt.plot(x,y,lw=3, label="{} ".format(i))
    plt.xlabel('x', size=14)
    plt.ylabel('$\psi$(x)',size=14)
plt.legend()
plt.title('Finite Difference Metod ile bulunan Normalize edilmiş Harmonik Osilatörün Dalga Fonksiyonu',size=14)
plt.show()

# %%
#print(V)
#print(T)
#print(H)
# %%
# %%
