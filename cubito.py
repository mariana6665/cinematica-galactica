###############################
#Como leer un cubo de datos
###############################


import astropy.io.fits as fits #modulo para leer archivos fits
import matplotlib.pyplot as plt #modulo para graficar
import numpy as np #este modulo es para trabajar con matrices como en matlab

#esta funcion les entrega un arreglo con los valores reales del eje correspondiente
#1 es velocidad, 2 es longirud galactica, 3 es latitud galactica
def values(h,j):
	N=h['NAXIS'+str(j)];
	val=np.zeros(N);
	for i in range(0,N):
		val[i] = (i+1-float(h['CRPIX'+str(j)]))*float(h['CDELT'+str(j)]) + float(h['CRVAL'+str(j)]);
	return val;

cubo	= fits.open("Cubo_de_datos.fits") #abrir objeto cubo de datos
data 	= cubo[0].data #extraer matriz de datos
header	= cubo[0].header #extraer el header del archivo fits
#Estos seran los tres arreglos con los valores reales de los tres ejes del cubo
velocidad=values(header,1)
longitud=values(header,2)
latitud=values(header,3)
j=1000

v = np.zeros(385)
velocidades = np.zeros(385)
velocidadesabs = np.zeros(385)
radios = np.zeros(385)
b = np.zeros(385)
c = np.zeros(33)
latitudes = np.zeros(385)
indice = 1000


for k in range(0,len(longitud)):
        print k
        for h in range(0,len(latitud)):
                j=1000
                for i,tem in enumerate(data[h][k][:]):
                        if i<305 and tem >= 0.5 and data[h][k][i+1]>=tem:
                                indice = i
                                break
                c[h]=i
        v[k]=min(c)
        for p,valor in enumerate (c):
                if valor==min(c):
                        break
        b[k]=p

for m in range (0,len(longitud)):
        print m
        velocidades[m]=velocidad[v[m]]

for q in range (0,len(longitud)):
        print q
        latitudes[q]=latitud[b[q]]
for vel in range (0,len(longitud)):
        velocidadesabs[vel]=abs(velocidades[vel])
for rad in range (0,len(longitud)):
        radios[rad]=-8.5*np.sin(3.14*longitud[rad]/180)
for rad in range (0,len(longitud)):
        radios[rad]=-8.5*np.sin(3.14*longitud[rad]/180)

vrot = velocidadesabs+radios*220/8.5
altura = abs(8.5*np.cos(longitud*3.14/180)*np.sin(latitudes*3.14/180))
distancia = abs(8.5*np.cos(longitud*3.14/180))
archivo = open('masa.dat','w')
import os
for i in range(len(vrot)):
	archivo.write(str(radios[i]) + ' ' + str(vrot[i]) +'\n')
archivo.close()

#Modelos de masa:

M1 = 2.95423*(10**21)
masa1 = np.sqrt(6.67384*(10**(-11))*M1/radios)/1000

p2 = 4.74717*(10**19)
masa2 = np.sqrt(6.67384*(10**(-11))*3.14*radios*p2)/1000

p3 = 6.83774*(10**17)
masa3 = np.sqrt(6.67384*(10**(-11))*4*3.14*(radios**2)*p3/3)/1000

M4 = 1.11313*(10**21)
p4 = 3.11928*(10**19)
masa4 = np.sqrt(6.67384*(10**(-11))*(M4+3.14*(radios**2)*p4)/radios)/1000

M5 = 1.61088*(10**21)
p5 = 3.23028*(10**18)
masa5 = np.sqrt(6.67384*(10**(-11))*(M5+4*3.14*(radios**3)*p5/3)/radios)/1000

#graficos
radioss=radios/8.5
xlargo=[0,8]
xancho=[220,220]
fig=plt.figure()
"""
plt.plot(radioss,vrot,'m.')
plt.xlabel("Radio/Ro")
plt.ylabel("Velocidad de rotacion [km/s]")
plt.title("Velocidad de rotacion en funcion del Radio")
plt.axis([0.1, 0.9, 180, 260])
plt.show()
fig.savefig('vel_rot.png')

plt.plot(radioss,velocidadesabs,'m.')
plt.xlabel("Radio/Ro")
plt.ylabel("V_LSR[km/s]")
plt.title("Velocidad LSR en funcion del Radio")
plt.show()
fig.savefig('vel_lsr.png')

plt.plot(distancia,altura,'m-')
plt.axis([4, 8.5, 0, 0.20])
plt.xlabel("Distancia del centro galactico [kpc]")
plt.ylabel("Z [kpc]")
plt.title("Corrugacion del plano")
plt.show()
fig.savefig('corrugacion.png')

plt.plot(radioss,vrot,'m.',radioss,masa1,'r-',radioss,masa2,'b-',radioss,masa3,'c-',radioss,masa4,'g-',radioss,masa5,'k-')
plt.title("Modelos de masa")
plt.ylabel("Velocidad de rotacion [km/s]")
plt.xlabel("Radio/Ro")

plt.axis([0.1, 0.9, 159, 300])
plt.show()



fig.savefig('masas.png')
"""
suma = 0
for i in range(len(vrot)):
	suma = suma + vrot[i]
	i = i + 1
print "velocidad de rotacion promedio:"
print suma/len(vrot)
