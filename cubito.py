import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np
import os


def values(h,j):
	N = h['NAXIS'+str(j)];
	val = np.zeros(N);
	for i in range(0,N):
		val[i] = (i + 1 - float(h['CRPIX' + str(j)]))*float(h['CDELT' + str(j)]) + float(h['CRVAL' + str(j)]);
	return val;

cubo = fits.open("Cubo_de_datos.fits")
data = cubo[0].data
header = cubo[0].header

#Arreglos con los valores de Velocidades, Longitudes y Latitudes
vel = values(header,1)
lon = values(header,2)
lat = values(header,3)
j = 1000

v = np.zeros(385)
radios = np.zeros(385)
velocidades = np.zeros(385)
Vabs = np.zeros(385)
b = np.zeros(385)
c = np.zeros(33)
latitudes = np.zeros(385)

indice = 1000


for k in range(0,len(lon)):
        print k
        for h in range(0,len(lat)):
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

for m in range (0,len(lon)):
        print m
        velocidades[m] = vel[v[m]]

for q in range (0,len(lon)):
        print q
        latitudes[q] = lat[b[q]]
for vel in range (0,len(lon)):
        Vabs[vel] = abs(velocidades[vel])
for rad in range (0,len(lon)):
        radios[rad] = -8.5*np.sin(np.pi*lon[rad]/180)
for rad in range (0,len(lon)):
        radios[rad] = -8.5*np.sin(np.pi*lon[rad]/180)

vrot = Vabs + radios*220/8.5
Z = abs(8.5*np.cos(lon*np.pi/180)*np.sin(latitudes*np.pi/180))
distancia = abs(8.5*np.cos(lon*np.pi/180))
archivo = open('masa.dat','w')

for i in range(len(vrot)):
	archivo.write(str(radios[i]) + ' ' + str(vrot[i]) +'\n')
archivo.close()

#Calcula la velocidad de rotacion promedio
suma = 0
for i in range(len(vrot)):
	suma = suma + vrot[i]
	i = i + 1
print "velocidad de rotacion promedio:"
print suma/len(vrot)


#Modelos de distribucion de masa
#modelo1: masa puntual
M1 = 2.95423*(10**21)
M1r = np.sqrt(6.67384*(10**(-11))*M1/radios)/1000

#modelo2: disco uniforme
sigma2 = 4.74717*(10**19)
M2r = np.sqrt(6.67384*(10**(-11))*np.pi*radios*sigma2)/1000

#modelo3: esfera uniforme
rho3 = 6.83774*(10**17)
M3r = np.sqrt(6.67384*(10**(-11))*4*np.pi*(radios**2)*rho3/3)/1000

#modelo4: masa puntual + disco uniforme
M4 = 1.11313*(10**21)
sigma4 = 3.11928*(10**19)
M4r = np.sqrt(6.67384*(10**(-11))*(M4 + np.pi*(radios**2)*sigma4)/radios)/1000

#masa puntual + esfera uniforme
M5 = 1.61088*(10**21)
rho5 = 3.23028*(10**18)
M5r = np.sqrt(6.67384*(10**(-11))*(M5 + 4*np.pi*(radios**3)*rho5/3)/radios)/1000


R = radios/8.5
xlargo = [0,8]
xancho = [220,220]
fig = plt.figure()

#Grafico velocidad de rotacion
plt.plot(R, vrot, 'm.')
plt.xlabel("Radio/Ro")
plt.ylabel("Velocidad de rotacion [km/s]")
plt.title("Velocidad de rotacion en funcion del Radio")
plt.axis([0.1, 0.9, 180, 260])
plt.show()
fig.savefig('vel_rot.png')

#Grafico Corrugacion del plano
plt.plot(distancia, Z, 'm-')
plt.axis([4, 8.5, 0, 0.20])
plt.xlabel("Distancia del centro galactico [kpc]")
plt.ylabel("Z [kpc]")
plt.title("Corrugacion del plano")
plt.show()
fig.savefig('corrugacion.png')

#Grafico Modelos de Distribucion de Masa
plt.plot(R, vrot, 'm.', R, M1r, 'r-', R, M2r, 'b-', R, M3r, 'c-', R, M4r, 'g-', R, M5r, 'k-')
plt.title("Modelos de masa")
plt.ylabel("Velocidad de rotacion [km/s]")
plt.xlabel("Radio/Ro")
plt.axis([0.1, 0.9, 159, 300])
plt.show()
fig.savefig('masas.png')
