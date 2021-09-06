
"""
import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
from scipy import signal
from scipy.fft import fft, ifft
from scipy.signal import lti

"""# PARTE 1: Parámetros del circuito RLC.

"""

# 1. Funcion de transferecia
# Declaramos las variables como simbolicas.
S=sym.Symbol('S')
t=sym.Symbol('t')
# Variables con los datos dados en el lab.
Resis=50
Lind=0.001
Capaci=0.0000001
# Funcion de transferencia.
ft=(1/((Lind*Capaci*(S*S)+(Resis*Capaci*S))+1))

# 2. Función para calcular frecuencia, factor de calidad y factor de amortiguamiento.
def fun(R,L,C,S):
    ft=(1/((L*C*(S*S)+(R*C*S))+1)) #Fun Transferencia
    f_res=1/(np.sqrt(L*C))# frecuencia de resonancia
    q=(f_res*(L/C)) # factor de calidad
    f_am=R/(2*L*f_res) #factor de amortiguamiento
    return(f_res, q, f_am)

# 3. Resultados de la funcion anteriror con los datos designados en este punto.
x,y,z=fun(Resis,Lind, Capaci,S)
print ("PARTE 1 DEL LAB:\n")
print('Frecuencia de resonancia : {} \nFactor de calidad:  {} \nfactor de amortiguamineto: {}'.format(x,y,z))
print (" ")

"""# PARTE 2: Convolución señal escalón, sinusoidal y tren de pulsos.


"""

# 4. Respuesta al impulso h(t)
sys1=signal.lti([1],[Lind*Capaci, Resis*Capaci, 1]) # ft
t, y1 = signal.impulse(sys1)

# Parametros para graficar.
plt.figure() 
plt.grid()
plt.plot(t, y1) 
plt.ylabel("Amplitud", fontsize = 20)
plt.xlabel("Tiempo (s)", fontsize = 20)
plt.title("4. Resp. Pulso H(t)", fontsize = 20)

# 5. Respuesta al escalon unitario
sys1=signal.lti([1],[Lind*Capaci, Resis*Capaci, 1]) # ft
t, y = signal.step2(sys1) # Respuesta a escalón unitario
# Parametros para graficar.
plt.figure()
plt.grid()
plt.plot(t, y) 
plt.ylabel("Amplitud", fontsize = 20)
plt.xlabel("Tiempo (s)", fontsize = 20)
plt.title("5. Resp. Escalon U", fontsize = 20)

# 6. Señal sinoidal con frec= frec_ resonancia
t=np.linspace(0,2,1000)
f0=99999.1
s=np.cos(2*np.pi*f0*t-np.pi/2)
conv=np.convolve(s,y1)   #Hallamos la convolucion entre la señal de exitacion y la respuesta al pulso unitario.

# Parametros para graficar.
plt.figure()
plt.grid()
plt.plot(conv) 
plt.xlim(0,10)
plt.ylabel("Amplitud", fontsize = 20)
plt.xlabel("Tiempo (s)", fontsize = 20)
plt.title("6. Resp. Sinusoidal.", fontsize = 20)

# 7. Diagrama de Bode
from scipy import signal
sys1=signal.lti([1],[Lind*Capaci, Resis*Capaci, 1]) # ft
w, mag, phase = signal.bode(sys1) # Diagrama de bode: frecuencias, magnitud y fase
plt.figure()
plt.semilogx(w, mag, label='Magnitud')    # Bode magnitude plot
plt.legend()
plt.title("7. Diagrama de bode.", fontsize = 20)
plt.figure()
plt.semilogx(w, phase, color='r', label='fase')  # Bode phase plot
plt.legend()
  
# 8 . Tren de pulsos

#junto con el dutty, define el ancho de la señal en alto
u = lambda t: np.piecewise(t,t>=0,[1,0])

def puerta(t,ancho):
    return u(t + ancho/2)-u(t-ancho/2)

T0=5  # Perido de la señal 
dutty=0.5 #Ciclo de dureza

t=np.arange(0, T0, T0/300)
tren=signal.square((99/10)*t,dutty)
# Parametros para graficar.
plt.figure()
plt.plot(t,tren)

plt.ylabel("Amplitud", fontsize = 20)
plt.xlabel("Tiempo (s)", fontsize = 20)
plt.title("8. Tren de pulsos.", fontsize = 20)

# Funcion trenp recibe el periodo de la señal
def trenp(T0):
    Npuntos=300  #Resolucion del muestreado
    Npulsos=80
    tt=np.arange(0, T0, T0/Npuntos)
    x=puerta(tt-T0/4,T0*dutty)
    base=np.ones([1,Npulsos])
    t=np.arange(0,Npulsos*T0,T0/Npuntos)
    pulsos=np.outer(base, x, out=None).reshape(1,Npuntos*Npulsos).T
    return t,pulsos

t,rr=trenp(T0)

t,rr=trenp(T0)

# Convolucion del tren de pulsos con el sistema
conv=np.convolve(np.ravel(rr),np.ravel(y1))
#resp=np.convolve(y1,tren)/fs

# Parametros para graficar.
plt.figure()
plt.plot(conv,color='r',label='convolucion TrenPulsos Vs Sistema')
plt.xlim(0,1200)
plt.legend()
plt.ylabel("Amplitud", fontsize = 20)
plt.xlabel("Tiempo (s)", fontsize = 20)

"""# PARTE 3: Propiedad de Fourier de la convolucion."""

# 9. Realizar la transformada de fourier de la señal sinusoidal del punto 6.

t=np.linspace(0,1,100)
f0=99999.1
s=np.cos(2*np.pi*f0*t-np.pi/2)
yt=fft(s)
sig=yt*y1
yti=ifft(sig)
# Parametros para graficar.
plt.figure()
plt.plot(yti)
plt.ylabel("Amplitud", fontsize = 20)
plt.title("9. transformada de Fourier.", fontsize = 20)
plt.xlabel("Tiempo (s)", fontsize = 20)
plt.xlim(0,100)

plt.show()

