#!/usr/bin/env python3

from scipy.special import legendre
import matplotlib.pyplot as plt
import numpy as np


def gaussxw(N):

    """
    La función `gaussxw` calcula los puntos y pesos de la cuadratura Gauss-Legendre para aproximar integrales.

    Utiliza el método de Newton para encontrar las raíces de los polinomios de Legendre y calcular los pesos correspondientes.
    Este método es eficiente para obtener una buena aproximación de la integral con pocos puntos de muestreo.

    Parámetros:
        N (int): Número de puntos de la cuadratura.

    Retorna:
        tuple: Una tupla `(x, w)` donde:
            - `x` (numpy.ndarray): Los puntos de muestreo obtenidos a partir de los polinomios de Legendre.
            - `w` (numpy.ndarray): Los pesos correspondientes para la cuadratura de Gauss-Legendre.

    Algoritmo:
        1. Se generan los puntos iniciales mediante un coseno modificado.
        2. Se emplea el método de Newton para refinar las raíces de los polinomios de Legendre.
        3. Se calculan los pesos con la fórmula:
        
            w = 2 / ((1 - x^2) * dp^2)

        donde `dp` es la derivada del polinomio de Legendre evaluada en `x`.

    Ejemplo:
        >>> x, w = gaussxw(5)
        >>> print(x)  # Muestra los puntos de muestreo
        array([-0.90617985, -0.53846931,  0.0,  0.53846931,  0.90617985])
        >>> print(w)  # Muestra los pesos correspondientes
        array([0.23692689, 0.47862867, 0.56888889, 0.47862867, 0.23692689])
    """


    # Aproximación inicial
    a = np.linspace(3, 4 * (N - 1), N) / ((4 * N) + 2)
    # Note código vectorial aquí
    x = np.cos(np.pi * a + 1 / (8 * N * N * np.tan(a)))

    # Ahora calculamos las raíces de los polinomios utilizando el método de Newton
    # Este es un tema que veremos la próxima semana!
    # De momento, puede ignorar el siguiente flujo de control con el "while" y saber que esto
    # devuelve los puntos de muestreo obtenidos con los polinomios de Legendre
    epsilon = 1e-15
    delta = 1.0
    while delta > epsilon:
        p0 = np.ones(N, dtype = float)
        # Deep copy
        p1 = np.copy(x)
        for k in range(1, N):
            p0, p1 = p1, ((2 * k + 1) * x * p1 - k * p0) / (k + 1)
        dp = (N + 1) * (p0 - x * p1) / (1 - x * x)
        dx = p1 / dp
        x -= dx
        delta = np.max(np.abs(dx))

    # Ahora calculamos los pesos
    w = 2 * (N + 1) * (N + 1)/(N * N * (1 - x * x) * dp * dp)

    # Note que la función devuelve un tuple
    return x,w



def gaussxwab(a, b, x, w):

    """
    La función `gaussxwab` escala los puntos de la cuadratura Gauss-Legendre al intervalo [a, b].

    Parámetros:
        a (float): Límite inferior del intervalo de integración.
        b (float): Límite superior del intervalo de integración.
        x (numpy.ndarray): Puntos de muestreo de la cuadratura.
        w (numpy.ndarray): Pesos de la cuadratura.

    Retorna:
        tuple: Una tupla `(x, w)` donde:
            - `x` (numpy.ndarray): Los puntos escalados para el intervalo [a, b].
            - `w` (numpy.ndarray): Los pesos escalados para el intervalo [a, b].

    Algoritmo:
        Los puntos `x` y los pesos `w` se escalan multiplicándolos por el factor `0.5 * (b - a)` y ajustando los puntos con una adición de `0.5 * (b + a)`.

    Ejemplo:
        >>> x = np.array([-0.90617985, -0.53846931, 0.0, 0.53846931, 0.90617985])
        >>> w = np.array([0.23692689, 0.47862867, 0.56888889, 0.47862867, 0.23692689])
        >>> x_scaled, w_scaled = gaussxwab(1, 3, x, w)
        >>> print(x_scaled)  # Muestra los puntos escalados al intervalo [1, 3]
        array([0.14636023, 1.17706031, 1.5       , 1.82293969, 2.85363977])
        >>> print(w_scaled)  # Muestra los pesos escalados
        array([0.23692689, 0.47862867, 0.56888889, 0.47862867, 0.23692689])
    """


    return 0.5 * (b - a) * x + 0.5 * (b + a), 0.5 * (b - a) * w






#A continuacion se definen los pesos y los puntos mediante la funcion que utilizamos para la aproximacion gaussiana

xN3,wN3 = gaussxw(3)

xN4,wN4 = gaussxw(4)

xN5,wN5 = gaussxw(5)

xN6,wN6 = gaussxw(6)

xN7,wN7 = gaussxw(7)

xN10,wN10 = gaussxw(10)


#Se utiliza la funcion que escala estos datos 

xN3escalado,wN3escalado = gaussxwab(1,3,xN3,wN3)

xN4escalado,wN4escalado = gaussxwab(1,3,xN4,wN4)

xN5escalado,wN5escalado = gaussxwab(1,3,xN5,wN5)

xN6escalado,wN6escalado = gaussxwab(1,3,xN6,wN6)

xN7escalado,wN7escalado = gaussxwab(1,3,xN7,wN7)

xN10escalado,wN10escalado = gaussxwab(1,3,xN10,wN10)



def poli6(x): #se define la funcion dada 
    return (x**6) - (x**2)*(np.sin(2*x))

    """
    La función `poli6` calcula el valor del polinomio dado en el punto `x`, definido como:

        f(x) = (x**6) - (x**2) * np.sin(2 * x)

    Parámetros:
        x (float o numpy.ndarray): La variable independiente o un arreglo de valores.

    Retorna:
        float o numpy.ndarray: El valor del polinomio evaluado en `x`.

    Algoritmo:
        Calcula el valor de la función `(x**6) - (x**2) * np.sin(2 * x)` para cada valor de `x`.

    Ejemplo:
        >>> poli6(1)
        1.0 - 0.9029742682

    """



#se aplica la sumatoria para obtener las areas 

N3 = np.sum(wN3escalado*poli6(xN3escalado))

N4 = np.sum(wN4escalado*poli6(xN4escalado))

N5 = np.sum(wN5escalado*poli6(xN5escalado))

N6 = np.sum(wN6escalado*poli6(xN6escalado))

N7 = np.sum(wN7escalado*poli6(xN7escalado))

N10 = np.sum(wN10escalado*poli6(xN10escalado))


#Se tomaron valores en N sabiendo que el resulatdo exacto se obtendra para 2N-1<= 6 como se nota en el N=3que no lo cumple y el valor se aleja bastante
#se imprimen los resultados para cada N. Se utilizaron valores que se pudieran acercar al valor de la integral resuelta de forma analitica que es aproximadamente 317 
#Note que los valores mas cercanos son los de N=4 y N=7, inclui uno mas para N=10 pues tal vez obedecia a lagun patron que saltaba de 3 en 3, pero no es asi. 



print(f'aproximacion N=3: {N3}')
print(f'aproximacion N=4: {N4}')
print(f'aproximacion N=5: {N5}')
print(f'aproximacion N=6: {N6}')
print(f'aproximacion N=7: {N7}')
print(f'aproximacion N=10: {N10}')
