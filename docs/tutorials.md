# Ejemplo de Uso
Aquí se muestra un ejemplo básico de cómo utilizar la implementación de la cuadratura gaussiana para calcular la integral de una función.


La resolucion analitica es la siguiente: 


Si utilizamos N=2 para aproximar la integral de $f(x) = e^x$ en $[-1,1]$. Los nodos y pesos de Gauss-Legendre son:

Para el ppolinomio de legendre de segundo grado se obtienen las raices de

$P_2= \frac{1}{2}(3x^2 -1)$ que son

 $x_1 = -\frac{1}{\sqrt{3}}, x_2 = \frac{1}{\sqrt{3}}$

al ser las raices de la funcion y los pesos de la funcion estan dados y son iguales para N=2

 $w_1 = w_2 = 1$. La aproximación de la integral de \(f(x) = e^x\) en \([-1,1]\) es:



Al aplicar
$$
I \approx \sum_{i=1}^{n} w_i f(x_i)
$$

Se obtiene que la integral de nuestra funcion es aproximadamente:

$$
I \approx e^{-1/\sqrt{3}} + e^{1/\sqrt{3}}.
$$



Para implementarla en python sera necesario algo del tipo: 

```python

from cuadrature.py import gaussxw

from cuadrature.py import gaussxwab



# Función a integrar

def f(x):
    return e**x

# Llamada a la función de cuadratura gaussiana y el escalado con 5 puntos

xN5,xN5=gaussxw(5)

xN5escalado,wN5escalado = gaussxwab(-1,1,xN5,wN5)

resultado = np.sum(wN5escalado*f(xN5escalado))




print(f"Resultado de la integral: {resultado}")

