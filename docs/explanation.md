# Explanation: Cuadratura Gaussiana con Polinomios de Legendre

Es un método numerico para la aproximacion de interales que utiliza puntos de cuadratura $x_i$ (raices de los polinomio s de Legendre) y pesos $w_i$, calculados mediante:

$$
w_i = \frac{2}{(1 - x_i^2) (P_n'(x_i))^2}
$$

La integral se aproxima como:

$$
I \approx \sum_{i=1}^{n} w_i f(x_i)
$$


Este metodo tiene una alta precision, en parte gracias a que utiliza puntos equidistantes, y es exacto para polinomios de gradomenos o igual a 2N-1 con N el numero de nodos.   

