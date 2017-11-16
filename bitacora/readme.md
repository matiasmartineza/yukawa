# Bitácora

## 16 / Noviembre / 2017

#### NOTA:
- Paralelización de M2P en el árbol de yukawa.
- Es necesario inicializar las variables privadas que vienen desde afuera dentro de un pragama omp for al iterar sobre un loop. Esto se hizo en las variables a_aux, ax_aux, ay_aux, az_aux.
- El pragama omp parallel for se ha introducido en la iteración sobre los targets “i” dentro de un twig determinado.

#### Experimento: 
Se hizo correr el código para 1000000 de targets y sources, en una versión paralelizada y en una versión sin paralelizar. Los resultados muestran un speedup de 3.9408 para el número de threads por defecto.

![alt text](https://github.com/matiasmartineza/yukawa/blob/master/bitacora/exp_1.png "Experimento 1")

#### Experimento: 
Se hizo correr el código para verificar la cantidad de tiempo tomado en el proceso de M2P según distintos número de threads. Los resultados se ven en la figura del Experimento 2.

![alt text](https://github.com/matiasmartineza/yukawa/blob/master/bitacora/exp_2.pdf "Experimento 2")

#### Siguientes pasos: 
- Experimentar para un número fijo de threads, cambiando el número de partículas.
- Paralelizar P2P.
