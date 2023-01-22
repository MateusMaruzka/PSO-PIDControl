# PSO-PIDControl
Optimization of a discrete PID control 

Trabalho de conclusão de curso - UTFPR

## PSO (Particle Swarm Optimization) Library
This code is an implementation of the Particle Swarm Optimization (PSO) algorithm for solving optimization problems. The PSO algorithm is a heuristic optimization method that simulates the behavior of a swarm of birds or insects searching for food. It is inspired by the social behavior of animals and is used to find the global minimum or maximum of a function.

## Getting Started
To use this code, you will need to have Python 3 and the following libraries installed:

numpy
matplotlib
You can install these libraries using pip by running the following command in your terminal:
```
pip install numpy matplotlib
```
## Usage
The main function of the code is the pso function, which takes the following arguments:

* <b>``fObj``</b>: a function that computes the fitness of a particle
* <b>``T_ENXAME``</b>: the number of particles in the swarm
* **``DIM``**: the dimension of the problem
* **``iterMax``**: the maximum number of iterations (default: 100)
* **``var``**: the initial range of the particles (default: 15)
* **``_Wmin``**: the minimum value of the inertia weight (default: 0.1)
* **``_Wmax``**: the maximum value of the inertia weight (default: 0.9)
* **``_c1``**: the cognitive coefficient (default: 2.05)
* **``_c2``**: the social coefficient (default: 2.05)

It also accepts any additional keyword arguments that will be passed to the fitness function.

Here is an example of how to use the pso function to minimize the sphere function:

```python
import numpy as np

def sphere(x, **kwargs):
    return np.sum(x**2, axis=1)

x_opt, f_opt = pso(sphere, T_ENXAME=100, DIM=5)
```
The function returns the optimal solution ``x_opt`` and the optimal fitness ``f_opt``.

This code is a basic implementation of the PSO algorithm that can be used to solve optimization problems. It can be easily modified to suit the specific requirements of your problem.

You can find more details in the code documentation and comments.

Authors
Mateus Maruzka Roncaglio
Prof Daniel Cavalcanti Jeronymo
Wesley Klewerton
