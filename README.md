## Non Linear Model Predictive Control

The work done here deals with 2 methods for determining Optimal Temperature Control profile for batch Cyrstallization process.

1. **Deterministic Optimal Control** aims at finding the an optimum temperature profile to maximise an objective function selected to achieve a desired volume of the product. Herein, the experimental kinetic parameters are employed to simulate a batch crystalllization process.
2. **Stochastic Optimal Control** undertakes the task of quantifying the uncertainites which creep in due to experimentation. It aims to achive a maximum expected value for the desired product, simultaneously incorporating randomness in the process parameters into the model. Namely, Two methods Ito Process and Polynomial Chaos Expansions are used.

### Deterministic Optimal Control 
This undertakes the use of [**Maximum Principle**](https://en.wikipedia.org/wiki/Pontryagin%27s_maximum_principle) using the Hamiltonian Derivative to move towards the _optimum_ value of Temperature at each time step.
