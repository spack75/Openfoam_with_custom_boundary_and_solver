# OpenFOAM projet with custom boundary and solver - Suspension Balance Model

## Problem definition 

We want to numerically solve a Couette flow with scalar transport. Thus, the equations are :

- Incompressible Navier-Stokes equations.

- Transport equation : $\frac{\partial C}{\partial t} + \vec{U}\cdot \nabla C = - \vec{\nabla \cdot \vec{j}}$

With $C$ the phase volume fraction, $\vec{U}$ the suspension mean velocity vector and $\vec{j}$ the particle migration flux.

## Solver programmation

The "starting point" solver is icoFoam. In the .C file we have to implement the transport equation, named CEqn.
It is also required to declare CEqn's variables in createField.H 

## Custom boundary condition programmation

We coded the condition $\vec{j}\cdot\vec{n}=0$ which is not provided by default by OpenFOAM.
