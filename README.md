# Gray-Scott Reaction Diffusion Model
Modeling Gray-Scott Reaction Diffusion Model using  MPI in C++

## C++ MPI Code

Basic [Code](https://github.com/nirajan-mandal/Gray-Scott-Reaction-Diffusion-Model/blob/main/grayscott_final.cpp) flow diagram:

![Code flow chart](https://github.com/nirajan-mandal/Gray-Scott-Reaction-Diffusion-Model/blob/main/Code_Flow_chart_2.jpg "Code flow chart")

## Example: simulation using three nodes. 

Each node has 100x100 grid. The boundry wraps such that the recation occours on a torus surface.

Parameters:

* Du=0.16 
* Dv=0.08 
* k=0.055 
* f=0.02

Last frame output, see vidoe simulation below, that matches leopard's skin pattern.
![Leopard pattern](https://github.com/nirajan-mandal/Gray-Scott-Reaction-Diffusion-Model/blob/main/7_graph102475.jpg "Leopard pattern")

Video Simulation
[![Watch the video](https://github.com/nirajan-mandal/Gray-Scott-Reaction-Diffusion-Model/blob/main/leopard.jpg)](https://github.com/nirajan-mandal/Gray-Scott-Reaction-Diffusion-Model/blob/main/exp7.mp4)



