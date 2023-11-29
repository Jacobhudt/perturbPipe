# perturbPipe
Perturbation application for pipes in OpenFOAMv2306 based on perturbU from Eugene de Villiers. 

## Installation and usage
Install with 

```sh
wmake all
```
Put the perturbPipeDict in the /constant directory and change parameters.


Correct formatting in "transportDict":
```c++

Ubar Ubar [0 1 -1 0 0 0 0] (0 0 0.795);

transportModel Newtonian;

nu nu [0 2 -1 0 0 0 0] 1.5e-05; 

Retau Retau [0 0 0 0 0 0 0] 180; 

```