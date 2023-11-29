# perturbPipe
Perturbation application for pipes in OpenFOAMv2306 based on perturbU from Eugene de Villiers. 

## Installation and usage
Install with 
'''shell
wmake all
'''
Put the perturbPipeDict in the /constant directory and change parameters.


Correct formatting in "transportDict":
'''shell
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2106                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      transportProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


Ubar Ubar [0 1 -1 0 0 0 0] (0 0 0.795);

transportModel Newtonian;

nu nu [0 2 -1 0 0 0 0] 1.5e-05; 

Retau Retau [0 0 0 0 0 0 0] 180; 



// ************************************************************************* //
'''