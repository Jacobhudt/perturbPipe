/*---------------------------------------------------------------------------*\
Confirmed to work with v2021 (.com) version.

Custom modification of perturbU utilities by Eugene de Villiers for perturbation in both stream, and spanwise direction in pipes.

Spanwise perturbations are introduced in a same manner as what is done in channels, with some small modifications.
Contrary to the original, this utility for turbulent pipes allows modifications of the paramateres post-compile via the perturbPipe.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "Random.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

    const vectorField centers(mesh.C());

    Info<< "Time = " << runTime.value() << endl;


    IOobject Uheader
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );


    //// Check U exists
    if (Uheader.typeHeaderOk<volVectorField>(true))
    {
        Info<< "    Reading U" << endl;
        volVectorField U(Uheader, mesh);

            IOdictionary perturbDict
            (
             IOobject
                (
                "perturbPipe",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
                )
             );


        IOdictionary transportProperties
        (
            IOobject
            (
                "transportProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )

        
        );
        dimensionedScalar nu
        (
            transportProperties.lookup("nu")
        );
        dimensionedVector Ubar
        (
            transportProperties.lookup("Ubar")
        );
 
        /* const vector Ubar(readLabel(transportProperties.lookup("Ubar")));

        const scalar nu(readLabel(transportProperties.lookup("nu")));
 */
        const direction streamDir(readLabel(perturbDict.lookup("streamwise")));
        const direction spanDir(readLabel(perturbDict.lookup("spanwise")));

        Info<< "Streamwise flow component = " << streamDir << nl
            << "Spanwise flow component   = " << spanDir << nl
            << endl;

        
        if (streamDir > 2 || spanDir > 2 || streamDir == spanDir)
        {
            FatalErrorIn(args.executable())
                << "Spanwise and streamwise components have to be 0,1 or 2 and"
                << " differ from one another." << nl
                << "streamDir:" << streamDir
                << " spanDir:" << spanDir
                << exit(FatalError);
        }

        // Get component normal to streamDir and spanDir. This is the height
        // component.
        direction heightDir = 0;
        if (streamDir == heightDir)
        {
            heightDir++;
        }
        if (spanDir == heightDir)
        {
            heightDir++;
        }

        
        
        //Read from perturbDict
        const scalar Retau(readScalar(perturbDict.lookup("Retau")));

        //Characteristic length
        const scalar d(readScalar(perturbDict.lookup("d")));
        
        const scalar utau = Retau*nu.value()/d;
        //wall normal circulation
        const scalar duplus = Ubar.value()[0]*0.5/utau;

        //spanwise wavenumber spacing
        const scalar bP(readScalar(perturbDict.lookup("betaPlus")));
        const scalar betaPlus = 2.0*constant::mathematical::pi*(1.0/bP);
        

        const scalar sigma(readScalar(perturbDict.lookup("sigma")));
        
        //streamwise wave number spacing
        const scalar aP(readScalar(perturbDict.lookup("alphaPlus")));
        const scalar alphaPlus = 2.0*constant::mathematical::pi*(1.0/aP);
        
        //pertubation amplitude
        //const scalar epsilon(readScalar(perturbDict.lookup("epsilon")));
        const scalar epsilonStream(readScalar(perturbDict.lookup("epsilonStream")));
        const scalar epsilonSpan(readScalar(perturbDict.lookup("epsilonSpan")));

        const scalar perturbationSeed(readScalar(perturbDict.lookup("perturbationSeed")));
        //Random perturbation(1234567);
        //Random perturbation(3456);
        Random perturbation(perturbationSeed);

        //Set velocity = 0 to avoid addition if stream/span direction is changed
        Info<< "Initializing U " << endl;
        forAll(centers, celli)
        {
            U[celli][0] = 0;
            U[celli][1] = 0;
            U[celli][2] = 0;
        }
        U.write();

        //Do pertubation
        Info << "Performing pertubation" << endl;

        forAll(centers, celli)
        {
            scalar deviation=1.0 + 0.2*perturbation.GaussNormal<scalar>();

            scalar& Uz(U[celli].z());
            vector cCenter = centers[celli];
            scalar r = ::sqrt(::sqr(cCenter.y()) + ::sqr(cCenter.x()));
            Uz = 2*mag(Ubar.value())*(1-::sqr(r/d));
            scalar y = d - r;
            r = r*Retau/d;
            y = y*Retau/d;

            scalar theta = ::atan(cCenter.y()/cCenter.x());
            scalar z = cCenter.z()*Retau/d;


            Uz = Uz + (utau*duplus/2.0)
                 *::cos(betaPlus*theta*r) *(y/30)
                 *::exp(-sigma*::sqr(y) + 0.5);
            
            scalar utheta = epsilonStream*::sin(alphaPlus*z)*y
                            *::exp(-sigma*::sqr(y));


            U[celli][streamDir] = U[celli][streamDir] + utheta*Uz;

            U[celli][spanDir] =   epsilonSpan
              * Foam::sin(alphaPlus*z)
              * y
              * Foam::exp(-sigma*Foam::sqr(y))
              * deviation;

        }

        U.write();
    }
    else
    {
        Info<< "    No U" << endl;
    }

    Info<< "Pertubation was successful. " << endl;


    return(0);
}


// ************************************************************************* //