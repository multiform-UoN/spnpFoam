/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    pnpFoam

Description
    Transient solver for multiple species transport with Nernst-Planck forcing
    and Poisson potential (electrostatic approximation) within multiple fluid/solid
    regions.
    Coupling is implemented using a segregated algorithm (PISO).
    Coupling between meshes using PIMPLE.

Authors:
    Federico Municchi, Nottingham (2019)
    Robert Barnett, Nottingham (2020)
\*---------------------------------------------------------------------------*/


#include "fvCFD.H"
//#include "constrainFluxes.H"
#include "regionProperties.H"
#include "fvOptions.H"
#include "coordinateSystem.H"
#include "pimpleMultiRegionControl.H"
#include "electrokineticConstants.H"
using namespace electrokineticConstants;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Define specie from dictionary
typedef dictionary specie;

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    // Construct region controls
    pimpleMultiRegionControl pimples(fluidRegions, solidRegions);

    #include "createTimeControls.H"
    #include "setAllPnpCourantNos.H"
    #include "setInitialDeltaT.H"

    Info<< "\nStarting time loop\n" << endl;

    while (pimples.run(runTime))
    {
        #include "readTimeControls.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- PIMPLE loop
        while(pimples.loop())
        {
            //- Solve for each fluid region
            forAll(fluidRegions,i)
            {
                Info<< "\nSolving for fluid region "
                    << fluidRegions[i].name() << endl;

                //- Set fields, mesh and fvOptions to new fluid region
                #include "setRegionFluidFields.H"

                //- Construct coupling controls for new mesh
                #include "solveFluid.H"
            }

            //- Solve for each solid region
            forAll(solidRegions,i)
            {
                Info<< "\nSolving for solid region "
                    << solidRegions[i].name() << endl;

                //- Set fields, mesh and fvOptions to new solid region
                #include "setRegionSolidFields.H"

                #include "solveSolid.H"
            }
            Info << endl;
        }

        #include "setAllPnpCourantNos.H"

        // Update maximum Co # of all regions and species
        CoNum = gMax(MaxCoFluid);
        // TODO COMPARE THIS WITH THE DIFFUSION NUMBERS AND CHOOSE THE MAXIMUM

        // Update global time step using algorithm
        #include "setDeltaT.H"

        #include "writeStats.H"

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

   return 0;
}


// ************************************************************************* //
