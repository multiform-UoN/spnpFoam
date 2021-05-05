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
and Poisson potential (electrostatic approximation).
Coupling is implemented using a segregated algorithm.

Authors:
Federico Municchi, Nottingham (2019)
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
//#include "constrainFluxes.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "electrokineticConstants.H"
using namespace electrokineticConstants;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Define specie from dictionary
typedef dictionary specie;

int main(int argc, char *argv[])
{
  #include "postProcess.H"

	#include "setRootCaseLists.H"
	#include "createTime.H"
  #include "createTimeControls.H"
	#include "createMesh.H"

	//Add pimple coupling controls
	pimpleControl pimple(mesh);

	#include "createFields.H"
	#include "initContinuityErrs.H"

	#include "setPnpCourantNos.H"
	#include "setInitialDeltaT.H"

	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	Info<< "\nStarting time loop\n" << endl;

	while (runTime.loop())
	{
		Info<< "Time = " << runTime.timeName() << nl << endl;

		// --- Coupling loop
		while (pimple.loop())
		{
      #include "readTimeControls.H"

			rho *= scalar(0); //- Simply reset charge density to zero

			if (!electroNeutrality)
			//-Update charge density
			{
				forAll(species,sp)
				{
					dimensionedScalar Z(species[sp].lookup("Z"));
										rho += Fy_*cPtrL[sp]*Z;
				}
			}

			//- Solve Stokes
			if(!pimple.frozenFlow())
			{
				#include "UEqn.H"

				while(pimple.correct())
				{
						#include "pEqn.H"
				}
			}

			scalar lambda(1.0);//lambda(0.01); // HARD-CODED body force switch

			//- Solve Poisson
			while (pimple.correctNonOrthogonal())
			{
				//V.storePrevIter();
				fvScalarMatrix VEqn
				(
					- fvm::laplacian(epsilon_,V) // Minus on laplacian side to get positive definite
					==
					lambda*rho
				);

				VEqn.relax();
        // Set the reference value for V to RefCell & RefValue, do nothing if not
				VEqn.setReference(VRefCell, VRefValue);

				VEqn.solve();
				V.relax();
			}

			//- Solve species
			forAll(species,sp)
			{

				//Gather pointers to specific ion species sp
				volScalarField&   C(cPtrL[sp]);
				dimensionedScalar D(species[sp].lookup("D"));
				dimensionedScalar Z(species[sp].lookup("Z"));

				//- Update Nerst-Planck flux
				phiNP = phi - fvc::flux(e_*fvc::grad(V)/(k_*T))*Z*D;

				//Check if electro-neutrality is wanted
				if (electroNeutrality && sp == species.size()-1)
				{

					//Look up the first species
					const dimensionedScalar& Z_ini(species[0].lookup("Z"));
					const volScalarField& C_ini(cPtrL[0]);

					//Set concentration C to first species
					C  = -Z_ini*C_ini/Z;

					//Loop through all species except final species
					for (int i=1;i <= species.size()-2;i++)
					{
						//Look up the next species
						const dimensionedScalar& Z_iter(species[i].lookup("Z"));
						const volScalarField& C_iter(cPtrL[i]);

						//Add next species to C
						C -=  Z_iter*C_iter/Z;

					}

					C.correctBoundaryConditions();

				}
				else
				{

					//- Non-orthogonal correction loop
					while (pimple.correctNonOrthogonal())
					{

						//C.storePrevIter();
						fvScalarMatrix CEqn
						(
							fvm::ddt(C)
							+ fvm::div(phiNP,C,"div(phi,C)")
							- fvm::laplacian(D,C,"laplacian(D,C)")
							==
							fvOptions(C)
						);

						//constrainFluxes(CEqn);

						CEqn.relax();
					  fvOptions.constrain(CEqn);
						CEqn.solve();
						fvOptions.correct(C);

						//- Regularise solution
						{
							//- Force solution lower bound
							//C = (mag(C) + C)/scalar(2.0); // ?? Unclear on what this is -R
							//C.correctBoundaryConditions();

							//- Relax
							C.relax();
						}
					}
				}
			}
		}


		#include "setPnpCourantNos.H"

		#include "setDeltaT.H"

		forAll(species,sp)
		{
			Info << "Average concentration specie " << sp << " = "
			<< Foam::gSum(cPtrL[sp]().field()*mesh.V())/Foam::gSum(mesh.V())
			<< endl;
		}

		runTime.write();

		Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
		<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
		<< nl << endl;
	}

	Info<< "End\n" << endl;

	return 0;
}


// ************************************************************************* //
