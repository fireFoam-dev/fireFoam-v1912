/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    fireFoam

Group
    grpCombustionSolvers

Description
    Transient solver for fires and turbulent diffusion flames with reacting
    particle clouds, surface film and pyrolysis modelling.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "turbulentFluidThermoModel.H"
#include "basicReactingCloud.H"
#include "surfaceFilmModel.H"
#include "pyrolysisModelCollection.H"
#include "radiationModel.H"
#include "SLGThermo.H"
#include "solidChemistryModel.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"

#include "singleStepCombustion.H"
#include "thermoPhysicsTypes.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCase.H"
    #include "printVersion.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "infoFieldsOutput.H" // fm
    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"
    #include "createRegionControls.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    Info<< "Start cpuTimeIncrement: "<<runTime.cpuTimeIncrement()<<endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "solidRegionDiffusionNo.H"
        #include "setMultiRegionDeltaT.H"

        if ( !updateMesh )
        {
           #include "setDeltaT.H"
        }
        else
        {
            #include "dySetDeltaT.H"
            #include "updateMesh.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        Info<< "timeSpent-in-others: " << setprecision(5) << runTime.cpuTimeIncrement() << " s" << endl;

        parcels.evolve();

        Info<< "timeSpent-in-spray: " << setprecision(5) << runTime.cpuTimeIncrement() << " s" << endl;

        surfaceFilm.evolve();

        Info<< "timeSpent-in-film: " << setprecision(5) << runTime.cpuTimeIncrement() << " s" << endl;

        if(solvePyrolysisRegion)
        {
            pyrolysis.evolve();
        }

        Info<< "timeSpent-in-pyrolysis: " << setprecision(5) << runTime.cpuTimeIncrement() << " s" << endl;

        if (solvePrimaryRegion)
        {
            #include "rhoEqn.H"

            // --- PIMPLE loop
            while (pimple.loop())
            {
                #include "UEqn.H"
                #include "YEEqn.H"

                // --- Pressure corrector loop
                while (pimple.correct())
                {
                    #include "pEqn.H"
                }

                if (pimple.turbCorr())
                {
                    turbulence->correct();
                }
            }

            rho = thermo.rho();

            #include "infoOutput.H" // fm

        }

        Info<< "timeSpent-in-gasphase: " << setprecision(5) << runTime.cpuTimeIncrement() << " s" << endl;

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"; // kvm

        if (runTime.writeTime()) // kvm
        {
            Info<< " +"; // kvm
        }
        Info<< nl << endl; // kvm

    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
