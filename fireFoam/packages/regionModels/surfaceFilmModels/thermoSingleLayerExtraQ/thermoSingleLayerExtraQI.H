/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "thermoSingleLayerExtraQ.H"
#include "heatTransferModel.H"
#include "filmRadiationModel.H"
#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const SLGThermo& thermoSingleLayerExtraQ::thermo() const
{
    return thermo_;
}


inline label thermoSingleLayerExtraQ::liquidId() const
{
    return liquidId_;
}


inline tmp<scalarField> thermoSingleLayerExtraQ::hs
(
    const scalarField& T,
    const label patchI
) const
{
    const scalarField& Cp = Cp_.boundaryField()[patchI];
    return Cp*(T - 298.15);
}


inline tmp<volScalarField> thermoSingleLayerExtraQ::hs
(
    const volScalarField& T
) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "hs(" + T.name() + ")",
                time().timeName(),
                regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            Cp_*(T - (dimensionedScalar("Tstd", dimTemperature, 298.15))),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}


inline tmp<volScalarField> thermoSingleLayerExtraQ::T
(
    const volScalarField& hs
) const
{
    tmp<volScalarField> tT
    (
        new volScalarField
        (
            IOobject
            (
                "T(" + hs.name() + ")",
                time().timeName(),
                regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            hs/Cp_ + dimensionedScalar("Tstd", dimTemperature, 298.15),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tT().min(Tmax_);
    tT().max(Tmin_);

    return tT;
}


inline const volScalarField& thermoSingleLayerExtraQ::hsSp() const
{
    return hsSp_;
}


inline const volScalarField& thermoSingleLayerExtraQ::hsSpPrimary() const
{
    return hsSpPrimary_;
}


inline const volScalarField& thermoSingleLayerExtraQ::TPrimary() const
{
    return TPrimary_;
}


inline const PtrList<volScalarField>& thermoSingleLayerExtraQ::YPrimary() const
{
    return YPrimary_;
}


inline const heatTransferModel& thermoSingleLayerExtraQ::htcs() const
{
    return htcs_();
}


inline const heatTransferModel& thermoSingleLayerExtraQ::htcw() const
{
    return htcw_();
}


inline const phaseChangeModel& thermoSingleLayerExtraQ::phaseChange() const
{
    return phaseChange_();
}


inline const filmRadiationModel& thermoSingleLayerExtraQ::radiation() const
{
    return radiation_();
}


inline tmp<scalarField> thermoSingleLayerExtraQ::Qconvw(const label patchI) const
{
    const scalarField htc(htcw_->h()().boundaryField()[patchI]);
    const scalarField& Tp = T_.boundaryField()[patchI];
    const scalarField& Twp = Tw_.boundaryField()[patchI];

    return htc*(Tp - Twp);
}


inline tmp<scalarField> thermoSingleLayerExtraQ::Qconvp(const label patchI) const
{
    const scalarField htc(htcs_->h()().boundaryField()[patchI]);
    const scalarField& Tp = T_.boundaryField()[patchI];
    const scalarField& Tpp = TPrimary_.boundaryField()[patchI];

    return htc*(Tp - Tpp);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
