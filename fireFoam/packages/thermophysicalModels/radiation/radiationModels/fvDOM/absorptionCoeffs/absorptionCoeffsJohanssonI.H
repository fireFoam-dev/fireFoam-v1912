/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

inline bool Foam::radiation::absorptionCoeffsJohanssen::invTemp() const
{
    return  invTemp_;
}


inline Foam::scalar Foam::radiation::absorptionCoeffsJohanssen::Tcommon() const
{
    return  Tcommon_;
}


inline Foam::scalar Foam::radiation::absorptionCoeffsJohanssen::Tlow() const
{
    return  Tlow_;
}


inline Foam::scalar Foam::radiation::absorptionCoeffsJohanssen::Thigh() const
{
    return  Thigh_;
}


inline const Foam::radiation::absorptionCoeffsJohanssen::coeffArray&
Foam::radiation::absorptionCoeffsJohanssen::highACoeffs() const
{
    return  highACoeffs_;
}


inline const Foam::radiation::absorptionCoeffsJohanssen::coeffArray&
Foam::radiation::absorptionCoeffsJohanssen::lowACoeffs() const
{
    return  lowACoeffs_;
}


// ************************************************************************* //
