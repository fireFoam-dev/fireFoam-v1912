/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "CloudFunctionObject.H"
#include "cloud.H" // kvm
#include "Pstream.H" // kvm

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //
#define DEBUG(x) {                                              \
            std::streamsize p = std::cout.precision();              \
            std::ios::fmtflags myFlags;                             \
            myFlags = cout.flags();                                 \
            std::cout.precision(10);                                \
            std::cout.setf(std::ios::fixed,std::ios::floatfield);   \
            std::cout << "["<< __FILE__ << ":" << __LINE__ << "] "; \
            std::cout << "p" << Pstream::myProcNo();                \
            std::cout << " " << #x " = " << x << std::endl;         \
            std::cout.precision(p);                                 \
            std::cout.flags(myFlags);                               \
        }
#define TRACE(s) {                                              \
            std::cout << "["<< __FILE__ << ":" << __LINE__ << "] "; \
            std::cout << "p" << Pstream::myProcNo();                \
            std::cout << " " << #s << std::endl;                    \
            s;                                                      \
        }


template<class CloudType>
void Foam::CloudFunctionObject<CloudType>::write()
{
    NotImplemented;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CloudFunctionObject<CloudType>::CloudFunctionObject(CloudType& owner)
:
    CloudSubModelBase<CloudType>(owner),
    outputDir_(),
    callWriteEveryTime_() // ankur
{}


template<class CloudType>
Foam::CloudFunctionObject<CloudType>::CloudFunctionObject
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName,
    const word& objectType
)
:
    CloudSubModelBase<CloudType>(modelName, owner, dict, typeName, objectType),
    outputDir_(),
    callWriteEveryTime_(this->coeffDict().lookupOrDefault("callWriteEveryTime",false)) // ankur
{
    // Put in undecomposed case
    // (Note: gives problems for distributed data running)

    outputDir_ =
    (
        owner.mesh().time().globalPath()
      / functionObject::outputPrefix
      / this->localPath()
    );

    outputDir_.clean();  // Remove unneeded ".."
}


template<class CloudType>
Foam::CloudFunctionObject<CloudType>::CloudFunctionObject
(
    const CloudFunctionObject<CloudType>& ppm
)
:
    CloudSubModelBase<CloudType>(ppm),
    outputDir_(ppm.outputDir_),
    callWriteEveryTime_(ppm.callWriteEveryTime_) // ankur
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CloudFunctionObject<CloudType>::~CloudFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::CloudFunctionObject<CloudType>::preEvolve()
{}


template<class CloudType>
void Foam::CloudFunctionObject<CloudType>::postEvolve()
{
    if (callWriteEveryTime_) // ankur
    {
        this->write(); // ankur
    }
    else if (this->owner().time().writeTime())
    {
        this->write();
    }
}


template<class CloudType>
void Foam::CloudFunctionObject<CloudType>::postMove
(
    typename CloudType::parcelType&,
    const scalar,
    const point&,
    bool&
)
{}


template<class CloudType>
void Foam::CloudFunctionObject<CloudType>::postPatch
(
    const typename CloudType::parcelType&,
    const polyPatch&,
    bool&
)
{}


template<class CloudType>
void Foam::CloudFunctionObject<CloudType>::postFace
(
    const typename CloudType::parcelType&,
    bool&
)
{}


template<class CloudType>
const Foam::fileName& Foam::CloudFunctionObject<CloudType>::outputDir() const
{
    return outputDir_;
}


template<class CloudType>
Foam::fileName Foam::CloudFunctionObject<CloudType>::writeTimeDir() const
{
    return outputDir_/this->owner().time().timeName();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CloudFunctionObjectNew.C"

// ************************************************************************* //
