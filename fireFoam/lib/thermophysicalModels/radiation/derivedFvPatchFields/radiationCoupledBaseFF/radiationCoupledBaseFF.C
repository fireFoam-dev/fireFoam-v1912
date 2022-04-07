/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "radiationCoupledBaseFF.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "fvPatchFieldMapper.H"
#include "radiationModel.H"
#include "absorptionEmissionModel.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::radiationCoupledBaseFF::emissivityMethodType,
        3
    >::names[] =
    {
        "solidRadiation",
        "pyrolysisModel",
        "lookup"
    };
}


const Foam::NamedEnum<Foam::radiationCoupledBaseFF::emissivityMethodType, 3>
    Foam::radiationCoupledBaseFF::emissivityMethodTypeNames_;


namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::radiationCoupledBaseFF::absorptivityMethodType,
        4
    >::names[] =
    {
        "emissivity",
        "pyrolysisModel",
        "lookup",
        "solidRadiation"
    };
}


const Foam::NamedEnum<Foam::radiationCoupledBaseFF::absorptivityMethodType, 4>
    Foam::radiationCoupledBaseFF::absorptivityMethodTypeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationCoupledBaseFF::radiationCoupledBaseFF
(
    const fvPatch& patch,
    const word& calculationTypeE,
    const word& calculationTypeA,
    const scalarField& emissivity,
    const scalarField& absorptivity
)
:
    patch_(patch),
    method_(emissivityMethodTypeNames_[calculationTypeE]),
    absMethod_(absorptivityMethodTypeNames_[calculationTypeA]),
    emissivity_(emissivity),
    absorptivity_(absorptivity)
{}


Foam::radiationCoupledBaseFF::radiationCoupledBaseFF
(
    const fvPatch& patch,
    const dictionary& dict
)
:
    patch_(patch),
    method_(emissivityMethodTypeNames_.read(dict.lookup("emissivityMode"))),
    absMethod_(absorptivityMethodTypeNames_[dict.lookupOrDefault<Foam::word>("absorptivityMode",word("emissivity"))])
    //absMethod_(absorptivityMethodTypeNames_.read(dict.lookupOrDefault<Foam::Istream>("absorptivityMode","emissivity")))
    //absMethod_(absorptivityMethodTypeNames_.read(string(dict.lookupOrDefault<Foam::word>("absorptivityMode",word("emissivity")))))
    //absMethod_(absorptivityMethodTypeNames_.read(refCast<Foam::Istream>(dict.lookupOrDefault<Foam::word>("absorptivityMode",word("emissivity")))))
    //absMethod_(absorptivityMethodTypeNames_.read(refCast<Foam::Istream&>(dict.lookupOrDefault<const word&>("absorptivityMode","emissivity"))))
{
    switch (method_)
    {
        case SOLIDRADIATION:
        {
            if (!isA<mappedPatchBase>(patch_.patch()))
            {
                FatalIOErrorIn
                (
                    "radiationCoupledBaseFF::radiationCoupledBaseFF\n"
                    "(\n"
                    "    const fvPatch& p,\n"
                    "    const dictionary& dict\n"
                    ")\n",
                    dict
                )   << "\n    patch type '" << patch_.type()
                    << "' not type '" << mappedPatchBase::typeName << "'"
                    << "\n    for patch " << patch_.name()
                    << exit(FatalIOError);
            }

            emissivity_ = scalarField(patch_.size(), 0.0);
        }
        break;

        case PYROLYSISMODELE:
        {
            if (!isA<mappedPatchBase>(patch_.patch()))
            {
                FatalIOErrorIn
                (
                    "radiationCoupledBaseFF::radiationCoupledBaseFF\n"
                    "(\n"
                    "    const fvPatch& p,\n"
                    "    const dictionary& dict\n"
                    ")\n",
                    dict
                )   << "\n    patch type '" << patch_.type()
                    << "' not type '" << mappedPatchBase::typeName << "'"
                    << "\n    for patch " << patch_.name()
                    << exit(FatalIOError);
            }

           // Get the coupling information from the mappedPatchBase
            const mappedPatchBase& mpp =
                refCast<const mappedPatchBase>(patch_.patch());

            const polyMesh& nbrMesh = mpp.sampleMesh();

            HashTable<const pyroModelType*> models =
                nbrMesh.time().lookupClass<pyroModelType>();

            bool regionMatch = false;
            forAllConstIter(HashTable<const pyroModelType*>, models, iter)
            {
                if (iter()->regionMesh().name() == nbrMesh.name())
                {
                    regionMatch = true; 
                }
            } 

//            const IOdictionary& pyroZones  =
//                db().lookupObject<IOdictionary>
//                (
//                    "pyrolysisZones"
//                );

//            const wordList regionNames(pyroZones.lookUp("regionName"));

//            setSize(regionNames.size());

 //           for (label i = 0; i < regionNames.size(); i++)
 //           {  
 //             if (regionNames[i] == nbrMesh.name())
 //             {
 //               regionMatch = true; 
 //             }  
 //           } 
           
            if (!regionMatch)
            {
                FatalIOErrorIn
                (
                    "radiationCoupledBaseFF::radiationCoupledBaseFF\n"
                    "(\n"
                    "    const fvPatch& p,\n"
                    "    const dictionary& dict\n"
                    ")\n",
                    dict
                )   << "\n    Patch '" << patch_.name() 
                    << "' not connected to a pyrolysis region. Use different emissivity mode."
                    << exit(FatalIOError);

            }

            emissivity_ = scalarField(patch_.size(), 0.0);
        }
        break;

        case LOOKUPE:
        {
            if (!dict.found("emissivity"))
            {
                FatalIOErrorIn
                (
                    "radiationCoupledBaseFF::radiationCoupledBaseFF\n"
                    "(\n"
                    "    const fvPatch& p,\n"
                    "    const dictionary& dict\n"
                    ")\n",
                    dict
                )   << "\n    emissivity key does not exist for patch "
                    << patch_.name()
                    << exit(FatalIOError);
            }
            else
            {
                emissivity_ = scalarField("emissivity", dict, patch_.size());
            }
        }
        break;
    }

    switch (absMethod_)
    {
        case EMISSIVITY:
        {
            absorptivity_ = scalarField(emissivity_);
        }
        break;

        case PYROLYSISMODELA:
        {
            if (emissivityMethod() != absorptivityMethod())
            {

                if (!isA<mappedPatchBase>(patch_.patch()))
                {
                  FatalIOErrorIn
                  (
                      "radiationCoupledBaseFF::radiationCoupledBaseFF\n"
                      "(\n"
                      "    const fvPatch& p,\n"
                      "    const dictionary& dict\n"
                      ")\n",
                      dict
                  )   << "\n    patch type '" << patch_.type()
                      << "' not type '" << mappedPatchBase::typeName << "'"
                      << "\n    for patch " << patch_.name()
                      << exit(FatalIOError);
                }

                // Get the coupling information from the mappedPatchBase
                const mappedPatchBase& mpp =
                    refCast<const mappedPatchBase>(patch_.patch());

                const polyMesh& nbrMesh = mpp.sampleMesh();

                HashTable<const pyroModelType*> models =
                    nbrMesh.time().lookupClass<pyroModelType>();

                bool regionMatch = false;
                forAllConstIter(HashTable<const pyroModelType*>, models, iter)
                {
                   if (iter()->regionMesh().name() == nbrMesh.name())
                   {
                    regionMatch = true; 
                   }
                } 

//            const IOdictionary& pyroZones  =
//                db().lookupObject<IOdictionary>
//                (
//                    "pyrolysisZones"
//                );

//            const wordList regionNames(pyroZones.lookUp("regionName"));

//            setSize(regionNames.size());

 //           for (label i = 0; i < regionNames.size(); i++)
 //           {  
 //             if (regionNames[i] == nbrMesh.name())
 //             {
 //               regionMatch = true; 
 //             }  
 //           } 
           
                if (!regionMatch)
                {
                  FatalIOErrorIn
                  (
                      "radiationCoupledBaseFF::radiationCoupledBaseFF\n"
                      "(\n"
                      "    const fvPatch& p,\n"
                      "    const dictionary& dict\n"
                      ")\n",
                      dict
                  )   << "\n    Patch '" << patch_.name() 
                      << "' not connected to a pyrolysis region. Use different emissivity mode."
                      << exit(FatalIOError);

                }
               
            }

            absorptivity_ = scalarField(patch_.size(), 0.0);
 
        }
        break;

        case LOOKUPA:
        {
            if (!dict.found("absorptivity"))
            {
                FatalIOErrorIn
                (
                    "radiationCoupledBaseFF::radiationCoupledBaseFF\n"
                    "(\n"
                    "    const fvPatch& p,\n"
                    "    const dictionary& dict\n"
                    ")\n",
                    dict
                )   << "\n    absorptivity key does not exist for patch "
                    << patch_.name()
                    << exit(FatalIOError);
            }
            else
            {
                absorptivity_ = scalarField("absorptivity", dict, patch_.size());
            }
        }
        break;

        case SOLIDRADIATIONA:
        {
            if (!isA<mappedPatchBase>(patch_.patch()))
            {
                FatalIOErrorIn
                (
                    "radiationCoupledBaseFF::radiationCoupledBaseFF\n"
                    "(\n"
                    "    const fvPatch& p,\n"
                    "    const dictionary& dict\n"
                    ")\n",
                    dict
                )   << "\n    patch type '" << patch_.type()
                    << "' not type '" << mappedPatchBase::typeName << "'"
                    << "\n    for patch " << patch_.name()
                    << exit(FatalIOError);
            }

            absorptivity_ = scalarField(patch_.size(), 0.0);
        }
        break;
    }

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::radiationCoupledBaseFF::
pyroModelType&
Foam::radiationCoupledBaseFF::
pyroModel(const fvMesh& meshRef,const word& regName) const
{
    HashTable<const pyroModelType*> models =
        meshRef.time().lookupClass<pyroModelType>();

    forAllConstIter(HashTable<const pyroModelType*>, models, iter)
    {
        if (iter()->regionMesh().name() == regName)
        {
            return *iter();
        }
    }

    FatalErrorIn
    (
        "const Foam::radiationCoupledBaseFF::"
        "pyroModelType& "
        "Foam::radiationCoupledBaseFF::"
        "pyroModel() const"
    )
        << "Unable to locate pyrolysis region " << regName
        << abort(FatalError);

    return **models.begin();
}



Foam::scalarField Foam::radiationCoupledBaseFF::emissivity()
{
    switch (method_)
    {
        case SOLIDRADIATION:
        {
            // Get the coupling information from the mappedPatchBase
            const mappedPatchBase& mpp =
                refCast<const mappedPatchBase>(patch_.patch());

            const polyMesh& nbrMesh = mpp.sampleMesh();

            const radiation::radiationModel& radiation =
                nbrMesh.lookupObject<radiation::radiationModel>
                (
                    "radiationProperties"
                );


            const fvMesh& nbrFvMesh = refCast<const fvMesh>(nbrMesh);

            const fvPatch& nbrPatch =
                nbrFvMesh.boundary()[mpp.samplePolyPatch().index()];


            scalarField emissivity
            (
                radiation.absorptionEmission().e()().boundaryField()
                [
                    nbrPatch.index()
                ]
            );
            mpp.distribute(emissivity);
    
            emissivity_ = emissivity;

            return emissivity;

        }
        break;

       case PYROLYSISMODELE:
        {
            // Get the coupling information from the mappedPatchBase
            const mappedPatchBase& mpp =
                refCast<const mappedPatchBase>(patch_.patch());

            const polyMesh& nbrMesh = mpp.sampleMesh();

            const fvMesh& nbrFvMesh = refCast<const fvMesh>(nbrMesh);

            const fvPatch& nbrPatch =
                nbrFvMesh.boundary()[mpp.samplePolyPatch().index()];

//            const pyroModelType& pyroModelRef = pyroModel(nbrFvMesh,nbrFvMesh.name()); 
//            scalarField emissivity
//            (
//                pyroModelRef.kappaRad()().boundaryField()
//                [
//                    nbrPatch.index()
//                ]
//            );

            volScalarField& emmVol = const_cast<volScalarField& >(nbrFvMesh.lookupObject<volScalarField>("emmBnd"));
             
            scalarField emissivity = emmVol.boundaryField()[nbrPatch.index()];

            mpp.distribute(emissivity);

            emissivity_ = emissivity;

            return emissivity;

        }
        break;

        case LOOKUPE:
        {
            // return local value
            return emissivity_;
        }
        break;

        default:
        {
            FatalErrorIn
            (
                "radiationCoupledBaseFF::emissivity(const scalarField&)"
            )   << "Unimplemented method " << method_ << endl
                << "Please set 'emissivity' to one of "
                << emissivityMethodTypeNames_.toc()
                << exit(FatalError);
        }
        break;
    }

    return scalarField(0);
}



Foam::scalarField Foam::radiationCoupledBaseFF::absorptivity()
{
    switch (absMethod_)
    {
        case EMISSIVITY:
        {    
             absorptivity_ = Foam::radiationCoupledBaseFF::emissivity();

             //return Foam::radiationCoupledBaseFF::emissivity();
             return absorptivity_;
        }
        break;

       case PYROLYSISMODELA:
        {
            // Get the coupling information from the mappedPatchBase
            const mappedPatchBase& mpp =
                refCast<const mappedPatchBase>(patch_.patch());

            const polyMesh& nbrMesh = mpp.sampleMesh();

            const fvMesh& nbrFvMesh = refCast<const fvMesh>(nbrMesh);

            const fvPatch& nbrPatch =
                nbrFvMesh.boundary()[mpp.samplePolyPatch().index()];

//            const pyroModelType& pyroModelRef = pyroModel(nbrFvMesh,nbrFvMesh.name()); 
//            scalarField absorptivity
//            (
//                pyroModelRef.kappaRad()().boundaryField()
//                [
//                    nbrPatch.index()
//               ]
//            );

            const volScalarField& absVol = nbrFvMesh.lookupObject<volScalarField>("absBnd");
             
            scalarField absorptivity = absVol.boundaryField()[nbrPatch.index()];
      
            mpp.distribute(absorptivity);

            absorptivity_ = absorptivity; 

            return absorptivity;

        }
        break;

        case LOOKUPA:
        {
            // return local value
            return absorptivity_;
        }
        break;

        case SOLIDRADIATIONA:
        {
            // Get the coupling information from the mappedPatchBase
            const mappedPatchBase& mpp =
                refCast<const mappedPatchBase>(patch_.patch());

            const polyMesh& nbrMesh = mpp.sampleMesh();

            const radiation::radiationModel& radiation =
                nbrMesh.lookupObject<radiation::radiationModel>
                (
                    "radiationProperties"
                );


            const fvMesh& nbrFvMesh = refCast<const fvMesh>(nbrMesh);

            const fvPatch& nbrPatch =
                nbrFvMesh.boundary()[mpp.samplePolyPatch().index()];


            scalarField absorptivity
            (
                radiation.absorptionEmission().a()().boundaryField()
                [
                    nbrPatch.index()
                ]
            );
            mpp.distribute(absorptivity);

            absorptivity_ = absorptivity; 

            return absorptivity;

        }
        break;

        default:
        {
            FatalErrorIn
            (
                "radiationCoupledBaseFF::absorptivity(const scalarField&)"
            )   << "Unimplemented method " << absMethod_ << endl
                << "Please set 'absorptivity' to one of "
                << absorptivityMethodTypeNames_.toc()
                << exit(FatalError);
        }
        break;

    }

    return scalarField(0);
}

Foam::scalarField Foam::radiationCoupledBaseFF::filmTransmissivity()
{
	// Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch_.patch());

	const polyMesh& nbrMesh = mpp.sampleMesh();

	const fvMesh& nbrFvMesh = refCast<const fvMesh>(nbrMesh);

    const fvPatch& nbrPatch =
        nbrFvMesh.boundary()[mpp.samplePolyPatch().index()];

    volScalarField& fTransVol = const_cast<volScalarField& >(nbrFvMesh.lookupObject<volScalarField>("filmTransmissivity"));

    scalarField transmissivity = fTransVol.boundaryField()[nbrPatch.index()];

    mpp.distribute(transmissivity);

    return transmissivity;
}



void Foam::radiationCoupledBaseFF::write(Ostream& os) const
{
    os.writeKeyword("emissivityMode") << emissivityMethodTypeNames_[method_]
        << token::END_STATEMENT << nl;
    emissivity_.writeEntry("emissivity", os);

        os.writeKeyword("absorptivityMode") << absorptivityMethodTypeNames_[absMethod_]
        << token::END_STATEMENT << nl;
    absorptivity_.writeEntry("absorptivity", os);
}


// ************************************************************************* //
