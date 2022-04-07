{
    static volScalarField rhos
        (
         IOobject
         (
          "rhos",
          time().timeName(),
          regionMesh(),
          IOobject::NO_READ,
          IOobject::NO_WRITE
         ),
         regionMesh(),
         lc.density,
        zeroGradientFvPatchScalarField::typeName
        );
    static volScalarField cps
        (
         IOobject
         (
          "cps",
          time().timeName(),
          regionMesh(),
          IOobject::NO_READ,
          IOobject::NO_WRITE
         ),
         regionMesh(),
         lc.cp1,
        zeroGradientFvPatchScalarField::typeName
        );
    forAll(cps,ic){
        scalar Tratio=(Tw_[ic]-lc.T1.value())/(lc.T2.value()-lc.T1.value());
        cps[ic]=Tratio*(lc.cp2.value()-lc.cp1.value())+lc.cp1.value();
    }
    static volScalarField ks
        (
         IOobject
         (
          "ks",
          time().timeName(),
          regionMesh(),
          IOobject::NO_READ,
          IOobject::NO_WRITE
         ),
         regionMesh(),
         lc.k1,
        zeroGradientFvPatchScalarField::typeName
        );
    forAll(ks,ic){
        scalar Tratio=(Tw_[ic]-lc.T1.value())/(lc.T2.value()-lc.T1.value());
        ks[ic]=Tratio*(lc.k2.value()-lc.k1.value())+lc.k1.value();
    }
    static volScalarField emissivityBack
        (
         IOobject
         (
          "emissivityBack",
          time().timeName(),
          regionMesh(),
          IOobject::NO_READ,
          IOobject::NO_WRITE
         ),
         regionMesh(),
         lc.emissivityBack,
        zeroGradientFvPatchScalarField::typeName
        );
    forAll(emissivityBack,ic){
        emissivityBack[ic]=0.7916+0.000275*Tw_[ic];
    }
    static volScalarField emissivityFront
        (
         IOobject
         (
          "emissivityFront",
          time().timeName(),
          regionMesh(),
          IOobject::NO_READ,
          IOobject::NO_WRITE
         ),
         regionMesh(),
         lc.emissivityFront,
        zeroGradientFvPatchScalarField::typeName
        );
    forAll(emissivityFront,ic){
        emissivityFront[ic]=0.7916+0.000275*Tw_[ic];
    }
    const dimensionedScalar sigma = constant::physicoChemical::sigma;
    tmp<fvScalarMatrix> tTEqn
        (
         fvm::ddt(rhos*cps, Tw_)
         - fvm::laplacian(ks, Tw_)
         == (
        		  qDotFilm_
        		+ qDotBackDry_
        		+ qDotFrontDry_
        		+ emissivityFront*qRadAugmentedDry_ // front side absorption
        		- (1.0-wettedFraction_)*emissivityFront*sigma*(pow(Tw_,4)-1.0*pow(lc.TsurrFront,4)) 	// front side emission (view factor?)
        		- emissivityBack*sigma*(pow(Tw_,4)-pow(lc.TsurrBack,4)) 	// back side emission
        	  )/lc.thickness
        );

    fvScalarMatrix& TEqn = tTEqn.ref();
    TEqn.relax();
    TEqn.solve();

}
