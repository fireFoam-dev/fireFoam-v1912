
if(0){
    volScalarField procId
        (
         IOobject
         (
          "procId",
          time().timeName(),
          regionMesh(),
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
         ),
         regionMesh(),
         dimensionedScalar("zero", dimless, 0.0),
         pSp_.boundaryField().types()
        );
    const volVectorField& cellCenters = regionMesh().C();
    forAll(procId,i){
        const scalar x=cellCenters[i][0];
//gambitMesh        const scalar y=cellCenters[i][1]-3.0;
        const scalar y=cellCenters[i][1];
//gambitMesh        const scalar z=cellCenters[i][2];
        const scalar z=cellCenters[i][2]-.8;
        const scalar zero=0.0;
//gambitMesh        const scalar halfBox=0.6096; //21in+3in (center of boxes in x/z dir)
        const scalar halfBox=0.377; //.25m+.254/2m (center of boxes in x/y dir)

        if(0){
            //16 procs
            scalar id=0;
            if((fabs(x)-halfBox)>0&&x>zero&&y>zero&&z>zero){
                procId[i]=id;
            }
            id++;
            if((fabs(x)-halfBox)>0&&x>zero&&y>zero&&z<zero){
                procId[i]=id;
            }
            id++;
            if((fabs(x)-halfBox)>0&&x>zero&&y<zero&&z>zero){
                procId[i]=id;
            }
            id++;
            if((fabs(x)-halfBox)>0&&x>zero&&y<zero&&z<zero){
                procId[i]=id;
            }
            id++;
            if((fabs(x)-halfBox)>0&&x<zero&&y>zero&&z>zero){
                procId[i]=id;
            }
            id++;
            if((fabs(x)-halfBox)>0&&x<zero&&y>zero&&z<zero){
                procId[i]=id;
            }
            id++;
            if((fabs(x)-halfBox)>0&&x<zero&&y<zero&&z>zero){
                procId[i]=id;
            }
            id++;
            if((fabs(x)-halfBox)>0&&x<zero&&y<zero&&z<zero){
                procId[i]=id;
            }
            id++;
            if((fabs(x)-halfBox)<0&&x>zero&&y>zero&&z>zero){
                procId[i]=id;
            }
            id++;
            if((fabs(x)-halfBox)<0&&x>zero&&y>zero&&z<zero){
                procId[i]=id;
            }
            id++;
            if((fabs(x)-halfBox)<0&&x>zero&&y<zero&&z>zero){
                procId[i]=id;
            }
            id++;
            if((fabs(x)-halfBox)<0&&x>zero&&y<zero&&z<zero){
                procId[i]=id;
            }
            id++;
            if((fabs(x)-halfBox)<0&&x<zero&&y>zero&&z>zero){
                procId[i]=id;
            }
            id++;
            if((fabs(x)-halfBox)<0&&x<zero&&y>zero&&z<zero){
                procId[i]=id;
            }
            id++;
            if((fabs(x)-halfBox)<0&&x<zero&&y<zero&&z>zero){
                procId[i]=id;
            }
            id++;
            if((fabs(x)-halfBox)<0&&x<zero&&y<zero&&z<zero){
                procId[i]=id;
            }
            id++;
        }
        else{
            //8 procs
            scalar id=0;
            if(x>zero&&y>zero&&z>zero){
                procId[i]=id;
            }
            id++;
            if(x>zero&&y>zero&&z<zero){
                        procId[i]=id;
                    }
                    id++;
                    if(x>zero&&y<zero&&z>zero){
                        procId[i]=id;
                    }
                    id++;
                    if(x>zero&&y<zero&&z<zero){
                        procId[i]=id;
                    }
                    id++;
                    if(x<zero&&y>zero&&z>zero){
                        procId[i]=id;
                    }
                    id++;
                    if(x<zero&&y>zero&&z<zero){
                        procId[i]=id;
                    }
                    id++;
                    if(x<zero&&y<zero&&z>zero){
                        procId[i]=id;
                    }
                    id++;
                    if(x<zero&&y<zero&&z<zero){
                        procId[i]=id;
                    }
                    id++;
                }
            }
            procId.write();
        }
