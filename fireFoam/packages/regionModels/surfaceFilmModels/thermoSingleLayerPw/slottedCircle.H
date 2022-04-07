void thermoSingleLayerPw::slottedCircle(){
    Info << "Slotted Circle Initialization\n";

    static const vector center(0,0,0);
    static const scalar circleRadius=.01;
    static const scalar slotWidth=circleRadius*.0;
    static const volVectorField& cellCentres = regionMesh().C();

    forAll(delta_,i){
        scalar x=cellCentres[i][0];
        scalar y=cellCentres[i][1];
        scalar radius=sqrt(pow(x-center.x(),2)+pow(y-center.y(),2));
        if(radius<circleRadius){
            delta_[i]=0.0003;
        }
        else{
            delta_[i]=0.0;
        }
        /*cut out slot*/
        if(y<0&&x>-slotWidth&&x<slotWidth){
            delta_[i]=0.0;
        }
    }

    //Info <<"SplashedMass " << time().value() << " " << cummulatedMass << endl;

    return;
}
