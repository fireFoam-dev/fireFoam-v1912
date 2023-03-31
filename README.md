# FireFOAM
This is the public-facing release of the Fire Modeling team at FM Global's customized version of the FireFOAM solver for CFD simulations of fires. It is based on the ESI/OpenFOAM-v1912 implementation. The key external dependencies are [ESI/OpenFOAM](https://develop.openfoam.com/Development/openfoam.git)-v1912 and [ESI/ThirdParty](https://sourceforge.net/projects/openfoam/files/)-v1912.

# Installation instructions
We will install ESI/ThirdParty, ESI/OpenFOAM and FM-public/FireFOAM in a case-sensitive drive on a Linux machine at a location `/path/to/FM-public-FF/`.

1. Create the containing directory and enter it.
```
mkdir -p /path/to/FM-public-FF
cd /path/to/FM-public-FF
```

2. Copy the setup script to this location and make it executable. Do not clone the repository yet.
```
wget https://raw.githubusercontent.com/fireFoam-dev/fireFoam-v1912/main/fireFoam_setup.bash
chmod +x fireFoam_setup.bash
```

3. Check the top of the file and adjust options as desired. For this project, the defaults should suffice. If debugging capability is desired (at the expense of computational speed), change the `BUILD_TYPE` setting to `Debug`.

4. Run the installation script. This can take up to three hours to complete. The `&` puts the job in the background. Foreground and background processes can be toggled using the commands `fg` and `bg`.
```
./fireFoam_setup.bash > install_log.txt &
```
- The downloading of ThirdParty from SourceForge sometimes fails. In this case, simply repeat step 4. 
- Occasionally there may be an issue midway through the installation process. The script can be re-run without repeating successfully completed sections by specifying the `START_STEP` parameter based on the steps in the `fireFoam_setup.bash` script. See comments in the script for details.

After following these steps, the directory should look as follows:
```
/path/to/FM-public-FF/
  fireFoam-v1912/
  OpenFOAM-v1912/
  ThirdParty-v1912/
  fireFoam_setup.bash
  install_log.txt
  wget-log
```

# Tutorials

## Running instructions
1. Navigate to the directory where you want to run the tutorials. It should not be inside `/path/to/FM-public-FF/`.
```
mkdir -p /path/to/cases/
cd /path/to/cases/
```

2. Copy a case to this location, e.g.,
```
cp -r /path/to/FM-public-FF/fireFoam-v1912/tutorials/poolfireMcCaffrey .
cd poolfireMcCaffrey
```

3. Prepare the environment.
```
source /path/to/FM-public-FF/OpenFOAM-v1912/etc/bashrc
```

4. Set the number of cores to be used in the simulation in `system/decomposeParDict`. For cases with film and pyrolysis regions, the same must be done for each region, by setting the number of cores in `system/filmRegion/decomposeParDict` and `system/fuelRegion/decomposeParDict`.

5. Prepare the case.
```
./mesh.sh
```

6. Either run the job locally, or submit the job to a cluster. For the latter, consider the example submission scripts named `run.sh`.

7. Use other scripts, e.g., `postproc.sh` and `reconstruct.sh`, to perform initial analysis of the data.

## Case: poolfireMcCaffrey
58kW methane gas pool fire from McCaffrey's experiments. Simuation described in:

[1] https://github.com/MaCFP/macfp-db/tree/master/Gaseous_Pool_Fires/McCaffrey_Flames/Computational_Results/2017/FMGlobal

[2] Y. Wang, P. Chatterjee, J.L. de Ris, Large Eddy Simuulation of Fire Plumes, Proc. Combust. Inst. (2011)
    DOI: 10.1016/j.proci.2010.07.031

## Case: burningBoxSuppression
Description:


    Two cardboard boxes floating above ground

    Each box is 20 inch cube

    Burner centered in ground and beneath box, ~110kW

    Injector centered above box, 0.1 kg/s

    Order of events:

    * START
    * 0  s    : burner on
    * 90 s    : injector on
    * END

Configuration: 

        
             o          <------ Sprinkler
            /|\
           / | \        <------ Spray
          /  |  \

          -------
         |       |
         |       |      <------ Cardboard box
         |       |
          -------



          -------
         |       |
         |       |      <------ Cardboard box
         |       |
        ^ ------- ^
        ^^^^^^^^^^^
           ^^^^^
            ^^^
            ^^^         <------ Burner flame
            ^^^
            ^^^
        ====^^^====     <------ Ground
        
             ^
             |
             |
           Burner  

   
Purpose of case:

    Demonstrate coupling between:

        - Pyrolysis
        - Combustion
        - Radiation
        - Spray
        - Film

# Issues
Issues with installation, tutorials, code behavior and documentation can be opened, tagged and tracked via "Issues" tab at the top of the page.
