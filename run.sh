#!/bin/bash
# run.sh ULTRA DELUXE SUPREME EXTREME UPGRADED with Allpix-squared simulation and multi-threaded simulation

# Function for single-thread macro creation in multi-thread simulation
function writeSingleMacro() {
    # Remember to add the path before the macro file name
    # Also, all path names should be followed by "/"
    macroFile=$1
    nThread=$2
    sourcePos=("$3")
    sourceAxisDir=("$4")
    outputPath=$5
    nEvents=$6
    limRegion=$7
    checkUntrig=$8
    # Write macro file
    {
        echo "# Multi-threaded simulation macro for TPX3 simulation";
        echo "# Thread ${nThread}";
        echo "";
        echo "/g4TPX/phys/addPhysics emlivermore";
        echo "";
        echo "/process/em/pixe true";
        echo "/process/em/augerCascade true";
        echo "/process/em/deexcitationIgnoreCut true";
        echo "/process/em/Polarisation true";
        echo "";
        echo "/random/setSeeds $((nThread + 1)) $((nThread + 1))";
        echo "";
        echo "/g4TPX/SensorPos 0 0 0 cm";
        echo "/g4TPX/shieldSize 10 10 10 cm";
        echo "/g4TPX/setColDiameter 25 mm";
        echo "/g4TPX/setScattererMaterial G4_Air";
        echo "/g4TPX/setScattererSize 0.01 0.01 0.01 mm";
        echo "/g4TPX/setScattererPos -22 0 0 cm";
        echo "/g4TPX/setSensorThickness 1 mm";
        echo "";
        echo "/g4TPX/loadEnergyEdgeFile ../SingleRun/source_data/calibrate_source/EnergyEdges_Co60.txt";
        echo "/g4TPX/loadEnergySpectrumFile ../SingleRun/source_data/calibrate_source/EnergySpectrum_Co60.txt";
        echo "/g4TPX/loadPolarizationDegreeFile ../SingleRun/source_data/calibrate_source/PolarizationDegree_Co60.txt";
        echo "";
        echo "/g4TPX/setSourceRadius 1 mm";
        echo "/g4TPX/setSourcePos ${sourcePos[*]} cm";
        echo "/g4TPX/setSourceAxisDir ${sourceAxisDir[*]}";
        echo "/g4TPX/setSourceType radioactive";
        if [ "$limRegion" -eq 1 ]; then
            echo "/g4TPX/setLimitRegion true";
        else
            echo "/g4TPX/setLimitRegion false";
        fi
        echo "";
        echo "/g4TPX/setAllpix false";
        echo "/g4TPX/calRecomb false";
        if [ "$checkUntrig" -eq 1 ]; then
            echo "/g4TPX/checkUntriggered true";
        else
            echo "/g4TPX/checkUntriggered false";
        fi
        echo "/g4TPX/setTriggerThreshold 4.542";
        echo "/g4TPX/use2mmParam false";
        echo "";
        echo "/g4TPX/setOutputDir ${outputPath}";
        echo "";
        echo "/control/verbose 0";
        echo "/control/saveHistory";
        echo "/run/verbose 0";
        echo "";
        echo "/run/initialize";
        echo "";
        echo "/gun/particle gamma";
        echo "/run/beamOn ${nEvents}";
    } > "${macroFile}"
}

# Function for single-thread macro creation in multi-thread poin simulation
function writeSingleMacroPion() {
    # Remember to add the path before the macro file name
    # Also, all path names should be followed by "/"
    macroFile=$1
    nThread=$2
    sourcePos=("$3")
    sourceAxisDir=("$4")
    outputPath=$5
    nEvents=$6
    # Write macro file
    {
        echo "# Multi-threaded simulation macro for TPX3 pion simulation";
        echo "# Thread ${nThread}";
        echo "";
        echo "/g4TPX/phys/addPhysics emstandard_opt0";
        echo "";
        echo "/process/em/pixe true";
        echo "/process/em/augerCascade true";
        echo "/process/em/deexcitationIgnoreCut true";
        echo "/process/em/Polarisation true";
        echo "";
        echo "/random/setSeeds $((nThread + 1)) $((nThread + 1))";
        echo "";
        echo "/g4TPX/SensorPos 0 0 0 cm";
        echo "/g4TPX/shieldSize 10 10 10 cm";
        echo "/g4TPX/setColDiameter 25 mm";
        echo "/g4TPX/setScattererMaterial G4_Air";
        echo "/g4TPX/setScattererSize 0.01 0.01 0.01 mm";
        echo "/g4TPX/setScattererPos -22 0 0 cm";
        echo "/g4TPX/setSensorThickness 2 mm";
        echo "";
        echo "/g4TPX/setSourceType pion";
        echo "/g4TPX/setPionMomentum 40";
        echo "";
        echo "/g4TPX/setSourceRadius 1 mm";
        echo "/g4TPX/setSourcePos ${sourcePos[*]} cm";
        echo "/g4TPX/setSourceAxisDir ${sourceAxisDir[*]}";
        echo "/g4TPX/setLimitRegion false";
        echo "";
        echo "/g4TPX/setAllpix false";
        echo "/g4TPX/calRecomb false";
        echo "/g4TPX/useTheoCCE true";
        echo "/g4TPX/setTriggerThreshold 2.52";
        echo "/g4TPX/use2mmParam true";
        echo "";
        echo "/g4TPX/setOutputDir ${outputPath}";
        echo "";
        echo "/control/verbose 0";
        echo "/control/saveHistory";
        echo "/run/verbose 0";
        echo "";
        echo "/run/initialize";
        echo "";
        echo "/gun/particle pi-";
        echo "/run/beamOn ${nEvents}";
    } > "${macroFile}"
}

# Options for TPX3 simulation & processing
process=0
quit=0
runSimulation=1
useallpix=0
runAllpix=1
runPion=0
specifyScript=0
script=""

# Process input options/parameters
while getopts "pmqnadt:j:e:clif:" arg
do
    case $arg in
        p)
            # Process the ssimulation result
            process=1
        ;;
        q)
            # Quit root (".q" command) after processing the result
            quit=1
        ;;
        n)
            # Don't run Geant4 simulation part (usually used to process previously-generated data)
            runSimulation=0
        ;;
        a)
            # Use Allpix for charge transport simulation
            useallpix=1
        ;;
        d)
            # Skip Allpix simulation part (usually used to process previously-generated data)
            runAllpix=0
        ;;
        i)
            # Run pion simulation
            runPion=1
        ;;
        f)
            # Run specific script
            specifyScript=1
            script=$OPTARG
        ;;
        ?)
            echo "Unknown argument $arg"
        ;;
    esac
done

cd build || exit
if [ $runPion -eq 1 ]; then
    # Run pion simulation
    if [ $runSimulation -eq 1 ]; then
        ./TPX3 ../SingleRun/macro/runPionTest.mac
    fi

    # Process the results
    if [ $process -eq 1 ]; then
        cd ../SingleRun/PostProcess/runScript || exit
        if [ $quit -eq 1 ]; then
            root -l -q -x PionProcess.C
        else
            root -l -x PionProcess.C
        fi
    fi
else
    # Run photon simulation
    # Run simulation by using Allpix for charge carrier transport simulation
    if [ $useallpix -eq 1 ]; then
        echo "WARNING: due to recently discovered defects in the drift model of Allpix-squared, simulation with Allpix-squared is not recommended!"
        # Run Geant4 simulation part
        if [ $runSimulation -eq 1 ]; then
            ./TPX3 ../SingleRun/macro/runAllpix.mac
        fi

        # Run Allpix simulation for charge carrier transport
        cd ../allpix || exit
        if [ $runAllpix -eq 1 ]; then
            python addEvents.py
            allpix -c TPX3_transport.conf -j 8
        fi
        
        # Process the results
        if [ $process -eq 1 ]; then
            cd PostProcess/runScript || exit
            if [ $quit -eq 1 ]; then
                root -l -q -x PostProcess.C
            else
                root -l -x PostProcess.C
            fi
        fi
    # Run simulation by using built-in charge carrier transport simulation
    else
        # Run Geant4 simulation part
        if [ $runSimulation -eq 1 ]; then
            if [ $specifyScript -eq 1 ]; then
                if [ -f "$script" ]; then
                    echo "Executing specified macro \"$script\""
                    ./TPX3 "$script"
                else
                    echo "Error: file \"$script\" does not exist"
                fi
            else
                echo "Executing default macro \"../SingleRun/macro/runTest.mac\""
                ./TPX3 ../SingleRun/macro/runTest.mac
            fi
        fi

        # Process the results
        if [ $process -eq 1 ]; then
            cd ../SingleRun/PostProcess/runScript || exit
            if [ $quit -eq 1 ]; then
                root -l -q -x PostProcess.C
            else
                root -l -x PostProcess.C
            fi
        fi
    fi
fi
