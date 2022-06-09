// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief A test problem for the coupled Stokes/Darcy problem (1p)
 */
#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/istl/io.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/discretization/method.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/loadsolution.hh>

#include "main_macro_dummy.hh"

#include "dumux-precice/couplingadapter.hh"

//HK
#include <map>


int main(int argc, char** argv)
{   
    """
    Dummy macro simulation which is coupled to a set of micro simulations via preCICE and the Micro Manager
    """

    // initialize MPI, finalize is done automatically on exit
    const auto &mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);
    
    // parse command line arguments and input file
    Parameters::init(argc, argv);

    int nv = 25;
    int n = 0;
    int n_checkpoint = 0;
    int t = 0;
    int t_checkpoint = 0;

    // Initialize preCICE.Tell preCICE about:
    // - Name of solver
    // - What rank of how many ranks this instance is
    // Configure preCICE. For now the config file is hardcoded.

    std::string preciceConfigFilename = "precice-config.xml";
    if (argc > 2)
        preciceConfigFilename = argv[argc - 1];

    auto &couplingInterface = Dumux::Precice::CouplingAdapter::getInstance();
    couplingInterface.announceSolver("Macro-dummy", preciceConfigFilename,
                                     mpiHelper.rank(), mpiHelper.size());

    //define coupling meshes
    std::string meshName = "macro-mesh";
    std::map<std::string, int> readDataNames = {{"micro-scalar-data", 0}, {"micro-vector-data", 1}};
    std::map<std::string, int> writeDataNames = {{"macro-scalar-data", 0},{"macro-vector-data", 1}};

    //Coupling mesh
    std::vector<double> coords; //( dim * vertexSize );
    for (int x = 0; x < nv; x++){
        for (int d = 0; d < couplingInterface.getDimensions(), d++){
            coords.push_back(x);
        }
    }

    //Define Gauss points on entire domain as coupling mesh
    auto numberOfPoints = coords.size()/dim;
    couplingInterface.setMesh(meshName, numberOfPoints, coords);

    std::map<std::string, int> readDataIDs;
    auto iter == readDataNames.begin();
    while (iter != readDataNames.end()){
        readDataIDs[iter->first] = couplingInterface.announceQuantity(iter->first) //getDataID
        ++iter;
    }

    std::map<std::string, int> writeDataIDs;
    for (auto iter = writeDataNames.begin(); iter != writeDataNames.end(); iter++){
        writeDataIDs[iter->first] = couplingInterface.announceQuantity(iter->first); //getDataID
        ++iter;
    }

    //initialize preCICE
    const double preciceDt = couplingInterface.initialize();

    std::vector<double> writeScalarData;
    std::vector<double> writeVectorData; 

    for (int i = 0; i < nv, i++){
        writeScalarData.push_back(i);
        for (int d = 0; d < couplingInterface.getDimensions(); d++){
            writeVectorData.pushback(i);
        }
    }

    if couplingInterface.hasToWriteInitialData(){
        for (auto iter = writeDataNames.begin(); iter != writeDataNames.end(); iter++){
            if (iter->second == 0){
                couplingInterface.writeScalarQuantityToOtherSolver();
            }
            else if (iter->second == 1){
                couplingInterface.writeScalarQuantityToOtherSolver(); //Scalar meh.. but dealing only in scalars rn
            }
        }
        couplingInterface.announceInitialDataWritten();
    }
    couplingInterface.initializeData();

    //time loop
    while (couplingInterface.isCouplingOngoing()) {
        // write checkpoint
        if (couplingInterface.hasToWriteIterationCheckpoint()) {
            std::cout << "Saving macro state.";
            t_checkpoint = t;
            n_checkpoint = n;
            couplingInterface.announceIterationCheckpointWritten();
        }
        // Read porosity and apply
        for (auto iter = readDataNames.begin(); iter != readDataNames.end(); iter++){
            if (iter->second == 0){
                couplingInterface.readScalarQuantityFromOtherSolver(iter->first);
            }
            else if (iter->second == 1){ 
                couplingInterface.readScalarQuantityFromOtherSolver(iter->first); //again: meh
            }
        }
    }

    
    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
} // end main