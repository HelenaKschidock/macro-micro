// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \brief Dummy macro 1pni simulation which is coupled to a set of micro simulations via preCICE and the Micro Manager
 */

#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>

#include "properties.hh"

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/discretization/method.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/loadsolution.hh>

#include <map>
#include <list>

/*!
 * \brief Provides an interface for customizing error messages associated with
 *        reading in parameters.
 *
 * \param progName  The name of the program, that was tried to be started.
 * \param errorMsg  The error message that was issued by the start function.
 *                  Comprises the thing that went wrong and a general help message.
 */
void usage(const char *progName, const std::string &errorMsg)
{
    if (errorMsg.size() > 0) {
        std::string errorMessageOut = "\nUsage: ";
                    errorMessageOut += progName;
                    errorMessageOut += " [options]\n";
                    errorMessageOut += errorMsg;
                    errorMessageOut += "\n\nThe list of mandatory options for this program is:\n"
                                        "\t-TimeManager.TEnd      End of the simulation [s] \n"
                                        "\t-TimeManager.DtInitial Initial timestep size [s] \n"
                                        "\t-Grid.LowerLeft                 Lower left corner coordinates\n"
                                        "\t-Grid.UpperRight                Upper right corner coordinates\n"
                                        "\t-Grid.Cells                     Number of cells in respective coordinate directions\n";
        std::cout << errorMessageOut
                  << "\n";
    }
}

int main(int argc, char** argv)
{
    using namespace Dumux;
    using namespace Dumux::Precice; //for QuantityType 

    // define the type tag for this problem
    using TypeTag = Properties::TTag::OnePNIConductionCCTpfa; //using TypeTag = Properties::TTag::TYPETAG; : TYPETAG is a CMakeLists.txt input

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv, usage);
    Parameters::print();

    // try to create a grid (from the given grid file or the input file)
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    // the problem (initial and boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    const std::string paramGroup = GridGeometry::discMethod == DiscretizationMethods::ccmpfa ? "MpfaTest" : ""; //checks if CCMpfa, here CCTpfa: paramGroup = None
    auto problem = std::make_shared<Problem>(gridGeometry, paramGroup);

    // Initialize preCICE.Tell preCICE about:
    // - Name of solver
    // - What rank of how many ranks this instance is
    // Configure preCICE. For now the config file is hardcoded.

    std::string preciceConfigFilename = "precice-config-heat.xml";
    if (argc > 2)
        preciceConfigFilename = argv[argc - 1];

    auto &couplingInterface = Dumux::Precice::CouplingAdapter::getInstance();
    couplingInterface.announceSolver("Macro-heat", preciceConfigFilename,
                                     mpiHelper.rank(), mpiHelper.size());

    //verify that dimensions match
    const int dim = couplingInterface.getDimensions();
    std::cout << "coupling Dims = " << dim << " , leafgrid dims = " << int(leafGridView.dimension) << std::endl;
    if (dim != int(leafGridView.dimension)){
        DUNE_THROW(Dune::InvalidStateException, "Dimensions do not match");
    }

    //get mesh coordinates
    std::string meshName = "macro-mesh";
    std::vector<double> coords;  //( dim * vertexSize );
    std::vector<int> coupledElementIdxs;
    for (const auto &element : elements(leafGridView)) {
        auto fvGeometry = localView(*gridGeometry); 
        fvGeometry.bindElement(element);
        for (const auto &scv : scvs(fvGeometry)){ //only one SCV per element for CCTpfa (but 4 scvfs)
            coupledElementIdxs.push_back(scv.elementIndex());
            const auto &pos = scv.center();
            for (const auto p : pos)
                coords.push_back(p);
        }
    }

    std::cout << "coupledElementIdxs.size() = " << coupledElementIdxs.size() << std::endl;
    
    //initialize preCICE
    auto numberOfPoints = coords.size()/dim;
    const double preciceDt = couplingInterface.setMeshAndInitialize(
        meshName, numberOfPoints, coords);
    couplingInterface.createIndexMapping(coupledElementIdxs); 

    //coupling data
    std::list<std::string> readDataNames = {"k_00", "k_01", "k_10", "k_11", "porosity"};
    std::map<std::string,int> readDataIDs;
    for (auto iter = readDataNames.begin(); iter != readDataNames.end(); iter++){
        readDataIDs[iter->c_str()] = couplingInterface.announceScalarQuantity(iter->c_str());
    }
    std::string writeDataName = "temperature";
    int temperatureID = couplingInterface.announceScalarQuantity(writeDataName);

    std::vector<double> poroData;
    std::vector<double> k_00;
    std::vector<double> k_01;
    std::vector<double> k_10;
    std::vector<double> k_11;
    std::vector<double> temperatures;
    std::map<std::string, std::vector<double>> conductivityData {{"k_00", k_00}, {"k_01", k_01}, {"k_10", k_10}, {"k_11", k_11}};

    problem->updatePreciceDataIds(readDataIDs, temperatureID);    
 
    // get some time loop parameters
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    auto dt = preciceDt;

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    std::cout << "x:" << x << std::endl;                    //TODO: what is happening here? initialization with 0?
    problem->applyInitialSolution(x); 
    auto xOld = x;

    // the grid variables                           //TODO what is happening here? does this need to be modified? we only want to write the temperature, porosity, conductivity
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // intialize the vtk output module
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.write(0.0); //restart time = 0
    
    //initialize coupling data TODO
    for (const auto &elementIdx : coupledElementIdxs){
        temperatures.push_back(couplingInterface.getScalarQuantityOnFace(temperatureID, elementIdx)); 
    }
    
    couplingInterface.writeQuantityVector(temperatureID, temperatures);
    if (couplingInterface.hasToWriteInitialData()){
        couplingInterface.writeQuantityToOtherSolver(temperatureID, QuantityType::Scalar);
        couplingInterface.announceInitialDataWritten();
    }
    couplingInterface.initializeData();

    // output every vtkOutputInterval time step
    const int vtkOutputInterval = getParam<int>("Problem.OutputInterval");

    // instantiate time loop
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd);

    // the assembler with time loop for instationary problem
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    // the linear solver
    using LinearSolver = ILU0BiCGSTABBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // time loop
    std::cout << "Time Loop starts" << std::endl;
    const auto outputTimeInterval = getParam<Scalar>("TimeLoop.TOutput");
    auto tOld = timeLoop->time();
    timeLoop->start(); do
    {   // write checkpoint
        if (couplingInterface.hasToWriteIterationCheckpoint()) {
            xOld = x;
            couplingInterface.announceIterationCheckpointWritten();
        }

        //Read porosity
        couplingInterface.readQuantityFromOtherSolver(readDataIDs["porosity"], QuantityType::Scalar);
        poroData = couplingInterface.getQuantityVector(readDataIDs["porosity"]);
        
        //Read conductivity and TODO apply write into dumux TODO check
        for (auto iter = conductivityData.begin(); iter != conductivityData.end(); iter++){
            couplingInterface.readQuantityFromOtherSolver(readDataIDs[iter->first], QuantityType::Scalar);
            conductivityData[iter->first] = couplingInterface.getQuantityVector(readDataIDs[iter->first]);
        }

        //write Data to couplingInterface Faces
        for (const auto &elementIdx : coupledElementIdxs){
                couplingInterface.writeScalarQuantityOnFace(readDataIDs["porosity"], elementIdx, poroData[elementIdx]);
                std::cout << "poroData:" << poroData[elementIdx] << std::endl;
                for (auto iter = conductivityData.begin(); iter != conductivityData.end(); iter++){
                    couplingInterface.writeScalarQuantityOnFace(readDataIDs[iter->first], elementIdx, conductivityData[iter->first][elementIdx]);
                }
        }
        
        std::cout << "Solver starts" << std::endl;
        // linearize & solve
        nonLinearSolver.solve(x, *timeLoop);

        // make the new solution the old solution
        if (couplingInterface.hasToReadIterationCheckpoint()) {
            //            //Read checkpoint
            //            freeFlowVtkWriter.write(vtkTime);
            //            vtkTime += 1.;
            xOld = x;
            gridVariables->update(x);
            gridVariables->advanceTimeStep();
            couplingInterface.announceIterationCheckpointRead();
        }

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // report statistics of this time step
        timeLoop->reportTimeStep();

        //nested forloops
        //temperatures = gridVariables->curGridVolVars()[temperatureID]; //TODO
        //std::cout << gridVariables->curGridVolVars() << std::endl; //trying to output the solution
        //temperatures = x[temperatureIDx];
        couplingInterface.writeQuantityVector(temperatureID, temperatures);
        couplingInterface.writeQuantityToOtherSolver(temperatureID, QuantityType::Scalar);

        //advance precice
        const double preciceDt = couplingInterface.advance(dt);
        dt = std::min(preciceDt, nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

        // set new dt as suggested by the newton solver or by precice
        timeLoop->setTimeStepSize(dt);
        
        //TODO output every 0.1. currently does not exactly hit this.
        if (timeLoop->timeStepIndex()==0 || int(timeLoop->time()/outputTimeInterval) > int(tOld/outputTimeInterval) || timeLoop->finished())
            vtkWriter.write(timeLoop->time());
        tOld = timeLoop->time();

    } while (!timeLoop->finished() && couplingInterface.isCouplingOngoing());

    timeLoop->finalize(leafGridView.comm());

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    couplingInterface.finalize();

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
} // end main
