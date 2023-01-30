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

// based on dumux-phasefield/examples/plainallencahn
#include <config.h>

#include <ctime>
#include <iostream>
#include <fstream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/linear/pdesolver.hh>        // for LinearPDESolver
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/discretization/method.hh>
#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/vtkoutputmodule.hh>

#include "properties_allencahn.hh"
#include "properties_cellproblem.hh"

void usage(const char *progName, const std::string &errorMsg)
{
    if (errorMsg.size() > 0) {
        std::string errorMessageOut = "\nUsage: ";
                    errorMessageOut += progName;
                    errorMessageOut += " [options]\n";
                    errorMessageOut += errorMsg;
                    errorMessageOut += "\n\nThe list of important arguments for this program is:\n"
                            "\t-TimeLoop.TEnd             end time of the simulation\n"
                            "\t-TimeLoop.DtInitial        initial time step size\n"
                            "\t-TimeLoop.MaxTimeStepSize  maximal time step size\n"
                            "\t-Grid.LowerLeft            lower left (front) corner of the domain\n"
                            "\t-Grid.UpperRight           upper right (back) corner of the domain\n"
                            "\t-Grid.Cells                grid resolution in each coordinate direction\n"
                            "\t-Problem.xi                phasefield parameter\n"
                            "\t-Problem.omega             phasefield diffusivity/surface tension parameter\n"
                            "\t-Problem.delta             regularization parameter\n"
                            "\t-Problem.OutputInterval    interval size for VTK output\n"
                            "\t-Problem.Name              base name for VTK output files\n";
        std::cout << errorMessageOut
                  << "\n";
    }
}

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using AllenCahnTypeTag = Properties::TTag::PlainAllenCahn;
    using CellProblemTypeTag = Properties::TTag::CellProblem;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // try to create a grid (from the given grid file or the input file)
    GridManager<GetPropType<AllenCahnTypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<AllenCahnTypeTag, Properties::GridGeometry>;   //TODO care that works, is same for CP
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
    gridGeometry->update(leafGridView);

    ////////////////////////////////////////////////////////////
    // Set up the Allen-Cahn Problem
    ////////////////////////////////////////////////////////////
    
    // the allen-cahn problem (initial and boundary conditions)
    using ACProblem = GetPropType<AllenCahnTypeTag, Properties::Problem>;
    auto acProblem = std::make_shared<ACProblem>(gridGeometry);

    // get some time loop parameters
    using Scalar = GetPropType<AllenCahnTypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // the solution vector
    using ACSolutionVector = GetPropType<AllenCahnTypeTag, Properties::SolutionVector>;
    auto phiPtr = std::make_shared<ACSolutionVector>();
    ACSolutionVector phi = *phiPtr;
    acProblem->applyInitialSolution(phi);
    auto phiOldPtr = std::make_shared<ACSolutionVector>();
    *phiOldPtr = phi;

    // the grid variables
    using ACGridVariables = GetPropType<AllenCahnTypeTag, Properties::GridVariables>;
    auto acGridVariables = std::make_shared<ACGridVariables>(acProblem, gridGeometry);
    acGridVariables->init(phi);

    // instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    bool enableCheckPoints = hasParam("Problem.OutputInterval");
    if (enableCheckPoints)
        timeLoop->setPeriodicCheckPoint(getParam<double>("Problem.OutputInterval"));

    // the assembler with time loop for instationary problem
    using ACAssembler = FVAssembler<AllenCahnTypeTag, DiffMethod::numeric>;
    auto acAssembler = std::make_shared<ACAssembler>(acProblem, gridGeometry, acGridVariables, timeLoop, *phiOldPtr);

    // the linear solver
    using LinearSolver = UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolver<ACAssembler, LinearSolver>;
    auto nonLinearSolver = std::make_shared<NewtonSolver>(acAssembler, linearSolver);

    //TODO: solve one step of AC to init cell problem

    ////////////////////////////////////////////////////////////
    // Set up the Cell Problem
    ////////////////////////////////////////////////////////////

    //setup the cell problem
    using CPProblem = GetPropType<CellProblemTypeTag, Properties::Problem>;
    auto cpProblem = std::make_shared<CPProblem>(gridGeometry);

    // the solution vector
    using CPSolutionVector = GetPropType<CellProblemTypeTag, Properties::SolutionVector>;
    CPSolutionVector psi(leafGridView.size(0));

    // the grid variables
    using CPGridVariables = GetPropType<CellProblemTypeTag, Properties::GridVariables>;
    auto cpGridVariables = std::make_shared<CPGridVariables>(cpProblem, gridGeometry);
    cpGridVariables->init(psi);

    // TODO the assembler with time loop for stationary problem
    // (using nonlinear solver for now, should be linear)
    using CPAssembler = FVAssembler<CellProblemTypeTag, DiffMethod::numeric>; //DiffMethod::analytic?
    auto cpAssembler = std::make_shared<CPAssembler>(cpProblem, gridGeometry, cpGridVariables);
    
    
    //LinearPDESolver<CPAssembler, LinearSolver> cpSolver(cpAssembler, linearSolver);
    
    // the non-linear solver REPLACE
    using CPNewtonSolver = Dumux::NewtonSolver<CPAssembler, LinearSolver>;
    auto cpNonLinearSolver = std::make_shared<CPNewtonSolver>(cpAssembler, linearSolver);

    ////////////////////////////////////////////////////////////
    // Set up the Output
    ////////////////////////////////////////////////////////////

    // intialize the vtk output module
    using ACIOFields = GetPropType<AllenCahnTypeTag, Properties::IOFields>;
    VtkOutputModule<ACGridVariables, ACSolutionVector> vtkWriter(*acGridVariables, phi, acProblem->name());
    ACIOFields::initOutputModule(vtkWriter); //!< Add model specific output fields
    vtkWriter.addField(acProblem->getPorosityAsField(phi), "porosity");
    vtkWriter.addField(cpProblem->getK00AsField(psi), "k00");
    vtkWriter.addField(cpProblem->getK10AsField(psi), "k10");
    vtkWriter.addField(cpProblem->getK01AsField(psi), "k01");
    vtkWriter.addField(cpProblem->getK11AsField(psi), "k11");
    vtkWriter.addField(psi, "psi");
    vtkWriter.write(0.0);

    // time loop
    timeLoop->start(); do
    {
        // linearize & solve the allen cahn problem
        nonLinearSolver->solve(phi, *timeLoop);

        //calculate porosity
        //acProblem->calculatePorosity(phi); part of output writer

        //update Phi in the cell problem
        cpProblem->spatialParams().updatePhi(phi);

        //solve the cell problem
        cpNonLinearSolver->solve(psi); 

        //calculate the conductivity tensor
        //cpProblem->calculateConductivityTensorComponent(psi, 0, 0); //etc. for other indices. part of the output writer

        // make the new solution the old solution
        *phiOldPtr = phi;
        acGridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by the newton solver
        timeLoop->setTimeStepSize(nonLinearSolver->suggestTimeStepSize(timeLoop->timeStepSize()));

        if (timeLoop->isCheckPoint() || timeLoop->finished() || !enableCheckPoints)
        {
            vtkWriter.write(timeLoop->time());
        }
        //TODO add cp output

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

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
}