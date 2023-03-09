#include <iostream>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h> // numpy arrays
#include <pybind11/stl.h> // std::vector conversion
#include <numbers> // std::numbers

#include <config.h>

#include <ctime>
#include <iostream>
#include <fstream>

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

namespace py = pybind11; 

class MicroSimulation
{   
    using AllenCahnTypeTag = Dumux::Properties::TTag::PlainAllenCahn;
    using CellProblemTypeTag = Dumux::Properties::TTag::CellProblem;
    using ACSolutionVector = Dumux::GetPropType<AllenCahnTypeTag, Dumux::Properties::SolutionVector>;
    using CPSolutionVector = Dumux::GetPropType<CellProblemTypeTag, Dumux::Properties::SolutionVector>;
    using ACProblem = Dumux::GetPropType<AllenCahnTypeTag, Dumux::Properties::Problem>;
    using CPProblem = Dumux::GetPropType<CellProblemTypeTag, Dumux::Properties::Problem>;
    using ACGridVariables = Dumux::GetPropType<AllenCahnTypeTag, Dumux::Properties::GridVariables>;
    using CPGridVariables = Dumux::GetPropType<CellProblemTypeTag, Dumux::Properties::GridVariables>;
    using ACAssembler = Dumux::FVAssembler<AllenCahnTypeTag, Dumux::DiffMethod::numeric>;
    using CPAssembler = Dumux::FVAssembler<CellProblemTypeTag, Dumux::DiffMethod::numeric>; 
    using LinearSolver = Dumux::UMFPackBackend;
    using CPLinearSolver = Dumux::ILUnBiCGSTABBackend;
    using ACNewtonSolver = Dumux::NewtonSolver<ACAssembler, LinearSolver>;
    using GridGeometry = Dumux::GetPropType<AllenCahnTypeTag, Dumux::Properties::GridGeometry>;   
    using Scalar = Dumux::GetPropType<AllenCahnTypeTag, Dumux::Properties::Scalar>;
    using GridManager = Dumux::GridManager<Dumux::GetPropType<AllenCahnTypeTag, Dumux::Properties::Grid>>;
    using JacobianMatrix = Dumux::GetPropType<CellProblemTypeTag, Dumux::Properties::JacobianMatrix>;
    using SolutionVector = Dumux::GetPropType<CellProblemTypeTag, Dumux::Properties::SolutionVector>;

public:
    MicroSimulation(int sim_id);
    void initialize();
    // solve takes python dict for macro_write data, dt, and returns python dict for macro_read data
    py::dict solve(py::dict macro_write_data, double dt);
    void save_checkpoint();
    void reload_checkpoint();

private:
    const double pi_ = 3.14159265358979323846;
    int _sim_id;
    double _k_00;
    double _k_01;
    double _k_10;
    double _k_11;
    double _porosity;

    std::shared_ptr<ACNewtonSolver> _acNonLinearSolver;
    std::shared_ptr<LinearSolver> _linearSolver;
    std::shared_ptr<CPLinearSolver> _cpLinearSolver;
    std::shared_ptr<ACAssembler> _acAssembler;
    std::shared_ptr<CPAssembler> _cpAssembler;
    std::shared_ptr<Dumux::CheckPointTimeLoop<double> > _timeLoop;
    std::shared_ptr<ACProblem> _acProblem;
    std::shared_ptr<CPProblem> _cpProblem;
    std::shared_ptr<CPGridVariables> _cpGridVariables;
    std::shared_ptr<ACGridVariables> _acGridVariables;
    std::shared_ptr<GridGeometry> _gridGeometry;
    ACSolutionVector _phi;
    ACSolutionVector _phiOld;
    CPSolutionVector _psi1;
    CPSolutionVector _psi2;
    GridManager _gridManager;
    std::shared_ptr<JacobianMatrix> _A;
    std::shared_ptr<CPSolutionVector> _r;

};

// Constructor
MicroSimulation::MicroSimulation(int sim_id) : _sim_id(sim_id), _k_00(0), _k_01(0), _k_10(0),_k_11(0),_porosity(0), _phiOld(0) {};

// Initialize
void MicroSimulation::initialize()
{   
    using namespace Dumux;

    std::cout << "Initialize micro problem (" << _sim_id << ")\n";
    _k_00 = 0;
    _k_01 = 0;
    _k_10 = 0;
    _k_11 = 0;

    // parse command line arguments and input file
    Parameters::init("params.input");//argc, argv); TODO 

    // try to create a grid (from the given grid file or the input file)
    _gridManager.init();

    // we compute on the leaf grid view
    const auto& leafGridView = _gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    _gridGeometry = std::make_shared<GridGeometry>(leafGridView);
    _gridGeometry->update(leafGridView);

    ////////////////////////////////////////////////////////////
    // Set up the Allen-Cahn Problem
    ////////////////////////////////////////////////////////////
    
    // the allen-cahn problem (initial and boundary conditions)
    _acProblem = std::make_shared<ACProblem>(_gridGeometry);

    // get some time loop parameters
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // the solution vector
    auto phiPtr = std::make_shared<ACSolutionVector>();
    _phi = *phiPtr;
    _acProblem->applyInitialSolution(_phi);
    _phiOld = _phi;

    // the grid variables
    _acGridVariables = std::make_shared<ACGridVariables>(_acProblem, _gridGeometry);
    _acGridVariables->init(_phi);

    // instantiate time loop
    _timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0.0, dt, tEnd);
    _timeLoop->setMaxTimeStepSize(maxDt);

    // the assembler with time loop for instationary problem
    _acAssembler = std::make_shared<ACAssembler>(_acProblem, _gridGeometry, _acGridVariables, _timeLoop, _phiOld);

    // the linear solver
    _linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    _acNonLinearSolver = std::make_shared<ACNewtonSolver>(_acAssembler, _linearSolver);

    ////////////////////////////////////////////////////////////
    // Set up the Cell Problem
    ////////////////////////////////////////////////////////////

    //setup the cell problem
    _cpProblem = std::make_shared<CPProblem>(_gridGeometry);

    // the solution vectors
    CPSolutionVector psi1(leafGridView.size(0));
    _psi1 = psi1;
    CPSolutionVector psi2(leafGridView.size(0));
    _psi2 = psi2;

    // The jacobian matrix (`A`), the solution vector (`psi`) and the residual (`r`) make up the linear system.
    _A = std::make_shared<JacobianMatrix>();
    _r = std::make_shared<SolutionVector>();
    _cpLinearSolver = std::make_shared<CPLinearSolver>();

    // the grid variables
    _cpGridVariables = std::make_shared<CPGridVariables>(_cpProblem, _gridGeometry);
    _cpGridVariables->init(_psi1);

    _cpAssembler = std::make_shared<CPAssembler>(_cpProblem, _gridGeometry, _cpGridVariables);
    //LinearPDESolver<CPAssembler, LinearSolver> cpSolver(cpAssembler, linearSolver);
    _cpAssembler->setLinearSystem(_A,_r);

    // intialize the vtk output module
    _timeLoop->start();
}

// Solve
py::dict MicroSimulation::solve(py::dict macro_write_data, double dt)
{   
    // set leafgridView
    // TODO: find way to set leafGridView to shared private variable (then redundant)
    // we compute on the leaf grid view
    const auto& leafGridView = _gridManager.grid().leafGridView();
    _gridGeometry->update(leafGridView);

    std::cout << "Solve timestep of micro problem (" << _sim_id << ")\n";

    // assert(dt != 0);
    if (dt == 0)
    {
        std::cout << "dt is zero\n";
        exit(1);
    }
    
    _timeLoop->setTimeStepSize(dt);

    //read concentration from preCICE
    double conc = macro_write_data["concentration"].cast<double>();

    //input macro concentration into allen-cahn problem
    _acProblem->updateConcentration(conc);

    // linearize & solve the allen cahn problem
    _acNonLinearSolver->solve(_phi, *_timeLoop);

    //calculate porosity 
    _porosity = _acProblem->calculatePorosity(_phi);

    //update Phi in the cell problem
    _cpProblem->spatialParams().updatePhi(_phi);

    //solve the cell problems 
    std::cout << "Solve Cell Problem 1" << std::endl;
    int psiIdx = 0;
    _cpProblem->spatialParams().updatePsiIndex(psiIdx);
    _cpGridVariables->update(_psi1);

    // We solve the linear system `A psi = r`.
    _cpAssembler->assembleJacobianAndResidual(_psi1); 
    _cpLinearSolver->solve(*_A, _psi1, *_r);

    std::cout << "Solve Psi Derivative" << std::endl;
    _cpProblem->computePsiDerivatives(*_cpProblem, *_cpAssembler, *_cpGridVariables, _psi1, psiIdx);

    std::cout << "Solve Cell Problem 2" << std::endl;
    psiIdx = 1;
    _cpProblem->spatialParams().updatePsiIndex(psiIdx);
    _cpGridVariables->update(_psi2);

    _cpAssembler->assembleJacobianAndResidual(_psi2); 
    _cpLinearSolver->solve(*_A, _psi2, *_r);

    std::cout << "Solve Psi Derivative" << std::endl;
    _cpProblem->computePsiDerivatives(*_cpProblem, *_cpAssembler, *_cpGridVariables, _psi2, psiIdx);

    //calculate the conductivity tensor 
    _k_00 = _cpProblem->calculateConductivityTensorComponent(0,0);
    _k_10 = _cpProblem->calculateConductivityTensorComponent(1,0);
    _k_01 = _cpProblem->calculateConductivityTensorComponent(0,1);
    _k_11 = _cpProblem->calculateConductivityTensorComponent(1,1);

    // create python dict for micro_write_data
    py::dict micro_write_data;

    // add micro_scalar_data and micro_vector_data to micro_write_data
    micro_write_data["k_00"] = _k_00;
    micro_write_data["k_10"] = _k_10;
    micro_write_data["k_01"] = _k_01;
    micro_write_data["k_11"] = _k_11;
    micro_write_data["porosity"] = _porosity;
    micro_write_data["grain_size"] = std::sqrt((1-_porosity)/pi_);
    
    // return micro_write_data
    return micro_write_data;
}
// Save Checkpoint
void MicroSimulation::save_checkpoint()
{
    std::cout << "Saving state of micro problem (" << _sim_id << ")\n";
    _phiOld = _phi; 
}

// Reload Checkpoint
void MicroSimulation::reload_checkpoint()
{   
    std::cout << "Reverting to old state of micro problem (" << _sim_id << ")\n";
    _phi = _phiOld;
    _acGridVariables->update(_phi);
    _acGridVariables->advanceTimeStep();
}

PYBIND11_MODULE(micro_sim, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    py::class_<MicroSimulation>(m, "MicroSimulation")
        .def(py::init<int>())
        .def("initialize", &MicroSimulation::initialize)
        .def("solve", &MicroSimulation::solve)
        .def("save_checkpoint", &MicroSimulation::save_checkpoint)
        .def("reload_checkpoint", &MicroSimulation::reload_checkpoint);
}