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
 * \brief Macro simulation solving the unsteady heat equation; coupled to a set of micro simulations via preCICE and the Micro Manager
 */

#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/istl/io.hh>

#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>

#include "main_macro_heat.hh"

#include <map>

/* NOTES: 
from nutils import mesh, function, solver, export, cli
import treelog
import numpy as np
import precice
*/

int main(int argc, char** argv)
try {
    /*
    2D unsteady heat equation on a unit square.
    The material consists of a mixture of two materials, the grain and sand
    */  
    using namespace Dumux;
    using namespace Dumux::Precice; //for QuantityType 

    // initialize MPI, finalize is done automatically on exit
    const auto &mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);
    
    // parse command line arguments and input file
    Parameters::init(argc, argv);
    Parameters::print();
    //Parameters::init(params.input);
    // Elements in one direction
    int nelems = 5;
    
    // TODO topo, geom = mesh.unitsquare(nelems, 'square')

    /* TODO
    ns = function.Namespace(fallback_length=2)
    ns.x = geom
    ns.basis = topo.basis('std', degree=1)
    ns.kbasis = topo.basis('std', degree=1).vector(topo.ndims).vector(topo.ndims)
    ns.u = 'basis_n ?solu_n'
    */
    
    //Coupling quantities
    /* TODO
    ns.phi = 'basis_n ?solphi_n'
    ns.k_ij = 'kbasis_nij ?solk_n'
    */

    //initial values
    
    float phi = 0.5; 
    float k = 1.0;

    /* TODO
    ns.rhos = 1.0
    ns.rhog = 2.0
    ns.dudt = 'basis_n (?solu_n - ?solu0_n) / ?dt'
    */

    //Dirichlet BCs temperatures
    /* TODO
    ns.ubottom = 273
    ns.utop = 400
    */

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

    //define coupling meshes 
    std::string meshName = "macro-mesh";
    //TODO: confirm that not needed: mesh_id = interface.get_mesh_id(mesh_name)

    //TODO: from macro-dummy.cc: Coupling mesh (but done above already from macro-dummy.py)
    
    /*TODO Define Gauss points on entire domain as coupling mesh
    couplingsample = topo.sample('gauss', degree=2)  # mesh located at Gauss points
    vertex_ids = interface.set_mesh_vertices(mesh_id, couplingsample.eval(ns.x))

    c++:
    auto numberOfPoints = coords.size()/couplingInterface.getDimensions();
    couplingInterface.setMesh(meshName, numberOfPoints, coords);
    
    print("Number of coupling vertices = {}".format(len(vertex_ids)))
    */

    /* TODO
    sqrphi = couplingsample.integral((ns.phi - phi) ** 2)
    solphi = solver.optimize('solphi', sqrphi, droptol=1E-12)

    sqrk = couplingsample.integral(((ns.k - k * np.eye(2)) * (ns.k - k * np.eye(2))).sum([0, 1]))
    solk = solver.optimize('solk', sqrk, droptol=1E-12)
    */

    //initialize preCICE
    const double preciceDt = couplingInterface.initialize();

    //coupling data
    std::list<std::string> readDataNames = {"k_00", "k_01", "k_10", "k_11", "porosity"};
    std::map<std::string,int> readDataIDs;
    //TODO: verify that all are scalar quantities
    for (auto iter = readDataNames.begin(); iter != readDataNames.end(); iter++){
        readDataIDs[iter->c_str()] = couplingInterface.announceScalarQuantity(iter->c_str());
    }
    std::string writeDataName = "temperature";
    int temperatureID = couplingInterface.announceScalarQuantity(writeDataName); 

    std::vector<double> temperatures;
    std::vector<double> poroData;
    std::vector<double> k_00;
    std::vector<double> k_01;
    std::vector<double> k_10;
    std::vector<double> k_11;
    std::map<std::string, std::vector<double>> conductivityData {{"k_00", k_00}, {"k_01", k_01}, {"k_10", k_10}, {"k_11", k_11}};

    //Time related variables
    //TODO ns.dt = preciceDt;
    int n = 0;
    int n_checkpoint = 0;
    int t = 0;
    int t_checkpoint = 0;
    float t_out = 0.1;
    int n_out = int(t_out/preciceDt);

    //define the weak form
    //TODO res = topo.integral('((rhos phi + (1 - phi) rhog) basis_n dudt + k_ij basis_n,i u_,j) d:x' @ ns, degree=2)

    /* TODO Set Dirichlet boundary conditions
    sqr = topo.boundary['bottom'].integral('(u - ubottom)^2 d:x' @ ns, degree=2)
    sqr += topo.boundary['top'].integral('(u - utop)^2 d:x' @ ns, degree=2)
    cons = solver.optimize('solu', sqr, droptol=1e-15)
    */

    /* TODO Set domain to initial condition
    sqr = topo.integral('(u - ubottom)^2' @ ns, degree=2)
    solu0 = solver.optimize('solu', sqr)
    temperatures = couplingsample.eval('u' @ ns, solu=solu0)
    */
    couplingInterface.writeQuantityVector(temperatureID, temperatures);
    if (couplingInterface.hasToWriteInitialData()){
        couplingInterface.writeQuantityToOtherSolver(temperatureID, QuantityType::Scalar);
        couplingInterface.announceInitialDataWritten();
    }
    couplingInterface.initializeData();

    //prepare the post processing sample
    //TODO bezier = topo.sample('bezier', 2)

    /* TODO VTK output of initial state
    x, phi, u = bezier.eval(['x_i', 'phi', 'u'] @ ns, solphi=solphi, solu=solu0)
    with treelog.add(treelog.DataLog()):
        export.vtk('macro-heat-initial', bezier.tri, x, T=u)
    */
    
    //time loop
    auto dt = preciceDt;

    while (couplingInterface.isCouplingOngoing()) {
        // write checkpoint
        if (couplingInterface.hasToWriteIterationCheckpoint()) {
            //TODO solu_checkpoint = solu0
            t_checkpoint = t;
            n_checkpoint = n;
            couplingInterface.announceIterationCheckpointWritten();
        }
        //Read porosity and apply
        couplingInterface.readQuantityFromOtherSolver(readDataIDs["porosity"], QuantityType::Scalar);
        poroData = couplingInterface.getQuantityVector(readDataIDs["porosity"]);
        /* TODO 
        poro_data = interface.read_block_scalar_data(poro_id, vertex_ids) //see above
        poro_coupledata = couplingsample.asfunction(poro_data)

        sqrphi = couplingsample.integral((ns.phi - poro_coupledata) ** 2)
        solphi = solver.optimize('solphi', sqrphi, droptol=1E-12)
        */

        //Read conductivity and apply
        int i = 0;
        for (auto iter = conductivityData.begin(); iter != conductivityData.end(); iter++){
            couplingInterface.readQuantityFromOtherSolver(readDataIDs[iter->first], QuantityType::Scalar);
            conductivityData[iter->first] = couplingInterface.getQuantityVector(readDataIDs[iter->first]);
        }
        /* TODO
        k_00_c = couplingsample.asfunction(k_00)
        k_01_c = couplingsample.asfunction(k_01)
        k_10_c = couplingsample.asfunction(k_10)
        k_11_c = couplingsample.asfunction(k_11)

        k_coupledata = function.asarray([[k_00_c, k_01_c], [k_10_c, k_11_c]])
        sqrk = couplingsample.integral(((ns.k - k_coupledata) * (ns.k - k_coupledata)).sum([0, 1]))
        solk = solver.optimize('solk', sqrk, droptol=1E-12)
        */

        /* TODO solve timestep
       solu = solver.solve_linear('solu', res, constrain=cons,
                                   arguments=dict(solu0=solu0, dt=dt, solphi=solphi, solk=solk))

        temperatures = couplingsample.eval('u' @ ns, solu=solu)
        */
        
        couplingInterface.writeQuantityVector(temperatureID, temperatures);
        couplingInterface.writeQuantityToOtherSolver(temperatureID, QuantityType::Scalar);
        
        //do the coupling
        const double preciceDt = couplingInterface.advance(dt);
        dt = std::min(preciceDt, dt);

        //advance variables
        n += 1;
        t += dt;
        //TODO solu0 = solu

        if (couplingInterface.hasToReadIterationCheckpoint()){
            //TODO solu0 = solu_checkpoint
            t = t_checkpoint;
            n = n_checkpoint;
            couplingInterface.announceIterationCheckpointRead();
        }
        else { //go to next timestep
            if (n % n_out == 0){
                /* TODO
                x, phi, u = bezier.eval(['x_i', 'phi', 'u'] @ ns, solphi=solphi, solu=solu)
                with treelog.add(treelog.DataLog()):
                    export.vtk('macro-heat-' + str(n), bezier.tri, x, T=u, phi=phi)
                */
            }
        }
    }
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

    return 0; //end main
} catch (Dumux::ParameterException &e) {
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
} catch (Dune::DGFException &e) {
    std::cerr << "DGF exception thrown (" << e
              << "). Most likely, the DGF file name is wrong "
                 "or the DGF file is corrupted, "
                 "e.g. missing hash at end of file or wrong number "
                 "(dimensions) of entries."
              << " ---> Abort!" << std::endl;
    return 2;
} catch (Dune::Exception &e) {
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
} catch (...) {
    std::cerr << "Unknown exception thrown! ---> Abort!" << std::endl;
    return 4;
}