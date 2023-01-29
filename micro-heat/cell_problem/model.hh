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
// adapted from <dumux/porousmediumflow/nonisothermal/model.hh>
// see instead: <dumux/porousmediumflow/1p/model.hh>
#ifndef CELL_PROBLEM_MODEL_HH
#define CELL_PROBLEM_MODEL_HH

#include <string>
#include "indices.hh"

namespace Dumux {

struct CellProblemModelTraits 
{
    //! We solve for one more equation, i.e. the energy balance
    static constexpr int numEq() { return 2; }
    static constexpr int numComponents() {return 2;}
    /*
    static constexpr bool enableAdvection() { return true; }
    static constexpr bool enableMolecularDiffusion() { return false; }
    static constexpr bool enableEnergyBalance() { return false; }
    static constexpr bool enableThermalDispersion() { return false; }
    */
};

} // end namespace Dumux

#endif