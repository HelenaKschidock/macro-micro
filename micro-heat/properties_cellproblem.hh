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

//adapted from dumux/examples/1ptracer/properties.hh

#ifndef DUMUX_CELL_PROBLEM_PROPERTIES_HH
#define DUMUX_CELL_PROBLEM_PROPERTIES_HH

#include <dumux/porousmediumflow/1p/model.hh>

#include <dune/grid/yaspgrid.hh>
#include <dumux/discretization/cctpfa.hh>
#include "cell_problem/localresidual.hh"
#include "problem_cellproblem.hh"
#include "spatialparams_cellproblem.hh"
#include "cell_problem/model.hh"
#include "cell_problem/volumevariables.hh"

namespace Dumux::Properties {

namespace TTag {
struct CellProblem { using InheritsFrom = std::tuple<OneP, CCTpfaModel>; };
}

template<class TypeTag>
struct Grid<TypeTag, TTag::CellProblem> { using type = Dune::YaspGrid<2>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::CellProblem> { using type = CellProblemProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::CellProblem>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = CellProblemSpatialParams<GridGeometry, Scalar>;
};

template<class PV, class MT>
struct CellProblemVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using ModelTraits = MT;
};

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::CellProblem> { using type = CellProblemModelTraits; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::CellProblem>
{   
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;

    using Traits = PhasefieldVolumeVariablesTraits<PV, MT>;
public:
    using type = CellProblemVolumeVariables<Traits>;
};

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::CellProblem> { using type = CellProblemLocalResidual<TypeTag>; };

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::CellProblem> { static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::CellProblem> { static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::CellProblem> { static constexpr bool value = true; };

} // end namespace Dumux::Properties

#endif
