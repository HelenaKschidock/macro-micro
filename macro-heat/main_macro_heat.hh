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
 * \brief Macro simulation solving the unsteady heat equation
 */

#ifndef DUMUX_MACRO_HEAT_HH
#define DUMUX_MACRO_HEAT_HH

#include <dumux-precice/couplingadapter.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/common/numeqvector.hh>

#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/box.hh>

#include <dune/grid/yaspgrid.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/simpleh2o.hh>

#include <dumux/common/boundarytypes.hh>

namespace Dumux{ 
template<class TypeTag>
class MacroHeatProblem; //"OnePNIConductionProblem" in dumux-heat

namespace Properties 
{
//create new type tags
namespace TTag
{
struct OnePNIConduction {using InheritsFrom = std::tuple<OnePNI>;};
struct OnePNIConductionBox { using InheritsFrom = std::tuple<OnePNIConduction, BoxModel>; };
}; //end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePNIConduction> { using type = Dune::YaspGrid<2>; }; //structured parallel 2D grid

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePNIConduction> {using type = MacroHeatProblem<TypeTag>;};

//Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePNIConduction>
{   
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Dumux::Components::SimpleH2O<Scalar>>;
};
//Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePNIConduction>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    //TODO using type = OnePNISpacialParams<GridGeometry, Scalar>; //TODO from spatialparams.hh in dumux-heat
}; 
}; //end namespace Properties
/*! 
 * \brief  The 1-phase non-isothermal (1pni) macro problem
 */
template <class TypeTag>
class MacroHeatProblem : public PorousMediumFlowProblem<TypeTag>
{   using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GridGeometry::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    
    using Indices =
        typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    /*TODO add if useful 
    enum {
        // indices of the primary variables
        pressureIdx = Indices::pressureIdx,
        temperatureIdx = Indices::temperatureIdx
    };
    */

    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ThermalConductivityModel = GetPropType<TypeTag, Properties::ThermalConductivityModel>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using IapwsH2O = Components::SimpleH2O<Scalar>; //?

    //TODO something fishy about this
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    MacroHeatProblem(std::shared_ptr<const GridGeometry> gridGeometry, const std::string& paramGroup)
    : ParentType(gridGeometry, paramGroup),
      couplingInterface_(Dumux::Precice::CouplingAdapter::getInstance())
    {
        //initialize the fluid system
        FluidSystem::init();
        name_ = getParam<std::string>("Problem.Name");
        //TODO read out params
    };
    //TODO temperature calculation
    //TODO setboundaryconditions
private:
    Dumux::Precice::CouplingAdapter &couplingInterface_;
    // parameters
    std::string name_;
};
} //end namespace Dumux

#endif //DUMUX_MACRO_HEAT_HH