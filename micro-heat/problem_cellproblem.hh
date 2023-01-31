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

// adapted from dumux/examples/1ptracer/problem_1p.hh
#ifndef DUMUX_CELL_PROBLEM_HH
#define DUMUX_CELL_PROBLEM_HH

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/common/boundarytypes.hh>

namespace Dumux {

template<class TypeTag>
class CellProblemProblem : public PorousMediumFlowProblem<TypeTag>
{
    // A few convenience aliases used throughout this class.
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using BoundaryTypes = Dumux::BoundaryTypes<PrimaryVariables::size()>;

    static constexpr int dimWorld = GridView::dimensionworld;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

public:
    CellProblemProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {}

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes bcTypes;

        bcTypes.setAllNeumann();

        return bcTypes;
    }

    Scalar calculateConductivityTensorComponent(SolutionVector &psi1, SolutionVector &psi2, int iIdx, int jIdx) const //TODO
    {   
        return 0.0;
    }

    //to make available to vtkOutput, porosity has to be converted to a Field
    const std::vector<Scalar>& getK00AsField(SolutionVector &psi1, SolutionVector &psi2)
    {   
        std::vector<Scalar> k00(psi1.size(), calculateConductivityTensorComponent(psi1, psi2, 0, 0));
        k00_ = k00;
        return k00_; 
    }

    //to make available to vtkOutput, porosity has to be converted to a Field
    const std::vector<Scalar>& getK10AsField(SolutionVector &psi1, SolutionVector &psi2)
    {   
        std::vector<Scalar> k10(psi1.size(), calculateConductivityTensorComponent(psi1, psi2, 1, 0));
        k10_ = k10;
        return k10_; 
    }

    //to make available to vtkOutput, porosity has to be converted to a Field
    const std::vector<Scalar>& getK01AsField(SolutionVector &psi1, SolutionVector &psi2)
    {   
        std::vector<Scalar> k01(psi1.size(), calculateConductivityTensorComponent(psi1, psi2, 0, 1));
        k01_ = k01;
        return k01_; 
    }

    //to make available to vtkOutput, porosity has to be converted to a Field
    const std::vector<Scalar>& getK11AsField(SolutionVector &psi1, SolutionVector &psi2)
    {   
        std::vector<Scalar> k11(psi1.size(), calculateConductivityTensorComponent(psi1, psi2, 1, 1));
        k11_ = k11;
        return k11_; 
    }
private:
    // components of the effective upscaled conductivity matrix K
    std::vector<Scalar> k00_;
    std::vector<Scalar> k10_;
    std::vector<Scalar> k01_;
    std::vector<Scalar> k11_;
    //field of components of K before integration
    SolutionVector kij_; //


};
} // end namespace Dumux
// [[/codeblock]]
// [[/content]]
#endif
