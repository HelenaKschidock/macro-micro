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
//adapted from dumux/examples/1ptracer/spatialparams.hh
#ifndef DUMUX_CELL_PROBLEM_SPATIAL_PARAMS_HH
#define DUMUX_CELL_PROBLEM_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

template<class GridGeometry, class Scalar>
class CellProblemSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar,
                             CellProblemSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar,
                                           CellProblemSpatialParams<GridGeometry, Scalar>>;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename SubControlVolume::GlobalPosition;

    using AllenCahnTypeTag = Properties::TTag::PlainAllenCahn;
    using ACSolutionVector = GetPropType<AllenCahnTypeTag, Properties::SolutionVector>;

public:
    using PermeabilityType = Scalar;

    CellProblemSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)//, phi_(gridGeometry->gridView().size(0), 0.0)
    {   
        ks_ = getParam<Scalar>("Problem.ks");
        kg_ = getParam<Scalar>("Problem.kg");  
    }

    void updatePhi(ACSolutionVector& phi){ 
        phi_ = phi;
    }

    Scalar phasefield(const Element& element,
                        const SubControlVolume& scv)
    {
        return phi_[scv.elementIndex()];
    }

    Scalar phi0delta_(const Element& element,
                        const SubControlVolume& scv)
    {
        return phasefield(element, scv)*ks_ + (1-phasefield(element, scv))*kg_;
    }
private:
    ACSolutionVector phi_;
    Scalar ks_;
    Scalar kg_;
};
} // end namespace Dumux

#endif