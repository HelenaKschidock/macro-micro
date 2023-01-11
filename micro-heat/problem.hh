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
#ifndef DUMUX_PLAINALLENCAHN_PROBLEM_HH
#define DUMUX_PLAINALLENCAHN_PROBLEM_HH

#include <dumux/common/fvproblemwithspatialparams.hh>
#include <dumux/common/timeloop.hh>
#include <dune/istl/matrix.hh>
#include <fstream>
#include <iostream>

namespace Dumux {

template <class TypeTag >
class PlainAllenCahnProblem : public FVProblemWithSpatialParams<TypeTag>
{
    using ParentType = FVProblemWithSpatialParams<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag,
          Properties::ModelTraits>::numEq()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using TimeLoopPtr = std::shared_ptr<CheckPointTimeLoop<Scalar>>;
    using Matrix = Dune::FieldMatrix<Scalar,3,3>;
    using Vector = Dune::FieldVector<Scalar,3>;

    static constexpr int phiIdx = Indices::phiIdx;

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    PlainAllenCahnProblem( std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        omega_ = getParam<Scalar>("Problem.omega");
        alpha_ = 1.0;
        if (getParam<Scalar>("Problem.UseAlpha"))
        {
            alpha_ = 1.0/omega_;
            omega_ = 1.0;
        }
        xi_ = getParam<Scalar>("Problem.xi");
        lam_= 3/getParam<std::array<int, 2>>("Grid.Cells", std::array<int, 2>{{1, 1}})[0];
    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;

        bcTypes.setAllNeumann();

        return bcTypes;
    }

    /*!
     * \brief Calculates the derivative of the double-well potential P = 8 \phi^2 (1-\phi)^2.
     */
    Scalar pPrime(Scalar phi) const
    {
        return 16.0 * phi * (1.0 - phi) * (1.0 - 2.0*phi);
    }

    NumEqVector source(const Element &element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume &scv) const
    {
        // declare source vector
        NumEqVector source;
        // extract priVars
        const auto& priVars = elemVolVars[scv].priVars();

        source[phiIdx] = -omega_ * pPrime(priVars[phiIdx]);
        //source += conservationSource(fvGeometry, elemVolVars);
        source += 4*xi_*priVars[phiIdx]*(1.0-priVars[phiIdx])
            * interfaceVelocity(element, fvGeometry, elemVolVars, scv);
        return source;
    }

    /*!
     * \brief Returns the source term to enforce conservation of phase-field function calculated
     *        by updateConservationSource
     */
    Scalar conservationSource(const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars) const
    {
        return conservationSource_;
    }

    /*!
     * \brief Returns the interfaceVelocity to use in the source term
     */
    Scalar interfaceVelocity(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const SubControlVolume& scv) const
    {
        return 0.0;
    }

    Scalar poreVolume() const
    {
        return poreVolume_;
    }

    /*!
     * \brief Computes the expensive source term to enforce conservation of phase-field function
     */
    template <class Assembler, class SolutionVector>
    void updateConservationSource(const Assembler& assembler,
                                  const SolutionVector& curSol)
    {
        const auto& gridGeometry = this->gridGeometry();
        auto elemGeometry = localView(gridGeometry);
        const auto gridVolVars = assembler.gridVariables().curGridVolVars();
        ElementVolumeVariables elemVolVars = localView(gridVolVars);
        Scalar volume = 0;
        Scalar integral = 0;
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            elemGeometry.bind(element);
            elemVolVars.bind(element, elemGeometry, curSol);
            for (const auto& scv : scvs(elemGeometry))
            {
                volume += scv.volume();
                integral += pPrime(elemVolVars[scv].priVar(phiIdx)) * scv.volume();
            }
        }
        poreVolume_ = volume;
        conservationSource_ = omega_ * integral/volume;
    }

    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;

        const static Scalar centerX = getParam<Scalar>("Problem.CenterX");
        const static Scalar centerY = getParam<Scalar>("Problem.CenterY");
        const static Scalar radius = getParam<Scalar>("Problem.Radius");
        const static Scalar factor = getParam<Scalar>("Problem.PhasefieldICScaling");
        Scalar s = std::sqrt((globalPos[0]-centerX)*(globalPos[0]-centerX)
            +(globalPos[1]-centerY)*(globalPos[1]-centerY))
            - radius;
        values[phiIdx] = 1.0/(1.0 + std::exp(-factor*s/lam_));
        return values;
    }

    Scalar getAlpha() const
    {
        return alpha_;
    }

    Scalar getOmega() const
    {
        return omega_;
    }

private:
    Scalar lam_;
    Scalar xi_;
    Scalar omega_;
    Scalar alpha_;
    Scalar poreVolume_, conservationSource_;
};

} //end namespace Dumux

#endif // DUMUX_PLAINALLENCAHN_PROBLEM_HH