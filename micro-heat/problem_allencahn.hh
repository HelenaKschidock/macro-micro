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
#include <dumux/common/integrate.hh>

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
    PlainAllenCahnProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {   
        omega_ = getParam<Scalar>("Problem.omega");
        alpha_ = 1.0;
        xi_ = getParam<Scalar>("Problem.xi");
        kt_ = getParam<Scalar>("Problem.kt");
        eqconc_ = getParam<Scalar>("Problem.eqconc");
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
        source += 4*xi_*priVars[phiIdx]*(1.0-priVars[phiIdx]) * reactionRate();
        return source;
    }

    void updateConcentration(Scalar conc)
    {
        conc_ = conc;
    }
    /*!
     * \brief Returns the interfaceVelocity to use in the source term (F(T))
     */
    Scalar reactionRate() const
    {   
        return - kt_ * ((concentration()/eqconc_)*(concentration()/eqconc_) - 1) ;
    }

    /*!
     * \brief Returns the macro temperature/concentration (TODO)
     */
    Scalar concentration() const
    {
        return conc_;
    }

    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;

        const static Scalar centerX = getParam<Scalar>("Grid.CenterX");
        const static Scalar centerY = getParam<Scalar>("Grid.CenterY");
        const static Scalar radius = getParam<Scalar>("Problem.Radius");
        const static Scalar factor = getParam<Scalar>("Problem.PhasefieldICScaling");
        Scalar s = std::sqrt((globalPos[0]-centerX)*(globalPos[0]-centerX)
            +(globalPos[1]-centerY)*(globalPos[1]-centerY))
            - radius;
        values[phiIdx] = 1.0/(1.0 + std::exp(-factor*s/xi_));
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

    Scalar calculatePorosity(SolutionVector &sol) const //todo check
    {   
        std::size_t order = 2; 
        return integrateGridFunction(this->gridGeometry(), sol, order);
    }

    //to make available to vtkOutput, porosity has to be converted to a Field
    const std::vector<Scalar>& getPorosityAsField(SolutionVector &sol)
    {   std::vector<Scalar> poro(sol.size(), calculatePorosity(sol));
        poro_ = poro; 
        return poro_; 
    }

private:
    Scalar xi_;
    Scalar omega_;
    Scalar alpha_;
    Scalar kt_;
    Scalar eqconc_;
    std::vector<Scalar> poro_;
    Scalar conc_;
};

} //end namespace Dumux

#endif // DUMUX_PLAINALLENCAHN_PROBLEM_HH