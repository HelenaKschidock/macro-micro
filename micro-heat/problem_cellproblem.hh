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
#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/common/fvproblemwithspatialparams.hh>

namespace Dumux {

template<class TypeTag>
class CellProblemProblem : public FVProblemWithSpatialParams<TypeTag> 
{
    using ParentType = FVProblemWithSpatialParams<TypeTag>; 
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
    using DimWorldVector = Dune::FieldVector<Scalar, GridView::dimensionworld>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using Vector = Dune::FieldVector<Scalar, 1>; 
    using Extrusion = Extrusion_t<GridGeometry>;

    enum {
        psiIdx = Indices::psiIdx
    };
    

public:
    CellProblemProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    { 
        kij_.resize(gridGeometry->numDofs());
        d0psi1_.resize(gridGeometry->numDofs());
        d1psi1_.resize(gridGeometry->numDofs());
        d0psi2_.resize(gridGeometry->numDofs());
        d1psi2_.resize(gridGeometry->numDofs());
        dPsi_.resize(gridGeometry->numDofs());
        delta_ij_.resize(gridGeometry->numDofs());
        d_0Psi_.resize(gridGeometry->numDofs());
        d_1Psi_.resize(gridGeometry->numDofs());
        d_.resize(gridGeometry->numDofs());

    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes bcTypes;

        bcTypes.setAllNeumann();

        return bcTypes;
    }

    Scalar calculateConductivityTensorComponent(int psiIdx, int derivIdx) //TODO
    {   
        std::size_t order = 4; 
        return integrateGridFunction(this->gridGeometry(), effectiveConductivityField(psiIdx, derivIdx), order);
    }

    std::vector<Scalar>& effectiveConductivityField(int psiIdx, int derivIdx){
        
        if (psiIdx == derivIdx){
            std::fill(delta_ij_.begin(), delta_ij_.end(), 1.0);
        }
        else
        {
            std::fill(delta_ij_.begin(), delta_ij_.end(), 0.0);
        }
        dPsi_ = partialDerivativePsi(psiIdx, derivIdx);

        for (int i = 0; i < dPsi_.size(); ++i)
        {
            d_[i] = this->spatialParams().phi0deltaIdx(i)*(delta_ij_[i] + dPsi_[i]);
        }
        return d_;
    }

    std::vector<Scalar>& partialDerivativePsi(int psiIdx, int derivIdx)
    {   assert((psiIdx==0)||(psiIdx ==1));
        assert((derivIdx==0)||(derivIdx == 1));
        if ((psiIdx == 0) && (derivIdx == 0))
        {
            return d0psi1_;
        }
        else if ((psiIdx == 0) && (derivIdx == 1))
        {
            return d1psi1_;
        }
        else if ((psiIdx == 1) && (derivIdx == 0))
        {
            return d0psi2_;
        }
        else
        {
            return d1psi2_;
        }
    }

    //see dumux-adapter/examples/ff-pm/flow-over-square-2d/main_ff.cc "setInterfaceVelocities"
    template<class Problem, class Assembler, class GridVariables, class SolutionVector>
    void computePsiDerivatives(const Problem &problem,
                                const Assembler& assembler,
                                const GridVariables &gridVars,
                                const SolutionVector &psi, int psiIdx)
    {   
        const auto &gridGeometry = this->gridGeometry();
        auto fvGeometry = localView(gridGeometry);
        
        const auto gridVolVars = assembler.gridVariables().curGridVolVars();
        auto elemVolVars = localView(gridVolVars);
       
        DimWorldVector cellDeriv(0.0); 

        for (const auto &element : elements(gridGeometry.gridView())) {
            
            fvGeometry.bindElement(element);
            elemVolVars.bindElement(element, fvGeometry, psi);
            Scalar scvVolume(0.0);
            
            for (const auto &scvf : scvfs(fvGeometry)) {
                if (!scvf.boundary()) 
                {
                    int k = 0; 
                    
                    // Get the inside and outside volume variables
                    const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
                    const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
                    
                    const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
                    const auto valInside = insideVolVars.priVar(k);
                    
                    const Scalar ti = computeTpfaTransmissibility(fvGeometry, scvf, insideScv, 
                                                      insideVolVars.phi0delta(problem, element, insideScv), 
                                                      insideVolVars.extrusionFactor());
                    // cf. dumux/phasefield/localresidual.hh
                    // faces might lie on the periodic boundary, requiring the matching scvf of the scv
                    // on the other side of the periodic boundary.
                    auto outsideFvGeometry = localView(gridGeometry);
                    const auto& periodicElement = gridGeometry.element(outsideScv.elementIndex());
                    outsideFvGeometry.bind(periodicElement);
                    auto outsideElemVolVars = localView(gridVolVars);
                    outsideElemVolVars.bindElement(periodicElement, outsideFvGeometry, psi);
                    
                    Scalar tij = 0.0;
                    Scalar valOutside = 0.0;
                    for (const auto& outsideScvf : scvfs(outsideFvGeometry))
                    {
                        if (outsideScvf.unitOuterNormal() * scvf.unitOuterNormal() < -1 + 1e-6)
                        {   
                            const auto& outsideVolVars = outsideElemVolVars[outsideScvf.insideScvIdx()]; //
                            
                            valOutside = outsideVolVars.priVar(k);
                            const Scalar tj = computeTpfaTransmissibility(fvGeometry, outsideScvf, outsideScv, outsideVolVars.phi0delta(problem, element, outsideScv), outsideVolVars.extrusionFactor());
                            tij = scvf.area()*(tj)/(ti + tj);
                            break;
                        }
                    }

                    cellDeriv += scvf.area()*(valOutside - valInside)*tij*scvf.unitOuterNormal();
                    scvVolume = insideScv.volume(); 
                }
                
            }
            if (scvVolume > 0.0){
                cellDeriv /= scvVolume;
            }
            const int eIdxGlobal = gridGeometry.elementMapper().index(element);
            d_0Psi_[eIdxGlobal] = cellDeriv[0];
            d_1Psi_[eIdxGlobal] = cellDeriv[1];
            
        }
        if(psiIdx==0)
        {
            d0psi1_ = d_0Psi_;
            d1psi1_ = d_1Psi_;
        }
        else 
        {
            d0psi2_ = d_0Psi_;
            d1psi2_ = d_1Psi_;
        }
    }

    //to make available to vtkOutput, conductivity has to be converted to a Field
    const std::vector<Scalar>& getKijAsField(int derivIdx, int psiIdx)
    {   
        std::vector<Scalar> kij(d0psi1_.size(), calculateConductivityTensorComponent(psiIdx, derivIdx));
        kij_ = kij;
        return kij_; 
    }

private:
    std::vector<Scalar> kij_;
    std::vector<Scalar> dPsi_;
    std::vector<Scalar> d_0Psi_;
    std::vector<Scalar> d_1Psi_;
    std::vector<Scalar> delta_ij_;
    std::vector<Scalar> d0psi1_;
    std::vector<Scalar> d1psi1_;
    std::vector<Scalar> d0psi2_;
    std::vector<Scalar> d1psi2_;
    std::vector<Scalar> d_;
};
} // end namespace Dumux

#endif
