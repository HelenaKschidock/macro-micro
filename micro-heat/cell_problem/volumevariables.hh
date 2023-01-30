// -*- mode: C++; tab-width: 3; indent-tabs-mode: nil; c-basic-offset: 4 -*-
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
// adapted from <dumux/phasefield/volumevariables.hh>
#ifndef CELL_PROBLEM_VOLUME_VARIABLES_HH
#define CELL_PROBLEM_VOLUME_VARIABLES_HH

namespace Dumux {

template<class Traits>
class CellProblemVolumeVariables
{
    using Scalar = typename Traits::PrimaryVariables::value_type;

public:
    
    using PrimaryVariables = typename Traits::PrimaryVariables;
    
    using Indices = typename Traits::ModelTraits::Indices;


    template<class ElemSol, class Problem, class Element, class Scv>
    void update(const ElemSol& elemSol,
                const Problem& problem,
                const Element& element,
                const Scv& scv)
    {
        priVars_ = elemSol[scv.localDofIndex()];
        extrusionFactor_ = problem.spatialParams().extrusionFactor(element, scv, elemSol); //default is 1.0
    }

    const PrimaryVariables &priVars() const
    { return priVars_; }

    Scalar priVar(const int pvIdx) const
    { return priVars_[pvIdx]; }

    Scalar extrusionFactor() const
    { return extrusionFactor_; }

    template<class Problem, class Element, class Scv>
    Scalar phi0delta(const Problem& problem, 
                        const Element& element,
                              const Scv& scv) const
    {
        return problem.spatialParams().phi0delta(element, scv);
    }

private:
    PrimaryVariables priVars_;
    Scalar extrusionFactor_;
};

} // end namespace Dumux

#endif
