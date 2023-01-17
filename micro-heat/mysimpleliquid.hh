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
#ifndef DUMUX_SIMPLE_H2O_HH
#define DUMUX_SIMPLE_H2O_HH
 
#include <dumux/common/parameters.hh>
#include <dumux/material/idealgas.hh>
 
#include <cmath>
 
#include <dumux/material/components/base.hh>
#include <dumux/material/components/liquid.hh>
 
namespace Dumux::Components {
 
template <class Scalar>
class MySimpleLiquid
: public Components::Base<Scalar, MySimpleLiquid<Scalar> >
, public Components::Liquid<Scalar, MySimpleLiquid<Scalar> >
, public Components::Gas<Scalar, MySimpleLiquid<Scalar> >
{
    using IdealGas = Dumux::IdealGas<Scalar>;
 
public:
    static std::string name()
    { return "MySimpleLiquid"; }
 
    static constexpr bool liquidIsCompressible()
    { return false; }
 
    static constexpr bool liquidViscosityIsConstant()
    { return true; }
 
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    { 
        return 1.0; 
    }
 
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {   
        return 1.0;
    }
};
 
template <class Scalar>
struct IsAqueous<MySimpleLiquid<Scalar>> : public std::true_type {};
 
} // end namespace Dumux::Components
 
#endif