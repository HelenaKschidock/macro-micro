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
 * \brief Relation for the effective thermal conductivity.
 */
#ifndef DUMUX_LECTURE_FUELCELL_THERMALCONDUCTIVITY_CONSTANT_HH
#define DUMUX_LECTURE_FUELCELL_THERMALCONDUCTIVITY_CONSTANT_HH

namespace Dumux {

/*!
 * \brief A constant dummy effective thermal conductivity law
 */
template<class Scalar>
class UpscaledConductivity
{
public:
    //! effective thermal conductivity \f$[W/(m K)]\f$
    template<class VolumeVariables>
    static Scalar effectiveThermalConductivity(const VolumeVariables& volVars)
    { return volVars.solidThermalConductivity(); }
};

} // end namespace Dumumx
#endif
