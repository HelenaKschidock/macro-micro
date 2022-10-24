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
#include <dumux/material/components/gas.hh>
 
namespace Dumux::Components {
 
template <class Scalar>
class MySimpleH2O
: public Components::Base<Scalar, MySimpleH2O<Scalar> >
, public Components::Liquid<Scalar, MySimpleH2O<Scalar> >
, public Components::Gas<Scalar, MySimpleH2O<Scalar> >
{
    using IdealGas = Dumux::IdealGas<Scalar>;
 
public:
    static std::string name()
    { return "MySimpleH2O"; }
 
    static constexpr Scalar molarMass()
    { return 18e-3; }
 
    static Scalar criticalTemperature()
    { return 647.096; /* [K] */ }
 
    static Scalar criticalPressure()
    { return 22.064e6; /* [N/m^2] */ }
 
    static Scalar tripleTemperature()
    { return 273.16; /* [K] */ }
 
    static Scalar triplePressure()
    { return 611.657; /* [N/m^2] */ }
 
    static Scalar vaporPressure(Scalar T)
    {
        if (T > criticalTemperature())
            return criticalPressure();
        if (T < tripleTemperature())
            return 0; // water is solid: We don't take sublimation into account
 
        constexpr Scalar n[10] = {
            0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2,
            0.12020824702470e5, -0.32325550322333e7, 0.14915108613530e2,
            -0.48232657361591e4, 0.40511340542057e6, -0.23855557567849,
            0.65017534844798e3
        };
 
        const Scalar sigma = T + n[8]/(T - n[9]);
 
        const Scalar A = (sigma + n[0])*sigma + n[1];
        const Scalar B = (n[2]*sigma + n[3])*sigma + n[4];
        const Scalar C = (n[5]*sigma + n[6])*sigma + n[7];
 
        using std::sqrt;
        const Scalar term = 2.0*C/(sqrt(B*B - 4.0*A*C) - B);
 
        return 1e6*term*term*term*term;
    }
 
    static const Scalar gasEnthalpy(Scalar temperature,
                                    Scalar pressure)
    {
        static const Scalar tRef = getParam<Scalar>("SimpleH2O.ReferenceTemperature", 293.15);
        return gasHeatCapacity(temperature, pressure)*(temperature - tRef) + vaporizationEnthalpy();
    }
 
    static const Scalar liquidEnthalpy(Scalar temperature,
                                       Scalar pressure)
    {
        static const Scalar tRef = getParam<Scalar>("SimpleH2O.ReferenceTemperature", 293.15);
        return liquidHeatCapacity(temperature, pressure)*(temperature - tRef)
                + pressure/liquidDensity(temperature, pressure);
    }
 
    static Scalar vaporizationEnthalpy()
    {
        constexpr Scalar A = 2500.304;
        constexpr Scalar B = -2.2521025;
        constexpr Scalar C = -0.021465847;
        constexpr Scalar D = 3.1750136e-4 ;
        constexpr Scalar E = -2.8607959e-5;
 
        //tRef in °C
        static const Scalar tRef = getParam<Scalar>("SimpleH2O.ReferenceTemperature", 293.15) - 273.15;
 
        using std::pow;
        static const Scalar vaporizationEnthalpy = A + B*tRef + C*(pow(tRef, 1.5)) + D*(pow(tRef, 2.5)) + E*(pow(tRef, 3));
        return vaporizationEnthalpy;
    }
 
 
    static const Scalar gasInternalEnergy(Scalar temperature,
                                          Scalar pressure)
    {
        // 1/molarMass: conversion from [J/(mol K)] to [J/(kg K)]
        // R*T/molarMass: pressure *spec. volume for an ideal gas
        return gasEnthalpy(temperature, pressure)
                - 1.0/molarMass()*IdealGas::R*temperature;
    }
 
    static const Scalar liquidInternalEnergy(Scalar temperature,
                                             Scalar pressure)
    {
        return liquidEnthalpy(temperature, pressure)
                - pressure/liquidDensity(temperature, pressure);
    }
 
    static constexpr bool gasIsCompressible()
    { return true; }
 
    static constexpr bool liquidIsCompressible()
    { return false; }
 
    static constexpr bool gasViscosityIsConstant()
    { return true; }
 
    static constexpr bool liquidViscosityIsConstant()
    { return true; }
 
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        // Assume an ideal gas
        return IdealGas::density(molarMass(), temperature, pressure);
    }
 
    static Scalar gasMolarDensity(Scalar temperature, Scalar pressure)
    { return IdealGas::molarDensity(temperature, pressure); }
 
    static constexpr bool gasIsIdeal()
    { return true; }
 
    static Scalar gasPressure(Scalar temperature, Scalar density)
    {
        // Assume an ideal gas
        return IdealGas::pressure(temperature, density/molarMass());
    }
 
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        return getParam<Scalar>("Component.LiquidDensity");
    }
 
    static Scalar liquidMolarDensity(Scalar temperature, Scalar pressure)
    { return liquidDensity(temperature, pressure)/molarMass(); }
 
    static Scalar liquidPressure(Scalar temperature, Scalar density)
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The liquid pressure is undefined for incompressible fluids");
    }
 
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    {
        return 1e-05;
    }
 
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        return 1e-03;
    }
 
    static Scalar liquidHeatCapacity(Scalar temperature, Scalar pressure)
    {
        return getParam<Scalar>("Component.LiquidHeatCapacity");
    }
 
    static Scalar liquidThermalConductivity(Scalar temperature, Scalar pressure)
    {
       return getParam<Scalar>("Component.LiquidThermalConductivity");
    }
 
    static Scalar gasThermalConductivity(Scalar temperature, Scalar pressure)
    {
       return 0.025;
    }
 
    static Scalar gasHeatCapacity(Scalar temperature, Scalar pressure)
    {
        return 2.08e3;
    }
 
};
 
template <class Scalar>
struct IsAqueous<MySimpleH2O<Scalar>> : public std::true_type {};
 
} // end namespace Dumux::Components
 
#endif