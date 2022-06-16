// EnergyPlus, Copyright (c) 1996-2022, The Board of Trustees of the University of Illinois,
// The Regents of the University of California, through Lawrence Berkeley National Laboratory
// (subject to receipt of any required approvals from the U.S. Dept. of Energy), Oak Ridge
// National Laboratory, managed by UT-Battelle, Alliance for Sustainable Energy, LLC, and other
// contributors. All rights reserved.
//
// NOTICE: This Software was developed under funding from the U.S. Department of Energy and the
// U.S. Government consequently retains certain rights. As such, the U.S. Government has been
// granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable,
// worldwide license in the Software to reproduce, distribute copies to the public, prepare
// derivative works, and perform publicly and display publicly, and to permit others to do so.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted
// provided that the following conditions are met:
//
// (1) Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//
// (2) Redistributions in binary form must reproduce the above copyright notice, this list of
//     conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//
// (3) Neither the name of the University of California, Lawrence Berkeley National Laboratory,
//     the University of Illinois, U.S. Dept. of Energy nor the names of its contributors may be
//     used to endorse or promote products derived from this software without specific prior
//     written permission.
//
// (4) Use of EnergyPlus(TM) Name. If Licensee (i) distributes the software in stand-alone form
//     without changes from the version obtained under this License, or (ii) Licensee makes a
//     reference solely to the software portion of its product, Licensee must refer to the
//     software as "EnergyPlus version X" software, where "X" is the version number Licensee
//     obtained under this License and may not use a different name for the software. Except as
//     specifically required in this Section (4), Licensee shall not use in a company name, a
//     product name, in advertising, publicity, or other promotional activities any name, trade
//     name, trademark, logo, or other designation of "EnergyPlus", "E+", "e+" or confusingly
//     similar designation, without the U.S. Department of Energy's prior written consent.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
// AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include<string>
#include<algorithm>

#define AFN_TO_KELVIN(T) (T + 273.15)
#define AFN_TO_CELSIUS(T) (T - 273.15)

namespace airflownetwork {

enum class Temperature : int
{
  Kelvin,
  Celsius
};

template <typename T, bool> struct select_kelvin {
  static T convert (T v) { return v + 273.15; }
};

template <typename T> struct select_kelvin<T, true> {
  static T convert (T v) { return v; }
};

template <typename T, bool> struct select_celsius {
  static T convert(T v) { return v - 273.15; }
};

template <typename T> struct select_celsius<T, true> {
  static T convert(T v) { return v; }
};

template <typename U, Temperature units> U air_dynamic_viscosity(U T)
{
  return 1.71432e-5 + 4.828e-8 * select_celsius<U, units == Temperature::Celsius>::convert(T);
}

template <typename T, Temperature units> T air_density(T const pb,      // barometric pressure (Pascals)
                                                       T const tdb,     // dry bulb temperature (Celsius)
                                                       T const dw=0.0   // humidity ratio (kgWater/kgDryAir)
)
{
  // FUNCTION INFORMATION:
  //       AUTHOR         G. S. Wright
  //       DATE WRITTEN   June 2, 1994
  //       MODIFIED       na
  //       RE-ENGINEERED  na

  // PURPOSE OF THIS FUNCTION:
  // This function provides density of air as a function of barometric
  // pressure, dry bulb temperature, and humidity ratio.

  // METHODOLOGY EMPLOYED:
  // ideal gas law
  //    universal gas const for air 287 J/(kg K)
  //    air/water molecular mass ratio 28.9645/18.01534

  // REFERENCES:
  // Wylan & Sontag, Fundamentals of Classical Thermodynamics.
  // ASHRAE handbook 1985 Fundamentals, Ch. 6, eqn. (6),(26)
  return (pb / (287.0 * select_kelvin<T, units == Temperature::Kelvin>::convert(tdb) * (1.0 + 1.6077687 * std::max(dw, 1.0e-5))));
}

} // namespace airflownetwork
