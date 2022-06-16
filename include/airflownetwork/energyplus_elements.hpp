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

#ifndef AIRFLOWNETWORK_ENERGYPLUS_ELEMENTS_HPP

namespace airflownetwork {

template <typename T, bool> struct select_ratio {
  static T ratio(T v1, T v2) { return v1/v2; }
};

template <typename T> struct select_ratio<T, true> {
  static T ratio(T v1, T v2) { return (v1 + 273.15)/(v2 + 273.15); }
};

#ifndef AIRFLOWNETWORK_ELEMENTS_HPP
enum class Temperature; // Forward declare the enum so VS will not complain so much
#endif

template <typename L, Temperature units> struct Crack : public Element<L> // Power law crack element
{
  using value_type = Element<L>::value_type;

  Crack(const std::string &name, value_type coefficient, value_type exponent, value_type reference_density, value_type reference_viscosity)
    : Element<L>(name), coefficient(coefficient), exponent(exponent), reference_density(reference_density), reference_viscosity(reference_viscosity)
  {}

  virtual ~Crack()
  {}

  virtual int calculate(bool const init,              // Initialization flag.If true, use linear relationship
                        const L& link,                // Linkage using this element
                        std::array<value_type, 2>& F, // Airflow through the component [kg/s]
                        std::array<value_type, 2>& DF // Partial derivative:  DF/DP
  ) override
  {
    // SUBROUTINE INFORMATION:
    //       AUTHOR         George Walton
    //       DATE WRITTEN   Extracted from AIRNET
    //       MODIFIED       Lixing Gu, 2/1/04
    //                      Revised the subroutine to meet E+ needs
    //       MODIFIED       Lixing Gu, 6/8/05
    //       RE-ENGINEERED  Jason DeGraw

    // PURPOSE OF THIS SUBROUTINE:
    // This subroutine solves airflow for a surface crack component

    auto& propN{ link.nodes[0] };
    auto& propM{ link.nodes[1] };
    auto upwind{ link.nodes[0] };

    value_type VisAve{ 0.5 * (propN->viscosity + propM->viscosity) };
    value_type Tave{ 0.5 * (propN->temperature + propM->temperature) };

    value_type sign{ 1.0 };
    value_type abs_pdrop = link.pressure_drop;

    if (link.pressure_drop < 0.0) {
      sign = -1.0;
      abs_pdrop = -link.pressure_drop;
      upwind = link.nodes[1];
    }

    value_type coef = coefficient * link.control * link.multiplier / upwind->sqrt_density;

    // Linear calculation
    value_type RhoCor = select_ratio<value_type, units == Temperature::Kelvin>::ratio(upwind->temperature, Tave);
    value_type Ctl{ std::pow(reference_density / upwind->density / RhoCor, exponent - 1.0) * std::pow(reference_viscosity / VisAve, 2.0 * exponent - 1.0) };
    value_type CDM{ coef * upwind->density / upwind->viscosity * Ctl };
    value_type FL{ CDM * link.pressure_drop };
    value_type abs_FT;

    if (init) {
      DF[0] = CDM;
      F[0] = FL;
    } else {
      // Nonlinear flow.
      if (exponent == 0.5) {
        abs_FT = coef * upwind->sqrt_density * std::sqrt(abs_pdrop) * Ctl;
      } else {
        abs_FT = coef * upwind->sqrt_density * std::pow(abs_pdrop, exponent) * Ctl;
      }
      // Select linear or nonlinear flow.
      if (std::abs(FL) <= abs_FT) {
        F[0] = FL;
        DF[0] = CDM;
      } else {
        F[0] = sign * abs_FT;
        DF[0] = F[0] * exponent / link.pressure_drop;
      }
    }
    return 1;
  }

  virtual int calculate(const L& link,                // Linkage using this element
                        std::array<value_type, 2>& F, // Airflow through the component [kg/s]
                        std::array<value_type, 2>& DF // Partial derivative:  DF/DP
  ) override
  {
    // SUBROUTINE INFORMATION:
    //       AUTHOR         George Walton
    //       DATE WRITTEN   Extracted from AIRNET
    //       MODIFIED       Lixing Gu, 2/1/04
    //                      Revised the subroutine to meet E+ needs
    //       MODIFIED       Lixing Gu, 6/8/05
    //       RE-ENGINEERED  Jason DeGraw

    // PURPOSE OF THIS SUBROUTINE:
    // This subroutine solves airflow for a surface crack component

    auto propN{ link.nodes[0] };
    auto propM{ link.nodes[1] };
    auto upwind{ link.nodes[0] };

    value_type VisAve{ 0.5 * (propN->viscosity + propM->viscosity) };
    value_type Tave{ 0.5 * (propN->temperature + propM->temperature) };

    value_type sign{ 1.0 };
    value_type abs_pdrop = link.pressure_drop;

    if (link.pressure_drop < 0.0) {
      sign = -1.0;
      upwind = link.nodes[1];
      abs_pdrop = -link.pressure_drop;
    }

    value_type coef = coefficient * link.control * link.multiplier / upwind->sqrt_density;

    // Linear calculation
    value_type RhoCor = select_ratio<value_type, units == Temperature::Kelvin>::ratio(upwind->temperature, Tave);
    value_type Ctl{ std::pow(reference_density / upwind->density / RhoCor, exponent - 1.0) * std::pow(reference_viscosity / VisAve, 2.0 * exponent - 1.0) };
    value_type CDM{ coef * upwind->density / upwind->viscosity * Ctl };
    value_type FL{ CDM * link.pressure_drop };
    value_type abs_FT;

    // Nonlinear flow.
    if (exponent == 0.5) {
      abs_FT = coef * upwind->sqrt_density * std::sqrt(abs_pdrop) * Ctl;
    } else {
      abs_FT = coef * upwind->sqrt_density * std::pow(abs_pdrop, exponent) * Ctl;
    }
    // Select linear or nonlinear flow.
    if (std::abs(FL) <= abs_FT) {
      F[0] = FL;
      DF[0] = CDM;
    } else {
      F[0] = sign * abs_FT;
      DF[0] = F[0] * exponent / link.pressure_drop;
    }

    return 1;
  }

  virtual Type type() const override
  {
    return Type::SCR;
  }

  value_type coefficient;         // Air mass flow coefficient [kg/s at 1Pa]
  value_type exponent;            // Air mass flow exponent [dimensionless]
  value_type reference_density;   // Reference density for crack data
  value_type reference_viscosity; // Reference viscosity for crack data
};


template <typename L, Temperature units> struct TwoWayOpening : public Crack<L, units> // Single opening two-way flow element
{
  using value_type = Element<L>::value_type;

  constexpr static value_type sqrt2 = 1.414213562;
  constexpr static value_type two_thirds = 2.0/3.0;
  constexpr static value_type gravity = 9.8;
  constexpr static value_type rhoair = 1.2041;             // density of standard air
  constexpr static value_type vsair = 1.81625e-5;          // dynamic viscosity of standard air

  TwoWayOpening(const std::string& name, value_type linear, value_type nonlinear, value_type exponent, value_type discharge_coefficient,
                value_type density_difference)
    : Crack<L, units>(name, linear, nonlinear, exponent, rhoair, vsair), height(height), width(width), discharge_coefficient(discharge_coefficient)
  {}

  virtual int calculate(bool const init,              // Initialization flag.If true, use linear relationship
                        const L& link,                // Linkage using this element
                        std::array<value_type, 2>& F, // Airflow through the component [kg/s]
                        std::array<value_type, 2>& DF // Partial derivative:  DF/DP
  ) override
  {
    auto& pn{ reference_to(link.nodes[0]) };
    auto& pm{ reference_to(link.nodes[1]) };
    int nf{ 1 };

    auto drho = pn.density - pm.density;
    auto abs_pdrop = std::abs(link.pressure_drop);
    if ((std::abs(drho) < density_difference) || init) {
      // Treat as orifice flow
      return Crack<L, units>::calculate(init, link, F, DF);
    }

    auto gdrho = gravity * drho;
    auto Y = link.pressure_drop / gdrho;
    link.Y = Y;
    // F0 = lower flow, FH = upper flow.
    auto C = sqrt2 * link.width * discharge_coefficient;
    auto B = C * std::sqrt(abs_pdrop);
    auto DF0 = B / abs(gdrho);
    // F0 = 0.666667d0 * C * SQRT(ABS(GDRHO * Y)) * ABS(Y)
    auto F0 = two_thirds * B * abs(Y);
    auto DFH = C * std::sqrt(std::abs((link.height - Y) / gdrho));
    // FH = 0.666667d0 * DFH * ABS(GDRHO * (Height - Y))
    auto FH = two_thirds * DFH * std::abs(gdrho * (link.height - Y));

    if (Y <= 0.0) {
      // One - way flow(negative).
      if (gdrho >= 0.0) {
        F[0] = -pm.sqrt_density * abs(FH - F0);
        DF[0] = pm.sqrt_density * abs(DFH - DF0);
      } else {
        F[0] = pn.sqrt_density * abs(FH - F0);
        DF[0] = pn.sqrt_density * abs(DFH - DF0);
      }
    } else if (Y >= link.height) {
      // One - way flow(positive).
      if (gdrho >= 0.0) {
        F[0] = pn.sqrt_density * abs(FH - F0);
        DF[0] = pn.sqrt_density * abs(DFH - DF0);
      } else {
        F[0] = -pm.sqrt_density * abs(FH - F0);
        DF[0] = pm.sqrt_density * abs(DFH - DF0);
      }
    } else {
      // Two - way flow.
      nf = 2;
      if (gdrho >= 0.0) {
        F[0] = -pm.sqrt_density * FH;
        DF[0] = pm.sqrt_density * DFH;
        F[1] = pn.sqrt_density * F0;
        DF[1] = pn.sqrt_density * DF0;
      } else {
        F[0] = pn.sqrt_density * FH;
        DF[0] = pn.sqrt_density * DFH;
        F[1] = -pm.sqrt_density * F0;
        DF[1] = pm.sqrt_density * DF0;
      }
    }
    return nf;
  }

  virtual int calculate(const L& link,                // Linkage using this element
                        std::array<value_type, 2>& F, // Airflow through the component [kg/s]
                        std::array<value_type, 2>& DF // Partial derivative:  DF/DP
  ) override
  {
    auto& pn{ reference_to(link.nodes[0]) };
    auto& pm{ reference_to(link.nodes[1]) };
    int nf{ 1 };

    auto drho = pn.density - pm.density;
    auto abs_pdrop = std::abs(link.pressure_drop);
    if ((std::abs(drho) < density_difference)) {
      // Treat as orifice flow
      return Crack<L, units>::calculate(link, F, DF);
    }


    auto gdrho = gravity * drho;
    auto Y = link.pressure_drop / gdrho;
    link.Y = Y;
    // F0 = lower flow, FH = upper flow.
    auto C = sqrt2 * width * discharge_coefficient;
    auto B = C * std::sqrt(abs_pdrop);
    auto DF0 = B / abs(gdrho);
    // F0 = 0.666667d0 * C * SQRT(ABS(GDRHO * Y)) * ABS(Y)
    auto F0 = two_thirds * B * abs(Y);
    auto DFH = C * std::sqrt(std::abs((height - Y) / gdrho));
    // FH = 0.666667d0 * DFH * ABS(GDRHO * (Height - Y))
    auto FH = two_thirds * DFH * std::abs(gdrho * (height - Y));

    if (Y <= 0.0) {
      // One - way flow(negative).
      if (gdrho >= 0.0) {
        F[0] = -pm.sqrt_density * abs(FH - F0);
        DF[0] = pm.sqrt_density * abs(DFH - DF0);
      } else {
        F[0] = pn.sqrt_density * abs(FH - F0);
        DF[0] = pn.sqrt_density * abs(DFH - DF0);
      }
    } else if (Y >= height) {
      // One - way flow(positive).
      if (gdrho >= 0.0) {
        F[0] = pn.sqrt_density * abs(FH - F0);
        DF[0] = pn.sqrt_density * abs(DFH - DF0);
      } else {
        F[0] = -pm.sqrt_density * abs(FH - F0);
        DF[0] = pm.sqrt_density * abs(DFH - DF0);
      }
    } else {
      // Two - way flow.
      nf = 2;
      if (gdrho >= 0.0) {
        F[0] = -pm.sqrt_density * FH;
        DF[0] = pm.sqrt_density * DFH;
        F[1] = pn.sqrt_density * F0;
        DF[1] = pn.sqrt_density * DF0;
      } else {
        F[0] = pn.sqrt_density * FH;
        DF[0] = pn.sqrt_density * DFH;
        F[1] = -pm.sqrt_density * F0;
        DF[1] = pm.sqrt_density * DF0;
      }
    }
    return nf;
  }

  virtual Type type() const override
  {
    return Type::SOP;
  }

  value_type height;
  value_type width;
  value_type discharge_coefficient;
  value_type density_difference;

};


template <typename L, Temperature units> struct AltTwoWayOpening : public Crack<L, units> // Single opening two-way flow element
{
  using value_type = Element<L>::value_type;

  constexpr static value_type sqrt2 = 1.414213562;
  constexpr static value_type two_thirds = 2.0/3.0;
  constexpr static value_type gravity = 9.8;
  constexpr static value_type rhoair = 1.2041;             // density of standard air
  constexpr static value_type vsair = 1.81625e-5;          // dynamic viscosity of standard air

  AltTwoWayOpening(const std::string& name, value_type linear, value_type nonlinear, value_type exponent, value_type height, value_type width,
                   value_type discharge_coefficient, value_type density_difference)
    : Crack<L, units>(name, linear, nonlinear, exponent, rhoair, vsair), height(height), width(width), discharge_coefficient(discharge_coefficient)
  {}

  virtual int calculate(bool const init,              // Initialization flag.If true, use linear relationship
                        const L& link,                // Linkage using this element
                        std::array<value_type, 2>& F, // Airflow through the component [kg/s]
                        std::array<value_type, 2>& DF // Partial derivative:  DF/DP
  ) override
  {
    auto& pn{ reference_to(link.nodes[0]) };
    auto& pm{ reference_to(link.nodes[1]) };
    int nf{ 1 };

    auto drho = pn.density - pm.density;
    auto abs_pdrop = std::abs(link.pressure_drop);
    if ((std::abs(drho) < density_difference) || init) {
      // Treat as orifice flow
      return Crack<L, units>::calculate(init, link, F, DF);
    }


    auto gdrho = gravity * drho;
    auto Y = link.pressure_drop / gdrho;
    link.Y = Y;
    // F0 = lower flow, FH = upper flow.
    auto C = sqrt2 * width * discharge_coefficient;
    auto B = C * std::sqrt(abs_pdrop);
    auto DF0 = B / abs(gdrho);
    // F0 = 0.666667d0 * C * SQRT(ABS(GDRHO * Y)) * ABS(Y)
    auto F0 = two_thirds * B * abs(Y);
    auto DFH = C * std::sqrt(std::abs((height - Y) / gdrho));
    // FH = 0.666667d0 * DFH * ABS(GDRHO * (Height - Y))
    auto FH = two_thirds * DFH * std::abs(gdrho * (height - Y));

    if (Y <= 0.0) {
      // One - way flow(negative).
      if (gdrho >= 0.0) {
        F[0] = -pm.sqrt_density * abs(FH - F0);
        DF[0] = pm.sqrt_density * abs(DFH - DF0);
      } else {
        F[0] = pn.sqrt_density * abs(FH - F0);
        DF[0] = pn.sqrt_density * abs(DFH - DF0);
      }
    } else if (Y >= height) {
      // One - way flow(positive).
      if (gdrho >= 0.0) {
        F[0] = pn.sqrt_density * abs(FH - F0);
        DF[0] = pn.sqrt_density * abs(DFH - DF0);
      } else {
        F[0] = -pm.sqrt_density * abs(FH - F0);
        DF[0] = pm.sqrt_density * abs(DFH - DF0);
      }
    } else {
      // Two - way flow.
      nf = 2;
      if (gdrho >= 0.0) {
        F[0] = -pm.sqrt_density * FH;
        DF[0] = pm.sqrt_density * DFH;
        F[1] = pn.sqrt_density * F0;
        DF[1] = pn.sqrt_density * DF0;
      } else {
        F[0] = pn.sqrt_density * FH;
        DF[0] = pn.sqrt_density * DFH;
        F[1] = -pm.sqrt_density * F0;
        DF[1] = pm.sqrt_density * DF0;
      }
    }
    return nf;
  }

  virtual int calculate(const L& link,                // Linkage using this element
                        std::array<value_type, 2>& F, // Airflow through the component [kg/s]
                        std::array<value_type, 2>& DF // Partial derivative:  DF/DP
  ) override
  {
    auto& pn{ reference_to(link.nodes[0]) };
    auto& pm{ reference_to(link.nodes[1]) };
    int nf{ 1 };

    auto drho = pn.density - pm.density;
    auto abs_pdrop = std::abs(link.pressure_drop);
    if ((std::abs(drho) < density_difference)) {
      // Treat as orifice flow
      return Crack<L, units>::calculate(link, F, DF);
    }

    auto gdrho = gravity * drho;
    auto Y = link.pressure_drop / gdrho;
    link.Y = Y;
    // F0 = lower flow, FH = upper flow.
    auto C = sqrt2 * width * discharge_coefficient;
    auto B = C * std::sqrt(abs_pdrop);
    auto DF0 = B / abs(gdrho);
    // F0 = 0.666667d0 * C * SQRT(ABS(GDRHO * Y)) * ABS(Y)
    auto F0 = two_thirds * B * abs(Y);
    auto DFH = C * std::sqrt(std::abs((height - Y) / gdrho));
    // FH = 0.666667d0 * DFH * ABS(GDRHO * (Height - Y))
    auto FH = two_thirds * DFH * std::abs(gdrho * (height - Y));

    if (Y <= 0.0) {
      // One - way flow(negative).
      if (gdrho >= 0.0) {
        F[0] = -pm.sqrt_density * abs(FH - F0);
        DF[0] = pm.sqrt_density * abs(DFH - DF0);
      } else {
        F[0] = pn.sqrt_density * abs(FH - F0);
        DF[0] = pn.sqrt_density * abs(DFH - DF0);
      }
    } else if (Y >= height) {
      // One - way flow(positive).
      if (gdrho >= 0.0) {
        F[0] = pn.sqrt_density * abs(FH - F0);
        DF[0] = pn.sqrt_density * abs(DFH - DF0);
      } else {
        F[0] = -pm.sqrt_density * abs(FH - F0);
        DF[0] = pm.sqrt_density * abs(DFH - DF0);
      }
    } else {
      // Two - way flow.
      nf = 2;
      if (gdrho >= 0.0) {
        F[0] = -pm.sqrt_density * FH;
        DF[0] = pm.sqrt_density * DFH;
        F[1] = pn.sqrt_density * F0;
        DF[1] = pn.sqrt_density * DF0;
      } else {
        F[0] = pn.sqrt_density * FH;
        DF[0] = pn.sqrt_density * DFH;
        F[1] = -pm.sqrt_density * F0;
        DF[1] = pm.sqrt_density * DF0;
      }
    }
    return nf;
  }

  virtual Type type() const override
  {
    return Type::SOP;
  }

  value_type height;
  value_type width;
  value_type discharge_coefficient;
  value_type density_difference;

};

template <typename L> struct SpecifiedMassFlow : public Element<L> // Specified mass flow element
{
  using value_type = Element<L>::value_type;

  // Default Constructor
  SpecifiedMassFlow(const std::string& name, value_type mass_flow = 0.0) : Element<L>(name), mass_flow(mass_flow)
  {}

  virtual int calculate([[maybe_unused]] bool const init, // Initialization flag.If true, use linear relationship
                        const L& link,                    // Linkage using this element
                        std::array<value_type, 2>& F,     // Airflow through the component [kg/s]
                        std::array<value_type, 2>& DF     // Partial derivative:  DF/DP
  ) override
  {
    // SUBROUTINE INFORMATION:
    //       AUTHOR         Jason DeGraw and Prateek Shrestha
    //       DATE WRITTEN   June 2021
    //       MODIFIED       na
    //       MODIFIED       na
    //       RE-ENGINEERED  na

    // PURPOSE OF THIS SUBROUTINE:
    // This subroutine computes airflow for a specified mass flow element.

    F[0] = mass_flow * link.control * link.multiplier;
    DF[0] = 0.0;
    F[1] = 0.0;
    DF[1] = 0.0;

    return 1;
  }

  virtual int calculate(const L& link,                    // Linkage using this element
                        std::array<value_type, 2>& F,     // Airflow through the component [kg/s]
                        std::array<value_type, 2>& DF     // Partial derivative:  DF/DP
  ) override
  {
    // SUBROUTINE INFORMATION:
    //       AUTHOR         Jason DeGraw and Prateek Shrestha
    //       DATE WRITTEN   June 2021
    //       MODIFIED       na
    //       MODIFIED       na
    //       RE-ENGINEERED  na

    // PURPOSE OF THIS SUBROUTINE:
    // This subroutine computes airflow for a specified mass flow element.

    F[0] = mass_flow * link.control * link.multiplier;
    DF[0] = 0.0;
    F[1] = 0.0;
    DF[1] = 0.0;

    return 1;
  }

  virtual Type type() const override
  {
    return Type::SMF;
  }

  value_type mass_flow; // Mass Flow [kg/s]
};

template <typename L> struct SpecifiedVolumeFlow : public Element<L> // Specified volume flow element
{
  using value_type = Element<L>::value_type;

  SpecifiedVolumeFlow(const std::string& name, value_type volume_flow=0.0) : Element<L>(name), volume_flow(volume_flow)
  {}

  virtual int calculate([[maybe_unused]] bool const init, // Initialization flag.If true, use linear relationship
                        const L& link,                    // Linkage using this element
                        std::array<value_type, 2>& F,     // Airflow through the component [kg/s]
                        std::array<value_type, 2>& DF     // Partial derivative:  DF/DP
  ) override
  {
    // SUBROUTINE INFORMATION:
    //       AUTHOR         Jason DeGraw and Prateek Shrestha
    //       DATE WRITTEN   June 2021
    //       MODIFIED       na
    //       MODIFIED       na
    //       RE-ENGINEERED  na

    // PURPOSE OF THIS SUBROUTINE:
    // This subroutine computes airflow for a specified volume flow element.

    auto upwind{ link.nodes[0] };

    if (volume_flow < 0.0) {
      upwind = link.nodes[1];
    }

    F[0] = upwind->density * volume_flow * link.control * link.multiplier;
    DF[0] = 0.0;
    F[1] = 0.0;
    DF[1] = 0.0;

    return 1;
  }

  virtual int calculate(const L& link,                    // Linkage using this element
                        std::array<value_type, 2>& F,     // Airflow through the component [kg/s]
                        std::array<value_type, 2>& DF     // Partial derivative:  DF/DP
  ) override
  {
    // SUBROUTINE INFORMATION:
    //       AUTHOR         Jason DeGraw and Prateek Shrestha
    //       DATE WRITTEN   June 2021
    //       MODIFIED       na
    //       MODIFIED       na
    //       RE-ENGINEERED  na

    // PURPOSE OF THIS SUBROUTINE:
    // This subroutine computes airflow for a specified volume flow element.

    auto upwind{ link.nodes[0] };

    if (volume_flow < 0.0) {
      upwind = link.nodes[1];
    }

    F[0] = upwind->density * volume_flow * link.control * link.multiplier;
    DF[0] = 0.0;
    F[1] = 0.0;
    DF[1] = 0.0;

    return 1;
  }


  virtual Type type() const override
  {
    return Type::SVF;
  }

  value_type volume_flow; // Volume Flow [m3/s]
};

} // namespace airflownetwork

#endif // AIRFLOWNETWORK_ENERGYPLUS_ELEMENTS_HPP
