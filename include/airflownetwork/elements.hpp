// Copyright (c) 2020-2022, Oak Ridge National Laboratory, managed by UT-Battelle.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted
// provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this list of
//    conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice, this list of
//    conditions and the following disclaimer in the documentation and/or other materials
//    provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors may be used to
//    endorse or promote products derived from this software without specific prior written
//    permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
// AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef AIRFLOWNETWORK_ELEMENTS_HPP

#include <string>
#include <cmath>
#include <type_traits>
#ifdef DONT_HAVE_CPP20_NUMBERS
#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif
#include <math.h>
namespace std {
namespace numbers {
template<typename T> inline constexpr T sqrt2_v = M_SQRT2;
}
}
#ifdef _MSC_VER
#undef _USE_MATH_DEFINES
#endif
#else
#include <numbers>
#endif

#include "airflownetwork/properties.hpp"

namespace airflownetwork {

enum class Type : int
{
  DOP = 1, // Detailed large opening component
  SOP,     // Simple opening component
  SCR,     // Surface crack component
  SEL,     // Surface effective leakage ratio component
  PLR,     // Distribution system crack component
  DWC,     // Distribution system duct component
  CVF,     // Distribution system constant volume fan component
  FAN,     // Distribution system detailed fan component
  MRR,     // Distribution system multiple curve fit power law resistant flow component
  DMP,     // Distribution system damper component
  ELR,     // Distribution system effective leakage ratio component
  CPD,     // Distribution system constant pressure drop component
  COI,     // Distribution system coil component
  TMU,     // Distribution system terminal unit component
  EXF,     // Zone exhaust fan
  HEX,     // Distribution system heat exchanger
  HOP,     // Horizontal opening component
  RVD,     // Reheat VAV terminal damper
  OAF,     // Distribution system OA
  REL,     // Distribution system relief air
  SMF,     // Specified mass flow component
  SVF,     // Specified volume flow component
  FCN,     // ContamX mass flow power law component
  QCN,     // ContamX volume flow power law component
  OPN,     // Opening with specified area
  VOP      // Variable area opening
};

template <typename L> struct Element
{
  using value_type = decltype(L::pressure_drop);

  Element()
  {}

  Element(const std::string& name) : name(name)
  {}

  virtual ~Element()
  {}

  virtual int calculate(bool const init,              // Initialization flag.If true, use linear relationship
                        const L& link,                // Linkage using this element
                        std::array<value_type, 2>& F, // Airflow through the component [kg/s]
                        std::array<value_type, 2>& DF // Partial derivative:  DF/DP
  ) = 0;

  virtual int calculate(const L& link,                // Linkage using this element
                        std::array<value_type, 2>& F, // Airflow through the component [kg/s]
                        std::array<value_type, 2>& DF // Partial derivative:  DF/DP
  ) = 0;

  virtual Type type() const = 0;

  std::string name; // Name of airflow element
};

template <typename L> struct PowerLaw : public Element<L> // Power law element
{
  using value_type = Element<L>::value_type;

  PowerLaw(const std::string& name, value_type linear, value_type nonlinear, value_type exponent,
           value_type reference_density, value_type reference_kinematic_viscosity)
    : Element<L>(name), linear(linear), nonlinear(nonlinear), exponent(exponent), reference_density(reference_density),
    reference_kinematic_viscosity(reference_kinematic_viscosity)
  {}

  PowerLaw(const std::string& name, value_type nonlinear, value_type exponent,
    value_type reference_density, value_type reference_kinematic_viscosity)
    : Element<L>(name), nonlinear(nonlinear), exponent(exponent), reference_density(reference_density),
    reference_kinematic_viscosity(reference_kinematic_viscosity)
  {
    linear = compute_linear_coefficient(30.0, reference_density, reference_kinematic_viscosity * reference_density);
  }

  virtual ~PowerLaw()
  {}

  virtual value_type compute_linear_coefficient(value_type transition_reynum, value_type reference_density, value_type reference_viscosity,
                                                value_type dp_min=1.0e-10) const
  {
    value_type A = this->nonlinear / (0.6 * std::numbers::sqrt2_v<value_type>);
    value_type F = reference_viscosity * transition_reynum * std::sqrt(A);
    value_type pdrop = std::max(std::pow(F / (this->nonlinear * std::sqrt(reference_density)), 1.0 / this->exponent), dp_min);
    return reference_viscosity * F / (reference_density * pdrop);
  }

  virtual value_type correction(decltype(L::nodes[0]) node) const
  {
    if (exponent == 0.5) {
      return node->sqrt_density;
    }
    return std::exp(std::log(reference_density / node->density) * (exponent - 0.5) +
                    std::log(reference_kinematic_viscosity * node->rdv) * (2 * exponent - 1.0)) * node->sqrt_density;
  }

  virtual int calculate(bool const init,              // Initialization flag.If true, use linear relationship
                        const L& link,                // Linkage using this element
                        std::array<value_type, 2>& F, // Airflow through the component [kg/s]
                        std::array<value_type, 2>& DF // Partial derivative:  DF/DP
  ) override
  {
    value_type sign{ 1.0 };
    value_type abs_pdrop = link.pressure_drop;
    auto upwind{ link.nodes[0] };
    if (link.pressure_drop < 0.0) {
      sign = -1.0;
      abs_pdrop = -link.pressure_drop;
      upwind = link.nodes[1];
    }

    auto ctrl = link.control * link.multiplier;
    auto cdm = linear * upwind->rdv;
    auto abs_Fl = cdm * abs_pdrop;

    if (init) {
      F[0] = sign * abs_Fl * ctrl;
      DF[0] = cdm * ctrl;
      return 1;
    }

    value_type abs_Ft;

    if (exponent == 0.5) {
      abs_Ft = this->correction(upwind) * nonlinear * std::sqrt(abs_pdrop);
    } else {
      abs_Ft = this->correction(upwind) * nonlinear * std::pow(abs_pdrop, exponent);
    }

    if (abs_Fl <= abs_Ft) {
      F[0] = sign * abs_Fl * ctrl;
      DF[0] = cdm * ctrl;
    } else {
      F[0] = sign * abs_Ft * ctrl;
      DF[0] = F[0] * exponent / link.pressure_drop;
    }
    
    return 1;
  }

  virtual int calculate(const L& link,                // Linkage using this element
                        std::array<value_type, 2>& F, // Airflow through the component [kg/s]
                        std::array<value_type, 2>& DF // Partial derivative:  DF/DP
  ) override
  {
    value_type sign{ 1.0 };
    value_type abs_pdrop = link.pressure_drop;
    auto upwind{ link.nodes[0] };
    if (link.pressure_drop < 0.0) {
      sign = -1.0;
      abs_pdrop = -link.pressure_drop;
      upwind = link.nodes[1];
    }

    auto ctrl = link.control * link.multiplier;
    auto cdm = linear * upwind->rdv;
    auto abs_Fl = cdm * abs_pdrop;

    value_type abs_Ft;

    if (exponent == 0.5) {
      abs_Ft = this->correction(upwind) * nonlinear * std::sqrt(abs_pdrop);
    } else {
      abs_Ft = this->correction(upwind) * nonlinear * std::pow(abs_pdrop, exponent);
    }

    if (abs_Fl <= abs_Ft) {
      F[0] = sign * abs_Fl * ctrl;
      DF[0] = cdm * ctrl;
    }
    else {
      F[0] = sign * abs_Ft * ctrl;
      DF[0] = F[0] * exponent / link.pressure_drop;
    }

    return 1;
  }

  virtual Type type() const override
  {
    return Type::PLR;
  }

  value_type linear;                        // Linear air mass flow coefficient [kg/s at 1Pa]
  value_type nonlinear;                     // Nonlinear air mass flow coefficient [kg/s at 1Pa]
  value_type exponent;                      // Air mass flow exponent [dimensionless]
  value_type reference_density;             // Reference density
  value_type reference_kinematic_viscosity; // Reference viscosity
};

template <typename L> struct CxPowerLaw : public PowerLaw<L> // Power law element
{
  using value_type = Element<L>::value_type;

  constexpr static value_type rhoair = 1.2041;             // density of standard air
  constexpr static value_type sqrt_rho = 1.097315;         // sqrt(RHOAIR)
  constexpr static value_type vsair = 1.81625e-5;          // dynamic viscosity of standard air
  constexpr static value_type nuair = 1.50839e-5;          // kinematic viscosity (mu / rho) of standard air
  constexpr static value_type transition_reynum = 30.0;    // transition Reynolds number
  constexpr static value_type min_pressure_drop = 1.0e-10; // smallest allowed transition pressure drop

  CxPowerLaw(const std::string& name, value_type linear, value_type nonlinear, value_type exponent)
    : PowerLaw<L>(name, linear, nonlinear, exponent, rhoair, nuair)
  {}

  CxPowerLaw(const std::string& name, value_type nonlinear, value_type exponent)
    : PowerLaw<L>(name, nonlinear, exponent, rhoair, nuair)
  {}

  virtual ~CxPowerLaw()
  {}

  virtual value_type compute_linear_coefficient(value_type transition_reynum, value_type, value_type, value_type dp_min = 1.0e-10) const override
  {
    value_type reference_density = rhoair;
    value_type reference_viscosity = vsair;
    value_type reference_sqrt_density = sqrt_rho;
    value_type A = this->nonlinear / (0.6 * std::numbers::sqrt2_v<value_type>);
    value_type F = reference_viscosity * transition_reynum * std::sqrt(A);
    value_type pdrop = std::max(std::pow(F / (this->nonlinear * reference_sqrt_density), 1.0 / this->exponent), dp_min);
    return reference_viscosity * F / (reference_density * pdrop);
  }

  using PowerLaw<L>::correction;

  virtual int calculate(bool const init,              // Initialization flag.If true, use linear relationship
                        const L& link,                // Linkage using this element
                        std::array<value_type, 2>& F, // Airflow through the component [kg/s]
                        std::array<value_type, 2>& DF // Partial derivative:  DF/DP
  ) override
  {
    auto upwind{ link.nodes[0] };
    if (link.pressure_drop < 0.0) {
      upwind = link.nodes[1];
    }

    auto ctrl = link.control * link.multiplier;

    auto cdm = this->linear * upwind->rdv;
    auto Fl = cdm * link.pressure_drop;

    if (init) {
      F[0] = Fl;
      DF[0] = cdm;
      return 1;
    }

    auto Tadj = this->correction(upwind);
    value_type Flow;
    value_type Ft;
    value_type dF;

    if (link.pressure_drop >= 0) {
      Ft = Tadj * this->nonlinear * std::pow(link.pressure_drop, this->exponent);
    } else {
      Ft = -Tadj * this->nonlinear * std::pow(-link.pressure_drop, this->exponent);
    }

    if (std::abs(Fl) <= std::abs(Ft)) {
      Flow = Fl;
      dF = cdm;
    } else {
      Flow = Ft;
      dF = Ft * this->exponent / link.pressure_drop;
    }

    F[0] = Flow * ctrl;
    DF[0] = dF * ctrl;

    return 1;
  }

  virtual int calculate(const L& link,                // Linkage using this element
                        std::array<value_type, 2>& F, // Airflow through the component [kg/s]
                        std::array<value_type, 2>& DF // Partial derivative:  DF/DP
  ) override
  {
    auto upwind{ link.nodes[0] };
    if (link.pressure_drop < 0.0) {
      upwind = link.nodes[1];
    }

    auto ctrl = link.control * link.multiplier;

    auto cdm = this->linear * upwind->rdv;
    auto Fl = cdm * link.pressure_drop;

    auto Tadj = this->correction(upwind);
    value_type Flow;
    value_type Ft;
    value_type dF;

    if (link.pressure_drop >= 0) {
      Ft = Tadj * this->nonlinear * std::pow(link.pressure_drop, this->exponent);
    }
    else {
      Ft = -Tadj * this->nonlinear * std::pow(-link.pressure_drop, this->exponent);
    }

    if (std::abs(Fl) <= std::abs(Ft)) {
      Flow = Fl;
      dF = cdm;
    }
    else {
      Flow = Ft;
      dF = Ft * this->exponent / link.pressure_drop;
    }

    F[0] = Flow * ctrl;
    DF[0] = dF * ctrl;

    return 1;
  }

  virtual Type type() const override
  {
    return Type::PLR;
  }

};

template <typename L> struct CxVolumeFlow : public CxPowerLaw<L> // Power law element
{
  using value_type = Element<L>::value_type;

  CxVolumeFlow(const std::string& name, value_type linear, value_type nonlinear, value_type exponent)
    : CxPowerLaw<L>(name, linear, nonlinear, exponent)
  {}

  CxVolumeFlow(const std::string& name, value_type nonlinear, value_type exponent)
    : CxPowerLaw<L>(name, nonlinear, exponent)
  {}

  virtual ~CxVolumeFlow()
  {}

  virtual value_type compute_linear_coefficient(value_type transition_reynum, value_type, value_type, value_type dp_min = 1.0e-10) const override
  {
    value_type reference_density = CxPowerLaw<L>::rhoair;
    value_type reference_viscosity = CxPowerLaw<L>::vsair;
    value_type reference_sqrt_density = CxPowerLaw<L>::sqrt_rho;
    value_type A = reference_sqrt_density * this->nonlinear / (0.6 * std::numbers::sqrt2_v<value_type>);
    value_type F = reference_viscosity * transition_reynum * std::sqrt(A);
    value_type pdrop = std::max(std::pow(F / (this->nonlinear * reference_density), 1.0 / this->exponent), dp_min);
    return reference_viscosity * F / (reference_density * pdrop);
  }

  virtual value_type correction(decltype(L::nodes[0]) node) const override
  {
    if (this->exponent == 0.5) {
      return node->density;
    }
    return std::exp(std::log(this->reference_density / node->density) * this->exponent +
                    std::log(this->reference_kinematic_viscosity * node->rdv) * (2 * this->exponent - 1.0)) * node->density;
  }

  virtual Type type() const override
  {
    return Type::QCN;
  }
};

template <typename L> struct CxMassFlow : public CxPowerLaw<L> // Power law element
{
  using value_type = Element<L>::value_type;

  CxMassFlow(const std::string& name, value_type linear, value_type nonlinear, value_type exponent)
    : CxPowerLaw<L>(name, linear, nonlinear, exponent)
  {}

  CxMassFlow(const std::string& name, value_type nonlinear, value_type exponent)
    : CxPowerLaw<L>(name, nonlinear, exponent)
  {}

  virtual ~CxMassFlow()
  {}

  virtual value_type compute_linear_coefficient(value_type transition_reynum, value_type, value_type, value_type dp_min = 1.0e-10) const override
  {
    value_type reference_density = CxPowerLaw<L>::rhoair;
    value_type reference_viscosity = CxPowerLaw<L>::vsair;
    value_type reference_sqrt_density = CxPowerLaw<L>::sqrt_rho;
    value_type A = reference_sqrt_density * this->nonlinear / (0.6 * std::numbers::sqrt2_v<value_type>);
    value_type F = reference_viscosity * transition_reynum * std::sqrt(A);
    value_type pdrop = std::max(std::pow(F / this->nonlinear, 1.0 / this->exponent), dp_min);
    return reference_viscosity * F / (reference_density * pdrop);
  }

  virtual value_type correction(decltype(L::nodes[0]) node) const override
  {
    if (this->exponent == 0.5) {
      return 1.0;
    }
    return std::exp(std::log(this->reference_density / node->density) * (this->exponent - 1.0) +
                    std::log(this->reference_kinematic_viscosity * node->rdv) * (2 * this->exponent - 1.0));
  }

  virtual Type type() const override
  {
    return Type::FCN;
  }
};

template <typename L> struct CxDoor : public CxPowerLaw<L> // Single opening two-way flow element
{
  using value_type = Element<L>::value_type;

  constexpr static value_type sqrt2 = 1.414213562;
  constexpr static value_type two_thirds = 0.666667;
  constexpr static value_type gravity = 9.8;

  CxDoor(const std::string& name, value_type linear, value_type nonlinear, value_type exponent, value_type height, value_type width,
         value_type discharge_coefficient)
    : CxPowerLaw<L>(name, linear, nonlinear, exponent), height(height), width(width), discharge_coefficient(discharge_coefficient)
  {}

  CxDoor(const std::string& name, value_type exponent, value_type height, value_type width, value_type discharge_coefficient)
    : CxPowerLaw<L>(name, 0, 0, exponent), height(height), width(width), discharge_coefficient(discharge_coefficient)
  {
    auto A = height * width;
    auto d_h = 2.0 * A / (height + width);
    this->nonlinear = sqrt2 * discharge_coefficient * A;
    auto F = CxPowerLaw<L>::vsair * CxPowerLaw<L>::transition_reynum * A / d_h;
    auto dPt = std::pow(F / (this->nonlinear * CxPowerLaw<L>::sqrt_rho), 1.0 / this->exponent);
    dPt = std::max(dPt, CxPowerLaw<L>::min_pressure_drop);
    this->linear = F * this->reference_kinematic_viscosity / dPt;
  }

  virtual int calculate(bool const init,              // Initialization flag.If true, use linear relationship
                        const L& link,                // Linkage using this element
                        std::array<value_type, 2>& F, // Airflow through the component [kg/s]
                        std::array<value_type, 2>& DF // Partial derivative:  DF/DP
  )
  {
    auto& pn{ reference_to(link.nodes[0]) };
    auto& pm{ reference_to(link.nodes[1]) };
    int nf{ 1 };

    auto drho = pn.density - pm.density;
    if ((std::abs(drho) <= 0.0001 * std::abs(link.pressure_drop)) || init) {
      // Treat as orifice flow
      return CxPowerLaw<L>::calculate(init, link, F, DF);
    }
    
    auto gdrho = gravity * drho;
    auto Y = link.pressure_drop / gdrho;
    link.Y = Y;
    auto G = two_thirds * link.control * width * discharge_coefficient * std::sqrt(2.0 * abs(gdrho));
    auto c = std::abs(0.5 * height - Y);
    auto ft = c * std::sqrt(c);
    c = std::abs(0.5 * height + Y);
    auto fb = c * std::sqrt(c);
    if (Y < -0.5 * height) { // Door entirely below neutral plane
      if (gdrho > 0) {  // One-way flow (m to n)
        F[0] = -pm.sqrt_density * G * std::abs(ft - fb);
      } else {
        F[0] =  pn.sqrt_density * G * std::abs(ft - fb);
      }
    } else if (Y > 0.5 * height) { // Door entirely above neutral plane
      if (gdrho > 0) { // One-way flow (n to m)
        link.flow[0] =  pn.sqrt_density * G * abs(ft - fb);
      } else {
        link.flow[0] = -pm.sqrt_density * G * abs(ft - fb);
      }
    } else {
      nf = 2;
      if (gdrho > 0.0) { // Two-way flow
        F[0] = -pm.sqrt_density * G * ft;
        F[1] =  pn.sqrt_density * G * fb;
      } else {
        F[0] =  pn.sqrt_density * G * ft;
        F[1] = -pm.sqrt_density * G * fb;
      }
    }
    return nf;
  }

  virtual int calculate(const L& link,                // Linkage using this element
                        std::array<value_type, 2>& F, // Airflow through the component [kg/s]
                        std::array<value_type, 2>& DF // Partial derivative:  DF/DP
  )
  {
    auto& pn{ reference_to(link.nodes[0]) };
    auto& pm{ reference_to(link.nodes[1]) };
    int nf{ 1 };

    auto drho = pn.density - pm.density;
    if ((std::abs(drho) <= 0.0001 * std::abs(link.pressure_drop))) {
      // Treat as orifice flow
      return CxPowerLaw<L>::calculate(link, F, DF);
    }
    
    auto gdrho = gravity * drho;
    auto Y = link.pressure_drop / gdrho;
    link.Y = Y;
    auto G = two_thirds * link.control * width * discharge_coefficient * std::sqrt(2.0 * abs(gdrho));
    auto c = std::abs(0.5 * height - Y);
    auto ft = c * std::sqrt(c);
    c = std::abs(0.5 * height + Y);
    auto fb = c * std::sqrt(c);
    if (Y < -0.5 * height) { // Door entirely below neutral plane
      if (gdrho > 0) {  // One-way flow (m to n)
        F[0] = -pm.sqrt_density * G * std::abs(ft - fb);
      } else {
        F[0] =  pn.sqrt_density * G * std::abs(ft - fb);
      }
    } else if (Y > 0.5 * height) { // Door entirely above neutral plane
      if (gdrho > 0) { // One-way flow (n to m)
        link.flow[0] =  pn.sqrt_density * G * abs(ft - fb);
      } else {
        link.flow[0] = -pm.sqrt_density * G * abs(ft - fb);
      }
    } else {
      nf = 2;
      if (gdrho > 0.0) { // Two-way flow
        F[0] = -pm.sqrt_density * G * ft;
        F[1] =  pn.sqrt_density * G * fb;
      } else {
        F[0] =  pn.sqrt_density * G * ft;
        F[1] = -pm.sqrt_density * G * fb;
      }
    }
    return nf;
  }

  value_type height;
  value_type width;
  value_type discharge_coefficient;

};

template <typename L> struct Opening : public PowerLaw<L> // Specified volume flow element
{
  using value_type = Element<L>::value_type;

  Opening(const std::string& name, value_type linear, value_type nonlinear, value_type exponent, value_type area,
         value_type discharge_coefficient, value_type reference_density, value_type reference_kinematic_viscosity)
    : PowerLaw<L>(name, linear, nonlinear, exponent, reference_density, reference_kinematic_viscosity), discharge_coefficient(discharge_coefficient), area(area)
  {}

  Opening(const std::string& name, value_type area, value_type discharge_coefficient, value_type reference_density, value_type reference_kinematic_viscosity,
          value_type transition_reynum=30.0, value_type minimum_transition_pressure_drop=1.0e-10)
    : PowerLaw<L>(name, 0.0, std::numbers::sqrt2_v<value_type>* discharge_coefficient* area, 0.5, reference_density, reference_kinematic_viscosity),
    discharge_coefficient(discharge_coefficient), area(area)
  {
    value_type F = this->reference_kinematic_viscosity * this->reference_density * transition_reynum * std::sqrt(area);
    value_type dPt = F / (this->nonlinear * std::sqrt(this->reference_density));
    dPt *= dPt;
    dPt = std::max(dPt, minimum_transition_pressure_drop);
    this->linear = F * reference_kinematic_viscosity / dPt;
  }

  virtual Type type() const override
  {
    return Type::OPN;
  }

  value_type discharge_coefficient;
  value_type area; // Volume Flow [m2]
};

template <typename L> struct VariableOpening : public PowerLaw<L> // Specified volume flow element
{
  using value_type = Element<L>::value_type;

  VariableOpening(const std::string& name, value_type discharge_coefficient, value_type reference_density, value_type reference_kinematic_viscosity,
                  value_type transition_reynum = 30.0, value_type minimum_transition_pressure_drop = 1.0e-10)
    : PowerLaw<L>(name, 0.0, 0.0, 0.5, reference_density, reference_kinematic_viscosity), discharge_coefficient(discharge_coefficient),
    transition_reynum(transition_reynum), minimum_transition_pressure_drop(minimum_transition_pressure_drop)
  {}

  virtual int calculate(bool const init,              // Initialization flag.If true, use linear relationship
                        const L& link,                // Linkage using this element
                        std::array<value_type, 2>& F, // Airflow through the component [kg/s]
                        std::array<value_type, 2>& DF // Partial derivative:  DF/DP
  ) override
  {
    this->nonlinear = std::numbers::sqrt2_v<value_type> * this->discharge_coefficient * link.area;
    value_type w = this->reference_kinematic_viscosity * this->reference_density * this->transition_reynum * std::sqrt(link.area);
    value_type dPt = w / (this->nonlinear * std::sqrt(this->reference_density));
    dPt *= dPt;
    dPt = std::max(dPt, minimum_transition_pressure_drop);
    this->linear = w * this->reference_kinematic_viscosity / dPt;

    return PowerLaw<L>::calculate(init, link, F, DF);
  }

  virtual int calculate(const L& link,                // Linkage using this element
                        std::array<value_type, 2>& F, // Airflow through the component [kg/s]
                        std::array<value_type, 2>& DF // Partial derivative:  DF/DP
  ) override
  {
    this->nonlinear = std::numbers::sqrt2_v<value_type> *this->discharge_coefficient * link.area;
    value_type w = this->reference_kinematic_viscosity * this->reference_density * this->transition_reynum * std::sqrt(link.area);
    value_type dPt = w / (this->nonlinear * std::sqrt(this->reference_density));
    dPt *= dPt;
    dPt = std::max(dPt, minimum_transition_pressure_drop);
    this->linear = w * this->reference_kinematic_viscosity / dPt;

    return PowerLaw<L>::calculate(link, F, DF);
  }

  virtual Type type() const override
  {
    return Type::VOP;
  }

  value_type discharge_coefficient;
  value_type minimum_transition_pressure_drop;
  value_type transition_reynum;

};

} // namespace airflownetwork

#include "energyplus_elements.hpp"

#endif // !AIRFLOWNETWORK_ELEMENTS_HPP
