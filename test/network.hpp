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

#ifndef NETWORK_HPP
#include <memory>

#define BARPRES 101325.0    // standard air pressure
#define STDTEMP 293.15      // standard air temperature
#define AIR_R 287.055       // specific gas constant for air
#define RHOAIR 1.20410      // density of standard air
#define SQRT_RHO 1.097315   // sqrt(RHOAIR)
#define VSAIR 1.81625e-5    // dynamic viscosity of standard air
#define NUAIR 1.50839e-5   // kinematic viscosity (mu / rho) of standard air

struct Node
{
  void update_temperature(double P, double T = STDTEMP)
  {
    pressure = P;
    temperature = T;
    density = (P + BARPRES) / (AIR_R * T);
    sqrt_density = std::sqrt(density);
    viscosity = 3.7143e-6 + 4.9286E-8 * T;
    rdv = density / viscosity;
  }

  void update_density(double P, double rho = RHOAIR)
  {
    pressure = P;
    temperature = (P + BARPRES) / (AIR_R * rho);
    density = rho;
    sqrt_density = std::sqrt(density);
    viscosity = 3.7143e-6 + 4.9286E-8 * temperature;
    rdv = density / viscosity;
  }

  double temperature{ 20.0 };
  double pressure{ 0.0 };
  double humidity_ratio{ 0.0 };
  double density{ RHOAIR };
  //double density{airflownetwork::air_density<double,airflownetwork::Temperature::Celsius>(101325.0, 20.0, 0.0)};
  double sqrt_density{ SQRT_RHO };
  //double sqrt_density{std::sqrt(airflownetwork::air_density<double,airflownetwork::Temperature::Celsius>(101325.0, 20.0, 0.0))};
  double viscosity{ VSAIR };
  //double viscosity{airflownetwork::air_dynamic_viscosity<double,airflownetwork::Temperature::Celsius>(20.0)};
  double rdv{ 1.0/NUAIR };
};

struct BaseLink
{
  double control = 1.0;
  double multiplier = 1.0;
  double pressure_drop = 0.0;
  double area = 0.0;
  std::optional<double> neutral_height;
};

struct PtrLink : BaseLink
{
  PtrLink(Node* n0, Node* n1) : nodes({ n0, n1 })
  {}
  std::array<Node*, 2> nodes = { {nullptr, nullptr} };
};

struct SharedLink : BaseLink
{
  SharedLink(std::shared_ptr<Node>& n0, std::shared_ptr<Node>& n1) : nodes({ n0, n1 })
  {}
  std::array<std::shared_ptr<Node>, 2> nodes = {{nullptr, nullptr}};
};

#endif
