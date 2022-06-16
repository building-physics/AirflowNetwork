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

#ifdef __cpp_modules
#ifdef _MSC_VER
#define BOOST_UT_DISABLE_MODULE
#include <boost/ut.hpp>
#else
import boost.ut; // Doesn't appear to work yet with MSVC/CMake
#endif
#else
#include <boost/ut.hpp>
#endif

#include <cmath>
#include <optional>
#include "airflownetwork/elements.hpp"

#include "network.hpp"

boost::ut::suite two_way_flow_elements = [] {
  using namespace boost::ut;

  "simple_opening"_test = [] {
    int NF;
    std::array<double, 2> F = { 0.0, 0.0 };
    std::array<double, 2> DF = { 0.0, 0.0 };

    double reference_density = airflownetwork::air_density<double, airflownetwork::Temperature::Celsius>(101325.0, 20.0, 0.0);
    double reference_viscosity = airflownetwork::air_dynamic_viscosity<double, airflownetwork::Temperature::Celsius>(20.0);

    airflownetwork::Crack<PtrLink, airflownetwork::Temperature::Celsius> element("element", 0.001, 0.65, reference_density, reference_viscosity);

    Node state0, state1;
    state0.density = state1.density = reference_density;
    double viscosity = state0.viscosity = state1.viscosity = reference_viscosity;
    double sqrt_density = state0.sqrt_density = state1.sqrt_density = std::sqrt(reference_density);

    PtrLink link(&state0, &state1);

    double dp{ 10.0 };
    double dp_tiny{ 1.0e-14 };

    // Linear
    link.pressure_drop = dp;
    NF = element.calculate(true, link, F, DF);
    expect(near(0.01 * sqrt_density / viscosity, F[0], 1.0e-12));
    expect(eq(0.0, F[1]));
    expect(near(0.001 * sqrt_density / viscosity, DF[0], 1.0e-12));
    expect(eq(0.0, DF[1]));
    expect(eq(1, NF));

    link.pressure_drop = -dp;
    NF = element.calculate(true, link, F, DF);
    expect(near(-0.01 * sqrt_density / viscosity, F[0], 1.0e-12));
    expect(eq(0.0, F[1]));
    expect(near(0.001 * sqrt_density / viscosity, DF[0], 1.0e-12));
    expect(eq(0.0, DF[1]));
    expect(eq(1, NF));

    // Nonlinear
    link.pressure_drop = dp;
    NF = element.calculate(false, link, F, DF);
    expect(eq(0.001 * std::pow(10.0, 0.65), F[0]));
    expect(eq(0.0, F[1]));
    expect(near(0.000065 * std::pow(10.0, 0.65), DF[0], 1.0e-12));
    expect(eq(0.0, DF[1]));
    expect(eq(1, NF));

    link.pressure_drop = -dp;
    NF = element.calculate(false, link, F, DF);
    expect(eq(-0.001 * std::pow(10.0, 0.65), F[0]));
    expect(eq(0.0, F[1]));
    expect(near(0.000065 * std::pow(10.0, 0.65), DF[0], 1.0e-12));
    expect(eq(0.0, DF[1]));
    expect(eq(1, NF));

    // Linear, alternate call
    link.pressure_drop = dp_tiny;
    NF = element.calculate(link, F, DF);
    expect(near(0.001 * dp_tiny * sqrt_density / viscosity, F[0], 1.0e-12));
    expect(eq(0.0, F[1]));
    expect(near(0.001 * sqrt_density / viscosity, DF[0], 1.0e-12));
    expect(eq(0.0, DF[1]));
    expect(eq(1, NF));

    link.pressure_drop = -dp_tiny;
    NF = element.calculate(link, F, DF);
    expect(near(-0.001 * dp_tiny * sqrt_density / viscosity, F[0], 1.0e-12));
    expect(eq(0.0, F[1]));
    expect(near(0.001 * sqrt_density / viscosity, DF[0], 1.0e-12));
    expect(eq(0.0, DF[1]));
    expect(eq(1, NF));

    // Nonlinear, alternate call
    link.pressure_drop = dp;
    NF = element.calculate(link, F, DF);
    expect(eq(0.001 * std::pow(10.0, 0.65), F[0]));
    expect(eq(0.0, F[1]));
    expect(near(0.000065 * std::pow(10.0, 0.65), DF[0], 1.0e-12));
    expect(eq(0.0, DF[1]));
    expect(eq(1, NF));

    link.pressure_drop = -dp;
    NF = element.calculate(link, F, DF);
    expect(eq(-0.001 * std::pow(10.0, 0.65), F[0]));
    expect(eq(0.0, F[1]));
    expect(near(0.000065 * std::pow(10.0, 0.65), DF[0], 1.0e-12));
    expect(eq(0.0, DF[1]));
    expect(eq(1, NF));
  };

  "CONTAM door"_test = [] {
    int NF;
    std::array<double, 2> F = { 0.0, 0.0 };
    std::array<double, 2> DF = { 0.0, 0.0 };

    auto node0 = std::make_shared< Node>();
    auto node1 = std::make_shared< Node>();

    SharedLink link(node0, node1);

    airflownetwork::PowerLaw<SharedLink> element("element", 3.07052e-05, 0.01, 0.5, RHOAIR, NUAIR);

    // Nonlinear flow
    link.pressure_drop = 5.0;
    NF = element.calculate(link, F, DF);
    expect(eq(1, NF));
    expect(near(0.0245366810898831, F[0], 1.0e-7)); // Contam 3.4.0.1
    expect(eq(0.02453670932730182, F[0])); // Python
    expect(eq(0.0, F[1]));

    link.pressure_drop = -5.0;
    NF = element.calculate(link, F, DF);
    expect(eq(1, NF));
    expect(near(-0.0245366810898831, F[0], 1.0e-7)); // Contam 3.4.0.1
    expect(eq(-0.02453670932730182, F[0])); // Python
    expect(eq(0.0, F[1]));

    // Linear flow
    link.pressure_drop = 0.00001;
    NF = element.calculate(link, F, DF);
    expect(eq(1, NF));
    expect(near(2.03562659843231e-05, F[0], 1.0e-9)); // Contam 3.4.0.1
    expect(near(2.035630079559532e-05, F[0], 1.0e-10)); // Python
    expect(eq(0.0, F[1]));

    link.pressure_drop = -0.00001;
    NF = element.calculate(link, F, DF);
    expect(eq(1, NF));
    expect(near(-2.03562659843231e-05, F[0], 1.0e-9)); // Contam 3.4.0.1
    expect(near(-2.035630079559532e-05, F[0], 1.0e-10)); // Python
    expect(eq(0.0, F[1]));

    // Nonlinear flow
    link.pressure_drop = 5.0;
    NF = element.calculate(false, link, F, DF);
    expect(eq(1, NF));
    expect(near(0.0245366810898831, F[0], 1.0e-7)); // Contam 3.4.0.1
    expect(eq(0.02453670932730182, F[0])); // Python
    expect(eq(0.0, F[1]));

    link.pressure_drop = -5.0;
    NF = element.calculate(false, link, F, DF);
    expect(eq(1, NF));
    expect(near(-0.0245366810898831, F[0], 1.0e-7)); // Contam 3.4.0.1
    expect(eq(-0.02453670932730182, F[0])); // Python
    expect(eq(0.0, F[1]));

    // Linear flow
    link.pressure_drop = 0.00001;
    NF = element.calculate(false, link, F, DF);
    expect(eq(1, NF));
    expect(near(2.03562659843231e-05, F[0], 1.0e-9)); // Contam 3.4.0.1
    expect(near(2.035630079559532e-05, F[0], 1.0e-10)); // Python
    expect(eq(0.0, F[1]));

    link.pressure_drop = -0.00001;
    NF = element.calculate(false, link, F, DF);
    expect(eq(1, NF));
    expect(near(-2.03562659843231e-05, F[0], 1.0e-9)); // Contam 3.4.0.1
    expect(near(-2.035630079559532e-05, F[0], 1.0e-10)); // Python
    expect(eq(0.0, F[1]));

    link.pressure_drop = 0.00001;
    NF = element.calculate(true, link, F, DF);
    expect(eq(1, NF));
    expect(near(2.03562659843231e-05, F[0], 1.0e-9)); // Contam 3.4.0.1
    expect(near(2.035630079559532e-05, F[0], 1.0e-10)); // Python
    expect(eq(0.0, F[1]));

    link.pressure_drop = -0.00001;
    NF = element.calculate(true, link, F, DF);
    expect(eq(1, NF));
    expect(near(-2.03562659843231e-05, F[0], 1.0e-9)); // Contam 3.4.0.1
    expect(near(-2.035630079559532e-05, F[0], 1.0e-10)); // Python
    expect(eq(0.0, F[1]));
  };

  "new simple opening"_test = [] {
    int NF;
    std::array<double, 2> F = { 0.0, 0.0 };
    std::array<double, 2> DF = { 0.0, 0.0 };

    auto node0 = std::make_shared< Node>();
    auto node1 = std::make_shared< Node>();

    SharedLink link(node0, node1);

    airflownetwork::CxPowerLaw<SharedLink> element("element", 3.07052e-05, 0.01, 0.5);

    // Nonlinear flow
    link.pressure_drop = 5.0;
    NF = element.calculate(link, F, DF);
    expect(eq(1, NF));
    expect(near(0.0245366810898831, F[0], 1.0e-7)); // Contam 3.4.0.1
    expect(eq(0.02453670932730182, F[0])); // Python
    expect(eq(0.0, F[1]));

    link.pressure_drop = -5.0;
    NF = element.calculate(link, F, DF);
    expect(eq(1, NF));
    expect(near(-0.0245366810898831, F[0], 1.0e-7)); // Contam 3.4.0.1
    expect(eq(-0.02453670932730182, F[0])); // Python
    expect(eq(0.0, F[1]));

    // Linear flow
    link.pressure_drop = 0.00001;
    NF = element.calculate(link, F, DF);
    expect(eq(1, NF));
    expect(near(2.03562659843231e-05, F[0], 1.0e-9));
    expect(eq(0.0, F[1]));

    link.pressure_drop = -0.00001;
    NF = element.calculate(link, F, DF);
    expect(eq(1, NF));
    expect(near(-2.03562659843231e-05, F[0], 1.0e-9));
    expect(eq(0.0, F[1]));

    // Nonlinear flow
    link.pressure_drop = 5.0;
    NF = element.calculate(false, link, F, DF);
    expect(eq(1, NF));
    expect(near(0.0245366810898831, F[0], 1.0e-7)); // Contam 3.4.0.1
    expect(eq(0.02453670932730182, F[0])); // Python
    expect(eq(0.0, F[1]));

    link.pressure_drop = -5.0;
    NF = element.calculate(false, link, F, DF);
    expect(eq(1, NF));
    expect(near(-0.0245366810898831, F[0], 1.0e-7)); // Contam 3.4.0.1
    expect(eq(-0.02453670932730182, F[0])); // Python
    expect(eq(0.0, F[1]));

    // Linear flow
    link.pressure_drop = 0.00001;
    NF = element.calculate(false, link, F, DF);
    expect(eq(1, NF));
    expect(near(2.03562659843231e-05, F[0], 1.0e-9));
    expect(eq(0.0, F[1]));

    link.pressure_drop = -0.00001;
    NF = element.calculate(false, link, F, DF);
    expect(eq(1, NF));
    expect(near(-2.03562659843231e-05, F[0], 1.0e-9));
    expect(eq(0.0, F[1]));

    link.pressure_drop = 0.00001;
    NF = element.calculate(true, link, F, DF);
    expect(eq(1, NF));
    expect(near(2.03562659843231e-05, F[0], 1.0e-9));
    expect(eq(0.0, F[1]));

    link.pressure_drop = -0.00001;
    NF = element.calculate(true, link, F, DF);
    expect(eq(1, NF));
    expect(near(-2.03562659843231e-05, F[0], 1.0e-9));
    expect(eq(0.0, F[1]));
  };
};
