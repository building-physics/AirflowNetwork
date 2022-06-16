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

boost::ut::suite specified_flow_elements = [] {
  using namespace boost::ut;

  "specified mass flow"_test = [] {
    int NF;
    std::array<double, 2> F = { 0.0, 0.0 };
    std::array<double, 2> DF = { 0.0, 0.0 };

    auto node0 = std::make_shared< Node>();
    auto node1 = std::make_shared< Node>();

    SharedLink link(node0, node1);

    double f = 0.1;
    double dp = 10.0;

    airflownetwork::SpecifiedMassFlow<SharedLink> element("element", f);

    link.pressure_drop = dp;
    NF = element.calculate(link, F, DF);
    expect(eq(f, F[0]));
    expect(eq(0.0, F[1]));
    expect(eq(0.0, DF[0]));
    expect(eq(0.0, DF[1]));
    expect(eq(1, NF));

    link.pressure_drop = -dp;
    NF = element.calculate(link, F, DF);
    expect(eq(f, F[0]));
    expect(eq(0.0, F[1]));
    expect(eq(0.0, DF[0]));
    expect(eq(0.0, DF[1]));
    expect(eq(1, NF));

    // Linear
    link.pressure_drop = dp;
    NF = element.calculate(true, link, F, DF);
    expect(eq(f, F[0]));
    expect(eq(0.0, F[1]));
    expect(eq(0.0, DF[0]));
    expect(eq(0.0, DF[1]));
    expect(eq(1, NF));

    link.pressure_drop = -dp;
    NF = element.calculate(true, link, F, DF);
    expect(eq(f, F[0]));
    expect(eq(0.0, F[1]));
    expect(eq(0.0, DF[0]));
    expect(eq(0.0, DF[1]));
    expect(eq(1, NF));

    // Nonlinear
    link.pressure_drop = dp;
    NF = element.calculate(false, link, F, DF);
    expect(eq(f, F[0]));
    expect(eq(0.0, F[1]));
    expect(eq(0.0, DF[0]));
    expect(eq(0.0, DF[1]));
    expect(eq(1, NF));

    link.pressure_drop = -dp;
    NF = element.calculate(false, link, F, DF);
    expect(eq(f, F[0]));
    expect(eq(0.0, F[1]));
    expect(eq(0.0, DF[0]));
    expect(eq(0.0, DF[1]));
    expect(eq(1, NF));

    link.pressure_drop = dp;
    NF = element.calculate(link, F, DF);
    expect(eq(f, F[0]));
    expect(eq(0.0, F[1]));
    expect(eq(0.0, DF[0]));
    expect(eq(0.0, DF[1]));
    expect(eq(1, NF));

    link.pressure_drop = -dp;
    NF = element.calculate(link, F, DF);
    expect(eq(f, F[0]));
    expect(eq(0.0, F[1]));
    expect(eq(0.0, DF[0]));
    expect(eq(0.0, DF[1]));
    expect(eq(1, NF));

    f *= -1;
    element.mass_flow = f;
    // Linear
    link.pressure_drop = dp;
    NF = element.calculate(true, link, F, DF);
    expect(eq(f, F[0]));
    expect(eq(0.0, F[1]));
    expect(eq(0.0, DF[0]));
    expect(eq(0.0, DF[1]));
    expect(eq(1, NF));

    link.pressure_drop = -dp;
    NF = element.calculate(true, link, F, DF);
    expect(eq(f, F[0]));
    expect(eq(0.0, F[1]));
    expect(eq(0.0, DF[0]));
    expect(eq(0.0, DF[1]));
    expect(eq(1, NF));

    // Nonlinear
    link.pressure_drop = dp;
    NF = element.calculate(false, link, F, DF);
    expect(eq(f, F[0]));
    expect(eq(0.0, F[1]));
    expect(eq(0.0, DF[0]));
    expect(eq(0.0, DF[1]));
    expect(eq(1, NF));

    link.pressure_drop = -dp;
    NF = element.calculate(false, link, F, DF);
    expect(eq(f, F[0]));
    expect(eq(0.0, F[1]));
    expect(eq(0.0, DF[0]));
    expect(eq(0.0, DF[1]));
    expect(eq(1, NF));
  };

  "specified volume flow"_test = [] {
    int NF;
    std::array<double, 2> F = { 0.0, 0.0 };
    std::array<double, 2> DF = { 0.0, 0.0 };

    auto node0 = std::make_shared< Node>();
    auto node1 = std::make_shared< Node>();

    SharedLink link(node0, node1);

    double f = 0.1 * node0->density;
    double dp = 10.0;

    airflownetwork::SpecifiedVolumeFlow<SharedLink> element("element", 0.1);

    link.pressure_drop = dp;
    NF = element.calculate(link, F, DF);
    expect(eq(f, F[0]));
    expect(eq(0.0, F[1]));
    expect(eq(0.0, DF[0]));
    expect(eq(0.0, DF[1]));
    expect(eq(1, NF));

    link.pressure_drop = -dp;
    NF = element.calculate(link, F, DF);
    expect(eq(f, F[0]));
    expect(eq(0.0, F[1]));
    expect(eq(0.0, DF[0]));
    expect(eq(0.0, DF[1]));
    expect(eq(1, NF));

    // Linear
    link.pressure_drop = dp;
    NF = element.calculate(true, link, F, DF);
    expect(eq(f, F[0]));
    expect(eq(0.0, F[1]));
    expect(eq(0.0, DF[0]));
    expect(eq(0.0, DF[1]));
    expect(eq(1, NF));

    link.pressure_drop = -dp;
    NF = element.calculate(true, link, F, DF);
    expect(eq(f, F[0]));
    expect(eq(0.0, F[1]));
    expect(eq(0.0, DF[0]));
    expect(eq(0.0, DF[1]));
    expect(eq(1, NF));

    // Nonlinear
    link.pressure_drop = dp;
    NF = element.calculate(false, link, F, DF);
    expect(eq(f, F[0]));
    expect(eq(0.0, F[1]));
    expect(eq(0.0, DF[0]));
    expect(eq(0.0, DF[1]));
    expect(eq(1, NF));

    link.pressure_drop = -dp;
    NF = element.calculate(false, link, F, DF);
    expect(eq(f, F[0]));
    expect(eq(0.0, F[1]));
    expect(eq(0.0, DF[0]));
    expect(eq(0.0, DF[1]));
    expect(eq(1, NF));

    link.pressure_drop = dp;
    NF = element.calculate(link, F, DF);
    expect(eq(f, F[0]));
    expect(eq(0.0, F[1]));
    expect(eq(0.0, DF[0]));
    expect(eq(0.0, DF[1]));
    expect(eq(1, NF));

    link.pressure_drop = -dp;
    NF = element.calculate(link, F, DF);
    expect(eq(f, F[0]));
    expect(eq(0.0, F[1]));
    expect(eq(0.0, DF[0]));
    expect(eq(0.0, DF[1]));
    expect(eq(1, NF));

    f *= -1;
    element.volume_flow = -0.1;
    // Linear
    link.pressure_drop = dp;
    NF = element.calculate(true, link, F, DF);
    expect(eq(f, F[0]));
    expect(eq(0.0, F[1]));
    expect(eq(0.0, DF[0]));
    expect(eq(0.0, DF[1]));
    expect(eq(1, NF));

    link.pressure_drop = -dp;
    NF = element.calculate(true, link, F, DF);
    expect(eq(f, F[0]));
    expect(eq(0.0, F[1]));
    expect(eq(0.0, DF[0]));
    expect(eq(0.0, DF[1]));
    expect(eq(1, NF));

    // Nonlinear
    link.pressure_drop = dp;
    NF = element.calculate(false, link, F, DF);
    expect(eq(f, F[0]));
    expect(eq(0.0, F[1]));
    expect(eq(0.0, DF[0]));
    expect(eq(0.0, DF[1]));
    expect(eq(1, NF));

    link.pressure_drop = -dp;
    NF = element.calculate(false, link, F, DF);
    expect(eq(f, F[0]));
    expect(eq(0.0, F[1]));
    expect(eq(0.0, DF[0]));
    expect(eq(0.0, DF[1]));
    expect(eq(1, NF));
  };

};
