# AirflowNetwork

AirflowNetwork is a C++ implementation of the pressure network approach to bulk air movement and contaminant transport, primarily in the context of buildings. It is a
descendant of NIST AIRNET (the precursor to NIST's CONTAM) and COMIS. This modeling technique divides the building's total volume up into a set of well-mixed zones
that communicate via airflows through openings that are part of the building. Parts of this code have been extracted from EnergyPlus (U.S. DOE's whole building energy
simulation program) and other parts are new development. At present, the code in this repository is focused on the components that go into the pressure network
calculation but does not provide the means to solve the resulting equations. This will change over time.

## Background

The pressure network approach is formulated by considering the transport problem as a network of nodes (that represent a volume of air in a building) connected by a
set of paths or linkages (representing cracks, windows, doors, and other features of buildings that air flows through). The governing principle is conservation of mass,
which in the steady state case is written as
$$\sum_j F_{i,j} = 0,$$
where $F_{i,j}$ is the mass flow between nodes $i$ and $j$ in the network. This is a system of equations when it is written for each node in the network. The typical
approach is to consider the flow as pressure driven, so
$$ F_{i,j} = f\left(\Delta P_{i,j}\right)$$
where $\Delta P_{i,j}$ is the pressure difference between nodes $i$ and $j$ and $f$ is a nonlinear function. In reality, temperature also plays a role (especially in
natural ventilation) and this is accounted for in the calculation of $\Delta P_{i,j}$ and $f$. Often, a multivariate Newton-Raphson iteration is used to solve the
system of equations, in which case the derivatives of $f$ must be calculated. A popular form of $f$ is the so-called power law formulation:
$$ F_{i,j} = C\sqrt{\rho}\Delta P_{i,j}^n.$$
The coefficient $C$ is a flow parameter and the exponent $n$ is typically in the neighborhood of 0.5 and 0.65. Calculation of the derivative for this form is fairly
easy (except near zero), and once the Jacobian matrix is formed, the solution can be iteratively solved for.

## Building

AirflowNetwork's build is CMake based. At present it only builds on Windows, which is another thing that will change.

## Testing

At present, the only thing that can be built are the unit test. These are implemented using the [boost::ut](https://github.com/boost-ext/ut) microframework.
