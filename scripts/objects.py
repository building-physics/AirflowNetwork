# Copyright (c) 2020-2022, Oak Ridge National Laboratory, managed by UT-Battelle.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted
# provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of
#    conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of
#    conditions and the following disclaimer in the documentation and/or other materials
#    provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to
#    endorse or promote products derived from this software without specific prior written
#    permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
# AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

import math

BARPRES = 101325.0  # standard air pressure
STDTEMP = 293.15    # standard air temperature
AIR_R = 287.055     # specific gas constant for air
#SQRT2 = 1.414213562 # sqrt(2.0)
RHOAIR = 1.20410    # density of standard air
SQRT_RHO = 1.097315 # sqrt(RHOAIR)
VSAIR = 1.81625e-5  # dynamic viscosity of standard air
DVSAIR = 1.50839e-5 # kinematic viscosity (mu/rho) of standard air
MUAIR = 1.81625e-5  # dynamic viscosity of standard air
RE_LT = 30.0        # Reynolds number at L-T transition
DPT_MIN = 1.0e-10   # Minimum transition pressure difference

#TWO_THIRDS = 2.0/3.0
#SQRT2 = math.sqrt(2.0)
GRAVITY = 9.8


class Node:
    def __init__(self, nr=0):
        self.nr = nr
        self.pressure = 0.0                    # pressure (gage)
        self.temperature = STDTEMP             # temperature
        self.density = RHOAIR                  # density
        self.sqrt_density = SQRT_RHO           # sqrt(density)
        self.viscosity = VSAIR                 # dynamic viscosity
        self.rdv = self.density/self.viscosity # 1/nu or rho/mu
    def update(self, P, T=STDTEMP, rho=None):
        self.pressure = P
        if rho is None:
            self.temperature = T
            self.density = (P + BARPRES) / (AIR_R * T)
            self.sqrt_density = math.sqrt(self.density)
            self.viscosity = 3.7143e-6 + 4.9286E-8 * T
            self.rdv = self.density / self.viscosity
        else:
            self.temperature = (P + BARPRES) / (AIR_R * rho)
            self.density = rho
            self.sqrt_density = math.sqrt(self.density)
            self.viscosity = 3.7143e-6 + 4.9286E-8 * self.temperature
            self.rdv = self.density / self.viscosity


class Path:
    def __init__(self, nr, pn, pm, pe):
        self.element = pe        # element
        self.pressure_drop = 0.0 # pressure difference
        self.control = 1.0       # control value
        self.multiplier = 1.0    # multiplier
        self.nr = nr             # path number
        self.pn = pn             # node_0
        self.pm = pm             # node_1
        self.flow = [0.0, 0.0]   # Flow
        self.dFdP = 0.0          # dFdP
        self.DF = [0.0, 0.0]     # dFdP
        self.Y = None            # neutral height, if computed


class CxPlr:
    def __init__(self, nr, lam, turb, expt=0.65):
        self.nr = nr           # element number
        self.linear = lam      # linear flow coeff
        self.exponent = expt   # exponent
        self.coefficient = turb # nonlinear flow coeff

    def correction(self, node):
        # (RHOAIR/dens)^(expt-0.5) * (DVSAIR*Dvisc)^(2*expt-1) * sqrt(dens)
        if self.exponent == 0.5:
            return node.sqrt_density
        return math.exp(math.log( RHOAIR / node.density ) * ( self.exponent - 0.5 ) + 
                        math.log( DVSAIR * node.rdv) * ( 2 * self.exponent - 1.0 )) * node.sqrt_density

    def calculate(self, path):
        if path.pressure_drop >= 0:
            upwind_node = path.pn
        else:
            upwind_node = path.pm

        ctrl = path.control * path.multiplier

        cdm = self.linear * upwind_node.rdv
        Fl = cdm * path.pressure_drop

        Tadj = self.correction(upwind_node)

        if path.pressure_drop >= 0:
            Ft = Tadj * self.coefficient * math.pow(path.pressure_drop, self.exponent)
        else:
            Ft = -Tadj * self.coefficent * math.pow(-path.pressure_drop, self.exponent)

        if abs(Fl) <= abs(Ft):
            F = Fl
            dF = cdm
        else:
            F = Ft
            dF = Ft * self.exponent / path.pressure_drop

        path.flow[0] = F * ctrl
        path.dFdP = dF * ctrl
        path.flow[1] = 0
        return 1


class PowerLaw(CxPlr):
    def __init__(self, nr, lam, turb, expt=0.65, reference_conditions=None):
        super().__init__(nr, lam, turb, expt)
        if reference_conditions is not None:
            self.reference_density = reference_conditions.density
            self.reference_dynamic_viscosity = reference_conditions.viscosity / reference_conditions.density
        else:
            self.reference_density = RHOAIR
            self.reference_dynamic_viscosity = DVSAIR

    def correction(self, node):
        # (RHOAIR/dens)^(expt-0.5) * (DVSAIR*Dvisc)^(2*expt-1) * sqrt(dens)
        if self.exponent == 0.5:
            return node.sqrt_density
        return math.exp(math.log( self.reference_density / node.density ) * ( self.exponent - 0.5 ) + 
                        math.log( self.reference_dynamic_viscosity* node.rdv) * ( 2 * self.exponent - 1.0 )) * node.sqrt_density

    def calculate(self, path):
        upwind_node = path.pn
        sign = 1.0
        abs_pdrop = path.pressure_drop
        if path.pressure_drop < 0:
            upwind_node = path.pm
            sign = -1.0
            abs_pdrop = -path.pressure_drop

        ctrl = path.control * path.multiplier

        cdm = self.linear * upwind_node.rdv
        abs_Fl = cdm * abs_pdrop

        Tadj = self.correction(upwind_node)

        abs_Ft = Tadj * self.coefficient * math.pow(abs_pdrop, self.exponent)

        if abs_Fl <= abs_Ft:
            path.flow[0] = sign * abs_Fl * ctrl
            path.dFdP = cdm * ctrl
        else:
            path.flow[0] = sign * abs_Ft * ctrl
            path.dFdP = path.flow[0] * self.exponent / path.pressure_drop

        path.flow[1] = 0
        return 1

    @classmethod
    def fit(cls, dp0, F0, dp1=None, F1=None, exponent=0.65, reference_conditions=None):
        dp0 = abs(dp0)
        F0 = abs(F0)
        if reference_conditions is None:
            reference_conditions = Node()
        if dp1 is None or F1 is None:
            # Single point
            C = F0 / (reference_conditions.sqrt_density * math.pow(dp0, exponent))
        else:
            # Two point
            dp1 = abs(dp1)
            F1 = abs(F1)
            exponent = (math.log(F0) - math.log(F1)) / (math.log(dp0) - math.log(dp1))
            C = F0 / (reference_conditions.sqrt_density * math.pow(dp0, exponent))


def orifice_correction(node, exponent, reference_density, reference_dynamic_viscosity):
        if exponent == 0.5:
            return node.sqrt_density
        return math.exp(math.log( reference_density / node.density ) * ( exponent - 0.5 ) + 
                        math.log( reference_dynamic_viscosity* node.rdv) * ( 2 * exponent - 1.0 )) * node.sqrt_density


def orifice_flow(linear, nonlinear, exponent, path, reference_density=RHOAIR, reference_dynamic_viscosity=DVSAIR):
    upwind_node = path.pn
    sign = 1.0
    abs_pdrop = path.pressure_drop
    if path.pressure_drop < 0:
        upwind_node = path.pm
        sign = -1.0
        abs_pdrop = -path.pressure_drop
    ctrl = path.control * path.multiplier
    cdm = linear * upwind_node.rdv
    abs_Fl = cdm * abs_pdrop
    abs_Ft = orifice_correction(upwind_node, exponent, reference_density, reference_dynamic_viscosity) * nonlinear * math.pow(abs_pdrop, exponent)
    if abs_Fl <= abs_Ft:
        path.flow[0] = sign * abs_Fl * ctrl
        path.dFdP = cdm * ctrl
    else:
        path.flow[0] = sign * abs_Ft * ctrl
        path.dFdP = path.flow[0] * exponent / path.pressure_drop
    path.flow[1] = 0
    return 1


def orifice_coefficients(ht=2.0, wd=0.8, cd=0.78, exponent=0.5, Re_t=30.0, reference_density=RHOAIR, reference_sqrt_density=SQRT_RHO,
                         reference_viscosity=MUAIR):
    A = ht * wd
    d_h = 2.0 * A / ( ht + wd )
    nonlinear = math.sqrt(2.0) * cd * A
    F = reference_viscosity * Re_t * A / d_h
    dPt = math.pow(F / (nonlinear * reference_sqrt_density), 1.0/exponent)
    dPt = max(dPt, DPT_MIN)
    linear = F * reference_viscosity / (reference_density * dPt)
    return linear, nonlinear


class CxFcn(CxPlr):
    def correction(self, node):
        # (RHOAIR/dens)^(expt-1) * (DVSAIR*Dvisc)^(2*expt-1)
        if self.exponent == 0.5:
            return 1.0
        return math.exp(math.log( RHOAIR / node.density ) * ( self.exponent - 1.0 ) + 
                        math.log( DVSAIR * node.rdv) * ( 2 * self.exponent - 1.0 ))


class CxQcn(CxPlr):
    def correction(self, node):
        # (RHOAIR/dens)^expt * (DVSAIR*dvisc)^(2*expt-1) * density
        if self.exponent == 0.5:
            return node.density
        return math.exp( math.log( RHOAIR / node.density ) * self.exponent +
                         math.log( DVSAIR * node.rdv ) * ( 2 * self.exponent - 1.0 ) ) * node.density

def pow1p5(c):
    return c * math.sqrt(c)


class CxDoor(CxPlr):
    SQRT2 = 1.414213562 # sqrt(2.0)
    TWO_THIRDS = 0.666667
    def __init__(self, ht=2.0, wd=0.8, cd=0.78, exponent=0.5):
        self.height = ht
        self.width = wd
        self.cd = cd
        self.linear, self.coefficient = type(self).compute_coefficients(ht=ht, wd=wd, cd=cd, exponent=exponent)
    def calculate(self, path):
        pn = path.pn
        pm = path.pm
        drho = pn.density - pm.density
        if abs(drho) <= 0.0001 * abs(path.pressure_drop):
            # Treat as orifice flow
            return CxPlr.calculate(path)
        gdrho = GRAVITY * drho
        ctrl = 1.0
        Y = path.pressure_drop / gdrho
        path.Y = Y
        G = self.TWO_THIRDS * ctrl * self.width * self.cd * math.sqrt(2.0*abs(gdrho))
        c = abs(0.5 * self.height - Y)
        ft = c * math.sqrt(c)
        c = abs(0.5 * self.height + Y)
        fb = c * math.sqrt(c)
        if Y < -0.5*self.height:  # Door entirely below neutral plane
            if gdrho > 0:  # One-way flow (m to n)
                path.flow[0] = -pm.sqrt_density * G * abs(ft-fb)
            else:
                path.flow[0] =  pn.sqrt_density * G * abs(ft-fb)
            path.flow[1] = 0.0
            nf = 1
        elif Y > 0.5*self.height:  # Door entirely above neutral plane
            if gdrho > 0:  # One-way flow (n to m)
                path.flow[0] =  pn.sqrt_density * G * abs(ft-fb)
            else:
                path.flow[0] = -pm.sqrt_density * G * abs(ft-fb)
            path.flow[1] = 0
            nf = 1
        else:
            nf = 2
            if gdrho > 0.0:  # Two-way flow
                path.flow[0] = -pm.sqrt_density * G * ft
                path.flow[1] =  pn.sqrt_density * G * fb
            else:
                path.flow[0] =  pn.sqrt_density * G * ft
                path.flow[1] = -pm.sqrt_density * G * fb
        return nf
    @classmethod
    def compute_coefficients(cls, ht=2.0, wd=0.8, cd=0.78, exponent=0.5, Re_t=30.0, reference_conditions=None):
        if reference_conditions is None:
            reference_conditions = Node()
        A = ht * wd
        d_h = 2.0 * A / ( ht + wd )
        nonlinear = cls.SQRT2 * cd * A
        F = reference_conditions.viscosity * Re_t * A / d_h
        dPt = math.pow(F / (nonlinear * reference_conditions.sqrt_density), 1.0/exponent)
        dPt = max(dPt, DPT_MIN)
        linear = F * reference_conditions.viscosity / (reference_conditions.density * dPt)
        return linear, nonlinear


class AltDoor(CxDoor):
    TWO_THIRDS = 2.0/3.0
    def calculate(self, path):
        pn = path.pn
        pm = path.pm
        drho = pn.density - pm.density
        if abs(drho) <= 0.0001 * abs(path.pressure_drop):
            # Treat as orifice flow
            pass
        gdrho = GRAVITY * drho
        ctrl = 1.0
        Y = path.pressure_drop / gdrho
        path.Y = Y
        G = self.TWO_THIRDS * ctrl * self.width * self.cd * math.sqrt(2.0*abs(gdrho))
        if Y > 0.5*self.height:
            nf = 1
            absF = G*(pow1p5(Y + 0.5*self.height) - pow1p5(Y - 0.5*self.height))
            if gdrho > 0.0:
                path.flow[0] =  pn.sqrt_density * absF
            else:
                path.flow[0] = -pm.sqrt_density * absF
        elif Y < -0.5*self.height:
            nf = 1
            absF = G*(pow1p5(0.5*self.height - Y) - pow1p5(-Y - 0.5*self.height))
            if gdrho > 0.0:
                path.flow[0] = -pm.sqrt_density * absF
            else:
                path.flow[0] =  pn.sqrt_density * absF
        else:  # Bidirectional flow
            nf = 2
            if gdrho > 0.0:
                path.flow[0] = -pm.sqrt_density * G * pow1p5(0.5*self.height - Y)
                path.flow[1] =  pn.sqrt_density * G * pow1p5(Y + 0.5*self.height)
            else:
                path.flow[0] =  pn.sqrt_density * G * pow1p5(0.5*self.height - Y)
                path.flow[1] = -pm.sqrt_density * G * pow1p5(Y + 0.5*self.height)
        return nf


class SimpleOpening(CxDoor):
    SQRT2 = math.sqrt(2.0)
    TWO_THIRDS = 2.0/3.0
    def calculate(self, path):
        pn = path.pn
        pm = path.pm
        drho = pn.density - pm.density
        abs_pdrop = abs(path.pressure_drop)
        if abs(drho) <= 0.0001 * abs_pdrop:
            # Treat as orifice flow
            pass
        gdrho = GRAVITY * drho
        Y = path.pressure_drop/gdrho
        path.Y = Y
        # F0 = lower flow, FH = upper flow.
        C = self.SQRT2 * self.width * self.cd
        B = C * math.sqrt(abs_pdrop)
        DF0 = B / abs(gdrho)
        # F0 = 0.666667d0*C*SQRT(ABS(GDRHO*Y))*ABS(Y)
        F0 = self.TWO_THIRDS * B * abs(Y)
        DFH = C * math.sqrt(abs((self.height - Y) / gdrho))
        # FH = 0.666667d0*DFH*ABS(GDRHO*(Height-Y))
        FH = self.TWO_THIRDS * DFH * abs(gdrho * (self.height - Y))
        nf = 1
        if Y <= 0.0:
            # One-way flow (negative).
            if gdrho >= 0.0:
                path.flow[0] = -pm.sqrt_density * abs(FH - F0)
                path.DF[0] = pm.sqrt_density * abs(DFH - DF0)
            else:
                path.flow[0] = pn.sqrt_density * abs(FH - F0)
                path.DF[0] = pn.sqrt_density * abs(DFH - DF0)
        elif Y >= self.height:
            # One-way flow (positive).
            if gdrho >= 0.0:
                path.flow[0] = pn.sqrt_density * abs(FH - F0)
                path.DF[0] = pn.sqrt_density * abs(DFH - DF0)
            else:
                path.flow[0] = -pm.sqrt_density * abs(FH - F0)
                path.DF[0] = pm.sqrt_density * abs(DFH - DF0)
        else:
            # Two-way flow.
            nf = 2
            if gdrho >= 0.0:
                path.flow[0] = -pm.sqrt_density * FH
                path.DF[0] = pm.sqrt_density * DFH
                path.flow[1] = pn.sqrt_density * F0
                path.DF[1] = pn.sqrt_density * DF0
            else:
                path.flow[0] = pn.sqrt_density * FH
                path.DF[0] = pn.sqrt_density * DFH
                path.flow[1] = -pm.sqrt_density * F0
                path.DF[1] = pm.sqrt_density * DF0
        return nf

class NewSimpleOpening(PowerLaw):
    SQRT2 = math.sqrt(2.0)
    TWO_THIRDS = 2.0/3.0
    def __init__(self, nr=0, ht=2.0, wd=0.8, cd=0.78, exponent=0.5, reference_conditions=None):
        self.height = ht
        self.width = wd
        self.cd = cd
        if reference_conditions is None:
            reference_conditions = Node()
        self.reference_density = reference_conditions.density
        self.reference_dynamic_viscosity = reference_conditions.viscosity / reference_conditions.density
        self.linear, self.coefficient = orifice_coefficients(ht=ht, wd=wd, cd=cd, exponent=exponent, 
                                                             reference_density=reference_conditions.density,
                                                             reference_sqrt_density=reference_conditions.sqrt_density,
                                                             reference_viscosity=reference_conditions.viscosity)
        self.exponent = exponent
    def calculate(self, path):
        pn = path.pn
        pm = path.pm
        drho = pn.density - pm.density
        abs_pdrop = abs(path.pressure_drop)
        if abs(drho) <= 0.0001 * abs_pdrop:
            #return orifice_flow(self.linear, self.coefficient, self.exponent, path)
            return super().calculate(path)
        gdrho = GRAVITY * drho
        Y = path.pressure_drop/gdrho
        path.Y = Y
        # F0 = lower flow, FH = upper flow.
        C = self.SQRT2 * self.width * self.cd
        B = C * math.sqrt(abs_pdrop)
        DF0 = B / abs(gdrho)
        # F0 = 0.666667d0*C*SQRT(ABS(GDRHO*Y))*ABS(Y)
        F0 = self.TWO_THIRDS * B * abs(Y)
        DFH = C * math.sqrt(abs((self.height - Y) / gdrho))
        # FH = 0.666667d0*DFH*ABS(GDRHO*(Height-Y))
        FH = self.TWO_THIRDS * DFH * abs(gdrho * (self.height - Y))
        nf = 1
        if Y <= 0.0:
            # One-way flow (negative).
            if gdrho >= 0.0:
                path.flow[0] = -pm.sqrt_density * abs(FH - F0)
                path.DF[0] = pm.sqrt_density * abs(DFH - DF0)
            else:
                path.flow[0] = pn.sqrt_density * abs(FH - F0)
                path.DF[0] = pn.sqrt_density * abs(DFH - DF0)
        elif Y >= self.height:
            # One-way flow (positive).
            if gdrho >= 0.0:
                path.flow[0] = pn.sqrt_density * abs(FH - F0)
                path.DF[0] = pn.sqrt_density * abs(DFH - DF0)
            else:
                path.flow[0] = -pm.sqrt_density * abs(FH - F0)
                path.DF[0] = pm.sqrt_density * abs(DFH - DF0)
        else:
            # Two-way flow.
            nf = 2
            if gdrho >= 0.0:
                path.flow[0] = -pm.sqrt_density * FH
                path.DF[0] = pm.sqrt_density * DFH
                path.flow[1] = pn.sqrt_density * F0
                path.DF[1] = pn.sqrt_density * DF0
            else:
                path.flow[0] = pn.sqrt_density * FH
                path.DF[0] = pn.sqrt_density * DFH
                path.flow[1] = -pm.sqrt_density * F0
                path.DF[1] = pm.sqrt_density * DF0
        return nf