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

from objects import *

def setup_case(n0, n1, dP, drho, ht):
    '''
    Set up the nodes for a given density difference and return the
    pressure difference at the base of the opening.
    '''
    rho = n0.density - drho
    n1.update(n1.pressure, rho=rho)
    return dP + 0.5 * GRAVITY * drho * ht

element0 = CxDoor(ht=3.0, wd=1.0, cd=0.1)
element1 = AltDoor(ht=3.0, wd=1.0, cd=0.1)
element2 = SimpleOpening(ht=3.0, wd=1.0, cd=0.1)
element3 = NewSimpleOpening(ht=3.0, wd=1.0, cd=0.1)

node0 = Node(0)
node1 = Node(1)

path0 = Path(0, node0, node1, element0)
path1 = Path(1, node0, node1, element1)
path2 = Path(2, node0, node1, element2)
path3 = Path(3, node0, node1, element3)

dP = 3.0
drho = 0.1
dP2 = setup_case(node0, node1, dP, drho, element0.height)
path0.pressure_drop = dP
path1.pressure_drop = dP
path2.pressure_drop = dP2
path3.pressure_drop = dP2

nf0 = element0.calculate(path0)
nf1 = element1.calculate(path1)
nf2 = element2.calculate(path2)
nf3 = element3.calculate(path3)

print(nf0, path0.flow[0], path0.Y)
print(nf1, path1.flow[0], path1.Y)
print(nf2, path2.flow[0], path2.Y)
print(nf3, path3.flow[0], path3.Y)
print()

dP = 3.0
drho = -0.1
dP2 = setup_case(node0, node1, dP, drho, element0.height)
path0.pressure_drop = dP
path1.pressure_drop = dP
path2.pressure_drop = dP2
path3.pressure_drop = dP2

nf0 = element0.calculate(path0)
nf1 = element1.calculate(path1)
nf2 = element2.calculate(path2)
nf3 = element3.calculate(path3)

print(nf0, path0.flow[0], path0.Y)
print(nf1, path1.flow[0], path1.Y)
print(nf2, path2.flow[0], path2.Y)
print(nf3, path3.flow[0], path3.Y)
print()

dP = -3.0
drho = -0.1
dP2 = setup_case(node0, node1, dP, drho, element0.height)
path0.pressure_drop = dP
path1.pressure_drop = dP
path2.pressure_drop = dP2
path3.pressure_drop = dP2

nf0 = element0.calculate(path0)
nf1 = element1.calculate(path1)
nf2 = element2.calculate(path2)
nf3 = element3.calculate(path3)

print(nf0, path0.flow[0], path0.Y)
print(nf1, path1.flow[0], path1.Y)
print(nf2, path2.flow[0], path2.Y)
print(nf3, path3.flow[0], path3.Y)
print()

dP = -3.0
drho = 0.1
dP2 = setup_case(node0, node1, dP, drho, element0.height)
path0.pressure_drop = dP
path1.pressure_drop = dP
path2.pressure_drop = dP2
path3.pressure_drop = dP2

nf0 = element0.calculate(path0)
nf1 = element1.calculate(path1)
nf2 = element2.calculate(path2)
nf3 = element3.calculate(path3)

print(nf0, path0.flow[0], path0.Y)
print(nf1, path1.flow[0], path1.Y)
print(nf2, path2.flow[0], path2.Y)
print(nf3, path3.flow[0], path3.Y)
print()

dP = -1.0
drho = 0.1
dP2 = setup_case(node0, node1, dP, drho, element0.height)
path0.pressure_drop = dP
path1.pressure_drop = dP
path2.pressure_drop = dP2
path3.pressure_drop = dP2

nf0 = element0.calculate(path0)
nf1 = element1.calculate(path1)
nf2 = element2.calculate(path2)
nf3 = element3.calculate(path3)

print(nf0, path0.flow[0], path0.flow[1], path0.Y)
print(nf1, path1.flow[0], path1.flow[1], path1.Y)
print(nf2, path2.flow[0], path2.flow[1], path2.Y)
print(nf3, path3.flow[0], path3.flow[1], path3.Y)
print()

dP = -1.0
drho = -0.1
dP2 = setup_case(node0, node1, dP, drho, element0.height)
path0.pressure_drop = dP
path1.pressure_drop = dP
path2.pressure_drop = dP2
path3.pressure_drop = dP2

nf0 = element0.calculate(path0)
nf1 = element1.calculate(path1)
nf2 = element2.calculate(path2)
nf3 = element3.calculate(path3)

print(nf0, path0.flow[0], path0.flow[1], path0.Y)
print(nf1, path1.flow[0], path1.flow[1], path1.Y)
print(nf2, path2.flow[0], path2.flow[1], path2.Y)
print(nf3, path3.flow[0], path3.flow[1], path3.Y)
print()
