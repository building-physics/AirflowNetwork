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

element0 = CxPlr(0, 3.07052e-05, 0.01, 0.5)
element1 = CxFcn(1, 2.43435e-05, 0.01, 0.5)
element2 = CxQcn(1, 3.52946e-05, 0.01, 0.5)

node0 = Node(0)
node1 = Node(1)

path0 = Path(0, node0, node1, element0)
path1 = Path(1, node0, node1, element1)
path2 = Path(2, node0, node1, element2)

dp = 0.00001
#dp = 5.0

path0.pressure_drop = dp
path1.pressure_drop = dp
path2.pressure_drop = dp

element0.calculate(path0)

#element1.calculate(path1)

#element2.calculate(path2)

print(path0.flow[0])
#print(node0.sqrt_density * path1.flow[0], path1.flow[0])
#print(path2.flow[0] / node0.sqrt_density, path2.flow[0])