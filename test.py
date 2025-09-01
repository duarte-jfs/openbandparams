#
#   Copyright (c) 2013-2014, Scott J Maddox
#   Copyright (c) 2025, Duarte Silva
#
#   This file is part of openbandparams.
#
#   openbandparams is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Affero General Public License as published
#   by the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   openbandparams is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Affero General Public License for more details.
#
#   You should have received a copy of the GNU Affero General Public License
#   along with openbandparams.  If not, see <http://www.gnu.org/licenses/>.
#
#############################################################################
# '''
# Find and run all unit tests in the project.
# '''
# Make sure we import the local package
# import os
# import sys
# sys.path.insert(0,
#     os.path.abspath(os.path.join(os.path.dirname(__file__), 'src')))
# import nose
# nose.main()

from openbandparams.iii_v_zinc_blende_quaternaries import AlPAsSb, AlGaInAs, AlInPSb
from openbandparams.iii_v_zinc_blende_ternaries import AlPAs, AlPSb, AlAsSb
from openbandparams.iii_v_zinc_blende_binaries import InP, AlSb, AlAs

T=300
x=0.07929292929292929
y=0.3887442324479361

alloy = AlPAsSb(x=x, a = InP.a(T=T))
Eg = alloy.Eg(T=T)
print(alloy.binaries)
print((1-x-y)/2, (2-x-2*y)/2, (2-2*x-y)/2)

print(AlAsSb(a=InP.a(T=T)).Eg(T=T))

print(AlSb.Eg(T=T))
print(AlAs.Eg(T=T))

import numpy as np

x=1
print(AlInPSb(x=x, y=(0.610-0.343*x)/(0.610+0.075*x)).Eg(T=T))