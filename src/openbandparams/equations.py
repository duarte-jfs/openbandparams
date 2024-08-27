#
#   Copyright (c) 2013-2014, Scott J Maddox
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

from numpy import sqrt

def varshni(Eg_0, alpha, beta, T):
    return Eg_0 - alpha * T ** 2 / (T + beta)

def refractive_index(E0, Eg, SO, A, B):
    x0 = E0/Eg
    x0s = E0/(Eg + SO)
    fx0 = x0**-2 * (2-sqrt(1+x0)-sqrt(1-x0))
    fx0s = x0s**-2 * (2-sqrt(1+x0s)-sqrt(1-x0s))

    return sqrt(A * (fx0 + 0.5*(Eg/(Eg+SO))**1.5*fx0s)+B)
