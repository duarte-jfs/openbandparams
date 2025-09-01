#
#   Copyright (c) 2013-2015, Scott J Maddox
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
__all__ = ['IIIVZincBlendeQuaternary']

from .iii_v_zinc_blende_mixed_alloy import IIIVZincBlendeMixedAlloy
from .algorithms import bisect

from typing import overload

class IIIVZincBlendeQuaternary(IIIVZincBlendeMixedAlloy):
    '''
    The base class for all III-V zinc blende quaternary alloys.
    '''
    def __init__(self, name, elements, ternaries, parameters=None,
                 x=None, y=None, z=None):
        if (len(ternaries) == 3 and
            ternaries[0].elements[0] == ternaries[1].elements[0] and # A
            ternaries[0].elements[0] == ternaries[2].elements[0] and # A
            ternaries[0].elements[1] == ternaries[1].elements[1] and # B
            ternaries[0].elements[2] == ternaries[2].elements[1] and # C
            ternaries[1].elements[2] == ternaries[2].elements[2]):   # D
            # Type 1: AB_{x}C_{y}D_{1-x-y}
            # binaries = (AB, AC, AD)
            # ternaries = (ABC, ABD ,ACD)
            self._type = 1
            self._element_w = ternaries[0].elements[0]
            self._element_x = ternaries[0].elements[1]
            self._element_y = ternaries[0].elements[2]
            self._element_z = ternaries[1].elements[2]
            self.binaries = (ternaries[0].binaries[0],
                             ternaries[0].binaries[1],
                             ternaries[1].binaries[1],)
            calc_elements = (ternaries[0].elements[0],
                             ternaries[0].elements[1],
                             ternaries[0].elements[2],
                             ternaries[1].elements[2])
        elif (len(ternaries) == 3 and
              ternaries[0].elements[0] == ternaries[1].elements[0] and # A
              ternaries[0].elements[1] == ternaries[2].elements[0] and # B
              ternaries[1].elements[1] == ternaries[2].elements[1] and # C
              ternaries[0].elements[2] == ternaries[1].elements[2] and # D
              ternaries[0].elements[2] == ternaries[2].elements[2]):   # D
            # Type 2: A_{x}B_{y}C_{1-x-y}D
            # binaries = (AD, BD, CD)
            # ternaries = (ABD, ACD, BCD)
            self._type = 2
            self._element_x = ternaries[0].elements[0]
            self._element_y = ternaries[2].elements[0]
            self._element_z = ternaries[2].elements[1]
            self._element_w = ternaries[2].elements[2]
            self.binaries = (ternaries[0].binaries[0],
                             ternaries[0].binaries[1],
                             ternaries[1].binaries[1],)
            calc_elements = (ternaries[0].elements[0],
                             ternaries[2].elements[0],
                             ternaries[2].elements[1],
                             ternaries[2].elements[2])
        elif (len(ternaries) == 4 and
              ternaries[0].elements[0] == ternaries[1].elements[0] and # A
              ternaries[0].elements[0] == ternaries[2].elements[0] and # A
              ternaries[0].elements[1] == ternaries[1].elements[1] and # B
              ternaries[0].elements[1] == ternaries[3].elements[0] and # B
              ternaries[0].elements[2] == ternaries[2].elements[1] and # C
              ternaries[0].elements[2] == ternaries[3].elements[1] and # C
              ternaries[1].elements[2] == ternaries[2].elements[2] and # D
              ternaries[1].elements[2] == ternaries[3].elements[2]):   # D
            # Type 3: A_{x}B_{1-x}C_{y}D_{1-y}
            # binaries = (AC, AD, BC, BD)
            # ternaries = (ABC, ABD, ACD, BCD)
            self._type = 3
            self._element_x = ternaries[0].elements[0]
            self._element_1mx = ternaries[0].elements[1]
            self._element_y = ternaries[2].elements[1]
            self._element_1my = ternaries[2].elements[2]
            self.binaries = (ternaries[2].binaries[0],
                             ternaries[2].binaries[1],
                             ternaries[3].binaries[0],
                             ternaries[3].binaries[1],)
            calc_elements = (ternaries[0].elements[0],
                             ternaries[0].elements[1],
                             ternaries[2].elements[1],
                             ternaries[2].elements[2])
        else:
            raise ValueError()
        assert elements == calc_elements
        super(IIIVZincBlendeQuaternary, self).__init__(name, elements,
                                                    parameters=parameters)
        self.ternaries = ternaries
        if x is not None or y is not None or z is not None:
            self._xyz = self._parse_xyz(x, y, z)
        else:
            self._xyz = None

    def __eq__(self, other):
        return all((type(self) == type(other) and
                self.name == other.name and
                self.elements == other.elements,
                self.ternaries == other.ternaries,
                self._parameters == other._parameters,
                self._xyz == other._xyz))
    
    def __deepcopy__(self, memo):
        from copy import deepcopy
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result

        for attr in [
            'name', 
            'elements', 
            '_parameters', 
            '_aliases', 
            '_xyz', 
            'binaries',
            'ternaries', 
            '_type',
            '_x',
            '_y',
            '_z',
            ]:
            setattr(result, attr, deepcopy(getattr(self, attr), memo))

        type = getattr(self, '_type')
        if type == 1:
            for attr in ['_element_w', '_element_x', '_element_y', '_element_z']:
                setattr(result, attr, deepcopy(getattr(self, attr), memo))
        elif type == 2:
            for attr in ['_element_x', '_element_y', '_element_z', '_element_w']:
                setattr(result, attr, deepcopy(getattr(self, attr), memo))
        elif type == 3:
            for attr in ['_element_x', '_element_1mx', '_element_y', '_element_1my']:
                setattr(result, attr, deepcopy(getattr(self, attr), memo))

        return result

    def _parse_xyz(self, x, y, z):
        if self._type == 1 or self._type == 2:
            # Type 1: AB_{x}C_{y}D_{1-x-y}
            # binaries = (AB, AC, AD)
            # ternaries = (ABC, ABD ,ACD)
            # Type 2: A_{x}B_{y}C_{1-x-y}D
            # binaries = (AD, BD, CD)
            # ternaries = (ABD, ACD, BCD)
            if x is not None and y is not None and z is None:
                x = round(float(x), 6)
                y = round(float(y), 6)
                z = round(1. - x - y, 6)
            elif x is not None and y is None and z is not None:
                x = round(float(x), 6)
                z = round(float(z), 6)
                y = round(1. - x - z, 6)
            elif x is None and y is not None and z is not None:
                y = round(float(y), 6)
                z = round(float(z), 6)
                x = round(1. - y - z, 6)
            else:
                raise ValueError()
        elif self._type == 3:
            # Type 3: A_{x}B_{1-x}C_{y}D_{1-y}
            # binaries = (AC, AD, BC, BD)
            # ternaries = (ABC, ABD, ACD, BCD)
            if x is not None and y is not None and z is None:
                x = round(float(x), 6)
                y = round(float(y), 6)
                z = None
            else:
                raise ValueError()
        else:
            raise RuntimeError()
        if (not (0. <= x <= 1.) or not (0. <= y <= 1.) or
            z is not None and not (0. <= z <= 1.)        ):
            raise ValueError('The alloy fractions must be between 0 and 1')
        
        self._x = x
        self._y = y
        self._z = z
        return x, y, z
    
    def _instance(self, x=None, y=None, z=None):
        instance = IIIVZincBlendeQuaternary(self.name, self.elements,
                                            self.ternaries, x=x, y=y, z=z)
        for parameter in self._parameters.values():
            instance.set_parameter(parameter)
        return instance

    def _has_x(self, kwargs):
        '''Returns True if x is explicitly defined in kwargs'''
        return (('x' in kwargs) or (self._element_x in kwargs) or
                (self._type == 3 and self._element_1mx in kwargs))

    def _get_x(self, kwargs):
        '''
        Returns x if it is explicitly defined in kwargs.
        Otherwise, raises TypeError.
        '''
        if 'x' in kwargs:
            return round(float(kwargs['x']), 6)
        elif self._element_x in kwargs:
            return round(float(kwargs[self._element_x]), 6)
        elif self._type == 3 and self._element_1mx in kwargs:
            return round(1. - float(kwargs[self._element_1mx]), 6)
        else:
            raise TypeError()

    def _has_y(self, kwargs):
        '''Returns True if y is explicitly defined in kwargs'''
        return (('y' in kwargs) or (self._element_y in kwargs) or
                (self._type == 3 and self._element_1my in kwargs))

    def _get_y(self, kwargs):
        '''
        Returns y if it is explicitly defined in kwargs.
        Otherwise, raises TypeError.
        '''
        if 'y' in kwargs:
            return round(float(kwargs['y']), 6)
        elif self._element_y in kwargs:
            return round(float(kwargs[self._element_y]), 6)
        elif self._type == 3 and self._element_1my in kwargs:
            return round(1. - float(kwargs[self._element_1my]), 6)
        else:
            raise TypeError()

    def _has_z(self, kwargs):
        '''
        Returns True if type is 1 or 2 and z is explicitly defined in kwargs.
        '''
        return ((self._type == 1 or self._type ==2) and
                (('z' in kwargs) or (self._element_z in kwargs)))

    def _get_z(self, kwargs):
        '''
        Returns z if type is 1 or 2 and z is explicitly defined in kwargs.
        Otherwise, raises TypeError.
        '''
        if self._type == 1 or self._type == 2:
            if 'z' in kwargs:
                return round(float(kwargs['z']), 6)
            elif self._element_z in kwargs:
                return round(float(kwargs[self._element_z]), 6)
        raise TypeError()

    def __call__(self, **kwargs):
        '''
        
        Args:
            x: The x composition.
            y: The y composition.
            z: The z composition.
            a: The lattice constant.

        returns:
            An instance of IIIVZincBlendeQuaternary with updated parameters.
        '''
        if self._has_x(kwargs) and self._has_y(kwargs):
            x = self._get_x(kwargs)
            y = self._get_y(kwargs)
            z = None
        elif self._has_x(kwargs) and self._has_z(kwargs):
            x = self._get_x(kwargs)
            y = None
            z = self._get_z(kwargs)
        elif self._has_y(kwargs) and self._has_z(kwargs):
            x = None
            y = self._get_y(kwargs)
            z = self._get_z(kwargs)
        elif 'a' in kwargs and (self._has_x(kwargs) or
                                self._has_y(kwargs) or
                                self._has_z(kwargs)   ):
            # lattice match to the given lattice constant
            a = kwargs['a']
            T = kwargs.get('T', 300.)
            # make sure the lattice constant is available
            if self._has_x(kwargs):
                x = self._get_x(kwargs)
                if self._type in [1, 2]:
                    ymax = round(1. - x, 6)
                elif self._type == 3:
                    ymax = 1.
                else:
                    raise RuntimeError()
                z = None
                a0 = self(x=x, y=0.).a(T=T)
                a1 = self(x=x, y=ymax).a(T=T)
                amin = min(a0, a1)
                amax = max(a0, a1)
                if not (amin <= a <= amax):
                    raise ValueError('a of {:g} out of range [{:g}, {:g}]'
                                     ''.format(a, amin, amax))
                # find the correct composition, x
                y = bisect(func=lambda y: self(x=x, y=y).a(T=T) - a,
                           a=0, b=ymax)
            elif self._has_y(kwargs):
                y = self._get_y(kwargs)
                if self._type in [1, 2]:
                    xmax = round(1. - y, 6)
                elif self._type == 3:
                    xmax = 1.
                else:
                    raise RuntimeError()
                z = None
                a0 = self(x=0., y=y).a(T=T)
                a1 = self(x=xmax, y=y).a(T=T)
                amin = min(a0, a1)
                amax = max(a0, a1)
                if not (amin <= a <= amax):
                    raise ValueError('a of {:g} out of range [{:g}, {:g}]'
                                     ''.format(a, amin, amax))
                # find the correct composition, x
                x = bisect(func=lambda x: self(x=x, y=y).a(T=T) - a,
                           a=0, b=xmax)
            elif self._has_z(kwargs):
                y = None
                z = self._get_z(kwargs)
                xmax = round(1. - z, 6)
                a0 = self(x=0., z=z).a(T=T)
                a1 = self(x=xmax, z=z).a(T=T)
                amin = min(a0, a1)
                amax = max(a0, a1)
                if not (amin <= a <= amax):
                    raise ValueError('a of {:g} out of range [{:g}, {:g}]'
                                     ''.format(a, amin, amax))
                # find the correct composition, x
                x = bisect(func=lambda x: self(x=x, z=z).a(T=T) - a,
                           a=0, b=xmax)
        else:
            raise TypeError(
                "Missing required key word argument.\n" + self._get_usage())
        return self._instance(x=x, y=y, z=z)

    def _get_usage(self):
        if self._type == 1 or self._type == 2:
            return ("The supported kwarg combinations are as follows:"
                    "\n    - ('x' or '{x}') and ('y' or '{y}')"
                    "\n    - ('x' or '{x}') and ('z' or '{z}')"
                    "\n    - ('y' or '{y}') and ('z' or '{z}')"
                    "\n    - 'a' [and 'T'] and ('x' or '{x}')"
                    "\n    - 'a' [and 'T'] and ('y' or '{y}')"
                    "\n    - 'a' [and 'T'] and ('z' or '{z}')"
                    "".format(x=self._element_x, y=self._element_y,
                              z=self._element_z))
        elif self._type == 3:
            return ("The supported kwargs combinations are as follows:"
                    "\n    - ('x', '{x}' or '{_1mx}') and "
                            "('y', '{y}' or '{_1my}')"
                    "\n    - 'a' [and 'T'] and ('x', '{x}' or '{_1mx}')"
                    "\n    - 'a' [and 'T'] and ('y', '{y}' or '{_1my}')"
                    "".format(x=self._element_x, _1mx=self._element_1mx,
                              y=self._element_y, _1my=self._element_1my))
        else:
            raise RuntimeError()

    def __repr__(self):
        if self._xyz is None:
            return '{}'.format(self.name)
        elif self._type in [1, 2, 3]:
            x, y, _ = self._xyz
            return '{}({}={}, {}={})'.format(self.name,
                                             self._element_x, x,
                                             self._element_y, y)
        else:
            raise RuntimeError()
     
    def latex(self):
        e = {'A':self.elements[0],
             'B':self.elements[1],
             'C':self.elements[2],
             'D':self.elements[3]}
        if self._type == 1:
            if self._xyz is None:
                return ("{A}{B}_{{x}}{C}_{{y}}{D}_{{1-x-y}}"
                       "".format(**e))
            else:
                x, y, z = self._xyz
                return ("{A}{B}_{{{:g}}}{C}_{{{:g}}}{D}_{{{:g}}}"
                        "".format(x, y, z, **e))
        elif self._type == 2:
            if self._xyz is None:
                return ("{A}_{{x}}{B}_{{y}}{C}_{{1-x-y}}{D}"
                        "".format(**e))
            else:
                x, y, z = self._xyz
                return ("{A}_{{{:g}}}{B}_{{{:g}}}{C}_{{{:g}}}{D}"
                        "".format(x, y, z, **e))
        elif self._type == 3:
            if self._xyz is None:
                return ("{A}_{{x}}{B}_{{1-x}}{C}_{{y}}{D}_{{1-y}}"
                        "".format(**e))
            else:
                x, y, z = self._xyz
                _1mx = round(1. - x, 6)
                _1my = round(1. - y, 6)
                return ("{A}_{{{:g}}}{B}_{{{:g}}}{C}_{{{:g}}}{D}_{{{:g}}}"
                        "".format(x, _1mx, y, _1my, **e))
        else:
            raise RuntimeError()
 
    def element_fraction(self, element):
        if self._xyz is None:
            raise TypeError('Alloy composition has not been specified.')
        if self._type == 1 or self._type == 2:
            # AB_{x}C_{y}D_{1-x-y}
            if element == self._element_w:
                return 1
            elif element == self._element_x:
                return self._x
            elif element == self._element_y:
                return self._y
            elif element == self._element_z:
                return self._z
            else:
                return 0
        elif self._type == 3:
            if element == self._element_x:
                return self._x
            elif element == self._element_1mx:
                return (1 - self._x)
            elif element == self._element_y:
                return self._y
            elif element == self._element_1my:
                return (1 - self._y)
            else:
                return 0
        else:
            raise RuntimeError()
     
    def _get_bowing(self, name, kwargs):
        p = self.get_parameter(name+'_bowing', default=None)
        if p is None:
            return None
        x, y, z = self._xyz
        return p(x=x, y=y, z=z, **kwargs)
     
    def _interpolate(self, name, kwargs):
        if self._xyz is None:
            raise TypeError('Alloy composition has not been specified.')
        if self._type == 1 or self._type == 2:
            return self._interpolate1or2(name, kwargs)
        elif self._type == 3:
            return self._interpolate3(name, kwargs)
        else:
            raise RuntimeError()
        
    def _interpolate1or2(self, name, kwargs):
        """
        Implements equation 11 from [1]

        [1] - Olesberg, J.T., 2024. Interpolation of compound semiconductor alloy parameters from those of their constituents. Journal of Applied Physics 136, 215105. https://doi.org/10.1063/5.0217016

        """
        x, y, z = self._xyz

        if self._type == 1:
            # Type 1: AB_{x}C_{y}D_{1-x-y}
            # binaries = (AB, AC, AD)
            # ternaries = (ABC, ABD ,ACD)

            a = x
            b = y
            c = 1-x-y
            p = 1

            b_AP = self.binaries[0].get_parameter(name)
            if b_AP is None:
                raise AttributeError('"{}" is missing a required parameter: "{}".'
                                     ''.format(self.binaries[0], name))
            b_AP = b_AP(**kwargs)

            b_BP = self.binaries[1].get_parameter(name)
            if b_BP is None:
                raise AttributeError('"{}" is missing a required parameter: "{}".'
                                     ''.format(self.binaries[1], name))
            b_BP = b_BP(**kwargs)

            b_CP = self.binaries[2].get_parameter(name)
            if b_CP is None:
                raise AttributeError('"{}" is missing a required parameter: "{}".'
                                     ''.format(self.binaries[2], name))
            b_CP = b_CP(**kwargs)

            bow_ABP = self.ternaries[1](x=0)._get_bowing(name, kwargs)
            bow_ACP = self.ternaries[2](x=0)._get_bowing(name, kwargs)
            bow_BCP = self.ternaries[0](x=0)._get_bowing(name, kwargs)
            bow_ABCP = self._get_bowing(name, kwargs)

            #In case the bow if None, put it to 0
            if bow_ABP is None:
                bow_ABP = 0
            if bow_ACP is None:
                bow_ACP = 0
            if bow_BCP is None:
                bow_BCP = 0
            if bow_ABCP is None:
                bow_ABCP = 0

            return (
                a*p*b_AP + 
                b*p*b_BP +
                c*p*b_CP -
                a*b*p*bow_ABP -
                a*c*p*bow_ACP -
                b*c*p*bow_BCP -
                a*b*c*p*bow_ABCP
            )
        elif self._type == 2:
            # Type 2: A_{x}B_{y}C_{1-x-y}D
            # binaries = (AD, BD, CD)
            # ternaries = (ABD, ACD, BCD)

            a = x
            b = y
            c = 1-x-y
            p = 1

            b_AP = self.binaries[0].get_parameter(name)
            if b_AP is None:
                raise AttributeError('"{}" is missing a required parameter: "{}".'
                                     ''.format(self.binaries[0], name))
            b_AP = b_AP(**kwargs)

            b_BP = self.binaries[1].get_parameter(name)
            if b_BP is None:
                raise AttributeError('"{}" is missing a required parameter: "{}".'
                                     ''.format(self.binaries[1], name))
            b_BP = b_BP(**kwargs)

            b_CP = self.binaries[2].get_parameter(name)
            if b_CP is None:
                raise AttributeError('"{}" is missing a required parameter: "{}".'
                                     ''.format(self.binaries[2], name))
            b_CP = b_CP(**kwargs)

            bow_ABP = self.ternaries[0](x=0)._get_bowing(name, kwargs)
            bow_ACP = self.ternaries[1](x=0)._get_bowing(name, kwargs)
            bow_BCP = self.ternaries[2](x=0)._get_bowing(name, kwargs)
            bow_ABCP = self._get_bowing(name, kwargs)

            #In case the bow if None, put it to 0
            if bow_ABP is None:
                bow_ABP = 0
            if bow_ACP is None:
                bow_ACP = 0
            if bow_BCP is None:
                bow_BCP = 0
            if bow_ABCP is None:
                bow_ABCP = 0

            return (
                a*p*b_AP + 
                b*p*b_BP +
                c*p*b_CP -
                a*b*p*bow_ABP -
                a*c*p*bow_ACP -
                b*c*p*bow_BCP -
                a*b*c*p*bow_ABCP
            )

    def _interpolate3(self, name, kwargs):
        """
        Implements equation 10 from [1]

        [1] - Olesberg, J.T., 2024. Interpolation of compound semiconductor alloy parameters from those of their constituents. Journal of Applied Physics 136, 215105. https://doi.org/10.1063/5.0217016

        """
        x, y, _ = self._xyz

        # Type 3: A_{x}B_{1-x}C_{y}D_{1-y}
        #binaries = (AC, AD, BC, BD)
        #ternaries = (ABC, ABD, ACD, BCD)
        # B1 = Q(0, 0) = BD = InSb
        # B2 = Q(1, 0) = BC = InAs
        # B3 = Q(1, 1) = AC = AlAs
        # B4 = Q(0, 1) = AD = AlSb
        # t12 = Q(x, 0) = ABD = AlInSb
        # t23 = Q(1, y) = ACD = AlAsSb
        # t43 = Q(x, 1) = ABC = AlInAs
        # t14 = Q(0, y) = BCD = InAsSb

        b_AP = self.binaries[0].get_parameter(name)
        if b_AP is None:
            raise AttributeError('"{}" is missing a required parameter: "{}".'
                                 ''.format(self.binaries[0], name))
        b_AP = b_AP(**kwargs)
        
        b_BP = self.binaries[2].get_parameter(name)
        if b_BP is None:
            raise AttributeError('"{}" is missing a required parameter: "{}".'
                                 ''.format(self.binaries[2], name))
        b_BP = b_BP(**kwargs)

        b_BQ = self.binaries[3].get_parameter(name)
        if b_BQ is None:
            raise AttributeError('"{}" is missing a required parameter: "{}".'
                                 ''.format(self.binaries[3], name))
        b_BQ = b_BQ(**kwargs)

        b_AQ = self.binaries[1].get_parameter(name)
        if b_AQ is None:
            raise AttributeError('"{}" is missing a required parameter: "{}".'
                                 ''.format(self.binaries[1], name))
        b_AQ = b_AQ(**kwargs)

        bow_ABP = self.ternaries[0](x=0)._get_bowing(name, kwargs)
        bow_ABQ = self.ternaries[1](x=0)._get_bowing(name, kwargs)
        bow_APQ = self.ternaries[2](x=0)._get_bowing(name, kwargs)
        bow_BPQ = self.ternaries[3](x=0)._get_bowing(name, kwargs)
        bow_ABPQ = self._get_bowing(name, kwargs)

        #In case any of the bowing parameters is None, put it to 0
        if bow_ABP is None:
            bow_ABP = 0
        if bow_ABQ is None:
            bow_ABQ = 0
        if bow_APQ is None:
            bow_APQ = 0
        if bow_BPQ is None:
            bow_BPQ = 0
        if bow_ABPQ is None:
            bow_ABPQ = 0

        a=x
        b=1-x
        p=y
        q=1-y

        return (
            a*p*b_AP +
            b*p*b_BP +
            a*q*b_AQ +
            b*q*b_BQ - 
            a*b*p*bow_ABP -
            a*b*q*bow_ABQ -
            a*p*q*bow_APQ -
            b*p*q*bow_BPQ -
            a*b*p*q*bow_ABPQ
        )