from .cpropep._cpropep import ffi, lib

__all__ = ['Propellant']

EL_SYMBOLS = [
  "H", "HE", "LI", "BE", "B", "C", "N", "O",
  "F", "NE", "NA", "MG", "AL", "SI", "P", "S", "CL", "AR", "K", "CA",
  "SC", "TI", "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "GE",
  "AS", "SE", "BR", "KR", "RB", "SR", "Y", "ZR", "NB", "MO", "TC", "RU",
  "RH", "PD", "AG", "CD", "IN", "SN", "SB", "TE", "I ", "XE", "CS", "BA",
  "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER",
  "TM", "YB", "LU", "HF", "TA", "W ", "RE", "OS", "IR", "PT", "AU", "HG", "TL",
  "PB", "BI", "PO", "AT", "RN", "FR", "RA", "AC", "TH", "PA", "U", "NP",
  "U6", "U5", "U1", "U2", "U3", "U4", "FM",
  "E", "D" ]


EL_MOLAR_MASS = [
  1.00794,   4.002602, 6.941,      9.012182, 10.811,    12.0107,
  14.00674,  15.9994,  18.9984032, 20.11797, 22.989770, 24.305,
  26.981538, 28.0855,  30.973761,  32.066,   35.4527,   39.948,
  39.0983,   40.078,   44.95591,   47.88,    50.9415,   51.996,
  54.938,    55.847,   58.9332,    58.6934,  63.546,    65.39,
  69.723,    72.61,    74.9216,    78.96,    79.904,    83.80,
  85.4678,   87.62,    88.9059,    91.224,   92.9064,   95.94,
  98.0,      101.07,   102.9055,   106.42,   107.868,   112.41,
  114.82,    118.71,   121.757,    127.60,   126.9045,  131.29,
  132.9054,  137.33,   138.9055,   140.12,   140.9077,  144.24,
  145.,      150.36,   151.965,    157.25,   158.9253,  162.50,
  164.9303,  167.26,   168.9342,   173.04,   174.967,   178.49,
  180.9479,  183.85,   186.207,    190.2,    192.22,    195.08,
  196.9665,  200.59,   204.383,    207.2,    208.9804,  209.,
  210.,      222.,     223.,       226.0254, 227.,      232.0381,
  231.0359,  238.029,  237.0482,   244.,     12.011,    9.01218,
  10.811,    24.305,   26.98154,   257.0,    0,         2]

class Propellant(dict):
    '''
    Thin shell around attrdict for propellant type.
    Purpose of a new class is so we can do instance check
    elsewhere.
    '''
    def __init__(self, *args, **kwargs):
        super(Propellant, self).__init__(*args, **kwargs)

    def formula(self, tex=False):
        elem = self['elem']
        coef = self['coef']
        code = ""
        for e,c in zip(elem, coef):
            if c > 0:
                if tex is True:
                    code += "{0}_{{{1}}}".format(EL_SYMBOLS[e], c)
                else:
                    code += "{0}{1}".format(EL_SYMBOLS[e], c)
        return code


    def atoms_of(self, element):
        '''
        Looks up and returns the number of atoms of
        element in the propellant
        '''
        elem = self['elem']
        coef = self['coef']
        if element in EL_SYMBOLS:
            el_ind = EL_SYMBOLS.index(element)
        else:
            print('blah')
            RuntimeWarning("Element {0} does not exist".format(element))
            return 0

        if el_ind in elem:
            el_ind = elem.index(el_ind)
            return coef[el_ind]

        return 0

    @property
    def mw(self):
        mw = 0.
        elem = self['elem']
        coef = self['coef']
        for e,c in zip(elem, coef):
            if c > 0:
                mw += c * EL_MOLAR_MASS[e]
        return mw

    def __str__(self):
        return "Propellant: {} - {} [{}]".format(self.formula(),
                                                 self['name'],
                                                 self['id'])

    def __repr__(self):
        return self.__str__()
