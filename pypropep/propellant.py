from cpropep._cpropep import ffi, lib

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
            print 'blah'
            RuntimeWarning("Element {0} does not exist".format(element))
            return 0

        if el_ind in elem:
            el_ind = elem.index(el_ind)
            return coef[el_ind]

        return 0

    def __str__(self):
        return "Propellant: {} - {} [{}]".format(self.formula(),
                                                 self['name'],
                                                 self['id'])

    def __repr__(self):
        return self.__str__()
