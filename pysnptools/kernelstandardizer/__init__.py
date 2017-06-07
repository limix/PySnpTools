
#!!!move these to their own files

import numpy as np
import scipy as sp
import logging

class KernelStandardizer(object):
    '''
    A KernelStandardizer is a class such as :class:`.DiagKtoN` and :class:`.Identity` to be used by the :meth:`.KernelData.standardize` to standardize Kernel data.
    It always works in-place *and* returns the :class:`.KernelData` on which it works.


    Details of Methods & Properties:
    '''
    def standardize(self, kerneldata, return_trained=False, force_python_only=False):
        '''
        Applies standardization, in place, to :class:`.KernelData`. For convenience also returns the :class:`KernelData`.

        :param snps: kernel values to standardize
        :type snps: :class:`.KernelData`

        :param return_trained: If true, returns a second value containing a constant :class:`.KernelStandardizer` trained on this data.
        :type return_trained: boolean

        :param force_python_only: optional -- If False (default), may use outside library code. If True, requests that the read
            be done without outside library code.
        :type force_python_only: bool

        :rtype: :class:`.KernelData`, (optional) constant :class:`.KernelStandardizer`

        '''
        raise NotImplementedError("subclass {0} needs to implement method '.standardize'".format(self.__class__.__name__))

class Identity(KernelStandardizer):
    '''
    A :class:`.KernelStandardizer` that does nothing to kernel data.

    See :class:`.KernelStandardizer` for more information about standardization.

    >>> from pysnptools.kernelstandardizer import Identity as KS_Identity
    >>> from pysnptools.kernelreader import KernelNpz
    '''

    def __init__(self):
        super(Identity, self).__init__()

    def standardize(self, kerneldata, return_trained=False, force_python_only=False):
        if return_trained:
            return kerneldata, self
        else:
            return kerneldata

    def __repr__(self):
        return "{0}()".format(self.__class__.__name__)

from pysnptools.standardizer import DiagKtoN #as SN_DiagKtoN
from pysnptools.standardizer import DiagKtoNTrained #as SN_DiagKtoNTrained



if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
