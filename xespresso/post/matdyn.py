
from xespresso.post.base import PostCalculation


class EspressoMatdyn(PostCalculation):

    package = 'matdyn'
    package_parameters = {'INPUT': ['flfrc', 'asr', 'dos', 'nk1', 'nk2', 'nk3', 'deltaE', 'ndos', 'fldos',
                                    'flfrq', 'flvec', 'fleig', 'fldvn', 'at', 'l1', 'l2', 'l3', 'ntyp',
                                    'amass', 'readtau', 'fltau', 'la2F', 'q_in_band_form', 'q_in_cryst_coord',
                                    'eigen_similarity', 'fd', 'na_ifc', 'nosym', 'loto_2d',
                                    'loto_disable']
                          }

    def __init__(self, parent_directory, prefix, queue=False, parallel='', **kwargs) -> None:
        super().__init__(parent_directory, prefix, queue=False, parallel='', **kwargs)
