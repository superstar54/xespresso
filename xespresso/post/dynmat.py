
from xespresso.post.base import PostCalculation


class EspressoDynmat(PostCalculation):

    package = 'dynmat'

    package_parameters = {'INPUT': ['fildyn', 'q', 'amass', 'asr', 'axis', 'lperm', 'lplasma', 'filout',
                                    'fileig', 'filmol', 'filxsf', 'loto_2d', 'el_ph_nsig', 'el_ph_sigma']
                          }

    def __init__(self, parent_directory, prefix, queue=False, parallel='', **kwargs) -> None:
        super().__init__(parent_directory, prefix, queue=False, parallel='', **kwargs)
