
from xespresso.post.base import PostCalculation


class EspressoBands(PostCalculation):

    package = 'bands'
    package_parameters = {'BANDS': ['prefix', 'outdir', 'filband', 'spin_component', 'lsigma',
                                    'lp', 'filp', 'lsym', 'no_overlap', 'plot_2d', 'firstk', 'lastk']
                          }

    def __init__(self, parent_directory, prefix, queue=False, parallel='', **kwargs) -> None:
        super().__init__(parent_directory, prefix, queue=queue, parallel=parallel, **kwargs)
