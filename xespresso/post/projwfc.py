
from xespresso.post.base import PostCalculation


class EspressoProjwfc(PostCalculation):

    package = 'projwfc'
    package_parameters = {'PROJWFC': ['prefix', 'outdir', 'ngauss', 'degauss',
                                      'Emin', 'Emax', 'DeltaE', 'lsym', 'pawproj', 'filpdos', 'filproj',
                                      'lwrite_overlaps', 'lbinary_data', 'kresolveddos', 'tdosinboxes',
                                      'n_proj_boxes', 'irmin(3,n_proj_boxes)', 'irmax(3,n_proj_boxes)',
                                      'plotboxes', ]
                          }

    def __init__(self, parent_directory, prefix, queue=False, parallel='', **kwargs) -> None:
        super().__init__(parent_directory, prefix, queue=queue, parallel=parallel, **kwargs)
