
from xespresso.post.base import PostCalculation

class EspressoPp(PostCalculation):
    
    package = 'pp'
    state_parameters = {'INPUTPP': ['prefix', 'outdir', 'filplot', 'plot_num', 
                       'spin_component', 'spin_component', 'emin', 'emax', 
                       'delta_e', 'degauss_ldos', 'sample_bias', 'kpoint', 
                       'kband', 'lsign', 'spin_component', 'emin', 'emax', 
                       'spin_component', 'spin_component', 'spin_component', 
                       'spin_component'], 
                'PLOT': ['nfile', 'filepp', 'weight', 'iflag', 'output_format', 
                    'fileout', 'interpolation', 'e1', 'x0', 'nx', 'e1', 'e2', 
                    'x0', 'nx', 'ny', 'e1', 'e2', 'e3', 'x0', 'nx', 'ny', 'nz', 
                    'radius', 'nx', 'ny']
                }

    def __init__(self, parent_directory, prefix, queue = False, parallel = '', **kwargs) -> None:
        super().__init__(parent_directory, prefix, queue = False, parallel = '', **kwargs)

