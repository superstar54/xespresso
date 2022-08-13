
from xespresso.post.base import PostCalculation


class EspressoHp(PostCalculation):

    package = 'hp'
    package_parameters = {'INPUTHP': ['prefix', 'outdir', 'iverbosity', 'max_seconds', 'nq1', 'nq2',
                                      'nq3', 'skip_equivalence_q', 'determine_num_pert_only',
                                      'find_atpert', 'docc_thr', 'skip_type', 'equiv_type',
                                      'perturb_only_atom', 'start_q', 'last_q', 'sum_pertq',
                                      'compute_hp', 'conv_thr_chi', 'thresh_init', 'ethr_nscf',
                                      'niter_max', 'alpha_mix(i)', 'nmix', 'num_neigh', 'lmin', 'rmax', ]
                          }

    def __init__(self, parent_directory, prefix, queue=False, parallel='', **kwargs) -> None:
        super().__init__(parent_directory, prefix, queue=queue, parallel=parallel, **kwargs)
