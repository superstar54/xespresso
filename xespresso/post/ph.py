
from xespresso.post.base import PostCalculation


class EspressoPh(PostCalculation):

    package = 'ph'
    package_parameters = {'INPUTPH': ['amass', 'outdir', 'prefix', 'niter_ph', 'tr2_ph',
                                      'alpha_mix(niter)', 'nmix_ph', 'verbosity', 'reduce_io',
                                      'max_seconds', 'fildyn', 'fildrho', 'fildvscf', 'epsil',
                                      'lrpa', 'lnoloc', 'trans', 'lraman', 'eth_rps', 'eth_ns',
                                      'dek', 'recover', 'low_directory_check', 'only_init', 'qplot',
                                      'q2d', 'q_in_band_form', 'electron_phonon', 'el_ph_nsigma',
                                      'el_ph_sigma', 'ahc_dir', 'ahc_nbnd', 'ahc_nbndskip',
                                      'skip_upperfan', 'lshift_q', 'zeu', 'zue', 'elop', 'fpol',
                                      'ldisp', 'nogg', 'asr', 'ldiag', 'lqdir', 'search_sym',
                                      'nq1', 'nq2', 'nq3', 'nk1', 'nk2', 'nk3', 'k1', 'k2', 'k3',
                                      'diagonalization', 'read_dns_bare', 'ldvscf_interpolate',
                                      'wpot_dir', 'do_long_range', 'do_charge_neutral', 'start_irr',
                                      'last_irr', 'nat_todo', 'modenum', 'start_q', 'last_q',
                                      'dvscf_star', 'drho_star', ],
                          'LINE': ['xq', 'atom'],
                          # 'qPointsSpecs': ['nqs', 'xq1', 'xq2', 'xq3', 'nq'],
                          }

    def __init__(self, parent_directory, prefix, queue=False, parallel='', **kwargs) -> None:
        super().__init__(parent_directory, prefix, queue=False, parallel='', **kwargs)
