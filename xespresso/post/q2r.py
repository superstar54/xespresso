from xespresso.post.base import PostCalculation


class EspressoQ2r(PostCalculation):

    package = 'q2r'
    package_parameters = {'INPUT': ['fildyn', 'flfrc', 'zasr', 'loto_2d'], }

    def __init__(self, parent_directory, prefix, queue=False, parallel='', **kwargs) -> None:
        super().__init__(parent_directory, prefix, queue=queue, parallel=parallel, **kwargs)
