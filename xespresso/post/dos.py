from xespresso.post.base import PostCalculation


class EspressoDos(PostCalculation):
    """_summary_

    Args:
        PostCalculation (_type_): _description_
    """

    package = "dos"
    package_parameters = {
        "DOS": [
            "prefix",
            "outdir",
            "bz_sum",
            "ngauss",
            "degauss",
            "Emin",
            "Emax",
            "DeltaE",
            "fildo",
        ]
    }

    def __init__(
        self, parent_directory, prefix, queue=False, parallel="", **kwargs
    ) -> None:
        super().__init__(
            parent_directory, prefix, queue=queue, parallel=parallel, **kwargs
        )
