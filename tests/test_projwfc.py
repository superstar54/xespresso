def test_projwfc():
    from xespresso.post.projwfc import EspressoProjwfc

    projwfc = EspressoProjwfc(
        parent_directory="calculations/scf/co", prefix="co", DeltaE=0.01
    )
    projwfc.run()


def test_pdos_analysis():
    from xespresso.dos import DOS
    import matplotlib.pyplot as plt

    # DOS analysis
    dos = DOS(label="calculations/scf/co", prefix="co")
    dos.read_pdos()
    dos.plot_pdos(Emin=-10, Emax=10, smearing=[0.02, 0.01], legend=True)
    plt.savefig("images/co-pdos.png")
