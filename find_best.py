def get_boltz_prob(seq, ss, scale=True):
    """viennaRNA boltzmann probability"""
    fc = RNA.fold_compound(seq)
    if scale:
        _, mfe = fc.mfe()
        fc.exp_params_rescale(mfe)
    fc.pf()
    pr = fc.pr_structure(ss)
    return pr