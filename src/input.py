import noddy2


def test():
    mynod = noddy2.Noddy()
    mynod.k_max = 78
    mynod.l_max = 76
    mynod.kc_max = 77
    mynod.lc_max = 75
    mynod.dt = 0.01
    mynod.t0 = 0
    mynod.time_stop = 0.03
    return mynod