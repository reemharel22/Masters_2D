import df


def test():
    data_f = df.Datafile()
    data_f.k_max = 78
    data_f.l_max = 76
    data_f.kc_max = 77
    data_f.lc_max = 75
    data_f.dt = 0.01
    data_f.t0 = 0
    data_f.time_diagnostic = 1
    data_f.time_stop = 0.03
    return data_f
