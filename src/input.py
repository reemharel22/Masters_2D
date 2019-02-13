import df


def test():
    #IN CGS
    data_f = df.Datafile()
    data_f.k_max = 78
    data_f.l_max = 76
    data_f.kc_max = 77
    data_f.lc_max = 75
    data_f.dt = 0.01
    data_f.t0 = 0
    data_f.time_diagnostic = 1
    data_f.time_stop = 0.03
    data_f.sigma_boltzman = 5.6704E-5
    data_f.c = 2.99792458E10
    data_f.a_rad = 4.0 * data_f.sigma_boltzman / data_f.c
    data_f.g = 1.0 / 9175.0
    data_f.f = 8.78E13
    data_f.alpha = 3.53
    data_f.beta = 1.1
    data_f.lambda1 = 0.75
    data_f.mu = 0.09
    data_f.T0 = 300 #kelvin
    data_f
    return data_f
