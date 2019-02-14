import df


def test():
    #IN CGS
    data_f = df.Datafile()
    data_f.k_max = 51
    data_f.l_max = 4
    data_f.kc_max = 50
    data_f.lc_max = 3
    data_f.dt = 1E-11
    data_f.dt_max = 1E-11
    data_f.dt_factor = 2
    data_f.t0 = 0
    data_f.time_diagnostic = 1E-11
    data_f.time_stop = 1.5E-9
    data_f.sigma_boltzman = 5.6704E-5
    data_f.c = 2.99792458E10
    data_f.a_rad = 4.0 * data_f.sigma_boltzman / data_f.c
    data_f.T0 = 300 #kelvin
    data_f.num_mats = 1
    data_f.bc_type = 1
    data_f.beta = 1.2
    return data_f

def material1():#silicon dioxcide
    data_f = df.Datafile()
    data_f.g = 1.0 / 9175.0
    data_f.f = 8.78E13
    data_f.alpha = 3.53
    data_f.beta = 1.1
    data_f.lambda1 = 0.75
    data_f.mu = 0.09
    data_f.i_start = 0
    data_f.i_end = 1001
    data_f.rho = 50E-3
    data_f.j_start = 0
    data_f.j_end = 4
    return data_f

def material2():
    data_f = df.Datafile()
    data_f.g = 1.0/7200.0
    data_f.f = 3.4E13
    data_f.alpha = 1.5
    data_f.beta = 1.6
    data_f.mu = 0.14
    data_f.lambda1 = 0.2