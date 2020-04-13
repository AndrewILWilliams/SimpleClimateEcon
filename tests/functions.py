import numpy as np
from scipy import stats
from scipy import integrate
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt

import fair
from fair.forward import fair_scm
from fair.RCPs import rcp3pd, rcp45, rcp85

## 10/02/2019 --- Useful functions, updated regularly.


def func1(x, a, b): # use for linear regressions in SciPy
    return a + b*x 

def func2(x, a, b, c): # use for linear regressions in SciPy
    return a + b*x + c*x**2 

def func3(x, a, b, c, d): # use for linear regressions in SciPy
    return a + b*x + c*x**2 + d*x**3

def running_mean(x, N):
    out = np.zeros_like(x, dtype=np.float64)
    dim_len = len(x)
    for i in range(dim_len):
        if N%2 == 0:
            a, b = i - (N-1)//2, i + (N-
                                      1)//2 + 2
        else:
            a, b = i - (N-1)//2, i + (N-1)//2 + 1

        #cap indices to min and max indices
        a = max(0, a)
        b = min(dim_len, b)
        out[i] = np.mean(x[a:b])
    return out

def fettweis(T): # fettweis et al 2013 parameterization of GrIS SMB as a function of globas TAS relative to 1980-1999
    return -71.5 * (T) - 20.4 * (T**2) - 2.8 * (T**3)

def differentiate(x, dt): # differentiate timeseries
    diff    = np.zeros(len(x))
    
    for i in range(len(x)):
        #if i ==0:
        #    diff[i] = x[0]     
        if i>0:
            diff[i] = np.divide((x[i] - x[i-1] ), dt)
    return diff

def thermosteric1(T, T_res, eps, kappa): # generate thermosteric sea level rise
    '''
    T: temperature pathway as a function of time (yearly FAIR output)
    T_res: time-resolution of temp series in seconds (use 3600*365*24 if using FAIR data)
    t: associated time-series (yearly FAIR output) (repressed for now)
    eps: expansion efficiency of heat (0.11 +- 0.01) m YJ^{-1}
    kappa: ocean heat uptake efficiency (~ 0.75, need to narrow this down)
    '''
    
    A_ocean    = 3.6704*10**14
    #A_world    = 5.10 * 10**14
    
    pre_factor = eps*kappa*A_ocean
    
    h          = integrate.cumtrapz(T, dx = T_res, initial = 0)
    
    return pre_factor*h

def thermosteric2(F, T, eps, alpha):
    '''
    F: Forcing timeseries (yearly FAIR output)
    T: Temperature timeseries (yearly FAIR output)
    eps: expansion efficiency of heat (0.11 +- 0.01 m YJ^{-1} )
    alpha: climate sensitivity (1.13 +- 0.325 Wm^-2 K^-1)
    '''
    
    A_ocean    = 3.6704*10**14

    N          = F - alpha*T # net energy flux into the climate system, assume goes entirely into ocean

    pre_factor = eps*A_ocean

    h          = integrate.cumtrapz(N, dx = 3600*24*365, initial = 0)

    return pre_factor*h

def glaciers(T, f, p):
    '''
    T: Temperature timeseries (yearly FAIR output)
    f, p: fitting parameters based of parameterization from de Vries et al (2014)
        : f has units of [mm K^-1 yr^-1]
        : p is dimensionless
    '''

    f = f*10**(-3)

    #dT = T - T[2006-1765]
    dT = T
    I = integrate.cumtrapz(dT, dx = 1, initial = 0)

    return f*(np.float_power(I,p)) + 9.5*10**(-3)

def AIS_smb(T, ref, P_amp, perc_inc): # sea level change due to AIS SMB (not partitioned)
    '''
    T: temperature time-series
    ref: reference SMB value, Church et al 2013 (1983+-122 Gt/yr)
    P_amp: South-polar amplification factor (1.1 +- 0.2)
    perc_inc = fractional increase in precipitation (~ 5.1% +- 1.5% /K)  
    '''
    #ref   = ref*10**12 # reference SMB value in kg

    delta = integrate.cumtrapz(T, dx = 1, initial = 0)

    ## TO-DO: Check if need to include the - 0.35 Normal_dist term from AR5 to match Luke's paper!!! (Maybe he just missed it out in the paper)

    return 0.825*ref*np.divide(1,361.8*10**3)*P_amp*(np.divide(perc_inc, 100))*delta

    ## WHY IS THERE A FACTOR OF 0.1 HERE?!? (Below, in the return.)
    #return 0.85*np.divide(ref, 3.67*10**17)*P_amp*(1+np.divide(perc_inc, 100))*delta

def GrIS_smb(T):
    
    delta = integrate.cumtrapz(fettweis(T)+420, dx = 1, initial = 0)

    return np.divide(1,361.8*10**3)*delta

def discharge(T, a, b):
    '''
    T: temperature time-series
    T_res: time-resolution of temp series in seconds (use 3600*365*24 if using FAIR data)
    a: coefficient
    b: coefficient
    '''

    delta = integrate.cumtrapz(a+b*T, dx = 1, initial = 0)

    return np.divide(1,361.8*10**3)*delta

def idealised_RCP(t_c, t_0, dE):
    '''
    Description:
    Outputs RCP2.6 emissions scenario but with the anthropogenic CO2 forcing changed.
    ---
    Inputs:
    t_c: Year of peak emissions (follow RCP8.5 until then)
    t_0: Year of net zero
    dE: magnitude of (constant) reduction in emissions (-ve number) -- for RCP2.6, dE ~1
    N.B. t_0 MUST be greater than t_c !!!
    ---
    Outputs:
    RCP2.6 emissions timeseries (multiple forcing agents), but with CO2_fossil edited
    '''

    copy   = np.multiply(rcp85.Emissions.co2_fossil, 1)
    output = rcp3pd.Emissions

    if t_0 < t_c:
        print('Error: Time to net zero emissions needs to be LATER than time to peak emissions.')
        exit()
    
    dE     = float(dE) # Magnitude of (constant) reduction in atmospheric carbon (-ve number) -- for RCP2.6, dE ~ -1
   
    t_c_i  = t_c - 1765 # year index for peak emissions
    t_0_i  = t_0 - 1765 # year index for net zero emissions

    # find year at which we reach emissions target, dE, from peak emissions at Tc

    t_e    = t_c + (t_0-t_c)*(1-np.divide(dE, rcp85.Emissions.co2_fossil[t_c_i]))
    t_e    = np.floor(t_e) 

    t_e_i  = t_e - 1765 # year index for reaching emissions target, dE.
    
    for j in np.linspace(t_c_i, t_e_i, t_c_i).astype(int):
        copy[j] = copy[t_c_i]*(1 - np.divide((j-t_c_i),(t_0_i-t_c_i)))

    copy[int(t_e_i):] = dE

    output.co2_fossil[:] = copy[:]
    
    return output











    


