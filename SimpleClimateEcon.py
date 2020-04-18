"""
Andrew Williams -- April 2020

Class-based implementation of the simple CBA 
model developed in the Jupyter Notebooks.

Currently no carbon-cycle feedbacks at play or sequestration tech,
so can -at best- only get stable temps and not declining temps.
"""

class CBA:
    def __init__(self, W0=75*10**12, D0=0.00267, pback=200, gamma=2,
                 gsig=-0.0125, TCRE=0.00054, ks=0.12, g=0.02, r=0.035,
                 E0=39, T0=0.9):
        # Model parameters
        self.W0    = W0
        self.T0    = T0
        self.E0    = E0
        self.sig0  = E0/W0
        self.D0    = D0
        self.pback = pback
        self.gamma = gamma
        self.gsig  = gsig
        self.TCRE  = TCRE
        self.ks    = ks
        self.kd    = r-g
        self.g     = g
        
    def increment_T(self, T, dTdt):
        return T+dTdt

    def dTdt(self, E):
        """
        TCRE relationship between instantaneous emissions 
        and annual rate of change of temperature.
        """
        return self.TCRE*E
    
    def consumption(self, year):
        return self.W0 * np.exp(self.g*(year-2005))

    def SCC(self, year, T, dTdt):
        """
        Social Cost of Carbon ($/tCO2)

        Derived in text.
        """
        prefactor=self.gamma*self.TCRE*self.D0*consumption(year)* T**(self.gamma-1)
        bracket  = np.divide(1, self.kd) - np.divide(1,self.ks+self.kd) 
        + np.divide((self.gamma-1)*dTdt, T)*((1/self.kd)**2 - (1/(self.ks+self.kd))**2)
        return prefactor*bracket*10**-9

    def abatement_rate(self, SCC):
        """
        Abatement rate, used in emissions equation.
        E = sigma * (1-abatement_rate) * W

        Abatement rate is a function of SCC. **NOT SURE WHAT STUART USES?**
        """
        test = 0.003*SCC
        if test >=1:
            return 1
        else:
            return test

    def emissions(self, year, abate, SCC):
        """
        E = sigma * (1-abatement_rate) * W

        sigma: ghg emissions/total output declines 
        exogeneously w time @ 1.25%/yr
        """

        sigma = self.sig0 * np.exp(self.gsig*(year-2005)) 
        W     = self.W_0  * np.exp(self.g*(year-2005))
        return sigma * (1-abate) * W

    def initialize(self, nyears):
        # Initialise arrays w initial values
        self.nyears = nyears
        self.years  = np.linspace(2005, 2005+self.nyears, self.nyears+1)
        self.temps  = np.zeros(self.nyears+2)
        self.dTdt   = np.zeros(self.nyears+2)
        self.abate  = np.zeros(self.nyears+2)
        self.emiss  = np.zeros(self.nyears+2)
        self.SCC    = np.zeros(self.nyears+2)
        
        self.temps[0] = self.T0
        self.emiss[0] = self.E0
            
    def forward_integration(self, nyears):
        
        self.initialize(nyears)
        
        """Run main loop"""
        for idx, year in enumerate(self.years):
            self.dTdt[idx]    = dTdt(self.emiss[idx])
            self.temps[idx+1] = increment_T(self.temps[idx], self.dTdt[idx])
            self.SCC[idx]     = SCC(year, self.temps[idx], self.dTdt[idx])
            self.abate[idx]   = abatement_rate(self.SCC[idx])
            self.emiss[idx+1] = emissions(year, self.abate[idx], self.SCC[idx])

            #if self.emiss[idx+1]==0:
            #    self.temps[idx+1:]=np.nan
            #    self.dTdt[idx+1:]=np.nan
            #    self.SCC[idx+1:]=np.nan
            #    self.abate[idx+1:]=np.nan
            #    break
                
        """Put outputs into nice dictionary"""
        self.outps = {'Temperature anomaly': self.temps, 
                      'dTdt': self.dTdt,
                      'Emissions': self.emiss, 
                      'Abatement rate': self.abate,
                      'Social Cost of Carbon': self.SCC}
    
    def plot(self, ax=None, fsize=(15,7), dpi_=300):
        
        label_dict = {'Temperature anomaly': r'$\Delta$ T [K]',
                      'dTdt': r'dTdt [K yr$^{-1}$]',
                      'Emissions': r'Emissions [GtCO$_{2}$ yr$^{-1}$]',
                      'Abatement rate': 'Abatement rate',
                      'Social Cost of Carbon': r'SCC [\$ /tCO$_{2}$]'}
        
        ylim_dict  = {'Temperature anomaly': (-0.1, 5),
                      'dTdt': (-0.001, 0.03),
                      'Emissions': (-5, 60),
                      'Abatement rate': (-0.1, 1.1),
                      'Social Cost of Carbon': (0, 500)}

        if ax is None:
            fig, axs = plt.subplots(nrows=2, ncols=3, dpi=dpi_, figsize=fsize)
                
        for idx, ax in enumerate(axs.flatten()):
            if idx == 5:
                ax.set_visible(False)
                continue
            else:
                ax.plot(self.years, list(self.outps.items())[idx][1][:-1])
                ax.set_title(list(self.outps.items())[idx][0], fontsize=15)
                ax.set_ylim(list(ylim_dict.values())[idx])
                ax.set_ylabel(list(label_dict.values())[idx], fontsize=12)
                ax.tick_params(axis='both', which='major', labelsize=10)
                
        fig.tight_layout()
        return axs # For external manipulation, if needed -- otherwise remove 