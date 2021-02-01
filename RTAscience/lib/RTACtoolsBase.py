from RTAscience.cfg.Config import Config

class RTACtoolsBase:

    def __init__(self):
        self.caldb = 'prod2'  # production name in calibration database ---!
        self.irf = 'South_0.5h'  # irf ID name ---!

        # data ---!
        self.tmax = 1800  # maximum exposure time needed (s) ---!
        self.e = [0.03, 150.0]  # energy range (TeV) ---!
        self.roi = 5  # region of indeterest (deg) ---!
    
    def configure(self, cfg: Config):
        self.caldb = cfg.get('caldb')
        self.irf = cfg.get('irf')
        self.tmax = cfg.get('tobs') - cfg.get('onset') # time of bkg only a.k.a. delayed onset of burst (s)
        self.e = [cfg.get('emin'), cfg.get('emax')] 
        self.roi = cfg.get('roi')