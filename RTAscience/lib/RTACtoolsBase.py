from RTAscience.cfg.Config import Config

class RTACtoolsBase:

    def __init__(self):
        self.caldb = 'prod2'  # production name in calibration database ---!
        self.irf = 'South_0.5h'  # irf ID name ---!
        # data ---!
        self.e = [0.03, 150.0]  # energy range (TeV) ---!
        self.roi = 5  # region of indeterest (deg) ---!
    
    def configure(self, cfg: Config):
        self.caldb = cfg.get('caldb')
        self.irf = cfg.get('irf')
        self.e = [cfg.get('emin'), cfg.get('emax')] 
        self.roi = cfg.get('roi')