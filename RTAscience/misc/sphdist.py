import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from os.path import join

path = '/home/ambra/Desktop/CTA/projects/rta-pipe/archive_tests/paper2021/newsim_data/'

total = join(path, 'deg_flux1-3_off2-4_del50-150.txt')    
data = pd.read_csv(total, sep=' ')
print(len(data))

# SphDistance calc ---!
compute = 'sph_dist' not in data.keys()
print(f"Check if sph. dist. is missing from DF: {compute}")
if compute:
    print('add sph. dist. to DF')
    trueRA = 33.057
    trueDEC = -51.841
    true_coord = SkyCoord(ra = trueRA*u.deg, dec = trueDEC*u.deg, frame='fk5')
    dist = []
    for i in range(len(data)):
        print(i+1)
        dist.append(float(true_coord.separation(SkyCoord(ra=data['ra'][i]*u.deg, dec=data['dec'][i]*u.deg, frame='fk5')).deg))
    print(len(dist))
    data['sph_dist'] = dist
    data.to_csv(join(path, 'deg_flux1-3_off2-4_del50-150_updated.txt'), sep=' ', index=False, header=True)

print('exit')

