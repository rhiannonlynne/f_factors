import numpy as np
import os
import pandas as pd
from scipy.interpolate import interp1d

from lsst.afw.cameraGeom import FIELD_ANGLE, PIXELS
from lsst.daf.persistence import Butler, NoResults

from rubin_sim.photUtils import PhotometricParameters, Bandpass, Sed, LSSTdefaults
from rubin_sim.utils import angularSeparation
import syseng_throughputs as st


# Read substitute readnoise files
#Run 13040 has 2.232 sec Readout
#Run 13057 has 2.090 sec Readout
#Run 13060 has 2.374 sec Readout

series = '13057'
series = '13040'
filename = f'readnoise_{series}.pkl'
readnoise_val = pd.read_pickle(filename)

raftNames = []

DATADIR = f"{os.environ['OBS_LSST_DIR']}/lsstcam/CALIB"
print(DATADIR)
butler = Butler(DATADIR)
cam = butler.get('camera')


detectors = []
rafts = []
chips = []
for det in cam:
    raftName, chipName = det.getName().split('_')
    # Ignore the corner rafts
    if raftName in ['R00', 'R04', 'R40', 'R44']:
        continue
    detectors.append(det.getName())
    rafts.append(raftName)
    chips.append(chipName)
# de-duplicate the raft/chip names
# ... we might not need these really
rafts = list(set(rafts))
chips = list(set(chips))

# this sets up a default QE curve to fall back to, as it maps raft -> vendor in dd.vendor[raftname]
# mapping based on 
#https://confluence.slac.stanford.edu/pages/viewpage.action?spaceKey=LSSTCAM&title=Raft+Delivery+and+Acceptance+Testing+Status
dd = pd.read_csv('raftInstall.csv',index_col=0)


# ### Default photometric parameters, as used in standard m5 calculations
# Note that effarea is not in this list here, because it varies with field.
# The read noise is not in this list either, because it varies by amp.

exptime=15
nexp=2
othernoise=0
darkcurrent=0.2
X=1.0

skymag_sim = {'u': 22.68,
           'g': 22.11,
           'r': 21.11,
           'i': 20.39,
           'z': 19.43,
           'y': 18.63}

seeing_sim = 0.72 # fwhm500 with standard seeing model
# also u = 30s visits
# also we add losses to all components

#we do not need this cell to calculate m5, but these FWHMeff are the default values we use actually
lsstDefaults = LSSTdefaults()

# ### Set up throughputs for hardware and atmosphere.
addLosses = True
defaultDirs = st.setDefaultDirs()

darksky = Sed()
darksky.readSED_flambda(os.path.join(defaultDirs['atmosphere'], 'darksky.dat'))
atmos = st.readAtmosphere(defaultDirs['atmosphere'], atmosFile='atmos_10_aerosol.dat')
mirror1 = st.buildMirror(defaultDirs['mirror1'], addLosses)
mirror2 = st.buildMirror(defaultDirs['mirror2'], addLosses)
mirror3 = st.buildMirror(defaultDirs['mirror3'], addLosses)
lens1 = st.buildLens(defaultDirs['lens1'], addLosses)
lens2 = st.buildLens(defaultDirs['lens2'], addLosses)
lens3 = st.buildLens(defaultDirs['lens3'], addLosses)
filters = st.buildFilters(defaultDirs['filters'], addLosses)
detLosses = Bandpass()
# Losses for e2v and ITL are currently modeled as the same, just use ITL
detLosses.readThroughput(os.path.join(defaultDirs['detector'], 'joint_losses/det_Losses.dat'))
# combine common elements:
hardware_nodet = {}
system_nodet = {}
wavelen = atmos.wavelen

for f in filters:
    sb = mirror1.sb * mirror2.sb * mirror3.sb
    sb *= lens1.sb * lens2.sb * lens3.sb * filters[f].sb
    if addLosses:
        sb *= detLosses.sb
    hardware_nodet[f] = Bandpass()
    hardware_nodet[f].setBandpass(wavelen, sb)
    system_nodet[f] = Bandpass()
    system_nodet[f].setBandpass(wavelen, sb * atmos.sb)

m5SRD = np.array([23.9, 25.0, 24.7, 24.0, 23.3, 22.1])
# Nv1 from SRD table 24
Nv1 = np.array([56, 80, 184, 184, 160, 160])
omega = Nv1/sum(Nv1)
fidx = 'ugrizy' #very important!


# ### Prepare vignetting function

# below we use v3.11 values
vfile = f"{os.environ['HOME']}/other_repos/f_factors/data/vignettingF.txt"
M1D = 8.36 #clear aperture as in Optical design
aa = np.loadtxt(vfile, skiprows=12)
vr = aa[:,0]
vv = aa[:,1]


# ### Create the dataframes for the camera

filterlist = ['u', 'g', 'r', 'i', 'z', 'y', 'fS', 'u_30', 'g_30', 'r_30', 'i_30', 'z_30', 'y_30']
alist = ('raDeg', 'decDeg', 'radDeg', 'effarea', 'readnoise', 'gain', 'saturation')

adf = pd.DataFrame(index=alist, columns=detectors, dtype=object)
m5df = pd.DataFrame(index=filterlist, columns=detectors, dtype=object)
Tdf = pd.DataFrame(index=filterlist[:6], columns=detectors, dtype=object)
Sdf = pd.DataFrame(index=filterlist[:6], columns=detectors, dtype=object)


# Read the data from the butler and camera model about location (thus vignetting) and readnoise

for det in cam:
    key = det.getName()
    raft, chip = key.split('_')
    if raft in ['R00', 'R04', 'R40', 'R44']:
        continue

    vendor = dd.vendor[raft]
    vendorDir = defaultDirs['detector'] + '/../' + vendor.lower()
    addLosses = False
    # Get the design QE
    detector0 = st.buildDetector(vendorDir, addLosses)  # design QE from this vendor
    detector = Bandpass()

    # these values are filled in per amp for each raft/chip combo, then overwritten
    raDeg = {}
    decDeg = {}
    readnoise = {}
    gain = {}
    saturation = {}
    for amp in det:
        i = amp.getName()
        amp_point = amp.getBBox().getCenter()
        raDec = det.transform(amp_point, PIXELS, FIELD_ANGLE)
        [raDeg[i], decDeg[i]] = np.degrees(raDec)
        gain[i] = amp.getGain()
        saturation[i] = amp.getSaturation()
        # print(key, i, raDec, amp.getGain(), amp.getSaturation())
        readnoise[i] = amp.getReadNoise()
    adf[key].loc['raDeg'] = list(raDeg.values())
    adf[key].loc['decDeg'] = list(decDeg.values())
    #adf[key].loc['readnoise'] = list(readnoise.values())
    adf[key].loc['gain'] = list(gain.values())
    adf[key].loc['saturation'] = list(saturation.values())

    # Substitute in the readnoise from the pickle file
    adf[key].loc['readnoise'] = readnoise_val[key]
    updated = (adf[key].loc['readnoise'] == readnoise_val[key])
    if not updated:
        print(key)

    # effective area
    radius = angularSeparation(0., 0., adf[key]['raDeg'], adf[key]['decDeg'])
    adf[key].loc['radDeg'] = list(radius)
    adf[key].loc['effarea'] = list(np.interp(radius, vr, vv) * np.pi * (M1D / 2) ** 2)


    ## Get QE from butler/amp info (have to interpolate to build)
    # this is needed if there are only 6 QE measurements per amp
    idx = np.where(detector0.sb>0.01)
    idx1=idx[0][0]-1
    idx2=idx[0][-1]+1

    x1 = detector0.wavelen[idx1]
    y1 = detector0.sb[idx1]
    x2 = detector0.wavelen[idx2]
    y2 = detector0.sb[idx2]

    vendor = det.getSerial()[:3].lower()
    vendor = dd.vendor[raft].lower()
    vendorDir = defaultDirs['detector']+'/../'+vendor
    print('Calculating m5 for %s_%s'%(raft,chip))

    ampList = list(raDeg.keys())
    for f in filterlist:
        m5df[key][f] = [-1.]*len(ampList)
    for f in filterlist[:6]:
        Tdf[key][f] = [-1.]*len(ampList)
        Sdf[key][f] = [-1.]*len(ampList)

    try:
        qe_curve = butler.get('qe_curve', raftName=raft, detectorName=chip,
                              calibDate='1970-01-01T00:00:00')
    except NoResults:
        print(f'AHHH {chip} {raft} failed to get butler qe data, using defaults')
        hardware = {}
        system = {}
        for f in filters:
            sb = detector0.sb * hardware_nodet[f].sb
            hardware[f] = Bandpass()
            hardware[f].setBandpass(wavelen, sb)
            sb = detector.sb * system_nodet[f].sb
            system[f] = Bandpass()
            system[f].setBandpass(wavelen, sb)
    for amp in det:
        wavelen = detector0.wavelen

        k = amp.getName()
        if len(qe_curve.data[k][0])>10:
            wlen = qe_curve.data[k][0]
            eff = qe_curve.data[k][1]
            f = interp1d(wlen.value, eff.value, fill_value=0, bounds_error=False, kind='quadratic')
        else:
            aa = np.append(x1, qe_curve.data[k][0].value)
            aa = np.append(aa, x2)
            wlen = aa * qe_curve.data[k][0].unit

            aa = np.append(y1, qe_curve.data[k][1].value)
            aa = np.append(aa, y2)
            eff = aa * qe_curve.data[k][1].unit
            f = interp1d(wlen.value, eff.value, fill_value=0, bounds_error=False, kind='slinear')#quadratic causes overshoot

        sb = f(wavelen)*0.01
        #alternatively we could do (only for >10 QE measurements)
        #amp_point = amp.getBBox().getCenter()
        #sb = qe_curve.evaluate(det, amp_point, wavelen* u.nm, kind='quadratic').value*.01 #unit was percent in CALIB data

        sb[np.isnan(sb)] = 0
        if np.max(sb)>1.5:
            print('These seem too LARGE ', k)
            print(np.max(sb))
            sb = detector0.sb
        if np.max(sb)<0.2: #3 dead channels, 1 out of each of R01, R10, and R30; see camera confluence page table
            print('dead channel: %s %s, max sb = %.2f'%(key, amp.getName(), np.max(sb)))
            sb = detector0.sb

        detector.setBandpass(wavelen, sb)

        #build hardware and system *with* detector
        hardware = {}
        system = {}
        for f in filters:
            sb = detector.sb * hardware_nodet[f].sb
            hardware[f] = Bandpass()
            hardware[f].setBandpass(wavelen, sb)
            sb = detector.sb * system_nodet[f].sb
            system[f] = Bandpass()
            system[f].setBandpass(wavelen, sb)
            
        #calculate m5      
        iamp = ampList.index(amp.getName())
        effarea = adf[key]['effarea'][iamp]*100**2 #convert to cm^2
        readnoise = adf[key]['readnoise'][iamp]

        othernoise = 0
        darkcurrent = 0.2

        #effarea = np.pi*(6.423/2*100)**2
        m5 = st.makeM5(hardware, system, darksky=darksky, sky_mags=skymag_sim,
                    exptime=15, nexp=2, readnoise=readnoise, othernoise=othernoise, darkcurrent=darkcurrent,
                    effarea=effarea, X=1.0, fwhm500=seeing_sim)
        for f in filters:
            m5df[key][f][iamp] = m5.m5[f]
            Tdf[key][f][iamp] = m5.Tb[f]
            Sdf[key][f][iamp] = m5.Sb[f]
        m5amp = np.array([m5.m5[f] for f in fidx])
        if np.all(m5amp>0):
            m5df[key]['fS'][iamp] = sum(omega*10**(0.8*(m5amp - m5SRD)))
        
        # 1x30s visits
        m5 = st.makeM5(hardware, system, darksky=darksky, sky_mags=skymag_sim,
                    exptime=30, nexp=1, readnoise=readnoise, othernoise=othernoise, darkcurrent=darkcurrent,
                    effarea=effarea, X=1.0, fwhm500=seeing_sim)
        for f in filters:
            m5df[key][f'{f}_30'][iamp] = m5.m5[f]


# write output
dfDir = os.path.join('m5_output')
suffix = 'sim' #_novignetting'
if not os.path.exists(dfDir):
    os.mkdir(dfDir)
dfPath = os.path.join(dfDir, f'adf_{series}_{suffix}.csv')
adf.to_csv(dfPath)
dfPath = os.path.join(dfDir, f'm5df_{series}_{suffix}.csv')
m5df.to_csv(dfPath)
dfPath = os.path.join(dfDir, f'Tdf_{series}_{suffix}.csv')
Tdf.to_csv(dfPath)
dfPath = os.path.join(dfDir, f'Sdf_{series}_{suffix}.csv')
Sdf.to_csv(dfPath)
