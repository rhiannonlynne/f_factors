#converted from m5_by_amp.ipynb

import numpy as np
import os
import matplotlib.pyplot as plt
from astropy import units as u
import pandas as pd
from collections import OrderedDict
from scipy.interpolate import interp1d

from lsst.geom import Point2D, Point2I
from lsst.afw.cameraGeom import FIELD_ANGLE, PIXELS
import lsst.afw.cameraGeom.utils as cgUtils
from lsst.daf.persistence import Butler, NoResults

import lsst.syseng.throughputs as st
from lsst.sims.photUtils import PhotometricParameters, Bandpass, LSSTdefaults
from lsst.sims.utils import angularSeparation

DATADIR = f"{os.environ['OBS_LSST_DIR']}/lsstcam/CALIB"
print(DATADIR)
butler = Butler(DATADIR)
cam = butler.get('camera')

vfile = f"{os.environ['HOME']}/notebooks/f_factors/data/vignettingF.txt"
M1D = 8.36 #clear aperture as in Optical design
aa = np.loadtxt(vfile, skiprows=12)
vr = aa[:,0]
vv = aa[:,1]

dd = pd.read_csv('raftInstall.csv',index_col=0)
rnames = []
for i in dd.index:
    try:
        if np.isnan(dd.rtm[i]):
            print(i)
    except TypeError: #when dd.rtm[i] is a str
        rnames.append(i)

exptime=15
nexp=2
othernoise=0
darkcurrent=0.2
X=1.0

lsstDefaults = LSSTdefaults()
addLosses = True
defaultDirs = st.setDefaultDirs()

atmos = st.readAtmosphere(defaultDirs['atmosphere'], atmosFile='atmos_10_aerosol.dat')
mirror1 = st.buildMirror(defaultDirs['mirror1'], addLosses)
mirror2 = st.buildMirror(defaultDirs['mirror2'], addLosses)
mirror3 = st.buildMirror(defaultDirs['mirror3'], addLosses)
lens1 = st.buildLens(defaultDirs['lens1'], addLosses)
lens2 = st.buildLens(defaultDirs['lens2'], addLosses)
lens3 = st.buildLens(defaultDirs['lens3'], addLosses)
filters = st.buildFilters(defaultDirs['filters'], addLosses)

for rname in rnames:
    vendor = dd.vendor[rname]
    vendorDir = defaultDirs['detector']+'/../'+vendor.lower()
    addLosses = False
    detector0 = st.buildDetector(vendorDir, addLosses) #design QE from this vendor
    detector = Bandpass()

    filterlist = tuple([s for s in filters] + ['fS'])
    alist = ('raDeg', 'decDeg', 'radDeg', 'effarea', 'readnoise')
    detectors = []
    for det in cam:
        rname1, dname = det.getName().split('_')
        if rname1 != rname:
            continue;
        detectors.append(det.getName())
    adf = pd.DataFrame(index=alist, columns=detectors, dtype=object)
    m5df = pd.DataFrame(index=filterlist, columns=detectors, dtype=object)
    Tdf = pd.DataFrame(index=filterlist[:6], columns=detectors, dtype=object)
    Sdf = pd.DataFrame(index=filterlist[:6], columns=detectors, dtype=object)

    m5SRD = np.array([23.9, 25.0, 24.7, 24.0, 23.3, 22.1])
    #m5SRDmin = []
    # Nv1 from SRD table 24
    Nv1 = np.array([56, 80, 184, 184, 160, 160])
    omega = Nv1/sum(Nv1)
    fidx = 'ugrizy' #very important!

    for det in cam:
        rname1, dname = det.getName().split('_')
        if rname1 != rname:
            continue;
        key = rname+'_'+dname
        raDeg = {}
        decDeg = {}
        readnoise = {}
        for amp in det:
            i = amp.getName()
            amp_point = amp.getBBox().getCenter()
            raDec = det.transform(amp_point, PIXELS, FIELD_ANGLE)
            [raDeg[i], decDeg[i]] = np.degrees(raDec)
            #print(key, i, raDec)

            readnoise[i] = amp.getReadNoise()
        adf[key].loc['raDeg'] = list(OrderedDict(sorted(raDeg.items())).values())
        adf[key].loc['decDeg'] = list(OrderedDict(sorted(decDeg.items())).values())
        adf[key].loc['readnoise'] = list(OrderedDict(sorted(readnoise.items())).values())

        #effetive area
        radius = angularSeparation(0., 0., adf[key]['raDeg'], adf[key]['decDeg'])
        adf[key].loc['radDeg'] = list(radius)
        adf[key].loc['effarea'] = list(np.interp(radius, vr, vv)*np.pi*(M1D/2)**2)

    ampList = list(OrderedDict(sorted(raDeg.items())).keys())

    # this is needed if there are only 6 QE measurements per amp
    idx = np.where(detector0.sb>0.01)
    idx1=idx[0][0]-1
    idx2=idx[0][-1]+1

    x1 = detector0.wavelen[idx1]
    y1 = detector0.sb[idx1]
    x2 = detector0.wavelen[idx2]
    y2 = detector0.sb[idx2]

    for det in cam:
        rname1, dname = det.getName().split('_')
        if rname1 != rname:
            continue;

        vendor = det.getSerial()[:3].lower()
        assert dd.vendor[rname].lower() == vendor
        vendorDir = defaultDirs['detector']+'/../'+vendor
        print('Calculating m5 for %s_%s'%(rname,dname))

        key = rname+'_'+dname
        for f in filterlist:
            m5df[key][f] = [-1.]*len(ampList)
        for f in filterlist[:6]:
            Tdf[key][f] = [-1.]*len(ampList)
            Sdf[key][f] = [-1.]*len(ampList)

        for amp in det:

            try:
                qe_curve = butler.get('qe_curve', raftName=rname, detectorName=dname, calibDate='1970-01-01T00:00:00')
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
                    sb = 0
                if np.max(sb)<0.01:
                    print('dead channel: %s, max sb = %.2f'%(key, np.max(sb)))
                    continue;

                detector.setBandpass(wavelen, sb)

                #detector losses
                #os.listdir(vendorDir)
                detLosses = Bandpass()
                detLosses.readThroughput(os.path.join(vendorDir, '%s_Losses/det_Losses.dat' % (vendor)))

                #build hardware and system
                hardware = {}
                system = {}
                for f in filters:
                    sb = mirror1.sb * mirror2.sb *mirror3.sb
                    sb *= lens1.sb * lens2.sb * lens3.sb * filters[f].sb
                    sb *= detector.sb * detLosses.sb

                    hardware[f] = Bandpass()
                    hardware[f].setBandpass(wavelen, sb)
                    system[f] = Bandpass()
                    system[f].setBandpass(wavelen, sb * atmos.sb)

            except NoResults:
                if not useDefault:
                    print('No results found for this detector')
                assert useDefault

                hardware, system = st.buildHardwareAndSystem(defaultDirs)

            #calculate m5
            iamp = ampList.index(amp.getName())
            effarea = adf[key]['effarea'][iamp]*100**2 #convert to cm^2
            readnoise = adf[key]['readnoise'][iamp]

            m5 = st.makeM5(hardware, system, darksky=None,
                        exptime=15, nexp=2, readnoise=readnoise, othernoise=0, darkcurrent=0.2,
                        effarea=effarea, X=1.0)
            for f in filters:
                m5df[key][f][iamp] = m5.m5[f]
                Tdf[key][f][iamp] = m5.Tb[f]
                Sdf[key][f][iamp] = m5.Sb[f]
            m5amp = np.array([m5.m5[f] for f in fidx])
            if np.all(m5amp>0):
                m5df[key]['fS'][iamp] = sum(omega*10**(0.8*(m5amp - m5SRD)))

    dfDir = os.path.join('m5_output', rname)
    if not os.path.exists(dfDir):
        os.mkdir(dfDir)
    dfPath = os.path.join(dfDir, 'adf_%s.csv'%rname)
    adf.to_csv(dfPath)
    dfPath = os.path.join(dfDir, 'm5df_%s.csv'%rname)
    m5df.to_csv(dfPath)
    dfPath = os.path.join(dfDir, 'Tdf_%s.csv'%rname)
    Tdf.to_csv(dfPath)
    dfPath = os.path.join(dfDir, 'Sdf_%s.csv'%rname)
    Sdf.to_csv(dfPath)
    