import os
import argparse
import copy
import numpy as np
import pandas as pd
from astropy import units as u
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

import syseng_throughputs as st
import rubin_sim.photUtils as photUtils
import rubin_sim.site_models as site_models

from lsst.daf.persistence import Butler, NoResults
from lsst.geom import Point2D, Point2I
from lsst.afw.cameraGeom import FIELD_ANGLE, PIXELS
from lsst.obs.lsst.lsstCamMapper import LsstCamMapper


def calcM5(hardware, system, skysed, skymags, fwhmEff, photParams):
    """Calculate m5 given the hardware, system and sky sed.
    Also, skymags (dictionary or dataframe), fwhmEff (dictionary or dataframe), and single photParams.

    Call with photParams infinity (no noise) to calculate m5_infinity.
    """
    sky_copy = copy.deepcopy(skysed)
    m5 = {}
    for i, f in enumerate(skymags):
        fluxNorm = sky_copy.calcFluxNorm(skymags[f], hardware[f])
        sky_copy.multiplyFluxNorm(fluxNorm)
        m5[f] = photUtils.calcM5(sky_copy, system[f], hardware[f],
                                 photParams, FWHMeff=fwhmEff[f])
    return m5


def readQE():

    detectors = {}
    vendors = {}

    dd = pd.read_csv('raftInstall.csv', index_col=0)
    defaultDirs = st.setDefaultDirs()

    vendorQE = {}
    vendorQE['e2v'] = st.buildDetector(defaultDirs['detector'] + '/../e2v', addLosses=False)
    vendorQE['ITL'] = st.buildDetector(defaultDirs['detector'] + '/../itl', addLosses=False)

    for det in cam:
        key = det.getName()
        raft, chip = key.split('_')
        if raft in ['R00', 'R04', 'R40', 'R44']:
            continue

        vendor = dd.vendor[raft]
        # Get the design QE
        detector0 = copy.deepcopy(vendorQE[vendor])  # design QE from this vendor

        detectors[f'{raft}_{chip}'] = {}

        ## Get QE from butler/amp info (have to interpolate to build)
        # this is needed if there are only 6 QE measurements per amp
        idx = np.where(detector0.sb > 0.01)
        idx1 = idx[0][0] - 1
        idx2 = idx[0][-1] + 1

        x1 = detector0.wavelen[idx1]
        y1 = detector0.sb[idx1]
        x2 = detector0.wavelen[idx2]
        y2 = detector0.sb[idx2]

        print('Reading QE for %s_%s' % (raft, chip))
        qe_curve = butler.get('qe_curve', raftName=raft, detectorName=chip,
                              calibDate='1970-01-01T00:00:00')
        for amp in det:
            wavelen = detector0.wavelen

            k = amp.getName()
            if len(qe_curve.data[k][0]) > 10:
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
                # quadratic causes overshoot here
                f = interp1d(wlen.value, eff.value, fill_value=0, bounds_error=False, kind='slinear')

            sb = f(wavelen) * 0.01
            sb[np.isnan(sb)] = 0
            if np.max(sb) > 1.5:
                print('These seem too LARGE ', k)
                print(np.max(sb))
                sb = 0
            if np.max(
                    sb) < 0.2:  # 3 dead channels, 1 out of each of R01, R10, and R30; see camera confluence page table
                print('dead channel: %s %s, max sb = %.2f' % (key, amp.getName(), np.max(sb)))
                continue;

            detectors[f'{raft}_{chip}'][k] = Bandpass()
            detectors[f'{raft}_{chip}'][k].setBandpass(wavelen, sb)
            vendors[f'{raft}_{chip}'] = vendor
            iamp = ampList.index(amp.getName())