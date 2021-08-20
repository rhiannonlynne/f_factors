The code you find here enables us to evaluate m5 values for individual amplifiers and make plots to visualize the variation on the focal plane.

In each optical band, the m5 variation with amplifier takes into account 3 factors
  1. Vignetting ratio
  2. Detector QE
  3. Detector read noise

These have been verified to run on the LSP (https://lsst-lsp-stable.ncsa.illinois.edu).

The notebooks are mostly for running things for individual rafts step by step, so that we can examine results along the way. The main notebook is `m5_by_amp.ipynb`, which reads out the QE curves and read noise from DM butler, and uses code from `syseng_throughputs` to calculate m5.

If you want to loop over all the rafts without examining the results for each, use the `.py` files.

### Here is my workflow on updating "m5 by amp" predictions

(If you want to simply reproduce the "m5 by amp" plots, go to step 6 below. If you need to update the input camera data, you need to start from step 1. In either case, please note that our workflow depends on how DM stores the data and tracks the history. If that has changed, the workflow would need to change accordingly.)

1. Camera Project Scientist provides or approves data 
  * What we need for updating "m5 by amp" are the QE curves and read noise. The data usually come as pickle files and fits files, with gain and fullwell measurements as well. The fullwell here is the flat field blooming fullwell.
  * So far the data have been provide via the following Github repos:
    - https://github.com/smr456/raftexplorer
    - https://github.com/smr456/raftResults2020jan
  * The camera team documents their raft acceptance testing status here
    - https://confluence.slac.stanford.edu/display/LSSTCAM/
2. Once we have the data, we run `rexplorer1.ipynb` to examine the data.
  * This notebook is a slight modified version of Steve Ritz's notebook https://github.com/smr456/raftexplorer/blob/master/rexplorer.ipynb
  * The plots it produces are saved to cam_plots/, as a record of how they look like in the original data provided.
  * If the plots look good, we copy the data to /project/shared/ (on https://lsst-lsp-stable.ncsa.illinois.edu)
3. import QE data into ts8 camera
  * make clones of `obs_lsst` and `obs_lsst_data` repos
  * create ticket branches for both repos
  * on command line, set up the LSST software stack and the ticket branches of `obs_lsst` and `obs_lsst_data`.
  * run the command below (with raft ID and pickle file modified according to which raft you are running on)
  ~~~~
  python ./obs_lsst/bin.src/rewrite_ts8_qe_files.py /project/shared/bxin/cam_as_built/R22/RaftRun11351.p --out_root=$OBS_LSST_DATA_DIR/ts8/qe_curve --valid_start 1970-01-01T00:00:00
  ~~~~
  * If updating multiple rafts at the same time, run
  ~~~~
  ls /project/shared/bxin/cam_as_built/*/*p | awk '//{print "python $OBS_LSST_DIR/bin.src/rewrite_ts8_qe_files.py ", $0, "--out_root=$OBS_LSST_DATA_DIR/ts8/qe_curve --valid_start 1970-01-01T00:00:00"}' > batch/injest_ts8.sh
  Source inject_ts8.sh
  ~~~~
4. copy `obs_lsst_data/ts8` QE data into `obs_lsst_data/lsstcam`, and make read noise etc available.
  * run the following on command line (with raft ID and pickle file modified according to which raft you are running on)
  ~~~~
  cd ~/lsst_stack/obs_lsst/policy/
  python update_raft.py /project/shared/bxin/cam_as_built/R11/RaftRun10669.p R11 E2V lsstCam/R11.yaml
  ~~~~
  * If updating multiple rafts at the same time, run
  ~~~~
  ls /project/shared/bxin/cam_as_built/*/*p | awk '//{i=index($1,"Raft")-4;r=substr($1,i,3);print "python update_raft.py",$0,r,"E2V","lsstCam/"r".yaml”}’ > injestRN.sh then modify E2V to ITL by hand.
  source injestRN.sh
  ~~~~
5. Use ticket branches for both `obs_lsst` and `obs_lsst_data`, do `scons lsstcam`, all the data I need for m5 should now be available.
  * If `lsstcam/CALIB/calibRegistry.sqlite3` already exists, `scons` will skip ingestion. That is the way `lsstcam/SConscript` is written right now. It doesn’t check if there are more fits files available. The way to solve that is `scons -c`.
  * `scons ts8` will do the `ts8` ingestion only.
  * Use `Comp_QE_curve.ipynb` to read out QE curves from the butler and compare to the original E2V and ITL design curves. It also plots the measured data points to let you examine how well the interpolation/extrapolation is working. The vendor design curves were made by Steve Ritz using vendor data https://github.com/smr456/QE-working-area/blob/master/QE_v3.ipynb
  * Run `m5_by_amp.ipynb`, which calculates the "m5 by amp" results. It dumps the results, together with all the intermediate results that I thought might be of interest, into csv files.
  * If you are only interested in updating the csv files, you can use `batch/make_csv_outputs.py` to run over all the rafts.
  * use `m5_hist.ipynb` to produce histograms of m5 values

6. visualize the focal plane using the heat map tool by Richard Dubois
  * Richard’s confluence page on his heat map tool. https://confluence.slac.stanford.edu/pages/viewpage.action?pageId=240276456
  * on command line, do
  ~~~~
  source startHeatmap.sh
  ~~~~
  * Have a browser window point to https://lsst-lsp-stable.ncsa.illinois.edu/nb/user/bxin/proxy/5006/serveRenderFP and you should see the Bokeh plots.
