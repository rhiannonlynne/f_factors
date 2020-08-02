#source /opt/lsst/software/stack/loadLSST.bash; setup lsst_distrib
export PYTHONPATH=$HOME/notebooks/f_factors/source/m5_by_amp/:$HOME/notebooks/datacat-utilities/python:$HOME/notebooks/eTraveler-clientAPI/python:$HOME/notebooks/EO-plotsDisplays/python
bokeh serve $HOME/notebooks/EO-plotsDisplays/python/serveRenderFP.py --allow-websocket-origin=lsst-lsp-stable.ncsa.illinois.edu --args -e SR_emulation_config.txt --hook m5_hook -t "User m5 u"
#bokeh serve $HOME/notebooks/EO-plotsDisplays/python/serveRenderFP.py --allow-websocket-origin=tucson-teststand.lsst.codes --args -e SR_emulation_config.txt --hook m5_hook -t "User m5 u"
#bokeh serve $HOME/notebooks/EO-plotsDisplays/python/serveRenderFP.py --allow-websocket-origin=tucson-teststand.lsst.codes --args -e SR_emulation_config.txt --hook t_user_hook -t User
