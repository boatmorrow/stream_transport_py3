StreamTran model for simulating 1d transport in streams subject to atmospheric equilibration, groundwater inflow, first-order reaction.

It is especially suited for transport of environmental tracers such as radon, noble gases, dissolved CFC's, SF6, etc.

Dependencies:

Python package requirements:
numpy
scipy
matplotlib
pandas
lmfit
fipy
tracer_tools [package info.](https://github.com/boatmorrow/tracer_tools_py3)

    -> The package installation should take care of everything except tracer_tools, which is not available on a package manger yet.  Install the tracer_tools package according to the information for that package.

An example and template input deck (moderately well commented) can be found in the examp directory.

INSTALLATION
To install StreamTran:

Navigate to a directory where you would like StreamTran to reside.

Clone the git repository.  For example, using git at the terminal

git clone https://github.com/boatmorrow/stream_transport_py3

This will add the stream_transport_py3 directory.

        cd /path/to/clone/location/stream_transport_py3

Install the package:

        pip install .

StreamTran and tracer_tools written by:

W. Payton Gardner
Dept. of Geosciences
University of Montana
