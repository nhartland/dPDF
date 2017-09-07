#!/bin/bash
python3 ./histogram.py $1/
python3 ./plot_data.py $1/
python3 ./plot_f2r.py  $1/
python3 ./plot_pdf.py  $1/
