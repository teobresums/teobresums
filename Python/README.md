### Python scripts


## PlotWave.py

To plot anything relevant related to the waveform, use PlotWave.py

# Usage documentation

Run
	python PlotWave.py 

without any arguments to display the usage documentation:

	usage: PlotWave.py [-h] -i FILE [FILE ...] [-o [OUTPUTDIR]] [-m MODE]
		           [-s SHOW] [-a YRANGE [YRANGE ...]]

	Simple plot of waveform files from TEOBResumS

	optional arguments:
	  -h, --help            show this help message and exit
	  -i FILE [FILE ...]    Input files
	  -o [OUTPUTDIR]        Output directory
	  -m MODE               Type of plot. A string containing one or more of the
		                chars 'apori' for amplitude, phase, frequency, real,
		                and imaginary
	  -s SHOW               Show plot
	  -a YRANGE [YRANGE ...]
		                Ranges of yaxis


# Example run

	python PlotWave.py -i path_to_file/hlm_insplunge_l2_m2.txt -m a -a {-0.5,0.5}

plots the amplitude of the (2,2)-mode up to merger within the vertical range of [-0.5, 0.5]


# Hints/Suggestions

- The argument -m a plots the amplitude vs. time, likewise -m p plots the phase, etc. Many arguments can be plotted together, e.g., -m ar. Default is set to "aor".

- The argument -s is set to "yes" by default so can be omitted.

- The argument -a {x, y} sets the vertical-axis range for the plot. The default is set to {-0,75, 0.75} which suffices for a, r, i, and even o. However the phase "p" will range from 0 to hundreds, even thousands.



