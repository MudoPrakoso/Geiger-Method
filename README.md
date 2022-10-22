This program contains a procedure to determine hypocentre and origin time of an earthquake by inversion method. The inversion method used is Geiger Method (Least Square Error Optimization) (Wysession and Seth, 2007). A modification, damping factor, is added to original solution of least square error. This damping factor is used to prevent the solution being divergent.

The arrival time data of some earthquake around East Java, Indonesia, is provided by Meteorological, Climatological, and Geophysical Agency. These data contain the arrival time (GMT+8) of seismic wave at a specified station. The data is named after date of event.

The station data (sta_loc) is compilation of coordinate of seismic station within and around Indonesia.

The velocity model used is IASP91, that is sufficient to determine a local source which doesn’t consist a deep phase of seismic wave.
The damping factor should be large enough at first iteration and decreases to small value when the solution is convergent (Transtum and Sethna, 2012). Some report (Lienert et al, 1986; Transtrum and Sethna, 2012) suggested some adjusting factor to damping factor. Yet, based on empirical trial and error, 1 – 100 initial damping factor with adjusting factor 2 – 10 would give minimum root mean square error.

This program will give you a report about the coordinate of a hypocentre (latitude(in degree), longitude(in degree), altitude(in km), as well as the uncertainty(in km)) and an origin time (in GMT+8 as well as the uncertainty(in second)) of an earthquake. Also, a graph depicting how the root mean square changes through iteration.

*) This program doesn’t include data weighting, that is very important if the source far enough from seismic station (laterally or vertically). Seismic wave may travel through some heterogeneity that doesn’t represent our velocity model. So it will give infinite loop or no logic solution (the altitude above earth surface). Try with another combination of damping factor and adjusting factor.

References:

Lienert, B., E. Berg and N. Frazer. 1986. HYPOCENTER: An Earthquake Location Method Using Centered, Scaled and Adaptively Damped Least Squares. Bulletin of the Seismological Society of America. 76: 771 – 783.

Stein, S. and M. Wysession. 2003. An Introduction to Seismology, Earthquake, and Earth Structure. Blackwell Publishing. Massachusestts.

Transtrum, M. and J. Sethna. 2012. Improvement to the Levenberg-Marquartdt Algorithm for Nonlinear Least-Squares Minimization.
