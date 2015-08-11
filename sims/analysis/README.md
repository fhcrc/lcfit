

## `lcfit` failure analysis

aggfit.Rmd - Display of aggregated failure modes of lcfit.

shiny/ - interactive display of lcfit failure modes.
Allows selection and browsing of failures.
User can adjust initial model parmeters and immediately see the effect on the lcfit curve.

## Running the shiny app

	cd lcfit/sims/analysis/shiny
	Rscript -e ' library(methods);  shiny::runApp()'

This should launch a browser on your desktop with the shiny app.
If you do not wish to launch a browser, you can use:

	Rscript -e ' library(methods);  shiny::runApp(".", launch.browser=FALSE, host="140.107.116.67")'

This will start the shiny server on a randomly assigned port,
listening on ip address `140.107.116.67` (stoat.fhcrc.org).
Once it has launched the server, Shiny will print the URL that you
should open in a browser, e.g.

	Listening on http://140.107.116.67:4040

