
#### Spectral processing ####
## --------------------------
setGeneric("smoothSignal", 
	function(object, ...) standardGeneric("smoothSignal"))
setGeneric("reduceBaseline", 
	function(object, ...) standardGeneric("reduceBaseline"))

#### Peak processing ####
## -----------------------
setGeneric("peakPick", 
	function(object, ...) standardGeneric("peakPick"))
setGeneric("peakGroup", 
	function(object, ...) standardGeneric("peakGroup"))
setGeneric("peakAlign", 
	function(object, ...) standardGeneric("peakAlign"))
setGeneric("peakMerge", 
	function(object, ...) standardGeneric("peakMerge"))

