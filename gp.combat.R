setup <-
function(libdir)
{
	source(paste(libdir, "ComBat.R", sep=''))
}

parseCmdLine <- function(...) {
	suppressMessages(.parseCmdLine(...))
}

.parseCmdLine <- function(...)
{
	args <- list(...)
	input.file.name <- ''
	sample.info.file.name <- ''
	output.file.name <- ''
	libdir <- ''

	for(i in 1:length(args))
	{
	    flag <- substring(args[[i]], 0, 2)
		value <- substring(args[[i]], 3, nchar(args[[i]]))
		if(flag=='-i')
		{
			input.file.name <- value
		}
		else if(flag=='-s')
		{
			sample.info.file.name <- value
   		}
		else if(flag=='-0')
		{
			output.file.name <- value
   		}
		else if(flag=='-l')
		{
			libdir <- value
   		}

   		else
   		{
			stop(paste("unknown flag ", flag, " value ", value, sep=""), .call=FALSE)
		}
    }

    gp.combat.R(input.file.name, libdir)

}

gp.combat.R <- function(input.file.name, libdir)
{
    setup(libdir)

    dataset <- read.dataset(dataset.file)
}