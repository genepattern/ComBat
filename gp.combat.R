setup <-
function(libdir)
{
	source(paste(libdir, "ComBat.R", sep=''))
	source(paste(libdir, "common.R", sep=''))
}

parseCmdLine <- function(...) {
	suppressMessages(.parseCmdLine(...))
}

.parseCmdLine <- function(...)
{
	args <- list(...)
	input.file.name <- ''
	sample.info.file.name <- ''
	covariates <- ''
	filter <- ''
	prior.plots <- ''
	par.prior <- ''
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
   		else if(flag=='-c')
		{
			covariates <- value
   		}
   		else if(flag=='-f')
		{
			filter <- value
   		}
   		else if(flag=='-p')
		{
			prior.plots <- as.logical(value)
   		}
   		else if(flag=='-m')
		{
			par.prior <- as.logical(value)
   		}
		else if(flag=='-o')
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

    gp.combat.R(input.file.name, sample.info.file.name, libdir, output.file.name, prior.plots, par.prior)

}

gp.combat.R <- function(input.file.name, sample.info.file.name, libdir, output.file.name, prior.plots, par.prior)
{
    setup(libdir)

    dataset <- read.dataset(input.file.name)
    data.matrix <- dataset$data

    combat.input.file.name <- "combat_input_file.txt"
    on.exit(unlink(combat.input.file.name))

    #check if it is a gct file
    if(is.null(dataset$calls))
    {
        column.names <- c("ProbeID", colnames(data.matrix))
        cat(column.names, file=combat.input.file.name, sep="\t")
        cat("\n", file = combat.input.file.name, append=TRUE)
        write.table(data.matrix, file=combat.input.file.name, sep="\t", col.names = FALSE, quote=FALSE, append=TRUE)
    }
    else
    {
        dataset$calls
    }

    pdf(file=paste(output.file.name, ".plot", sep=''), height = 12, width = 15)
    combat.result <- ComBat(combat.input.file.name, sample.info.file.name, skip = 1, write = F, prior.plots = prior.plots, par.prior = par.prior)
    dev.off();

    combat.result <- subset(combat.result, select=colnames(data.matrix))
    row.names(combat.result) <- row.names(data.matrix)


    dataset$data <- combat.result
    write.gct(dataset, output.file.name)
}