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
			filter <- as.numeric(value)
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

    gp.combat.R(input.file.name = input.file.name, sample.info.file.name =sample.info.file.name, libdir = libdir, output.file.name = output.file.name, prior.plots = prior.plots , par.prior = par.prior, filter = filter)

}

gp.combat.R <- function(input.file.name, sample.info.file.name, libdir, output.file.name, prior.plots, par.prior, filter)
{
    setup(libdir)

    dataset <- read.dataset(input.file.name)
    data.matrix <- dataset$data

    combat.input.file.name <- "combat_input_file.txt"
    on.exit(unlink(combat.input.file.name))

    column.names <- c("ProbeID", colnames(data.matrix))
    #check if it is a res file
    if(!is.null(dataset$calls))
    {
        column.names <- NULL
        input.table <- NULL
        for(i in 1:ncol(data.matrix))
        {
            col.name <- colnames(data.matrix)[i]
            calls.column.name <- paste(col.name, "call", sep=' ')
            column.names <- c(column.names, col.name, calls.column.name)
		    input.table <- cbind(input.table, data.matrix[,i])
		    input.table <- cbind(input.table, as.character(dataset$calls[,i]))
	    }

        row.names(input.table) <- row.names(data.matrix)
        column.names <- c("ProbeID", column.names)
    }

    cat(column.names, file=combat.input.file.name, sep="\t")
    cat("\n", file = combat.input.file.name, append=TRUE)
    write.table(input.table, file=combat.input.file.name, sep="\t", col.names = FALSE, quote=FALSE, append=TRUE)

    if(prior.plots)
        pdf(file=paste(output.file.name, ".plot", sep=''), height = 12, width = 15)

    combat.result <- ComBat(combat.input.file.name, sample.info.file.name, filter = filter, skip = 1, write = F, prior.plots = prior.plots, par.prior = par.prior)

    if(prior.plots)
        dev.off();

    gene.info <- combat.result[,1]
    combat.result <- subset(combat.result, select = colnames(data.matrix))

    row.names(combat.result) <- gene.info
    dataset$data <- combat.result

    if(!is.null(dataset$row.descriptions))
    {
        descriptions.row <- rbind(row.names(data.matrix), dataset$row.descriptions)
        descriptions.row <- descriptions.row[gene.info]
        dataset$row.descriptions <- descriptions.row
    }


    if(!is.null(dataset$column.descriptions))
    {
        descriptions.column <- rbind(colnames(data.matrix), dataset$column.descriptions)
        descriptions.column <- descriptions.column[colnames(combat.result)]
        dataset$column.descriptions <- descriptions.column
    }


    if(!is.null(dataset$calls))
    {
        write.res(dataset, output.file.name)
    }
    else
    {
        write.gct(dataset, output.file.name)
    }
}