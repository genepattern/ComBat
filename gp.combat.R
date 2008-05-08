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
   		else if(flag == '-c')
		{
			covariates <- c(3, value)
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

    sample.info <- read.table(sample.info.file.name, header=T, sep='\t',comment.char='')

    if(length(colnames(sample.info)) < 3)
    {
        stop("The sample info file must have the 3 columns: Array, Sample, Batch.")
    }


    if(is.null(filter) || filter == '')
    {
        if(is.null(dataset$calls))
        {
            filter <- F
        }
        else
        {
           filter <- 1
        }
    }
    else if(!is.null(dataset$calls) && (filter <= 0 || filter >= 1))
    {
        stop("Absent calls filter must be greater than 0 and less than 1")
    }


    gct.ext <- regexpr(paste(".gct","$",sep=""), tolower(output.file.name))
    res.ext <- regexpr(paste(".res","$",sep=""), tolower(output.file.name))
	if(gct.ext[[1]] != -1 || res.ext[[1]] != -1)
	{
	    output.file.name <- substr(output.file.name, 0, (nchar(output.file.name)-4))
	}
		
	combat.input.file.name <- paste(output.file.name, ".txt")

    on.exit(unlink(combat.input.file.name))

    column.names <- c("ProbeID", colnames(data.matrix))
    input.table <- data.matrix
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
    {
        if (capabilities("jpeg"))
        {
            jpeg(filename = paste(output.file.name, ".plot.jpeg", sep=''), width = 800, height = 720) 
        }
        else
        {
           pdf(file = paste(output.file.name, ".plot.pdf", sep=''), height = 12, width = 15)
        }
    }

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
        colnames(descriptions.row) <- row.names(data.matrix)
        descriptions.row <- subset(descriptions.row, select = row.names(combat.result))
        dataset$row.descriptions <- descriptions.row[2,]
    }

    if(!is.null(dataset$column.descriptions))
    {
        dataset$column.descriptions <- dataset$column.descriptions[ dataset$column.descriptions!=""]
        descriptions.column <- rbind(colnames(data.matrix), dataset$column.descriptions)
        column.names <- colnames(combat.result)

        colnames(descriptions.column) <- descriptions.column[1,] 
        descriptions.column <- subset(descriptions.column, select=column.names)
        descriptions.column <- descriptions.column[2,]

        dataset$column.descriptions <- descriptions.column        
    }

    if(!is.null(dataset$calls))
    {
        rownames <- row.names(dataset$calls)
        calls <- t(dataset$calls)
        colnames(calls) <- rownames

        calls <- subset(calls, select = row.names(combat.result))
        calls <- t(calls)
        row.names(calls) <- gene.info
        colnames(calls) <- colnames(dataset$calls)


        dataset$calls <- calls
        write.res(dataset, output.file.name)
    }
    else
    {
        write.gct(dataset, output.file.name)
    }
}