compute.q.values <- function(p.value.file='p.values.txt', indices.file='indices.txt', previous.task.values, output.file.name='r_output.txt') {
	#on.exit(unlink(c(p.value.file, indices.file, previous.task.values)))
	library(qvalue)
	
	p.values <-scan(p.value.file)
	result <- qvalue(p.values)
	qvalues <- result$qvalues
	
	previous.values <- read.table(previous.task.values, sep="\t", row.names=1)
	
	m <- matrix(nrow=NROW(previous.values), ncol=NCOL(previous.values)+1)
	
	row.names(m) <- row.names(previous.values)
	for(i in 1:NCOL(previous.values)) {
		m[, i] <- previous.values[, i]
	}

	indices <- scan(indices.file)
	last.column <- NCOL(previous.values)+1
	for(i in 1:NROW(previous.values)) {
		m[i, last.column] <- qvalues[indices[[i]]+1]
	}

	write.table(m, output.file.name, sep="\t", col.names=FALSE, quote=FALSE)
	
	return(output.file.name)
}