compute.q.values <- function(previous.task.values, output.file.name='r_output.txt') {
	library(qvalue)
	values <- read.table(previous.task.values, sep="\t", row.names=1)
	result <- qvalue(values[,3])
	qvalues <- result$qvalues
	values <- cbind(values, qvalues)
	write.table(values, output.file.name, sep="\t", col.names=FALSE, quote=FALSE)
	return(output.file.name)
}