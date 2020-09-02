# This function is used to run jobs in parallel on SGE cluster
clsapply.plus <- function(dat.list, func.list, df.list=NULL, tempdir.list="~/project/tempZ", resdir,
		queque="1-day", timespan="540", mem="8G", req.pack='stats',
		tidy=FALSE) {
	nList <- length(df.list) 
	list.names <- names(df.list)
	if (is.null(list.names)) {
		list.names <- paste0('Sim', 1 : nList)
	}
	prefix.list <- list()
	for (i in 1 : nList) {
		cat('Job list ', i, '...\n')
		df <- df.list[[i]]
		tempdir <- tempdir.list[[i]]
		
		
		if (length(dat.list) == 1) {
			dat <- dat.list[[1]]
		} else {
			dat <- dat.list[[i]]
		}
		
		if (length(func.list) == 1) {
			func <- func.list[[1]]
		} else {
			func <- func.list[[i]]
		}
		
		if (!file.exists(tempdir)) {
			dir.create(tempdir)
		}
		setwd(tempdir)
		
		save(dat, func, df, file=("dat.RData"))
		
		nJob <- length(dat)
		# Create R script
		txt <- paste('
						args=(commandArgs(TRUE))
						if(length(args)==0){
						print("No arguments supplied.")
						##supply default values
						} else{
						for(i in 1:length(args)){
						eval(parse(text=args[[i]]))
						}
						}',
				paste(paste0('\nrequire(',req.pack, ')'), collapse="\n"),    
				'\n date()
						load("dat.RData")	
						if (is.list(dat)) {
						item <- dat[[part]]
						} 
						if (is.vector(dat)) {
						item <- dat[part]
						}           
						res0 <- func(item, df)
						save(res0, file=paste(part, ".res", sep=""))
						date()
						')
		writeLines(txt, "script.R")
		
		# Submit job
		cat("Submit jobs ...")
		prefix <- paste(c("J", sample(LETTERS, 4, repl=TRUE)), collapse="")
		prefix.list[[i]] <- prefix
		for (part in 1:nJob) {
			rfile <- "script.R"
			rout <- paste(part, ".Rout", sep="")
			
			sh <- paste(
					"qsub",
					paste0("-N ", prefix, part),
					"-j y",
					"-cwd",
					#			"-t", timespan,
					"-q", queque,
					"-m abe",
					paste0("-l h_vmem=", mem),
					"-V",
					paste("-o ",  part, ".out", sep=""),
					paste("-e ",  part, ".err", sep=""),
					"-b y",
					paste("\'R CMD BATCH --no-save --no-restore \"--args  part=", part, "\" ", 
							rfile, " ", rout, "\'", sep="")
			)
			print(sh)
			system(sh)
		}
	}
	

	cat("Please wait ...\n")
	
	# To see if all the jobs have been finished
	remaining.ids <- 1 : nList
	while (TRUE) {
		Sys.sleep(20)
		if (length(remaining.ids) == 0) {
			break
		}
		
		for (i in remaining.ids) {
			prefix <- prefix.list[[i]]
			tempdir <- tempdir.list[[i]]
			setwd(tempdir)
			
			output1 <- system("qstat ", intern=TRUE)
			output2 <- system("qstat ", intern=TRUE)
			missingfile <- NULL
			if (length(grep(prefix, output1)) == 0 & length(grep(prefix, output2)) == 0) {	
				cat("Job ", i, ": all runs finished and begin to combine\n")
				res <- list()
				for (part in 1:nJob) {
					resfile <- paste(part, ".res", sep="")
					if (file.exists(resfile)) {
						load(resfile)
						res[[paste(part)]] <-  res0
						rm(res0)
						cat(".")
					} else {
						# cat("\n", resfile, " not found!\n")
						missingfile <- c(missingfile, resfile)
					}
				}
				cat("\n")
				cat("Missing files:\n")
				if (length(missingfile) == 0) {
					cat("No jobs dropped!\n")
				} else{
					cat(paste(missingfile, collapse="\n"))
				} 
				# Clean
				if (tidy == TRUE) {
					cat("\nBegin to tidy up ...")
					system("rm *")
				}
				
				setwd(resdir)
				save(res, file=file.path(resdir, paste0(list.names[i], "_res.Rdata")))
				rm(res)
				remaining.ids <- setdiff(remaining.ids, i)
				cat("\nDone!\n")
			}
			
		}
		gc()
	}
	
}

