pp.kSamples <- function(x){
if(class(x) != "kSamples") stop('input not of class "kSamples"')
if(x$test.name == "Anderson-Darling"){
	version <- "0"
	while(substr(version,1,1) != "1" & substr(version,1,1) != "2"){
		version <- readline("choose version of test statistic: 1 or 2\n")
	} 
	if(names(x)[[2]]=="k") {
		dist <- x[[8+as.numeric(version)]]}else{
		dist <- x[[12+as.numeric(version)]]
	}
	if(!is.null(dist)){
		xy <- table(dist)
		N <- length(dist)
		if(names(x)[[2]]=="k") {
			xx <- (as.numeric(names(xy))-(x$k-1))/x$sig
			px <- ad.pval(xx,x$k-1,version)}else{
			xx <- (as.numeric(names(xy))-(x$mu.c))/x$sig.c
			px <- ad.pval(xx,x$mu.c,version)
		}
		px.exact <- rev(cumsum(rev(as.numeric(xy)))/N)
		m <- min(px,px.exact)
		M <- max(px,px.exact)
		par(pty="s")
		if(names(x)[[2]]=="k") titlex <- paste(x$test.name,"k-sample test, version", version)
		if(names(x)[[2]]=="M") titlex <- paste("combined", x$test.name,"k-sample tests, version", version)
		plot(px.exact,px,log="xy",xlim=c(m,M), ylim=c(m,M), main= titlex,
		xlab=paste(x$method,"significance probability"),ylab="asymptotic significance probability",pch=16,cex=.5)
		text(.055,m,".05",adj=0)
		abline(0,1,lty=2)
		abline(h=.05,lty=2)
		abline(v=.05,lty=2)
		text(.012,.5,".01",adj=0)
		abline(h=.01,lty=2)
		abline(v=.01,lty=2)
		if(x$method=="simulated"){ 
			text(m,M,paste("Nsim =",x$Nsim),adj=0)}else{
		        text(m,M,paste("Ncomb =",length(dist)),adj=0)
		}
	}else{
		cat("no distribution data\n")
	}
}
if(x$test.name == "Kruskal-Wallis" |
   	x$test.name == "van der Waerden scores" |
   	x$test.name == "normal scores"){
	dist <- x$null.dist

	if(!is.null(dist)){
		xy <- table(dist)
		N <- length(dist)
		xx <- as.numeric(names(xy))
		if(names(x)[[2]] == "k"){
			px <- 1-pchisq(xx,x$k - 1)}else{
			df <- sum(sapply(x$n.samples,length))-length(x$n.samples)
			px <- 1-pchisq(xx,df)
		}
		px.exact <- rev(cumsum(rev(as.numeric(xy)))/N)
		m <- min(px,px.exact)
		M <- max(px,px.exact)
		par(pty="s")
		if(names(x)[[2]]=="k") titlex <- paste(x$test.name,"k-sample test")
		if(names(x)[[2]]=="M") titlex <- paste("combined", x$test.name,"k-sample tests")
		plot(px.exact,px,log="xy",xlim=c(m,M), ylim=c(m,M), main= titlex,
			xlab=paste(x$method,"significance probability"),
			ylab="asymptotic significance probability",pch=16,cex=.5)
		text(.055,m,".05",adj=0)
		abline(0,1,lty=2)
		abline(h=.05,lty=2)
		abline(v=.05,lty=2)
		text(.012,.5,".01",adj=0)
		abline(h=.01,lty=2)
		abline(v=.01,lty=2)
		if(x$method=="simulated"){ 
			text(m,M,paste("Nsim =",x$Nsim),adj=0)}else{
	    	    text(m,M,paste("Ncomb =",length(dist)),adj=0)
		}
	}else{
		cat("no distribution data\n")
	}
}

if(x$test.name == "2 x t Contingency Table" | 
   	x$test.name == "Combined 2 x t Contingency Tables"){
	dist <- x$null.dist
	if(!is.null(dist)){
		xx <- dist[,1]
	        px.exact <- rev(cumsum(rev(dist[,2])))
		N <- length(dist[,1])
		if(names(x)[[3]] == "M"){
			M <- x$M
	    		tvec <- x$t
			df <- sum(tvec)-M 
		}else{
			df <- x$t -1
		}
		px <- pchisq(xx,df,lower.tail=FALSE)
		m <- min(px,px.exact)
		M <- max(px,px.exact)
		plot(px.exact,px,log="xy",xlim=c(m,M), ylim=c(m,M), main= x$test.name,
		xlab=paste(x$method,"significance probability"),ylab="asymptotic significance probability",pch=16,cex=.5)
		text(.055,m,".05",adj=0)
		abline(0,1,lty=2)
		abline(h=.05,lty=2)
		abline(v=.05,lty=2)
		text(.008,.5,".01",adj=1)
		abline(h=.01,lty=2)
		abline(v=.01,lty=2)
		if(x$method=="simulated"){ 
			text(m,M,paste("Nsim =",x$Nsim),adj=0)}else{
		        text(m,M,paste("Nsupp =",length(dist[,1])),adj=0)
		}
	}else{
		cat("no distribution data\n")
	}
}
}
