options(encoding="utf-8")
options("scipen"=999,"digits"=4)

REx_ANA_PLOT <- function(w=500,h=500) {
  ## save plot as temp file
  REx.plot.tempfn <<- paste(tempdir(), "\\REx_temp.png", sep="")
  png(filename=REx.plot.tempfn, width=w, height=h)
}

REx_ANA_PLOT_OFF <- function(caption) {
  dev.off()
  load.pkg("markdown")
  ## read temp file as a binary string
  img <- paste(markdown:::.b64EncodeFile(REx.plot.tempfn))
  R2HTML::HTML(paste("<p align=left><img src='", img, "' /><br /><font class=caption>", caption, "</font></p>", sep=""),file=stdout())
}

local.install.packages <- function(pkgs, repos, type='win.binary') {
	rversion = paste(R.version$major, strsplit(R.version$minor, '.', fixed=TRUE)[[1]][1], sep='.')
	local.contriburl = contrib.url(paste('file:///', repos, sep=''), type=type)
	local.pkgs = available.packages(local.contriburl)

	dir = '/src/contrib/'
	extension = '.tar.gz'
	if(type == 'mac.binary.leopard' | type == 'mac.binary') {
		dir = paste('/bin/macosx/leopard/contrib/', rversion, '/', sep='')
		extension = '.tgz'
	} else if(type == 'win.binary') {
		dir = paste('/bin/windows/contrib/', rversion, '/', sep='')
		extension = '.zip'
	}
	installed = installed.packages()
	#package.dependencies(local.pkgs[local.pkgs[,'Package'] %in% pkgs,])
	
	d = c("Depends", "Imports", "LinkingTo", "Suggests", "Enhances")
	toinstall = utils:::getDependencies(pkgs, d, local.pkgs)
	ininstall2 = utils:::getDependencies(toinstall, d, local.pkgs)
	while(length(toinstall) != length(ininstall2)) {
		toinstall = ininstall2
		ininstall2 = utils:::getDependencies(toinstall, d, local.pkgs)		
	}
	
	for(i in 1:length(toinstall)) {
		toinstall[i] = paste(repos, dir,
				 toinstall[i], '_', 
				 local.pkgs[local.pkgs[,'Package'] == toinstall[i],'Version'], 
				 extension,
				 sep='')
	}
	print('Install the following packages:')
	print(toinstall)
	install.packages(toinstall, repos=NULL, contrib.url=local.contriburl, available=local.pkgs, type=type)
}

inst.pkg <- function(pkg, repos="http://healthstat.snu.ac.kr/CRAN") {
	if (!suppressWarnings(require(pkg, quietly=TRUE, character.only=TRUE))) {
		ret <- try(local.install.packages(pkg, repos=repos), T)
		if (class(ret) == 'try-error') 2
		else 1
	}
	0
}

load.pkg <- function(pkgs) {
  for (i in pkgs) {
    ret <- inst.pkg(i)
	if (ret != 0) {
	  if (ret == 2) stop(paste0("패키지 [", i, "] 설치에 실패했습니다. REx 개발팀에게 문의해 주세요!"))
	  else if (ret == 1) {
	    if (class(require(i, quietly=TRUE, character.only=TRUE)) == 'try-error')
		  stop(paste0("패키지 [", i, "] 가 설치되었지만 로드에 실패했습니다. REx 개발팀에게 문의해 주세요!"))
	  }
	}
  }
}

# Plot function for distribution
plotDistr <- function(x, p, discrete=FALSE, cdf=FALSE, regions = NULL, 
                          col = "gray", legend = TRUE, legend.pos = "topright",
                          main, xlab, ylab, ...){
  library(ggplot2)
  dat <- data.frame(x,p)
  
  theme_opt <- theme_bw() + 
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) + 
    theme(plot.title = element_text(hjust = 0.5))
  
  ggp <- ggplot(aes(x=x, y=p), data=dat)
  if(discrete) ggp <- ggp + geom_point(size=2) + geom_segment(aes(xend=x, yend=0),lineend="butt")
  else ggp <- ggp + geom_line()
  if(cdf) ggp <- ggp + geom_hline(yintercept = 0:1, color="grey")
  else ggp <- ggp + geom_hline(yintercept = 0, color="grey")
  
  ggp <- ggp + theme_opt + labs(title=main, x=xlab, y=ylab)
  return(ggp)
}

# Index plot

REx_indexplot <- function(varname, dataset, type="spike", id=TRUE, idnum=2, 
                          title=NULL, xlab="Observation Index", ylab=varname,
                          color="black"){
  load.pkg("ggplot2")
  # browser()
  suppressWarnings(datavar <- as.numeric(as.character(dataset[,varname])))
  if(sum(is.na(datavar))==length(datavar)) stop("There is no numeric values.")
  dataset <- data.frame(datavar)
  colnames(dataset) <- varname
  
  idn <- sort(order(datavar, decreasing = T)[1:idnum])
  datalab <- data.frame(x=idn, y=datavar[idn])
  
  if(type=="spike"){
    spike_index1 <- ggplot(aes(x=1:nrow(dataset), eval(parse(text=paste0("y=",varname)))), data=dataset) + 
      geom_segment(aes(xend=1:nrow(dataset), yend=0),size=1,lineend="butt",color=color, na.rm=T) + 
      coord_cartesian(ylim = range(datavar, na.rm=T)) + 
      labs(list(title=title, x=xlab, y=ylab)) +
      theme_bw()
    if(id==FALSE) return(spike_index1 + theme(plot.title = element_text(hjust = 0.5)))
    if(id==TRUE) {
      spike_index2 <- spike_index1 +
        geom_text(aes(x=x, y=y, label=x), data=datalab,
                  hjust = 1.2) + theme(plot.title = element_text(hjust = 0.5))
      return(spike_index2)
    }
  }
  if(type=="dot"){
    dot_index1 <- ggplot(aes(x=1:nrow(dataset), eval(parse(text=paste0("y=",varname)))), data=dataset) + 
      geom_point(color=color, na.rm=T) +
      labs(list(title=title, x=xlab, y=ylab)) +
      theme_bw()
    if(id==FALSE) return(dot_index1 + theme(plot.title = element_text(hjust = 0.5)))
    if(id==TRUE) {
      dot_index2 <- dot_index1 +
        geom_text(aes(x=x, y=y, label=x), data=datalab,
                  hjust = 1.2) + theme(plot.title = element_text(hjust = 0.5))
      return(dot_index2)
    }
  }
}


##### Dot plot

REx_dotplot <- function(varname, dataset, by=NULL, bin=NULL, type="stack", lgd.pos="right", 
                        title=NULL, xlab=varname,
                        color="black", color.group=NULL){
  load.pkg("ggplot2")
  # browser()
  suppressWarnings(datavar <- as.numeric(as.character(dataset[,varname])))
  if(sum(is.na(datavar))==length(datavar)) stop("There is no numeric values.")
  dataset2 <- data.frame(datavar)
  colnames(dataset2) <- "x"

  if(is.null(by)){
    res <- ggplot(aes(x=x), data=dataset2)
    if(is.null(bin)) res <- res + geom_dotplot(method="histodot", color=color, fill=color, na.rm=T)
    else res <- res + geom_dotplot(method="histodot", binwidth=diff(range(datavar, na.rm=T))/bin, dotsize=bin/30, color=color, fill=color, na.rm=T)
    res <- res + 
      scale_y_continuous(NULL, breaks = NULL) +    
      labs(list(title=title, x=xlab)) +
      theme_bw() +
      theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
      theme(plot.title = element_text(hjust = 0.5))
    return(res)
  } else {
    dataset2$group <- factor(dataset[,by])
    if(sum(is.na(dataset2$group))>0) dataset2 <- dataset2[-which(is.na(dataset2$group)),]
    if(type=="stack"){
      if(is.null(bin)){
        group_dot1 <- ggplot(aes(x=x, fill=group, color=group), data=dataset2) + 
          geom_dotplot(method="histodot", stackgroups = T, na.rm=T) + 
          scale_y_continuous(NULL, breaks = NULL) +
          labs(list(title=title, x=xlab)) +
          guides(color = guide_legend(by), fill = guide_legend(by)) +
          theme_bw() + 
          theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
          theme(plot.title = element_text(hjust = 0.5), legend.position = lgd.pos)
        if(is.null(color.group)) return(group_dot1)
        else {
          group_dot2 <- group_dot1 + 
            scale_fill_manual(values=color.group) + 
            scale_colour_manual(values=color.group)
          return(group_dot2)
        }
      } else {
        group_dot3 <- ggplot(aes(x=x, fill=group, color=group), data=dataset2) + 
          geom_dotplot(method="histodot", stackgroups = T, binwidth=diff(range(datavar, na.rm=T))/bin, na.rm=T) + 
          scale_y_continuous(NULL, breaks = NULL) +
          labs(list(title=title, x=xlab)) +
          guides(color = guide_legend(by), fill = guide_legend(by)) +
          theme_bw() + 
          theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
          theme(plot.title = element_text(hjust = 0.5), legend.position = lgd.pos)
        if(is.null(color.group)) return(group_dot3)
        else {
          group_dot4 <- group_dot3 + 
            scale_fill_manual(values=color.group) + 
            scale_colour_manual(values=color.group)
          return(group_dot4)
        }
      }
    }
    
    if(type=="facet"){
      if(is.null(color.group)){
        if(is.null(bin)) bin <- 30
        group_dot5 <- ggplot(aes(x=x), data=dataset2) + 
          geom_dotplot(method="histodot", binwidth=diff(range(datavar, na.rm=T))/bin, na.rm=T) + 
          scale_y_continuous(NULL, breaks = NULL) +
          facet_wrap(~group, ncol=1, scales="free", strip.position = "top") +
          coord_cartesian(xlim = range(datavar, na.rm=T)) +
          labs(list(title=title, x=xlab)) +
          theme_bw() + 
          theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
          theme(strip.background = element_blank(), strip.placement = "outside") +
          theme(plot.title = element_text(hjust = 0.5))
        return(group_dot5)
      } else {
        if(is.null(bin)) bin <- 30
        group_dot6 <- ggplot(aes(x=x, fill=group, color=group), data=dataset2) + 
          geom_dotplot(method="histodot", binwidth=diff(range(datavar, na.rm=T))/bin, na.rm=T) + 
          scale_y_continuous(NULL, breaks = NULL) +
          scale_fill_manual(values=color.group) + 
          scale_colour_manual(values=color.group) +
          facet_wrap(~group, ncol=1, scales="free", strip.position = "top") +
          coord_cartesian(xlim = range(datavar, na.rm=T)) +
          labs(list(title=title, x=xlab)) +
          theme_bw() + 
          theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
          theme(strip.background = element_blank(), strip.placement = "outside") +
          theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
        return(group_dot6)
      }
    }
  }
}

##### Histogram

REx_histogram <- function(varname, dataset, by=NULL, bin=30, scale="freq", type="stack", lgd.pos="right", 
                          title=NULL, xlab=varname, ylab=NULL,
                          color="grey", color.group=NULL){
  load.pkg("ggplot2")
  suppressWarnings(datavar <- as.numeric(as.character(dataset[,varname])))
  if(sum(is.na(datavar))==length(datavar)) stop("There is no numeric values.")
  dataset2 <- data.frame(datavar)
  colnames(dataset2) <- "x"
  
  theme_opt <- theme_bw() + 
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    theme(strip.background = element_blank(), strip.placement = "outside") +
    theme(plot.title = element_text(hjust = 0.5))
  
  if(is.null(by)){
    if(scale=="freq"){
      if(is.null(ylab)) ylab <- "Frequency"
      histEXR1 <- ggplot(aes(x), data=dataset2) + 
        geom_histogram(bins=bin, color="black", fill=color, na.rm=T) + 
        labs(list(title=title, x=xlab, y=ylab))
      return(histEXR1 + theme_opt)
    }
    if(scale=="percent"){
      if(is.null(ylab)) ylab <- "Percent"
      histEXR2 <- ggplot(aes(x), data=dataset2) + 
        geom_histogram(aes(y = (..count..)/sum(..count..)), bins=bin, color="black", fill=color, na.rm=T) + 
        scale_y_continuous(labels = scales::percent) +
        labs(list(title=title, x=xlab, y=ylab))
      return(histEXR2 + theme_opt)
    }
  } else {
    dataset2$group <- factor(dataset[,by])
    if(sum(is.na(dataset2$group))>0) dataset2 <- dataset2[-which(is.na(dataset2$group)),]
    if(type=="stack"){
      theme_opt1 <- theme_opt + theme(legend.position = lgd.pos)
      if(scale=="freq"){
        if(is.null(ylab)) ylab <- "Frequency"
        histEXR3 <- ggplot(aes(x, fill=group), data=dataset2) + 
          geom_histogram(bins=bin, col="black", na.rm=T) + 
          guides(color = guide_legend(by), fill = guide_legend(by)) +
          labs(list(title=title, x=xlab, y=ylab))
        if(is.null(color.group)) return(histEXR3 + theme_opt1)
        else {
          histEXR4 <- histEXR3 + 
            scale_fill_manual(values=color.group) + 
            scale_colour_manual(values=color.group)
          return(histEXR4 + theme_opt1)
        }
      }
      if(scale=="percent"){
        if(is.null(ylab)) ylab <- "Percent"
        histEXR5 <- ggplot(aes(x, fill=group), data=dataset2) + 
          geom_histogram(aes(y = (..count..)/sum(..count..)), bins=bin, col="black", na.rm=T) + 
          scale_y_continuous(labels = scales::percent) +
          guides(color = guide_legend(by), fill = guide_legend(by)) +
          labs(list(title=title, x=xlab, y=ylab))
        if(is.null(color.group)) return(histEXR5 + theme_opt1)
        else {
          histEXR6 <- histEXR5 + 
            scale_fill_manual(values=color.group) + 
            scale_colour_manual(values=color.group)
          return(histEXR6 + theme_opt1)
        }
      }
    }
    if(type=="par"){
      theme_opt1 <- theme_opt + theme(legend.position = lgd.pos)
      if(scale=="freq"){
        if(is.null(ylab)) ylab <- "Frequency"
        histEXR3 <- ggplot(aes(x, fill=group), data=dataset2) + 
          geom_histogram(bins=bin, col="black", na.rm=T, position=position_dodge()) + 
          guides(color = guide_legend(by), fill = guide_legend(by)) +
          labs(list(title=title, x=xlab, y=ylab))
        if(is.null(color.group)) return(histEXR3 + theme_opt1)
        else {
          histEXR4 <- histEXR3 + 
            scale_fill_manual(values=color.group) + 
            scale_colour_manual(values=color.group)
          return(histEXR4 + theme_opt1)
        }
      }
      if(scale=="percent"){
        if(is.null(ylab)) ylab <- "Percent"
        histEXR5 <- ggplot(aes(x, fill=group), data=dataset2) + 
          geom_histogram(aes(y = (..count..)/sum(..count..)), bins=bin, col="black", na.rm=T, position=position_dodge()) + 
          scale_y_continuous(labels = scales::percent) +
          guides(color = guide_legend(by), fill = guide_legend(by)) +
          labs(list(title=title, x=xlab, y=ylab))
        if(is.null(color.group)) return(histEXR5 + theme_opt1)
        else {
          histEXR6 <- histEXR5 + 
            scale_fill_manual(values=color.group) + 
            scale_colour_manual(values=color.group)
          return(histEXR6 + theme_opt1)
        }
      }
    }
    if(type=="facet"){
      theme_opt1 <- theme_opt + theme(legend.position = "none")
      binwid <- diff(range(datavar, na.rm=T))/bin
      if(is.null(color.group)){
        if(scale=="freq"){
          if(is.null(ylab)) ylab <- "Frequency"
          histEXR7 <- ggplot(aes(x), data=dataset2) + 
            geom_histogram(binwidth=binwid, col="black", na.rm=T) + 
            facet_wrap(~group, ncol=1, scales="free_x") +
            coord_cartesian(xlim = range(datavar, na.rm=T)) +
            labs(list(title=title, x=xlab, y=ylab))
          return(histEXR7 + theme_opt1)
        }
        if(scale=="percent"){
          if(is.null(ylab)) ylab <- "Percent"
          histEXR8 <- ggplot(aes(x), data=dataset2) + 
            geom_histogram(aes(y = (..count..)/sum(..count..)), binwidth=binwid, col="black", na.rm=T) + 
            scale_y_continuous(labels = scales::percent) +
            facet_wrap(~group, ncol=1, scales="free_x") +
            coord_cartesian(xlim = range(datavar, na.rm=T)) +
            labs(list(title=title, x=xlab, y=ylab))
          return(histEXR8 + theme_opt1)
        }
      } else {
        if(scale=="freq"){
          if(is.null(ylab)) ylab <- "Frequency"
          histEXR9 <- ggplot(aes(x, fill=group), data=dataset2) + 
            geom_histogram(binwidth=binwid, col="black", na.rm=T) + 
            scale_fill_manual(values=color.group) + 
            facet_wrap(~group, ncol=1, scales="free_x") +
            coord_cartesian(xlim = range(datavar, na.rm=T)) +
            labs(list(title=title, x=xlab, y=ylab))
          return(histEXR9 + theme_opt1)
        }
        if(scale=="percent"){
          if(is.null(ylab)) ylab <- "Percent"
          histEXR10 <- ggplot(aes(x, fill=group), data=dataset2) + 
            geom_histogram(aes(y = (..count..)/sum(..count..)), binwidth=binwid, col="black", na.rm=T) + 
            scale_y_continuous(labels = scales::percent) +
            scale_fill_manual(values=color.group) + 
            facet_wrap(~group, ncol=1, scales="free_x") +
            coord_cartesian(xlim = range(datavar, na.rm=T)) +
            labs(list(title=title, x=xlab, y=ylab))
          return(histEXR10 + theme_opt1)
        }
      }
    }
  }
}

##### Density

REx_densityplot <- function(varname, dataset, by=NULL, kernel="gaussian", binwid="nrd0", adj=1, densthres=NULL,
                            scale="density", type="stack", lgd.pos="right", 
                            title=NULL, xlab=varname, ylab="Density",
                            color="black", color.group=NULL, alpha=0){
  load.pkg("ggplot2")
  # browser()
  suppressWarnings(datavar <- as.numeric(as.character(dataset[,varname])))
  if(sum(is.na(datavar))==length(datavar)) stop("There is no numeric values.")
  dataset2 <- data.frame(datavar)
  colnames(dataset2) <- varname
  if(!is.null(by)){
    dataset2$group <- factor(dataset[,by])
    if(sum(is.na(dataset2$group))>0) dataset2 <- dataset2[-which(is.na(dataset2$group)),]
    colnames(dataset2)[2] <- by
  }

  if(scale=="density") count <- ""
  if(scale=="count"){
    count <- "..count.., "
    ylab <- "Count"
  }
   
  theme_opt <- theme_bw() + 
    theme(panel.grid.minor=element_blank()) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = lgd.pos)
  densitybase <- ggplot(eval(parse(text=paste0("aes(",varname,", ",count,"fill=",by,", color=",by,")"))), data=dataset2)
  
  if(is.null(by)){
    density1 <- densitybase + 
      stat_density(kernel=kernel, bw=binwid, adjust=adj, position=type, alpha=alpha, color=color, fill=color, na.rm=T) +
      labs(list(title=title, x=xlab, y=ylab)) + theme_opt
  } else {
    if(is.null(color.group)){
      density1 <- densitybase +
        stat_density(kernel=kernel, bw=binwid, adjust=adj, position=type, alpha=alpha, na.rm=T) +
        labs(list(title=title, x=xlab, y=ylab)) + theme_opt
    } else {
      density1 <- densitybase +
        stat_density(kernel=kernel, bw=binwid, adjust=adj, position=type, alpha=alpha, na.rm=T) +
        scale_fill_manual(values=color.group) + 
        scale_colour_manual(values=color.group) + 
        labs(list(title=title, x=xlab, y=ylab)) + theme_opt
    }
  }
  if(!is.null(densthres)){
    if(densthres > max(ggplot_build(density1)$data[[1]][,"ymax"])) stop("밀도 하한을 더 낮게 설정해야 합니다.")
    xrange <- range(ggplot_build(density1)$data[[1]][which(ggplot_build(density1)$data[[1]][,"ymax"] >= densthres),"x"])
    density1 <- density1 + xlim(xrange)
  }
  return(density1)
}

####### Box

REx_boxplot <- function(varname, dataset, by=NULL, flip=FALSE, lgd.pos="none", 
                        title=NULL, xlab=NULL, ylab=varname,
                        color="grey", color.group=NULL){
  load.pkg("ggplot2")
  
  suppressWarnings(datavar <- as.numeric(as.character(dataset[,varname])))
  if(sum(is.na(datavar))==length(datavar)) stop("There is no numeric values.")
  dataset2 <- data.frame(datavar)
  colnames(dataset2) <- "x"
  
  theme_opt <- theme_bw() + 
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = lgd.pos)
  if(is.null(by)){
    boxbase <- ggplot(aes(x='', y=x), data=dataset2) +
      geom_boxplot(fill=color, na.rm=T) +
      labs(list(title=title, x=xlab, y=ylab))
  } else {
    dataset2$group <- factor(dataset[,by])
    if(flip==TRUE) dataset2$group <- factor(dataset2$group, levels=rev(levels(dataset2$group)))
    if(sum(is.na(dataset2$group))>0) dataset2 <- dataset2[-which(is.na(dataset2$group)),]
    if(is.null(xlab)) xlab <- by
    boxbase <- ggplot(aes(x=group, y=x, fill=group), data=dataset2) +
      geom_boxplot(na.rm=T) +
      guides(color = guide_legend(by), fill = guide_legend(by)) +
      labs(list(title=title, x=xlab, y=ylab))
    if(!is.null(color.group)){
      if(flip==TRUE) color.group <- rev(color.group)
      boxbase <- boxbase + 
        scale_fill_manual(values=color.group) + 
        scale_colour_manual(values=color.group)
    }
  }
  if(flip==TRUE) boxbase <- boxbase + coord_flip()
  return(boxbase + theme_opt)
}

#####Q-Q plot

REx_qqplot <- function(varname, dataset, dist="norm", distopt=c(0,1), line=TRUE, 
                       title=NULL, xlab=NULL, ylab=varname,
                       color.dot="black", color.line="red"){
  
  load.pkg("ggplot2")

  suppressWarnings(datavar <- as.numeric(as.character(dataset[,varname])))
  if(sum(is.na(datavar))==length(datavar)) stop("There is no numeric values.")
  dataset2 <- data.frame(datavar)
  colnames(dataset2) <- "x"

  theme_opt <- theme_bw() + 
    theme(panel.grid.minor=element_blank(), plot.title = element_text(hjust = 0.5))
  
  get_qqline <- function(vec, method, distopt, xlab){
    # browser()
    y <- quantile(vec[!is.na(vec)], c(0.25, 0.75))
    n <- length(vec)
    ord <- order(vec)
    ord.x <- vec[ord]
    P <- ppoints(length(vec))
    switch(method,
           norm={
             mean1 <- distopt[1]
             sd1 <- distopt[2]
             if(sd1<0) stop("표준편차는 0 이상이어야 합니다.")
             x <- qnorm(c(0.25, 0.75), mean1, sd1)
             z <- qnorm(P, mean1, sd1)
             dv <- dnorm(z, mean1, sd1)
             if(is.null(xlab)) xlab <- substitute(paste("Quantiles of Normal Distribution (",mu,"=",mean1,", ",sigma^2,"=",sd1^2,")"))
           },
           t={
             df1 <- distopt
             if(df1<0) stop("자유도는 0보다 커야 합니다.")
             x <- qt(c(0.25, 0.75),df1)
             z <- qt(P,df1)
             dv <- dt(z,df1)
             if(is.null(xlab)) xlab <- substitute(paste("Quantiles of t Distribution (df=",df1,")"))
           },
           chisq={
             df1 <- distopt
             if(df1<0) stop("자유도는 0보다 커야 합니다.")
             x <- qchisq(c(0.25, 0.75),df1)
             z <- qchisq(P,df1)
             dv <- dchisq(z,df1)
             if(is.null(xlab)) xlab <- substitute(paste("Quantiles of Chi-squared Distribution (df=",df1,")"))
           },
           f={
             df1 <- distopt[1]
             df2 <- distopt[2]
             if(df1<0|df2<0) stop("자유도는 모두 0보다 커야 합니다.")
             x <- qf(c(0.25, 0.75),df1,df2)
             z <- qf(P,df1,df2)
             dv <- df(z,df1,df2)
             if(is.null(xlab)) xlab <- substitute(paste("Quantiles of F Distribution (df1=",df1,", df2=",df2,")"))
           },
           unif={
             min1 <- distopt[1]
             max1 <- distopt[2]
             if(min1>max1) stop("최댓값은 최솟값보다 크거나 같아야 합니다.")
             x <- qunif(c(0.25,0.75), min1, max1)
             z <- qunif(P, min1, max1)
             dv <- dunif(z, min1, max1)
             if(is.null(xlab)) xlab <- substitute(paste("Quantiles of Uniform Distribution (min=",min1,", max=",max1,")"))
           }
    )
    slope <- diff(y)/diff(x)
    int <- y[1L] - slope * x[1L]
    zz <- qnorm(1 - (1 - 0.95)/2) # confidence level = 0.95
    SE <- (slope/dv) * sqrt(P * (1 - P)/n)
    fit.value <- int + slope * z
    upper <- fit.value + zz * SE
    lower <- fit.value - zz * SE
    df <- data.frame(x=z,y=ord.x,upper=upper,lower=lower)
    return(list(intercept=int,slope=slope,df=df,xlab=xlab))
  }
  
  element_qqline <- get_qqline(na.exclude(datavar), dist, distopt, xlab)
  
  ggp <- ggplot(element_qqline$df,aes(x=x,y=y)) +
    geom_point(color=color.dot)
  
  ggp2 <- ggp + 
    geom_abline(color = color.line, slope = element_qqline$slope, intercept = element_qqline$intercept) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill=color.dot, alpha=0.2)
  
  if(line==TRUE){
    res <- ggp2 + labs(list(title=title, x=element_qqline$xlab, y=ylab)) +
      theme_opt
  } else {
    res <- ggp + labs(list(title=title, x=element_qqline$xlab, y=ylab)) +
      theme_opt
  }
  return(res)
}

#####Scatterplot

REx_scatterplot <- function(varname1, varname2, dataset, by=NULL, jitter.x=0, jitter.y=0, 
                            marginal="none", lineby=FALSE, LeastSq=FALSE, Smooth=FALSE, 
                            Sm.span=1, Sm.spread=FALSE, Ellipse=FALSE, Index=FALSE, lgd.pos="bottom",
                            ext.dot=NULL, ext.abline=NULL, ext.hline=NULL, ext.vline=NULL,
                            title=NULL, xlab=varname1, ylab=varname2,
                            color.dot="black", color.line="red", color.group=NULL, 
                            color.extdot="blue", color.extline="orange",
                            marg.hist.bin=30, marg.dens.adj=1, marg.box.size=1, color.marg="grey"){
  
  load.pkg(c("ggExtra","ggplot2"))

  suppressWarnings(x <- as.numeric(as.character(dataset[,varname1])))
  suppressWarnings(y <- as.numeric(as.character(dataset[,varname2])))
  if(sum(is.na(x))==length(x)) stop("There is no numeric values.")
  if(sum(is.na(y))==length(y)) stop("There is no numeric values.")
  dat <- data.frame(x,y)

  jitter.x <- diff(range(x, na.rm=T))*jitter.x*0.1
  jitter.y <- diff(range(y, na.rm=T))*jitter.y*0.1
  
  theme_opt <- theme_bw() + 
    theme(panel.grid.major=element_blank()) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = lgd.pos)
  
  if(is.null(by)){
    ggp <- ggplot(aes(x=x, y=y), data=dat) + geom_jitter(color=color.dot, width = jitter.x, height = jitter.y, na.rm=T)
    if(LeastSq==TRUE) ggp <- ggp + geom_smooth(method = "lm", se=FALSE, color=color.line, show.legend=F, na.rm=T)
    if(Smooth==TRUE) ggp <- ggp + geom_smooth(method="loess", span=Sm.span, se=Sm.spread, color=color.line, show.legend=F, na.rm=T)
    if(Ellipse==TRUE) ggp <- ggp + stat_ellipse(color=color.line, show.legend=F, na.rm=T)
    if(Index==TRUE) ggp <- ggp + geom_path(color=color.line, show.legend=F, na.rm=T)
  } else {
    group <- factor(dataset[,by])
    dat <- cbind(dat,group)
    if(sum(is.na(dat$group))>0) dat <- dat[-which(is.na(dat$group)),]
    ggp <- ggplot(aes(x=x, y=y, color=group, fill=group), data=dat) + geom_jitter(width = jitter.x, height = jitter.y, na.rm=T)
    if(lineby==TRUE){
      if(LeastSq==TRUE) ggp <- ggp + geom_smooth(method = "lm", se=FALSE, size=0.7, show.legend=F, na.rm=T)
      if(Smooth==TRUE) ggp <- ggp + geom_smooth(method="loess", span=Sm.span, se=Sm.spread, show.legend=F, na.rm=T)
      if(Ellipse==TRUE) ggp <- ggp + stat_ellipse(show.legend=F, na.rm=T)
      if(Index==TRUE) ggp <- ggp + geom_path(show.legend=F, na.rm=T)
    } else {
      if(LeastSq==TRUE) ggp <- ggp + geom_smooth(aes(x=x, y=y), method = "lm", se=FALSE, color=color.line, fill=color.line, size=0.7, show.legend=F, na.rm=T)
      if(Smooth==TRUE) ggp <- ggp + geom_smooth(aes(x=x, y=y), method="loess", span=Sm.span, se=Sm.spread, color=color.line, fill=color.line, show.legend=F, na.rm=T)
      if(Ellipse==TRUE) ggp <- ggp + stat_ellipse(aes(x=x, y=y), color=color.line, fill=color.line, show.legend=F, na.rm=T)
      if(Index==TRUE) ggp <- ggp + geom_path(aes(x=x, y=y), color=color.line, fill=color.line, show.legend=F, na.rm=T)
    }
    if(!is.null(color.group)) ggp <- ggp + scale_colour_manual(values=color.group) + scale_fill_manual(values=color.group)
  }
  if(!is.null(ext.dot)){
    ext.dot.dat <- grep("[0-9.]+[ ]*,[ ]*[0-9.]+",unlist(strsplit(ext.dot,"[()]")), value=T)
    ext.dot.dat <- do.call(rbind, strsplit(ext.dot.dat,","))
    ext.dot.dat <- data.frame(x=as.numeric(ext.dot.dat[,1]), y=as.numeric(ext.dot.dat[,2]))
    ggp <- ggp + geom_point(mapping=aes(x=x,y=y), data=ext.dot.dat, color=color.extdot, fill=color.extdot)
  } 
  if(!is.null(ext.abline)){
    ext.abline.dat <- grep("[0-9.]+[ ]*,[ ]*[0-9.]+",unlist(strsplit(ext.abline,"[()]")), value=T)
    ext.abline.dat <- do.call(rbind, strsplit(ext.abline.dat,","))
    ext.abline.dat <- data.frame(intercept=as.numeric(ext.abline.dat[,1]), slope=as.numeric(ext.abline.dat[,2]))
    ggp <- ggp + geom_abline(aes(intercept=intercept, slope=slope), ext.abline.dat, color=color.extline)
  }
  if(!is.null(ext.hline)) ggp <- ggp + geom_hline(yintercept = ext.hline, color=color.extline)
  if(!is.null(ext.vline)) ggp <- ggp + geom_vline(xintercept = ext.vline, color=color.extline)
  
  ggp <- ggp + theme_opt + 
    labs(list(title=title, x=xlab, y=ylab))
  if(!is.null(by)) ggp <- ggp + guides(color = guide_legend(by), fill = guide_legend(by))
  
  if(is.null(ext.dot) & is.null(ext.abline) & is.null(ext.hline) & is.null(ext.vline)) {
    if(marginal=="histogram") ggp <- ggMarginal(ggp, type = "histogram", bins=marg.hist.bin, fill=color.marg, na.rm=T)
    if(marginal=="boxplot") ggp <- ggMarginal(ggp, type = "boxplot", fill=color.marg, size=10/marg.box.size, na.rm=T)
    if(marginal=="density") ggp <- ggMarginal(ggp, type = "density", color=color.marg, adjust=marg.dens.adj, na.rm=T)
  }
  
  return(ggp)
}

#####Scatter Matrix

REx_scattermatrix <- function(varname, dataset, by=NULL, jitter.x=0, jitter.y=0,
                              diagonal="none", lineby=FALSE, LeastSq=FALSE, Smooth=FALSE, 
                              Sm.span=1, Sm.spread=FALSE, Ellipse=FALSE, Index=FALSE, lgd.pos="right",
                              title=NULL, xlab=NULL, ylab=NULL, 
                              color.dot="black", color.line="red", color.group=NULL,
                              diag.hist.bin=30, diag.dens.adj=1, color.diag="grey"){
  
  load.pkg(c("GGally","ggplot2"))
 
  dat <- dataset[,varname]
  for(i in 1:ncol(dat)){
    suppressWarnings(dat[,i] <- as.numeric(as.character(dat[,i])))
    if(sum(is.na(dat[,i]))==nrow(dat)) stop("There is no numeric values.")
  }
  n <- length(varname)
  
  theme_opt <- theme_bw() + 
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = lgd.pos)
  theme_set(theme_opt)
  
  point_jitter <- function(data, mapping, ...){
    x <- data[,as.character(mapping$x)]
    y <- data[,as.character(mapping$y)]
    jitter.x <- diff(range(x, na.rm=T))*jitter.x*0.1
    jitter.y <- diff(range(y, na.rm=T))*jitter.y*0.1
    ggplot(data = data, mapping=mapping) +
      geom_jitter(color=color.dot, width = jitter.x, height = jitter.y, na.rm=T) +
      guides(color = guide_legend(by), fill = guide_legend(by))
  }
  
  point_jitter_group <- function(data, mapping, ...){
    x <- data[,as.character(mapping$x)]
    y <- data[,as.character(mapping$y)]
    jitter.x <- diff(range(x, na.rm=T))*jitter.x*0.1
    jitter.y <- diff(range(y, na.rm=T))*jitter.y*0.1
    ggplot(data = data, mapping=mapping) +
      geom_jitter(width = jitter.x, height = jitter.y, na.rm=T) +
      guides(color = guide_legend(by), fill = guide_legend(by))
  }
  
  qqDiag <- function(data, mapping, ...){
    get_qqline <- function(vec){
      y <- quantile(vec[!is.na(vec)], c(0.25, 0.75))
      n <- length(vec)
      ord <- order(vec)
      ord.x <- vec[ord]
      P <- ppoints(length(vec))
      x <- qnorm(c(0.25, 0.75))
      z <- qnorm(P)
      dv <- dnorm(z)
      slope <- diff(y)/diff(x)
      int <- y[1L] - slope * x[1L]
      zz <- qnorm(1 - (1 - 0.95)/2) # confidence level = 0.95
      SE <- (slope/dv) * sqrt(P * (1 - P)/n)
      fit.value <- int + slope * z
      upper <- fit.value + zz * SE
      lower <- fit.value - zz * SE
      df <- data.frame(x=z,y=ord.x,upper=upper,lower=lower)
      return(list(intercept=int,slope=slope,df=df))
    }
    element_qqline <- get_qqline(na.exclude(data[,as.character(mapping$x)]))
    ggplot(element_qqline$df,aes(x=x,y=y)) +
      geom_point(color=color.diag) +
      geom_abline(slope = element_qqline$slope, intercept = element_qqline$intercept) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.2)
  }
  
  boxDiag <- function(data, mapping, ...){
    ggplot(data = data, mapping=aes(x="",y=data[,as.character(mapping$x)])) +
      geom_boxplot(fill=color.diag, na.rm=T)
  }
  
  histoDiag <- function(data, mapping, ...){
    hist <- ggplot(data = data, mapping=mapping) +
      geom_histogram(color="black", bins=diag.hist.bin, na.rm=T)
    if(!is.null(color.group)) hist <- hist + scale_fill_manual(values=color.group)
    return(hist)
  }
  
  densDiag <- function(data, mapping, ...){
    dens <- ggplot(data = data, mapping=mapping) +
      geom_density(adjust=diag.dens.adj, position="stack", na.rm=T)
    if(!is.null(color.group)) dens <- dens + scale_fill_manual(values=color.group) + scale_colour_manual(values=color.group)
    return(dens)
  }
  
  if(is.null(by)){
    if(diagonal=="none") ggp <- ggpairs(data=dat, lower=list(continuous=point_jitter), upper=list(continuous=point_jitter), diag=NULL)
    if(diagonal=="density") ggp <- ggpairs(data=dat, lower=list(continuous=point_jitter), upper=list(continuous=point_jitter), diag = list(continuous = densDiag))
    if(diagonal=="histogram") ggp <- ggpairs(data=dat, lower=list(continuous=point_jitter), upper=list(continuous=point_jitter), diag = list(continuous = histoDiag))
    if(diagonal=="qqplot") ggp <- ggpairs(data=dat, lower=list(continuous=point_jitter), upper=list(continuous=point_jitter), diag=list(continuous = qqDiag))
    if(diagonal=="boxplot") ggp <- ggpairs(data=dat, lower=list(continuous=point_jitter), upper=list(continuous=point_jitter), diag=list(continuous = boxDiag))
    for(i in 1:ggp$nrow){
      for(j in 1:ggp$ncol){
        if(i == j) next else {
          if(LeastSq==TRUE) ggp[i,j] <- ggp[i,j] + geom_smooth(method = "lm", se=FALSE, color=color.line, na.rm=T)
          if(Smooth==TRUE) ggp[i,j] <- ggp[i,j] + geom_smooth(method="loess", span=Sm.span, se=Sm.spread, color=color.line, na.rm=T)
          if(Ellipse==TRUE) ggp[i,j] <- ggp[i,j] + stat_ellipse(color=color.line, na.rm=T)
          if(Index==TRUE) ggp[i,j] <- ggp[i,j] + geom_path(color=color.line, na.rm=T)
        }
      }
    }
  } else {
    group <- factor(dataset[,by])
    dat <- cbind(dat,group)
    if(sum(is.na(dat$group))>0) dat <- dat[-which(is.na(dat$group)),]
    
    if(lgd.pos=="none") lgd.site <- NULL
    else lgd.site <- c(n,1)
    
    if(diagonal=="none") ggp <- ggpairs(data=dat, columns=1:n, ggplot2::aes(colour=group), lower=list(continuous=point_jitter_group), upper=list(continuous=point_jitter_group), diag=NULL, legend=lgd.site)
    if(diagonal=="density") ggp <- ggpairs(data=dat, columns=1:n, ggplot2::aes(colour=group, fill=group), lower=list(continuous=point_jitter_group), upper=list(continuous=point_jitter_group), diag = list(continuous = densDiag), legend=lgd.site)
    if(diagonal=="histogram") ggp <- ggpairs(data=dat, columns=1:n, ggplot2::aes(colour=group, fill=group), lower=list(continuous=point_jitter_group), upper=list(continuous=point_jitter_group), diag = list(continuous = histoDiag), legend=lgd.site)
    if(diagonal=="qqplot") ggp <- ggpairs(data=dat, columns=1:n, ggplot2::aes(colour=group), lower=list(continuous=point_jitter_group), upper=list(continuous=point_jitter_group), diag=list(continuous = qqDiag), legend=lgd.site)
    if(diagonal=="boxplot") ggp <- ggpairs(data=dat, columns=1:n, ggplot2::aes(colour=group), lower=list(continuous=point_jitter_group), upper=list(continuous=point_jitter_group), diag=list(continuous = boxDiag), legend=lgd.site)
    for(i in 1:ggp$nrow){
      for(j in 1:ggp$ncol){
        if(i == j) next else {
          if(lineby==TRUE){
            if(LeastSq==TRUE) ggp[i,j] <- ggp[i,j] + geom_smooth(method = "lm", se=FALSE, show.legend=FALSE, size=0.7, na.rm=T)
            if(Smooth==TRUE) ggp[i,j] <- ggp[i,j] + geom_smooth(method="loess", span=Sm.span, se=Sm.spread, show.legend=FALSE, na.rm=T)
            if(Ellipse==TRUE) ggp[i,j] <- ggp[i,j] + stat_ellipse(show.legend=FALSE, na.rm=T)
            if(Index==TRUE) ggp[i,j] <- ggp[i,j] + geom_path(show.legend=FALSE, na.rm=T)
          } else {
            if(LeastSq==TRUE) ggp[i,j] <- ggp[i,j] + geom_smooth(aes(color=NULL, fill=NULL), method = "lm", color=color.line, se=FALSE, show.legend=FALSE, size=0.7, na.rm=T)
            if(Smooth==TRUE) ggp[i,j] <- ggp[i,j] + geom_smooth(aes(color=NULL, fill=NULL), method="loess", color=color.line, span=Sm.span, se=Sm.spread, show.legend=FALSE, na.rm=T)
            if(Ellipse==TRUE) ggp[i,j] <- ggp[i,j] + stat_ellipse(aes(color=NULL, fill=NULL), color=color.line, show.legend=FALSE, na.rm=T)
            if(Index==TRUE) ggp[i,j] <- ggp[i,j] + geom_path(aes(color=NULL, fill=NULL), color=color.line, show.legend=FALSE, na.rm=T)
          }
          if(!is.null(color.group)) ggp[i,j] <- ggp[i,j] + scale_fill_manual(values=color.group) + scale_colour_manual(values=color.group)
        }
      }
    }
  }
  ggp <- ggp + labs(x=xlab, y=ylab, title=title)
  return(ggp)
}

##### XY condition plot
  
REx_xyplot <- function(varname1, varname2, dataset, by=NULL, jitter.x=0, jitter.y=0,
                       on_one_panel=TRUE, scales="fixed", direction_by=NULL,
                       title=NULL, xlab=paste0(varname1, collapse = "+"), ylab=paste0(varname2, collapse = "+"),
                       lgd.pos="top", color.dot="black", color.group=NULL){
  
  load.pkg("ggplot2")

  # browser()
  dat <- dataset[,c(varname1, varname2)]
  for(i in 1:ncol(dat)){
    suppressWarnings(dat[,i] <- as.numeric(as.character(dat[,i])))
    if(sum(is.na(dat[,i]))==nrow(dat)) stop("There is no numeric values.")
  }
  # n <- length(varname)
  xname <- paste0(varname1, collapse = ".")
  yname <- paste0(varname2, collapse = ".")
  if(xname==yname) xname <- paste0(xname, "_")
  
  theme_opt <- theme_bw() + 
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = lgd.pos)
  theme_set(theme_opt)
  
  nd2 <- 0
  dat2 <- list()
  for(j1 in 1:(length(varname1))){
    for(j2 in 1:(length(varname2))){
      nd2 <- nd2+1
      dat2[[nd2]] <- cbind(dat[,c(varname1[j1], varname2[j2])], varname1[j1], varname2[j2])
      colnames(dat2[[nd2]]) <- c("x", "y", xname, yname)
    }
  }
  dat2 <- do.call(rbind,dat2)
  
  jitter.x <- diff(range(dat2$x, na.rm=T))*jitter.x*0.1
  jitter.y <- diff(range(dat2$y, na.rm=T))*jitter.y*0.1
  
  if(is.null(by)){
    ggp <- ggplot(aes(x=x, y=y), data=dat2) + geom_jitter(color=color.dot, width = jitter.x, height = jitter.y, na.rm=T)
    ggp <- ggp + facet_grid(as.formula(paste0(yname, " ~ ",xname)), scales=scales)

  } else {
    dat3 <- dataset[by]
    for(i2 in 1:ncol(dat3)) dat3[,i2] <- factor(dat3[,i2])
    dat2 <- cbind(dat2, dat3)
    for(i3 in by) if(sum(is.na(dat2[,i3]))>0) dat2 <- dat2[-which(is.na(dat2[,i3])),]
    dat2$group <- interaction(dat2[,by])
    
    ggp <- ggplot(aes(x=x, y=y, color=group), data=dat2) + geom_jitter(width = jitter.x, height = jitter.y, na.rm=T) +
      guides(color = guide_legend(paste0(by, collapse=".")), fill = guide_legend(paste0(by, collapse=".")))
    
    if(on_one_panel==T) ggp <- ggp + facet_grid(as.formula(paste0(yname, " ~ ",xname)), scales=scales)
    else {
      if(length(by[direction_by=="v"])!=0) {
        xname <- paste0(xname, " + ", paste0(by[direction_by=="v"], collapse=" + "))
        xsp <- rep(c(rep(0,length(levels(interaction(dat2[,by[direction_by=="v"]])))-1),0.5), length(varname1))
        xsp <- xsp[-length(xsp)]
      } else {
        xsp <- 0
      }
      if(length(by[direction_by=="h"])!=0) {
        yname <- paste0(yname, " + ", paste0(by[direction_by=="h"], collapse=" + "))
        ysp <- rep(c(rep(0,length(levels(interaction(dat2[,by[direction_by=="h"]])))-1),0.5), length(varname2))
        ysp <- ysp[-length(ysp)]
      } else {
        ysp <- 0
      }
      # browser()
      
      ggp <- ggp + facet_grid(as.formula(paste0(yname, " ~ ",xname)), scales=scales) + 
        theme(panel.spacing.x=unit(xsp, "lines"), panel.spacing.y=unit(ysp, "lines"))
    }
    if(!is.null(color.group)) ggp <- ggp + scale_colour_manual(values=color.group)
  }
  ggp <- ggp + labs(x=xlab, y=ylab, title=title)
  return(ggp)
}

#####Mean plot

REx_meanplot <- function(varname, dataset, by1, by2=NULL, 
                         errbar="se", conflevel=0.95, line=TRUE, lgd.pos="right",
                         title=NULL, xlab=by1, ylab=varname,
                         color.dot="black", color.line="red", color.group=NULL){
  
  load.pkg("ggplot2")

  theme_opt <- theme_bw() + 
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = lgd.pos)
  pd <- position_dodge(0.1)
  
  suppressWarnings(datavar <- as.numeric(as.character(dataset[,varname])))
  if(sum(is.na(datavar))==length(datavar)) stop("There is no numeric values.")
  # dataset2 <- data.frame(datavar)
  # colnames(dataset2) <- "x"

  databy1 <- factor(dataset[,by1])
  gp1 <- levels(databy1)
  
  if(is.null(by2)){
    gpmean <- upper <- lower <- c()
    for(i in 1:length(gp1)){
      datapiece <- datavar[databy1==gp1[i]]
      gpmean[i] <- mean(datapiece, na.rm=TRUE)
      gpsd <- sd(datapiece, na.rm=TRUE)
      gpse <- gpsd/sqrt(length(datapiece))
      switch(errbar, 
             se={
               upper[i] <- gpmean[i] + gpse
               lower[i] <- gpmean[i] - gpse
             },
             sd={
               upper[i] <- gpmean[i] + gpsd
               lower[i] <- gpmean[i] - gpsd
             },
             conf={
               int <- qnorm((1-conflevel)/2, lower.tail = F)
               upper[i] <- gpmean[i] + gpse*int
               lower[i] <- gpmean[i] - gpse*int
             },
             none={
               upper[i] <- gpmean[i]
               lower[i] <- gpmean[i]
             }
      )
    } 
    meandat <- data.frame(x = gp1, y = gpmean, upper, lower)
    ggp <- ggplot(aes(x=x, y=y), data=meandat) + geom_point(color=color.dot, na.rm=TRUE)
    if(line==TRUE) ggp <- ggp + geom_line(aes(group=1), color=color.line, na.rm=TRUE)
  } else {
    databy2 <- factor(dataset[,by2])
    gp2 <- levels(databy2)
    gp2rep <- rep(gp2, each=length(gp1))
    gpmean <- upper <- lower <- c()
    for(i in 1:length(gp2)){
      for(j in 1:length(gp1)){
        ij <- (i-1)*length(gp1)+j
        datapiece <- datavar[databy1==gp1[j] & databy2==gp2[i]]
        gpmean[ij] <- mean(datapiece, na.rm=TRUE)
        gpsd <- sd(datapiece, na.rm=TRUE)
        gpse <- gpsd/sqrt(length(datapiece))
        switch(errbar, 
               se={
                 upper[ij] <- gpmean[ij] + gpse
                 lower[ij] <- gpmean[ij] - gpse
               },
               sd={
                 upper[ij] <- gpmean[ij] + gpsd
                 lower[ij] <- gpmean[ij] - gpsd
               },
               conf={
                 int <- qnorm((1-conflevel)/2, lower.tail = F)
                 upper[ij] <- gpmean[ij] + gpse*int
                 lower[ij] <- gpmean[ij] - gpse*int
               },
               none={
                 upper[i] <- gpmean[i]
                 lower[i] <- gpmean[i]
               }
        )
      }
    }
    
    meandat <- data.frame(x = gp1, y = gpmean, group = gp2rep, upper, lower)
    ggp <- ggplot(aes(x=x, y=y, color=group),data=meandat) + geom_point(position=pd, na.rm=TRUE) + 
      guides(color = guide_legend(by2), fill = guide_legend(by2))
    if(line==TRUE) ggp <- ggp + geom_line(aes(group=group),position=pd, na.rm=TRUE)
    if(!is.null(color.group)) ggp <- ggp + scale_colour_manual(values=color.group)
  }
  if(errbar!="none") ggp <- ggp + geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position=pd, na.rm=TRUE)
  return(ggp + labs(list(title=title, x=xlab, y=ylab)) + theme_opt)
}

#####Bar plot

REx_barplot <- function(varname, dataset, by=NULL, flip=FALSE, type="stack", 
                        title=NULL, xlab=varname, ylab="Frequency",
                        lgd.pos="right",color="grey", color.group=NULL){
  load.pkg("ggplot2")
  
  theme_opt <- theme_bw() + 
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    theme(strip.background = element_blank(), strip.placement = "outside") +
    theme(plot.title = element_text(hjust = 0.5), legend.position = lgd.pos)
  theme_set(theme_opt)
  # browser()
  x <- factor(dataset[,varname])
  if(flip==TRUE) x <- factor(x, levels=rev(levels(x)))
  dat <- data.frame(x)
  if(!is.null(by)) dat$group <- factor(dataset[,by])
  dat <- na.exclude(dat)
  
  if(is.null(by)){
    ggp <- ggplot(aes(x), data=dat) +
      geom_bar(fill=color, color="black", na.rm=T) +
      labs(list(title=title, x=xlab, y=ylab))
  } else {

    
    ggp <- ggplot(aes(x, fill=group), data=dat)
    if(type=="stack"){
      ggp <- ggp + geom_bar(color="black", na.rm=T)
    }
    if(type=="par"){
      ggp <- ggp + geom_bar(color="black", na.rm=T, position=position_dodge())
    }
    if(type=="facet"){
      ggp <- ggp + geom_bar(color="black", na.rm=T) +
        facet_wrap(~group, ncol=1, scales="fixed") +
        theme(legend.position = "none")
    }
    ggp <- ggp + guides(color = guide_legend(by), fill = guide_legend(by)) + labs(list(title=title, x=xlab, y=ylab))
    if(!is.null(color.group)) ggp <- ggp + scale_fill_manual(values=color.group)
  }
  if(flip==TRUE) ggp <- ggp + coord_flip()
  return(ggp)
}

#####Circle plot

REx_circleplot <- function(varname, dataset, showCount=F, title=varname, color.group=NULL, lgd.pos="right"){
  load.pkg("ggplot2")

  theme_opt <- theme_void() + 
    theme(plot.title = element_text(hjust = 0.5, size=14, face="bold"), legend.position = lgd.pos)

  x <- factor(dataset[,varname])
  dat <- data.frame(table(x))
  if(sum(is.na(x))>0) dat <- rbind(dat, data.frame(x="NA", Freq=sum(is.na(x))))
  
  pie <- ggplot(dat, aes(x = factor(1), y=Freq, fill = x)) + geom_bar(color="black", width = 1, stat = "identity")
  pie <- pie + coord_polar(theta = "y", start=pi/2) + labs(title=title) + guides(color = guide_legend(""), fill = guide_legend("")) + theme_opt
  
  if(!is.null(color.group)) pie <- pie + scale_fill_manual(values=color.group)
  if(showCount) pie <- pie + geom_text(aes(y = cumsum(Freq)-Freq/2, label=Freq), check_overlap = TRUE)

  return(pie)
}

### Get comb of group variables from XY Plot
xyplot_color.group <- function(dataset, by){
  dat3 <- dataset[by]
  for(i2 in 1:ncol(dat3)){
    dat3[,i2] <- factor(dat3[,i2])
    if(sum(is.na(dat3[,i2]))>0) dat3 <- dat3[-which(is.na(dat3[,i2])),]
  }
  group <- levels(interaction(dat3[,by]))
  return(paste(group, collapse=","))
}

### Get comb of var levels from Circle Plot
circleplot_color.group <- function(dataset, varname){
  x <- factor(dataset[,varname])
  lev <- levels(x)
  if(sum(is.na(x))>0) lev <- c(lev, "NA")
  return(paste(lev, collapse=","))
}

### Default group color
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  paste(hcl(h = hues, l = 65, c = 100)[1:n], collapse=",")
}

### Get Factor Levels
manual_color.group <- function(dataset, by){
  dat3 <- factor(dataset[,by])
  return(paste(levels(dat3), collapse=","))
}

### Testfunc
testfunc <- function() {
	list(a=rnorm(100), b=c("aa","bbb","ccccc",
	"dddd"))
}

##기술통계
REx_DESCSTAT <- function(dataset,quan_var=NULL,qual_var=NULL,group_var=NULL) {
	#### 변수 설명 ####
	## dataset : 분석에 이용할 데이터셋
	## quan_var : 데이터셋 내 분석할 양적 변수 (UI차원에서 type을 컨트롤 할 필요 없음). 필수아님. 여러변수가 선택되었을 시 c('','','')의 형식으로 입력.
	## qual_var : 데이터셋 내 분석할 질적 변수 (UI차원에서 type을 컨트롤 할 필요 없음). 필수아님. 여러변수가 선택되었을 시 c('','','')의 형식으로 입력.
	## quan_var와 qual_var 중 적어도 하나의 변수는 선택되어야 분석가능
	## group_var : 그룹변수. 최대 한개. 필수 아님.

	#############################################
	#### Required packages : R2HTML, moments ####
	#############################################
	load.pkg(c("R2HTML", "moments"))

	######################################################################################
	############################### Descriptive Statistics ###############################
	######################################################################################

	html.output <- capture.output({
		#분석 이름
		R2HTML::HTML(R2HTML::as.title("Descriptive Statistics"),HR=1,file=stdout(),append=FALSE) 
		
		# 변수 type 변환
		if(length(quan_var)!=0){
			isnom <- !sapply(dataset[,quan_var,drop=F],is.numeric)
			if(sum(isnom)>0){
				isare <- ifelse(sum(isnom)==1,'is','are')
				waswere <- ifelse(sum(isnom)==1,'was','were')
				ItThey <- ifelse(sum(isnom)==1,'It','They')
				warn.msg1 <- paste0("\a Warning : '",paste(names(isnom)[isnom],collapse=", "),"' ",isare," not continuous but ",waswere," chosen as continuous variable. ",ItThey," ",waswere," excluded from the analysis.")
				quan_var <- quan_var[!isnom]
				if(length(quan_var)==0) quan_var <- NULL
			}
		}
		if(length(qual_var)!=0){
			isnum <- sapply(dataset[,qual_var,drop=F],is.numeric)
			if(sum(isnum)>0){
				isare <- ifelse(sum(isnum)==1,'is','are')
				waswere <- ifelse(sum(isnum)==1,'was','were')
				ItThey <- ifelse(sum(isnum)==1,'It','They')
				warn.msg2 <- paste0("\a Warning : '",paste(names(isnum)[isnum],collapse=", "),"' ",isare," continuous but ",waswere," chosen as categorical variable. ",ItThey," ",waswere," coerced into character.")
				for(i in qual_var[isnum]) dataset[,i] <- as.factor(dataset[,i])
			}
		}
		if(!is.null(group_var)) dataset[,group_var] <- as.factor(dataset[,group_var])
		if(is.null(c(quan_var,qual_var))) warn.msg3 <- warn.msg2 <- "\a Error: At least one variable should be selected. Analysis has been stopped."

		if(exists('warn.msg3')){
			HTML(as.title("Warnings"),HR=2,file=stdout())
			if(exists('warn.msg1')) HTML(warn.msg1,file=stdout())
			HTML(warn.msg3,file=stdout())
		} else {
			## Data Structure
			HTML(as.title("Data structure"),HR=2,file=stdout())
			total.var <- ncol(dataset)
			used.var <- length(c(quan_var,qual_var))
			DS <- matrix(c('Number of observations','Number of total variables','Number of used variables',nrow(dataset),total.var,used.var),ncol=2)
			HTML(DS,file=stdout(), innerBorder = 1,align="left")

			#선택변수 리스트
			R2HTML::HTML(R2HTML::as.title("Variable List"),HR=2,file=stdout()) 
			VL <- matrix(0,nrow=1,ncol=2)
			if(length(quan_var)!=0) VL <- rbind(VL,c('Continuous variable',paste(quan_var,collapse=', ')))
			if(length(qual_var)!=0) VL <- rbind(VL,c('Categorical variable',paste(qual_var,collapse=', ')))
			if(length(group_var)!=0) VL <- rbind(VL,c('Group variable',group_var))
			VL <- VL[-1,,drop=F]
			HTML(VL,file=stdout(), innerBorder = 1,align="left")
			if(exists('warn.msg1')) R2HTML::HTML(warn.msg1,file=stdout())
			if(exists('warn.msg2')) R2HTML::HTML(warn.msg2,file=stdout())

			#분석결과	
			R2HTML::HTML(R2HTML::as.title("Results of Descriptive Statistics"),HR=2,file=stdout())
			if(length(group_var)==0){
				if(length(quan_var)!=0){
					R2HTML::HTML(R2HTML::as.title("Continuous variable"),HR=3,file=stdout())
					result <- list()
					resulttable <- data.frame()
					for (i in 1:length(quan_var)) {
						data <- dataset[,quan_var[i]] 
						n <- length(data[!is.na(data)]) #N
						min <- min(data[!is.na(data)]) #최소값
						max <- max(data[!is.na(data)]) #최대값
						range <- max-min #범위
						sum <- sum(data[!is.na(data)]) #합계
						mean <- mean(data[!is.na(data)]) #평균 
						median <- median(data[!is.na(data)]) #평균 
						mean_se <- mean/sqrt(n) #평균의 표준오차
						sd <- sd(data[!is.na(data)]) #표준편차
						var <- var(data[!is.na(data)]) #분산
						skew <- skewness(data[!is.na(data)]) #왜도
						ifelse(n>3, skew_se <- sqrt(6*n*(n-1)/((n-2)*(n+1)*(n+3))),NA) #왜도의 표준오차
						kurt <- kurtosis(data[!is.na(data)]) #첨도
						ifelse(n>3, kurt_se <- 2*skew_se*sqrt((n^2-1)/((n-3)*(n+5))),NA) #첨도의 표준오차
						num_miss <- sum(is.na(data))
						
						#결과 저장
						result[[i]] <- data.frame("N"=n,"Missing"=num_miss,"Min"=min,"Max"=max,"Range"=range,"Sum"=sum,"Mean"=mean,"se(Mean)"=mean_se,
										"Median"=median,"Standard Deviation"=sd,
										"Variance"=var,"Skewness"=skew,"se(Skewness)"=skew_se,
										"Kurtosis"=kurt,"se(Kurtosis)"=kurt_se)
						resulttable <- rbind(resulttable,result[[i]])
					}
					colnames(resulttable) <- c("N.observed","N.missing","Min","Max","Range","Sum","Mean","se(Mean)","Median","SD","Variance","Skewness","se(Skewness)","Kurtosis","se(Kurtosis)")
					rownames(resulttable) <- quan_var
					R2HTML::HTML(round(resulttable,4),file=stdout(),innerBorder = 1,align="left") #결과 출력
				}
				if(length(qual_var)!=0){
					R2HTML::HTML(R2HTML::as.title("Categorical variable"),HR=3,file=stdout())
					for(ii in qual_var){
						R2HTML::HTML(R2HTML::as.title(ii),HR=4,file=stdout())
						Nom <- dataset[,ii,drop=F]
						tab <- table(Nom)
						lev <- c(names(tab),'Missing','Total')
						freq <- c(as.vector(tab),sum(is.na(Nom)),nrow(Nom))
						rela.freq <- round(freq/sum(nrow(Nom)),4)
						fin.tab <- data.frame(Level=lev,Freq=freq,Rela_Freq=rela.freq)
						rownames(fin.tab) <- NULL
						R2HTML::HTML(fin.tab,file=stdout(), innerBorder = 1,align="left",row.names=F)
					}
				}
			} else {
				G <- sort(unique(dataset[,group_var]))
				if(length(quan_var)!=0){
					R2HTML::HTML(R2HTML::as.title("Continuous variable"),HR=3,file=stdout())
					for(g in G){
						R2HTML::HTML(R2HTML::as.title(paste0(group_var,' - ',g)),HR=4,file=stdout())
						ind <- which(dataset[,group_var]==g)
						result <- list()
						resulttable <- data.frame()
						for (i in 1:length(quan_var)) {
							data <- dataset[ind,quan_var[i]] 
							n <- length(data[!is.na(data)]) #N
							min <- min(data[!is.na(data)]) #최소값
							max <- max(data[!is.na(data)]) #최대값
							range <- max-min #범위
							sum <- sum(data[!is.na(data)]) #합계
							mean <- mean(data[!is.na(data)]) #평균 
							median <- median(data[!is.na(data)]) #평균 
							mean_se <- mean/sqrt(n) #평균의 표준오차
							sd <- sd(data[!is.na(data)]) #표준편차
							var <- var(data[!is.na(data)]) #분산
							skew <- skewness(data[!is.na(data)]) #왜도
							ifelse(n>3, skew_se <- sqrt(6*n*(n-1)/((n-2)*(n+1)*(n+3))),NA) #왜도의 표준오차
							kurt <- kurtosis(data[!is.na(data)]) #첨도
							ifelse(n>3, kurt_se <- 2*skew_se*sqrt((n^2-1)/((n-3)*(n+5))),NA) #첨도의 표준오차
							num_miss <- sum(is.na(data))
							
							#결과 저장
							result[[i]] <- data.frame("N"=n,"Missing"=num_miss,"Min"=min,"Max"=max,"Range"=range,"Sum"=sum,"Mean"=mean,"se(Mean)"=mean_se,
											"Median"=median,"Standard Deviation"=sd,
											"Variance"=var,"Skewness"=skew,"se(Skewness)"=skew_se,
											"Kurtosis"=kurt,"se(Kurtosis)"=kurt_se)
							resulttable <- rbind(resulttable,result[[i]])
						}
						colnames(resulttable) <- c("N.observed","N.missing","Min","Max","Range","Sum","Mean","se(Mean)","Median","SD","Variance","Skewness","se(Skewness)","Kurtosis","se(Kurtosis)")
						rownames(resulttable) <- quan_var
						R2HTML::HTML(round(resulttable,4),file=stdout(),innerBorder = 1,align="left") #결과 출력
					}
				}
				if(length(qual_var)!=0){
					R2HTML::HTML(R2HTML::as.title("Categorical variable"),HR=3,file=stdout())
					for(ii in qual_var){
						R2HTML::HTML(R2HTML::as.title(ii),HR=4,file=stdout())
						Nom <- dataset[,c(ii,group_var),drop=F]
						tab <- table(Nom)
						lev <- c(rownames(tab),'Missing','Total')
						freq <- rbind(matrix(tab,nrow=nrow(tab)),as.vector(table(Nom[is.na(Nom[,1]),2])))
						freq <- rbind(freq,colSums(freq))
						fin.tab <- data.frame(Level=lev,as.data.frame(matrix(NA,ncol=length(levels(Nom[,2]))*2,nrow=length(levels(Nom[,1]))+2)))
						fin.tab[,seq(2,length(levels(Nom[,2]))*2,2)] <- freq
						for(jj in 1:length(levels(Nom[,2]))) fin.tab[,jj*2+1] <- round(freq[,jj]/freq[nrow(freq),jj],4)
						rownames(fin.tab) <- NULL
						colnames(fin.tab) <- c('Level',paste0(rep(c('Freq(grp=','Rela_Freq(grp='),length(levels(Nom[,2]))),rep(levels(Nom[,2]),each=2),')'))
						R2HTML::HTML(fin.tab,file=stdout(), innerBorder = 1,align="left",row.names=F)
					}
				}
			}
		}

		#분석종료시간
		R2HTML::HTMLhr(file=stdout())
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Descriptive Statistics",sep=""),file=stdout())
		R2HTML::HTMLhr(file=stdout())
	})
	return(html.output)
}



REx_DataExplorer <- function(dataset,variables=NULL,group=NULL,allvar=FALSE){
	#### 변수 설명 ####
	## variables : UI의 b항목인 분석변수목록. c("","","")와 같은 형태로 변수 이름이 표시되어야함. ex) iris에서 c("Sepal.Length","Species")
	## group : UI의 d항목인 집단변수. nominal variable만 가능하고 최대 1개 선택가능. variables와 중복선택 될 수 없음
	## allvar : UI의 c항목을 선택하면 TRUE. 이 값이 TRUE이면 variable과 group은 선택될 수 없고 NULL을 가짐.
	###################
	#### Required packages : R2HTML ####
	###################
	load.pkg("R2HTML")
	
	# Summary function for numerical variable
	numsum <- function(Num,groups=NULL,G=NULL,Numeric){
		R2HTML::HTML(R2HTML::as.title("Numerical Variable"),HR=3,file=stdout())
		# Numerical
		if(length(groups)==0){
			summary.num <- summary(Num)
			summary.num <- gsub("^.+:","",summary.num)
			summary.num <- gsub("[ ]+","",summary.num)
			if(any(colSums(is.na(Num))>0)){
				summary.num[is.na(summary.num)] <- 0
			} else {
				summary.num <- rbind(summary.num,0)
			}
			summary.num <- t(summary.num)
			colnames(summary.num) <- c("Min","Q1","Median","Mean","Q3","Max","NA's")
			R2HTML::HTML(summary.num,file=stdout(), innerBorder = 1,align="left")
		} else {
			for(i in 1:length(groups)) {
				R2HTML::HTML(R2HTML::as.title(paste0("Group ",i,") ",groups[i])),HR=4,file=stdout())
				summary.num <- summary(Num[which(G==groups[i]),,drop=F])
				summary.num <- gsub("^.+:","",summary.num)
				summary.num <- gsub("[ ]+","",summary.num)
				if(any(colSums(is.na(Num))>0)){
					summary.num[is.na(summary.num)] <- 0
				} else {
					summary.num <- rbind(summary.num,0)
				}
				summary.num <- t(summary.num)
				colnames(summary.num) <- c("Min","Q1","Median","Mean","Q3","Max","NA's")
				R2HTML::HTML(summary.num,file=stdout(), innerBorder = 1,align="left")
			}
		}
	}

	# Summary function for nominal variable
	nomsum <- function(Nom,newdat,groups=NULL,G=NULL,Numeric){
		R2HTML::HTML(R2HTML::as.title("Nominal Variable"),HR=3,file=stdout())
		# Nominal
		if(length(groups)==0){
			summary.nom <- summary(Nom)
			summary.nom[is.na(summary.nom)] <- ""
			rownames(summary.nom) <- NULL
			R2HTML::HTML(summary.nom,file=stdout(), innerBorder = 1,align="left")
		} else {
			for(i in 1:length(groups)) {
				R2HTML::HTML(R2HTML::as.title(paste0("Group ",i,") ",groups[i])),HR=4,file=stdout())
				summary.nom <- summary(Nom[which(G==groups[i]),,drop=F])
				summary.nom[is.na(summary.nom)] <- ""
				rownames(summary.nom) <- NULL
				R2HTML::HTML(summary.nom,file=stdout(), innerBorder = 1,align="left")
			}
		}
	}
	
	html.output <- capture.output({
		## re-denote dataset in order to contain list of variables
		if(allvar) {
			newdat <- dataset
		} else {
			newdat <- dataset[,variables,drop=F]
		}
		
		### HTML Output
		## Title
		R2HTML::HTML(R2HTML::as.title("Data Exploration"),HR=1,file=stdout(),append=FALSE)

		## Data structure
		R2HTML::HTML(R2HTML::as.title("Data Structure"),HR=2,file=stdout())
		if(length(group)==0){
			DS <- matrix(c('Number of observations','Number of variables',nrow(newdat),ncol(newdat)),nrow=2)
		} else {
			DS <- matrix(c('Number of observations','Number of variables',nrow(newdat),paste(ncol(newdat)+1,"(including 1 group variable)")),nrow=2)
		}	
		R2HTML::HTML(DS,file=stdout(), innerBorder = 1,align="left")

		## Variable List
		R2HTML::HTML(R2HTML::as.title("Variable List"),HR=2,file=stdout())
		Numeric <- which(sapply(newdat,class)=="numeric")
		CN <- colnames(newdat)
		if(length(Numeric)==0){
			varlist <- matrix(c('Nominal Variable',paste(CN,collapse=", ")),ncol=2)
		} else if(length(Numeric)==length(CN)){
			varlist <- matrix(c('Numerical Variable',paste(CN,collapse=", ")),ncol=2)
		} else {
			num <- CN[Numeric]
			nom <- CN[-Numeric]
			varlist <- matrix(c('Numerical Variable','Nominal variable',paste(num,collapse=", "),paste(nom,collapse=", ")),ncol=2)
		}
		if(!length(group)==0) varlist <- rbind(varlist,c('Group Variable',group))

		R2HTML::HTML(varlist,file=stdout(), innerBorder = 1,align="left")

		## Variable summary
		R2HTML::HTML(R2HTML::as.title("Variable Summary"),HR=2,file=stdout())
		if(length(group)==0){
			if(length(Numeric)==0){
				nomsum(Nom=newdat,newdat=newdat,Numeric=Numeric)
			} else if(length(Numeric)==length(CN)){
				numsum(Num=newdat,Numeric=Numeric)
			} else {
				Num <- newdat[,num,drop=F]
				Nom <- newdat[,nom,drop=F]
				numsum(Num=Num,Numeric=Numeric)
				nomsum(Nom=Nom,newdat=newdat,Numeric=Numeric)
			}
		} else {
			G <- as.character(dataset[,group])
			groups <- unique(G)
			if(length(Numeric)==0){
				nomsum(Nom=newdat,newdat=newdat,groups=groups,G=G,Numeric=Numeric)
			} else if(length(Numeric)==length(CN)){
				numsum(Num=newdat,groups=groups,G=G,Numeric=Numeric)
			} else {
				Num <- newdat[,num,drop=F]
				Nom <- newdat[,nom,drop=F]
				numsum(Num=Num,groups=groups,G=G,Numeric=Numeric)
				nomsum(Nom=Nom,newdat=newdat,groups=groups,G=G,Numeric=Numeric)
			}
		}

		# Analysis end time
		R2HTML::HTMLhr(file=stdout())
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Data Exploration.",sep=""),file=stdout())
		R2HTML::HTMLhr(file=stdout())

	})
	return(html.output)
}

REx_VAR_TEST <- function(dataset,res_var,group_var,methods='Bartlett',centering='mean',alternative_hypothesis='two.sided',CI=FALSE,confi.level=0.95,digits=3){
  #### 변수 설명 ####
  ## res_var: 반응변수 (변수선택-반응변수) 1개만 선택가능. 반드시 1개 선택
  ## group_var: 집단변수 (변수선택-집단변수). 반드시 1개 선택. 여러개 선택시 + 로 묶어서 들어감.
  ## methods: 분석방법선택(Bartlett, Levene, F-test 가능)
  ## centering: methods에서 Levene이 선택되었을 때 활성화되는 중앙(센터)에 해당. ('median','mean' - 각각 중앙값과 평균)
  ## alternative: 대립가설 형식('two.sided','greater','less'), default는 two.sided, methods에서 F-test를 선택하였을때만 활성화.
  ## CI : 신뢰구간을 선택하면 TRUE(default는 False). methods에서 F-test를 선택하였을때만 활성화.
  ## conf.level: 신뢰수준 (신뢰수준을 선택하면 활성화된다. default=0.95)
  ## digits : 반올림 소수점 자리수 (출력되는 소수점 자릿수-몇번째 소수점까지 출력할지)
  ###################
  #### Required packages : R2HTML, car
  ###################
  
  load.pkg(c("car","plyr"))

  html.output <- capture.output({
    R2HTML::HTML(R2HTML::as.title("Test for Homogeneity of Variances"),HR=1,file='./VarTest.html',append=FALSE)
    
    # Warning
    if(class(dataset[,res_var])!="numeric") {
      warn.msg1 <- paste0("\a Error : Response variable is numeric but selected as the qualitative variable. It was coreced into character.")
    }
    
    #group_var
    groupvar_list<-unlist(strsplit(group_var,"\\+"))
    var_num<-which(colnames(dataset) %in% groupvar_list)
    dataset$group_var <- factor(do.call("paste", c(dataset[var_num], sep=" ")))
    is.num <- sapply(groupvar_list,function(i) is.numeric(dataset[,i]))
    if(any(is.num)) warn.msg2 <- paste0("\a Warning : The type of variable '",paste(groupvar_list[is.num],collapse=', '),"' is numeric but selected as the group variable. It was coreced into character.")
    
    # Warning
    if(sum(table(unlist(strsplit(group_var,"\\+"))) <2)>=1) {
      warn.msg2 <- paste0("\a Error : Group variable must be at least 2 observations in each group.")
    }
    
    if(exists('warn.msg1')){
      R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file='./VarTest.html')
      R2HTML::HTML(warn.msg1,file='./VarTest.html')
    } else if(exists('warn.msg2')) {
      R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file='./VarTest.html')
      R2HTML::HTML(warn.msg2,file='./VarTest.html')
    } else{
      ## Data structure
      R2HTML::HTML(R2HTML::as.title("Data Structure"),HR=2,file='./VarTest.html')
      DS <- matrix(c('# of obs','# of total variables',nrow(dataset),length(strsplit(group_var,'\\+')[[1]])+1),nrow=2)
      R2HTML::HTML(DS,file='./VarTest.html', innerBorder = 1,align="left")
      
      ## Analysis description
      R2HTML::HTML(R2HTML::as.title("Analysis Description"),HR=2,file='./VarTest.html') 
      AD <- matrix(c('Method',methods,'Response variable',res_var,'Group variable',group_var),byrow=T,ncol=2)
      if(exists('warn.msg2')) R2HTML::HTML(warn.msg2,file='./VarTest.html')
      
      if(methods=='Bartlett'){
        AD1 <- matrix(c('Null Hypothesis (H0)','Variances in each of the groups are the same.','Alternative Hypothesis (H1)','not H0.'),byrow=T,ncol=2)
        AD <- rbind(AD,AD1)
        R2HTML::HTML(AD,file='./VarTest.html', innerBorder = 1,align="left")
      } else if(methods=='Levene') {
        AD1 <- matrix(c('Centering Method',centering,'Null Hypothesis (H0)','Variances in each of the groups are the same.','Alternative Hypothesis (H1)','not H0.'),byrow=T,ncol=2)
        AD <- rbind(AD,AD1)
        R2HTML::HTML(AD,file='./VarTest.html', innerBorder = 1,align="left")
      } else {
        H1 <- ifelse(alternative_hypothesis=="two.sided","Variance ratios between groups are not equal to 1",ifelse(alternative_hypothesis=="less","Variance ratios between groups are less than 1","Variance ratios between groups are greater than 1"))
        AD1 <- matrix(c('Null Hypothesis (H0)','Variances in each of the groups are the same.','Alternative Hypothesis (H1)',H1,'Note','F-test is conducted in pairs of groups.'),byrow=T,ncol=2)
        AD <- rbind(AD,AD1)
        R2HTML::HTML(AD,file='./VarTest.html', innerBorder = 1,align="left")
      }
      
      #분석결과
      if(methods=='Bartlett'){
        # Bartlett methods 
        R2HTML::HTML(R2HTML::as.title("Bartlett Test of homogeneity of variances"),HR=2,file='./VarTest.html') 
        
        if(length(strsplit(group_var,split='\\+')[[1]])==1){
          result<-bartlett.test(formula(paste0(res_var, '~ ' ,group_var)), data=dataset)
          results <- matrix(c("Bartlett's K-squared",round(result$statistic,digits),"degree of freedom",result$parameter,"P-value",round(result$p.value,digits)),byrow=T,ncol=2)
        } else{
          interaction.term<-gsub(pattern="\\+", replacement=", ", x=group_var)
          result<- bartlett.test(formula(paste0(res_var, '~interaction(' , interaction.term,')')), data=dataset)
          results <- matrix(c("Bartlett's K-squared",round(result$statistic,digits),"degree of freedom",result$parameter,"P-value",round(result$p.value,digits)),byrow=T,ncol=2)
        }
        R2HTML::HTML(results,file='./VarTest.html', innerBorder = 1,align="left")
      } else if(methods=="Levene"){
        # Levene methods 
        R2HTML::HTML(R2HTML::as.title("Levene's Test for homogeneity of variance"),HR=2,file='./VarTest.html') 
        
        if(length(strsplit(group_var,split='\\+')[[1]])==1){
          result<-leveneTest(formula(paste0(res_var ,'~ as.factor(' ,group_var,')')), data=dataset,center=centering)
          results <- matrix(c("F value",round(result$'F value'[1],digits),"degree of freedom(numerator)",result$Df[1],"degree of freedom(denominator)",result$Df[2],'p-value',round(result$'Pr(>F)'[1],digits)),byrow=T,ncol=2)
        } else {
          result<-leveneTest(formula(paste0(res_var, '~ as.factor(' , gsub(pattern="\\+", replacement="*", x=group_var,')'))), data=dataset,center=centering)
          results <- matrix(c("F value",round(result$'F value'[1],digits),"degree of freedom(numerator)",result$Df[1],"degree of freedom(denominator)",result$Df[2],'p-value',round(result$'Pr(>F)'[1],digits)),byrow=T,ncol=2)
        }
        R2HTML::HTML(results,file='./VarTest.html', innerBorder = 1,align="left")
      } else if(methods=="F-test") {
        # F-test
        R2HTML::HTML(R2HTML::as.title("F Test for homogeneity of variance"),HR=2,file='./VarTest.html') 
        
        grp.var0 <- strsplit(group_var,"\\+")[[1]]
        grp.var <- apply(dataset[,grp.var0,drop=F],1,function(x) paste(x,collapse='.'))
        dataset$grp.var <- grp.var
        
        a <- combn(unique(grp.var),2)
        for(i in 1:ncol(a)){
          dataset1 <- dataset[dataset$grp.var%in%a[,i],]
          result<-var.test(formula(paste0(res_var,"~ grp.var")),alternative=alternative_hypothesis,conf.level=confi.level,data=dataset1) 
          H1 <- ifelse(alternative_hypothesis=="two.sided","Variance ratio is not equal to 1",ifelse(alternative_hypothesis=="less","Variance ratio is less than 1","Variance ratio is greater than 1"))
          if(CI){
            Fresult <- matrix(c("Comparison groups",paste0(a[1,i],' & ',a[2,i]),"F statistic",round(result$statistic,digits),"degree of freedom(numerator)",result$parameter[1],"degree of freedom(denominator)",result$parameter[2],"P-value",round(result$p.value,digits),"Ratio Estimate",round(result$estimate,digits),paste0(confi.level*100,"% CI"),paste0("(",round(result$conf.int[1],digits),",",round(result$conf.int[2],digits),")")),byrow=T,ncol=2)
          } else {
            Fresult <- matrix(c("Comparison groups",paste0(a[1,i],' & ',a[2,i]),"F statistic",round(result$statistic,digits),"degree of freedom(numerator)",result$parameter[1],"degree of freedom(denominator)",result$parameter[2],"P-value",round(result$p.value,digits),"Ratio Estimate",round(result$estimate,digits)),byrow=T,ncol=2)
          }
          cap <- "Ratio Estimate : Estimate of variance ratio between groups"
          R2HTML::HTML(Fresult,file='./VarTest.html',caption=cap, innerBorder = 1,align="left") # 분산 비 추정값
        }
        
      }
    }
    #분석종료시간
    R2HTML::HTMLhr(file='./VarTest.html')
    R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Test for Homogeneity of Variances",sep=""),file=stdout())
    R2HTML::HTMLhr(file='./VarTest.html')   
  })
  
  return(html.output)       
}


REx_TSLS <- function(dataset,res_var,quan_var=NULL,qual_var=NULL,exp_var=NULL,ins_var=NULL,noint=FALSE,CI=TRUE,confint.level=0.95,ANOVA=FALSE,Predict=FALSE,Resid=FALSE,stdResid=FALSE){
	#### 변수 설명 ####
	## res_var : 종속변수 (필수 1개 선택되어야 함. 숫자형 변수만 가능)
	## qual_var : 질적변수 (2개 이상의 변수가 있는 경우 c(,,)로 구분, 필수사항아님)
	## quan_var : 양적변수 (2개 이상의 변수가 있는 경우 c(,,)로 구분, 필수사항아님)
	## exp_var : 설명변수 (2개 이상의 변수가 있는 경우 c(,,)로 구분, 필수사항아님)
	## ins_var : 도구변수 (2개 이상의 변수가 있는 경우 c(,,)로 구분, 필수사항아님)
	## noint : 모델설정 탭의 '상수항 포함하지 않음' 옵션, 이 옵션이 활성화되면 양적변수, 질적변수, 도구변수는 모두 선택되어야 함.
	## CI : 신뢰구간
	## confint.level : 신뢰수준
	## ANOVA : 출력옵션 탭의 'ANOVA Table' 옵션
	## Predict : 출력옵션 탭의 '예측값' 옵션
	## Resid : 출력옵션 탭의 '잔차' 옵션
	## stdResid: 출력옵션 탭의 '표준화잔차' 옵션
	###################
	#### Required packages : AER, R2HTML ####
	###################
	load.pkg("AER")
	options("scipen"=999,"digits"=4)

	html.output <- capture.output({
		# Title
		R2HTML::HTML(R2HTML::as.title("Two-Stage Least Square Regression"),HR=1,file=stdout(),append=FALSE)
	
		## Warnings
		# Response variable type
		if(!is.numeric(dataset[,res_var])) warn.msg1 <- '\a Error : Response variable should be numeric. Analysis has been stopped.'

		vars <- c(exp_var,ins_var)
		if(is.null(vars)) {
			Vars <- c()
		} else {
			Vars <-  unique(unlist(strsplit(vars,":")))
		}
		qual_var <- qual_var[qual_var%in%Vars]
		quan_var <- c(quan_var[quan_var%in%Vars],res_var)

		# explanatory & instrumental variable type
		if(!is.null(quan_var)) {
			is.nom <- sapply(quan_var,function(i) !is.numeric(dataset[,i]))
			if(any(is.nom)){
				warn.msg2 <- paste0("\a Warning : The type of variable '",paste(quan_var[is.nom],collapse=', '),"' is not continuous but selected as the quantitative variable. It was excluded from the analysis.")
				ec <- quan_var[is.nom]
				# explanatory variables
				exp_var <- exp_var[!exp_var%in%ec]
				if(length(exp_var)==0) exp_var <- NULL
				# instrumental variables
				ins_var <- ins_var[!ins_var%in%ec]
				if(length(ins_var)==0) ins_var <- NULL
			}
		}
		if(!is.null(qual_var)) {
			is.num <- sapply(qual_var,function(i) is.numeric(dataset[,i]))
			if(any(is.num)) warn.msg3 <- paste0("\a Warning : The type of variable '",paste(qual_var[is.num],collapse=', '),"' is continuous but selected as the qualitative variable. It was coreced into character.")
			for(i in qual_var) dataset[,i] <- factor(dataset[,i])
		}

		# About explanatory & instrumental variables
		dim.exp <- ifelse(is.null(exp_var),0,ncol(model.matrix(formula(paste0(res_var,"~",paste(exp_var,collapse="+"))),data=dataset)))
		dim.ins <- ifelse(is.null(ins_var),0,ncol(model.matrix(formula(paste0(res_var,"~",paste(ins_var,collapse="+"))),data=dataset)))
		if(dim.exp>dim.ins) warn.msg4 <- paste("\aError : The number of explanatory variables should be less than that of instruments. Your analysis has ",dim.exp," explanatory variables and ",dim.ins," instruments. (All qualitative variables are transformed to dummy variables.)",sep="")

		# no intercept & no explanatory variable : error
		if(noint & is.null(exp_var)) warn.msg5 <- "\aError : With no intercept, at least 1 explanatory variable should be selected. Analysis has been stopped."

		# warnings
		if(exists('warn.msg1')|exists('warn.msg4')|exists('warn.msg5')){
			R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
			if(exists('warn.msg1')) R2HTML::HTML(warn.msg1,file=stdout())
			if(exists('warn.msg5')){
				if(exists('warn.msg2')) R2HTML::HTML(warn.msg2,file=stdout())
				R2HTML::HTML(warn.msg5,file=stdout())
			}
			if(exists('warn.msg4')) R2HTML::HTML(warn.msg4,file=stdout())
		} else {
			### Model formula
			if(is.null(exp_var)){
				if(is.null(ins_var)){
					form.0 <- paste0(res_var,' ~ 1 | 1')
				} else {
					form.0 <- paste0(res_var,' ~ 1 | ',paste(ins_var,collapse=' + '))
				}
			} else {
				form.0 <- paste0(res_var,' ~ ',paste(exp_var,collapse=' + '),' | ',paste(ins_var,collapse=' + '))
			}
			if(noint) form.0 <- paste(gsub('\\|','-1 \\|',form.0),-1)
			form.1 <- as.formula(form.0)

			### Default output
			res_TSLS <- AER::ivreg(form.1,x=TRUE,data=dataset)

			## Data Structure
			R2HTML::HTML(R2HTML::as.title("Data Structure"),HR=2,file=stdout())
			total.var <- ncol(dataset)
			used.var <- ifelse(is.null(c(exp_var,ins_var)),0,length(unique(unlist(strsplit(c(exp_var,ins_var),":")))))
			DS <- matrix(c('Number of observations','Number of total variables','Number of used variables',nrow(dataset),total.var,used.var+1),ncol=2)
			R2HTML::HTML(DS,file=stdout(), innerBorder = 1,align="left")

			## Varibale list
			R2HTML::HTML(R2HTML::as.title("Variable List"),HR=2,file=stdout())
			if(is.null(c(exp_var,ins_var))) {
				Vars <- c()
			} else {
				Vars <-  unique(unlist(strsplit(c(exp_var,ins_var),":")))
			}
			qual <- qual_var[qual_var%in%Vars]
			quan <- c(quan_var[quan_var%in%Vars],res_var)
			varlist <- matrix(c('Quantitative variable',paste(quan,collapse=', ')),ncol=2,byrow=T)
			if(length(qual)!=0) varlist <- rbind(varlist,c('Qualitative variable',paste(qual,collapse=', ')))
			R2HTML::HTML(varlist,file=stdout(), innerBorder = 1,align="left")
			if(exists('warn.msg2')) R2HTML::HTML(warn.msg2,file=stdout())
			if(exists('warn.msg3')) R2HTML::HTML(warn.msg3,file=stdout())

			## Analysis Description
			R2HTML::HTML(R2HTML::as.title("Analysis Description"),HR=2,file=stdout())
			AD <- matrix(c('Response variable',res_var),ncol=2,byrow=T)
			if(!is.null(exp_var)) AD <- rbind(AD,c('Explanatory variable',paste(exp_var,collapse=', ')))
			if(!is.null(ins_var)) AD <- rbind(AD,c('Instrumental variable',paste(ins_var,collapse=', ')))
			AD <- rbind(AD,matrix(c('Intercept included',!noint),ncol=2,byrow=T))
			R2HTML::HTML(AD,file=stdout(), innerBorder = 1,align="left")

			## Analysis Results
			R2HTML::HTML(R2HTML::as.title("Results of Two-Stage Least Square Regression"),HR=2,file=stdout())

			# Coefficient Estimate
			R2HTML::HTML(R2HTML::as.title("Coefficients"),HR=3,file=stdout())
			if(!is.null(c(exp_var,ins_var))){
				CE <- summary(res_TSLS)$coef
				if(CI){
					CItab <- confint(res_TSLS, level=confint.level)
					CItab <- CItab[match(rownames(CE),rownames(CItab)),]
					CE <- cbind(CE,CItab)
					colnames(CE)[(ncol(CE)-1):ncol(CE)] <- paste0(c('Lower','Upper'),' bound of ',confint.level*100,'% CI')
				}
				R2HTML::HTML(round(CE,4),file=stdout(), innerBorder = 1,align="left")
			} else {
				R2HTML::HTML(res_TSLS$coef,file=stdout(), innerBorder = 1,align="left")
				R2HTML::HTML("\aWarning : Only estimate for intercept is provided for for intercept only model.",file=stdout())
			}
		
			# Anova table
			if(ANOVA){
				R2HTML::HTML(R2HTML::as.title("ANOVA Table"),HR=3,file=stdout())
				if(!noint){
					if(!is.null(c(exp_var,ins_var))){
						cc <- rownames(model.matrix(res_TSLS))		# Complete cases
						temp.dataset <- dataset[rownames(dataset)%in%cc,]
						Anova <- as.matrix(anova(AER::ivreg(formula(paste(res_var,'~1|1',sep='')),x=TRUE,data=temp.dataset),res_TSLS))
						A1 <- round(rbind(Anova[2,4:3], Anova[2:1,c(2,1)]),5)
						MS <- c(round(A1[1:2,1]/A1[1:2,2],5)," ")
						A2 <- data.frame(round(A1,5),MS=MS,F=c(round(Anova[2,5],5)," "," "),Pvalue=c(format(Anova[2,6],scientific=T)," "," "))
						colnames(A2)[1] <- 'SS'
						rownames(A2) <- c('SSR','SSE','SST')
						R2HTML::HTML(as.matrix(A2),file=stdout(), innerBorder = 1,align="left")
					} else {
						warn.msg6 <- "\a Warning : ANOVA is not supported for intercept only model."
						R2HTML::HTML(warn.msg6,file=stdout())
					}
				} else {
					warn.msg7 <- "\a Warning : ANOVA is not supported for no intercept model."
					R2HTML::HTML(warn.msg7,file=stdout())
				}
			}

			# 다음의 값들은 엑셀 시트에 저장
			if(Predict) {
				temp.0 <- data.frame(Fitted_TSLS=fitted(res_TSLS))
				temp <- data.frame(Fitted_TSLS=rep(NA,nrow(dataset)))
				row.dataset <- rownames(dataset)
				row.output <- rownames(temp.0)
				for(i in row.dataset) if(i%in%row.output) temp[row.dataset==i,1] <- temp.0[row.output==i,1]
				O <- temp
			}

			if(Resid) {
				temp.0 <- data.frame(Resid_TSLS=resid(res_TSLS))
				temp <- data.frame(Resid_TSLS=rep(NA,nrow(dataset)))
				row.dataset <- rownames(dataset)
				row.output <- rownames(temp.0)
				for(i in row.dataset) if(i%in%row.output) temp[row.dataset==i,1] <- temp.0[row.output==i,1]
				if(exists('O')){
					O <- cbind(O,temp)
				} else {
					O <- temp
				}
			}

			if(stdResid) {
				X <- as.matrix(res_TSLS$x$regressor)
				is.singular <- class(try(solve(t(X)%*%X),silent=T))=="try-error"
				if(is.singular){
					warn.msg6 <- "\a Error : Cannnot calculate standardized residuals because of the (near) linear dependences in explatory variables. Please check multicollinearity problem."
					R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
					R2HTML::HTML(warn.msg6,file=stdout())
				} else {
					H <- X%*%solve(t(X)%*%X)%*%t(X)
					oo <- res_TSLS$resid/(res_TSLS$sigma*sqrt(1-diag(H)))
					temp.0 <- data.frame(stdResid_TSLS=oo)
					temp <- data.frame(stdResid_TSLS=rep(NA,nrow(dataset)))
					row.dataset <- rownames(dataset)
					row.output <- rownames(temp.0)
					for(i in row.dataset) if(i%in%row.output) temp[row.dataset==i,1] <- temp.0[row.output==i,1]
					if(exists('O')){
						O <- cbind(O,temp)
					} else {
						O <- temp
					}
				}
			}
		}

		# Analysis end time
		R2HTML::HTMLhr(file=stdout())
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Two-Stage Least Square Regression",sep=""),file=stdout())
		R2HTML::HTMLhr(file=stdout())
	})

	if(exists('O')){
		return(list(html=html.output,Output=O))
	} else {
		return(html.output)
	}
}

REx_OSPT <- function(exact=FALSE,dataset,res_var,event,x,n,p,alternative,conf.level,correct){
	#### 변수 설명 ####
	## dataset : 데이터셋이름
	## res_var : 2개의 카테고리로 구성된 변수. 필수로 입력되어야 하며 최대 1개 입력가능
	## event : dataset의 2개의 값 중 사건으로 정의하는 값 지정, dataset이 있을 때만 활성화
	## x : 사건횟수, dataset이 들어오면 직접 입력을 받지 않음
	## n : 시행횟수, dataset이 들어오면 직접 입력을 받지 않음
	## alternative : "two.sided", "less", "greater"
	## conf.level : confidence level, default:0.95
	## correct : exact=FALSE일 때만 활성화됨, TRUE, FALSE
	###################
	#### Required packages : R2HTML, testit ####
	###################
	load.pkg(c("R2HTML","testit"))

	global_op<-options()	#수정
	options("scipen"=0,"digits"=7)	#수정

	html.output <- capture.output({
		if(class(try(nrow(dataset),s=T))!='try-error'){
			dataset <- dataset[,res_var,drop=F]
			x<-length(which(dataset==event))
			n<-length(which(!is.na(dataset)))	
			data_desc<-c("Data type : Original data",
					paste("Dependent variable : ",colnames(dataset),sep=''),	#수정
					paste("Event : ",event,"-Success, ",unique(dataset[,1])[which(!unique(dataset[,1])%in% c(event,NA))],"-Failure",sep=''))	
		}

		if(class(try(nrow(dataset),s=T))=='try-error'){
			data_desc<-c("Data type : Summarized data",
					paste("Number of trials : ",n,sep=''),
					paste("Number of Successes : ",x,sep=''))
		}


		if(!exact){
			temp <- capture.output(prop.test(x=x,n=n,p=p,alternative=alternative,conf.level=conf.level,correct=correct))
			ftemp<-temp
			ftemp[2]<-paste('Method: Normal approximation',gsub('\t1-sample proportions test','',temp[2]),sep='')
			ftemp[5]<-strsplit(temp[5],',')[[1]][3]
			ftemp[7]<-paste(temp[9],gsub(' ','',temp[11]),sep=' ')
			ftemp[9]<-paste(temp[7],"(",strsplit(temp[8],' ')[[1]][2],",",strsplit(temp[8],' ')[[1]][3],")",sep=' ')
			warn<-'Warning : Approximation may be incorrect.'		#수정

				if(!testit::has_warning(prop.test(x=x,n=n,p=p,alternative=alternative,conf.level=conf.level,correct=correct))){
					res_OSPT<-gsub('<',': <',c(data_desc,ftemp[c(2,6,5,7,9)]))	
				}
				else{
				      res_OSPT<-c(gsub('<',': <',c(data_desc,ftemp[c(2,6,5,7,9)])),warn)	
				}

		} else { 
			temp <- capture.output(binom.test(x=x,n=n,p=p,alternative=alternative,conf.level=conf.level))
			ftemp<-temp
			ftemp[2]<-paste('Method: ',gsub('\t','',temp[2]),sep='')
			ftemp[6]<-gsub('probability of success','p',temp[6])
			ftemp[5]<-strsplit(temp[5],',')[[1]][3]
			ftemp[7]<-paste(temp[9],gsub(' ','',temp[11]),sep=' ')
			ftemp[9]<-paste(temp[7],"(",strsplit(temp[8],' ')[[1]][2],",",strsplit(temp[8],' ')[[1]][3],")",sep=' ')


			res_OSPT<-gsub('<',': <',c(data_desc,ftemp[c(2,6,5,7,9)]))	

		}

		f_res_OSPT<-matrix(unlist(lapply(res_OSPT,FUN=strsplit,split=': | = ')),ncol=2,byrow=T)

		R2HTML::HTML(R2HTML::as.title("One sample proportion test"),HR=1,file=stdout(),append=FALSE) 
		
		R2HTML::HTML(R2HTML::as.title("Analysis description"),HR=2,file=stdout())
		R2HTML::HTML(f_res_OSPT[1:4,],file=stdout(), innerBorder = 1,align="left")

		R2HTML::HTML(R2HTML::as.title("Result"),HR=2,file=stdout())
		R2HTML::HTML(f_res_OSPT[5:nrow(f_res_OSPT),],file=stdout(), innerBorder = 1,align="left")

		R2HTML::HTMLhr(file=stdout())
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : One sample proportion test",sep=""),file=stdout())
		R2HTML::HTMLhr(file=stdout())
	})
	options(global_op)	#수정
	return(html.output)
}

REx_LifeTable <- function(dataset, time, event, confi=0.95, survival_plot=FALSE, a1_survival_plot=FALSE, log_survival_plot=FALSE, CHazard_plot=FALSE, Density_plot=FALSE){
	#### 변수 설명 ####
	## dataset : 데이터셋 이름
	## time : 변수선택 탭의 '시간변수'
	## event : 변수선택 탭의 '상태변수'
	## survival_plot : 출력옵션 탭의 '생존함수' 옵션
	## a1_survival_plot : 출력옵션 탭의 '1-생존함수' 옵션
	## Hazard_plot : 출력옵션 탭의 '위험함수' 옵션 
	## Density_plot : 출력옵션 탭의 '밀도함수' 옵션 
	## log_survival_plot : 출력옵션 탭의 '로그생존함수' 옵션
	###################
	#### Required packages : R2HTML, fishmethods
	###################
	load.pkg(c("fishmethods", "survival"))
  
	html.output <- capture.output({
		### Title
		R2HTML::HTML(R2HTML::as.title('Life Table'),HR=1,file=stdout(),append=FALSE) 

		### Warnings
		if(class(dataset[,time])!='numeric') warn.msg1 <- 'Time variable should be numeric.'
		if(!all(unique(dataset[,event])%in%c(0,1))) warn.msg2 <- 'Invalid status values are observed (Expected value : 0 - censoring, 1 - event).'
		if(exists("warn.msg1")|exists("warn.msg2")) {
			R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
			if(exists("warn.msg1")) R2HTML::HTML(warn.msg1,file=stdout())
			if(exists("warn.msg2")) R2HTML::HTML(warn.msg2,file=stdout())

		} else {

		ex <- cbind(dataset[,time], dataset[,event]) # time and event column
		ex2 <- ex[which(ex[,2]==1),] # only event==1 data 
		TIME <- as.numeric(paste(names(table(ex2[,1])))) # time variable
		EVENT <- as.vector(table(ex2[,1])) # event variable 

		ev1 <- cbind(TIME, EVENT) # only event==1 data 
		ev0 <- ex[(which(ex[,2]==0)),] # only event==0 data 
		ev0_1 <- ev0[!duplicated(ev0),]
		ev0_2 <- ev0_1[-which(ev0_1[,1] %in% ev1[,1]),]

		final <- rbind(ev1, ev0_2) 
		final_sort <- final[order(final[,1]),] # sort by TIME 

		LT <- lifetable(final_sort[,1], final_sort[,2])		
		LT2 <- data.frame(cbind(LT$age, LT$dx, LT$nx, LT$lx, LT$qx, LT$px, LT$ex), row.names=NULL)
		names(LT2) <- c("time","dx","nx","lx","qx","px","ex")

		
		my.surv = Surv(dataset[,time], dataset[,event]) ;
		my.fit = survfit(formula = my.surv ~ 1, data=dataset, conf.int=confi) ;
		

		### Description 
		R2HTML::HTML(R2HTML::as.title("Data structure"),HR=2,file=stdout()) 
		R2HTML::HTML(matrix(c("Number of observations",paste(dim(ex)[1], " (Event:", paste0(sum(final_sort[,2]), ","), " Censoring:", paste0(dim(ex)[1]-sum(LT$dx),")",sep=""))),ncol=2), file=stdout(), innerBorder = 1,align="left")

		R2HTML::HTML(R2HTML::as.title("Variable List"),HR=2,file=stdout(),append=TRUE) 
		R2HTML::HTML(matrix(c("Time variable", "Censoring variable", time, event),ncol=2),HR=2,file=stdout(),append=TRUE,innerBorder = 1,align="left") 

		### Default output  
		R2HTML::HTML(R2HTML::as.title("Life Table"),HR=2,file=stdout()) 	
		R2HTML::HTML(LT2,file=stdout(), row.names=F, caption="<div style='text-align:left'> dx=number of event <br> nx=number of censoring <br> lx=(cumulative) survivorship <br> qx=finite rate of mortality <br> px=finite rate of survival <br> ex=mean expectation of life", innerBorder = 1,align="left")
		

		### Option_plot : Survival function
			if(survival_plot) {
				R2HTML::HTML(R2HTML::as.title("Survival function"),HR=2,file=stdout(),append=TRUE)
				REx_ANA_PLOT()
				plot(my.fit, xlab="Time", ylab="Survival probability")
				REx_ANA_PLOT_OFF(paste(paste0(confi*100,"%"), "confidence interval bounded"))		 	
			}

		### Option_plot : a1-Survival function
			if(a1_survival_plot) {
				R2HTML::HTML(R2HTML::as.title("1-Survival function"),HR=2,file=stdout(),append=TRUE)
				REx_ANA_PLOT()
				plot(my.fit, xlab="Time", ylab="1-Survival probability", fun="event", ylim=c(0,1))
				REx_ANA_PLOT_OFF(paste(paste0(confi*100,"%"), "confidence interval bounded"))		 	
			}

				
		### log_survival_plot : Log Survival function
			if(log_survival_plot) {
				R2HTML::HTML(R2HTML::as.title("Log Survival function"),HR=2,file=stdout(),append=TRUE)
				REx_ANA_PLOT()
				plot(my.fit, xlab="Time", ylab="Log survival probability", fun="log")
				REx_ANA_PLOT_OFF(paste(paste0(confi*100,"%"), "confidence interval bounded"))
			}

		### CHazard_plot : Cumulative Hazard function
			if(CHazard_plot) {
				R2HTML::HTML(R2HTML::as.title("Cumulative Hazard function"),HR=2,file=stdout(),append=TRUE)
				REx_ANA_PLOT()
				H.hat <- -log(my.fit$surv); H.hat <- c(H.hat, H.hat[length(H.hat)])
				H.hat[!is.finite(H.hat)] <- NA
				h.sort.of <- my.fit$n.event / my.fit$n.risk
				H.tilde <- vector()
				for(i in 1:length(h.sort.of)) H.tilde[i] <- sum(h.sort.of[1:i])
				H.tilde <- c(H.tilde, H.tilde[length(H.tilde)])
				plot(my.fit, xlab="Time", ylab="Cumulative hazard", fun="cumhaz", ylim=range(c(na.omit(H.hat), na.omit(H.tilde))))
				REx_ANA_PLOT_OFF(paste(paste0(confi*100,"%"), "confidence interval bounded"))
			}		

		### Density_plot : Density function 
			if(Density_plot) {
				R2HTML::HTML(R2HTML::as.title("Density function"),HR=2,file=stdout(),append=TRUE)
				REx_ANA_PLOT()
				plot(LT2$time, LT2$qx, xlab="Times", ylab="Density")
				REx_ANA_PLOT_OFF("")
			}

		### Analysis time 
		R2HTML::HTMLhr(file=stdout())
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Life Table",sep=""),file=stdout())
		R2HTML::HTMLhr(file=stdout())
		}
	})
	return(html.output)
}

REx_PoisReg <- function(dataset, res_var, quan_var=NULL, qual_var=NULL, vars=NULL, offset=NULL, link=c("log", "I", "sqrt"), over_disp=FALSE, 
			ss_type=c(2, "2", "II", 1, "1", "I", 3, "3", "III"), chisq_mode=c("likelihood", "wald", "both"), 
			CI_mode=c("likelihood", "wald", "both"), CI_level=0.95, fit_info=TRUE, summary_stat=TRUE, fitness=TRUE,
			estim=TRUE, exp_estim=FALSE, testing=TRUE, covar_mat=FALSE, corr_mat=FALSE,
			fitted_value=FALSE, resid=FALSE, resid_mode=NULL, linear_pred=FALSE, 
			hat_value=FALSE, cook_distance=FALSE, num_digit=3){
	#### Arguments ####
	## rep_var : GUI의 '반응변수', 필수로 1개 필요
	## quan_var : 양적변수(2개 이상의 변수가 있는 경우 +로 구분, 필수아님. 아무변수도 선택되지 않으면 1로 코딩, 질적변수와 중복되어 사용될 수 없음)
	## qual_var : 질적변수(2개 이상의 변수가 있는 경우 +로 구분, 필수아님. 아무변수도 선택되지 않으면 1로 코딩, 양적변수와 중복되어 사용될 수 없음)
	## offset: offset for the Poisson regression model /  GUI의 '오프셋'과 대응
	## link: link function, must be one of "log", "identity", and "sqrt" / GUI의 '연결함수'와 대응
	## over_disp: overdispersion, if TRUE, use quasi-Poisson model instead / GUI의 '과분산 적합'과 대응
	## ss_type: type of chi-squared statistics / GUI의 '통계량 유형'과 대응
	## chisq_mode, CI_mode: Chi-squared statistics & confidence intervals mode / GUI의 '카이제곱 통계량', '신뢰구간'에 대응
	## CI_level: confidence level of the confidence interval / GUI의 '신뢰수준'에 대응
	## fit_info: fitting information, formula, distribution, link function, etc... / GUI의 '모형 정보'와 대응
	## summary_stat: summary statistics of variables / GUI의 '기술통계량'과 대응
	## fitness: model fitness measures such as deviance / GUI의 '적합성 측도'와 대응
	## estim, exp_estim: parameter estimates & exponential mode / GUI의 '모수 추정값' 및 '지수 모수 추정값'에 대응
	## testing: testing whether each parameter is equal to 0 / GUI의 '가설 검정'에 대응
	## covar_mat: corr_mat: covariance(correlation) matrix of parameter estimates / GUI의 '추정값 공분산 행렬' 및 '추정값 상관행렬'에 대응
	## fiited_value: print/save fitted values if TRUE / GUI의 '적합값'에 대응
	## resid: print/save residuals if TRUE / GUI의 '잔차'에 대응
	## resid_mode: modes of residuals, can designate multiple ("original", "pearson", "pearson_std", "deviance", "deviance_std", "likelihood")
	## 	// GUI의 '원래 잔차', 'Pearson 잔차', '표준화된 Pearson 잔차', 'Deviance 잔차', '표준화된 Deviance 잔차', '가능도 잔차'에 대응, resid=TRUE이면 반드시 하나 이상 선택되어야 함.
	## linear_pred: print/save linear predictors for each observations / GUI의 '선형 예측값'에 대응
	## hat_value: diagonals of hat matrix / GUI의 'Hat 행렬 대각값'에 대응
	## cook_distance: Cook'd distance / GUI의 'Cook 거리'에 대응
	## num_digit: number of output digits
	###################

	###################
	#### Required packages : R2HTML, car, AICcmodavg ####
	###################
	load.pkg(c("R2HTML", "car", "AICcmodavg"))

	###################
	html.output <- capture.output({
		ss_type	<- as.character(ss_type) ;
		inter_var	<- vars[sapply(vars, function(s) grepl(":", s))] ;
		
		## Warnings
		# Response variable type
		if(!is.numeric(dataset[,res_var])) {
			warn.msg1 <- '\a Error : Response variable should be numeric. Analysis has been stopped.'
		} else if (any(dataset[,res_var] < 0)) {
			warn.msg1 <- "\a Error : Response variable should be non-negative numeric. Analysis has been stopped."
		}
		# explanatory variable type
		if(!is.null(quan_var)) {
			is.nom <- sapply(quan_var, function(i) !is.numeric(dataset[,i]))
			if(any(is.nom)){
				warn.msg2	<- paste0("\a Warning : The type of variable '", paste(quan_var[is.nom],collapse=', '), "' is not numeric but selected as the quantitative variable. It was excluded from the analysis.")
				ec		<- quan_var[is.nom]
				vars	<- vars[-grep(paste0('[',paste(ec,collapse=','),']'), vars)]
				if(length(vars)==0)	vars <- NULL
			}
		}
		if(!is.null(qual_var)) {
			is.num <- sapply(qual_var, function(i) is.numeric(dataset[,i]))
			if(any(is.num)) warn.msg3 <- paste0("\a Warning : The type of variable '", paste(qual_var[is.num],collapse=', '), "' is numeric but selected as the qualitative variable. It was coreced into character.")
			for(i in qual_var) dataset[,i] <- factor(dataset[,i])
		}

		link	<- match.arg(link) ;
		ss_type	<- match.arg(ss_type) ;
		chisq_mode	<- match.arg(chisq_mode) ;
		CI_mode	<- match.arg(CI_mode) ;

		form.0 <- ifelse(is.null(c(quan_var,qual_var)), paste0(res_var,' ~ 1'),  
			paste0(res_var, ' ~ ', paste(c(vars),collapse=' + '), 
					ifelse(is.null(offset), "", paste0(" + offset(", paste(offset, collapse=" + "), ")")))) ;
		form.1 <- as.formula(form.0)

		link_actual	<- ifelse(link == "I", "identity", link) ;
		dist	<- ifelse(over_disp, quasipoisson, poisson) ;

		fit_glm	<- glm(form.1, data=dataset, family=dist(link=link_actual)) ;

		R2HTML::HTML(R2HTML::as.title("Poisson Regression"), HR=1, file=stdout(), append=FALSE) ;
		
		# warnings
		if(exists('warn.msg1')|exists('warn.msg4')){
			R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
			if(exists('warn.msg1')) R2HTML::HTML(warn.msg1,file=stdout())
			if(exists('warn.msg4')){
				if(exists('warn.msg2')) R2HTML::HTML(warn.msg2,file=stdout())
				R2HTML::HTML(warn.msg4,file=stdout())
			}
		} else {
			variables	<- as.character(form.1) ;
			var_split	<- unlist(strsplit(variables[3], c("offset"))) ;
			
			var_list	<- unlist(strsplit(var_split[1], c(" \\* | \\+ |:"))) ;
			var_list	<- unique(var_list) ;
			if (length(var_split) > 1) {
				off_vars	<- gsub("\\(|\\)", "", var_split[2]) ;
				off_list	<- unlist(strsplit(off_vars, c(" \\* | \\+ |:"))) ;
				off_list	<- unique(off_list) ;
			}
			if(variables[3] == "1"){
				summ_var <- 1
				var_lists <- variables[2]
			} else {
				var_list <- var_list[var_list != "1"]
				var_lists <- c(variables[2], var_list)
				summ_var <- length(var_lists) + ifelse(length(var_split) > 1, length(off_list), 0) ;
			}
			data_str	<- matrix(c("Number of observations", "Number of variables", nobs(fit_glm), summ_var), ncol=2) ;

			R2HTML::HTML(R2HTML::as.title("Data Structure"), HR=2, file=stdout()) ;
			R2HTML::HTML(data_str, file=stdout(), row.names=FALSE, col.names=TRUE, innerBorder = 1, align="left") ;


			if (fit_info) {
				link_print	<- switch(link, 
							log = "Log",
							I = "Identity",
							sqrt = "Squareroot") ;

				fit_print	<- matrix(c(variables[2], 
							paste0(link, "(E(", variables[2], ")) = ", variables[3]),
							ifelse(over_disp, "Quasi-Poisson", "Poisson"),
							link_print)) ;
				fit_print	<- cbind(c("Response variable", "Model", "Distribution", "Link function"), fit_print) ;
				colnames(fit_print)	<- c("Category", "Value") ;
				
				R2HTML::HTML(R2HTML::as.title("Model Information"), HR=2, file=stdout()) ;
				R2HTML::HTML(fit_print, file=stdout(), row.names=FALSE, col.names=TRUE, innerBorder = 1, align="left") ;
				if(exists('warn.msg2')) R2HTML::HTML(warn.msg2, file=stdout(), row.names=FALSE)
				if(exists('warn.msg3')) R2HTML::HTML(warn.msg3, file=stdout(), row.names=FALSE)
			}

			if (summary_stat) {
				summ_print <- summary(dataset[,var_lists,drop=F]) ;
				summ_print <- data.frame(unclass(summ_print), stringsAsFactors=FALSE) ;
				colnames(summ_print)[1] <- variables[2] ; colnames(summ_print)[-1] <- var_list ;
				colnames(summ_print) <- gsub(" ", "", colnames(summ_print), fixed=TRUE) ;

				R2HTML::HTML(R2HTML::as.title("Summary Statistics"), HR=2, file=stdout()) ;
				R2HTML::HTML(summ_print, file=stdout(),row.names=FALSE,na="", innerBorder = 1, align="left") ;
			}

			if (fitness) {
				fitness_print	<- data.frame(numeric(0), numeric(0)) ;

				fitness_print	<- rbind(fitness_print, c(deviance(fit_glm), fit_glm$df.residual)) ;

				chisq_pearson	<- sum(residuals(fit_glm, type="pearson")^2) ;
				fitness_print	<- rbind(fitness_print, c(chisq_pearson, summary(fit_glm)$df.residual)) ;

				fitness_print	<- rbind(fitness_print, c(as.numeric(logLik(fit_glm)), attr(logLik(fit_glm), "df"))) ;
				
				fitness_print	<- rbind(fitness_print, c(AIC(fit_glm), NA)) ;
				fitness_print	<- rbind(fitness_print, c(AICc(fit_glm), NA)) ;
				fitness_print	<- rbind(fitness_print, c(BIC(fit_glm), NA)) ;
				names(fitness_print)	<- c("Value", "DF") ;
				row.names(fitness_print)	<- c("Deviance", "Pearson's chi-square", "log-likelihood", "AIC", "AICc", "BIC") ;
				fitness_print[1]	<- format(round(fitness_print[1], num_digit), nsmall=num_digit) ;

				R2HTML::HTML(R2HTML::as.title("Model Fitness Measurements"), HR=2, file=stdout()) ;
				R2HTML::HTML(fitness_print, file=stdout(), na="", innerBorder = 1, align="left") ;
			}

			if (estim) {
				coef_print	<- data.frame(coef(summary(fit_glm))) ;
				colnames(coef_print)	<- c("Estimate","Std.Error","Z-value","P-value") ;

				if (CI_mode %in% c("wald", "both")) {
					coef_print	<- cbind(coef_print, confint.default(fit_glm, level=CI_level)) ;
					names(coef_print)[ncol(coef_print)-1:0]	<- paste0(100*CI_level, "% CI ", c("lower", "upper"), " (Wald)") ;
				}
				if (CI_mode %in% c("likelihood", "both")) {
					lCI <- matrix(confint(fit_glm, level=CI_level),ncol=2)
					rownames(lCI) <- rownames(coef_print)
					coef_print	<- cbind(coef_print,lCI) ;
					names(coef_print)[ncol(coef_print)-1:0]	<- paste0(100*CI_level, "% CI ", c("lower", "upper"), " (Profile likelihood)") ;
				}
				if(!variables[3]=="1"){
					row.names(coef_print)	<- c("(Intercept)", colnames(model.matrix(fit_glm))[-1]) ;
				}
				R2HTML::HTML(R2HTML::as.title("Regression Coefficient Estimates and Confidence Intervals"), HR=2, file=stdout()) ;
				R2HTML::HTML(format(round(coef_print, num_digit), nsmall=num_digit), file=stdout(), na="", innerBorder = 1, align="left") ;

				if (estim & exp_estim) {
					coef_print	<- data.frame(coef(summary(fit_glm)));
					colnames(coef_print)	<- c("Estimate","Std.Error","Z-value","P-value") ;
					exp_Est <- exp(coef_print[,1])
					coef_exp_print <- data.frame(exp_Est,coef_print)

					if (CI_mode %in% c("wald", "both")) {
						coef_exp_print	<- cbind(coef_exp_print, exp(confint.default(fit_glm, level=CI_level))) ;
						names(coef_exp_print)[ncol(coef_exp_print)-1:0]	<- paste0(100*CI_level, "% CI ", c("lower", "upper"), " (Wald)") ;
					}
					if (CI_mode %in% c("likelihood", "both")) {
						lCI <- matrix(exp(confint(fit_glm, level=CI_level)),ncol=2)
						rownames(lCI) <- rownames(coef_exp_print)
						coef_exp_print	<- cbind(coef_exp_print,lCI) ;
						names(coef_exp_print)[ncol(coef_exp_print)-1:0]	<- paste0(100*CI_level, "% CI ", c("lower", "upper"), " (Profile likelihood)") ;
					}
					if(!variables[3]=="1"){
						row.names(coef_exp_print)	<- c("(Intercept)", colnames(model.matrix(fit_glm))[-1]) ;
					} 
					R2HTML::HTML(R2HTML::as.title("Exponential Regression Coefficient Estimates and Confidence Intervals"), HR=2, file=stdout()) ;
					R2HTML::HTML(format(round(coef_exp_print, num_digit), nsmall=num_digit), file=stdout(), na="", innerBorder = 1, align="left") ;
				}
			}
			
			nowald_flag	<- FALSE ;
			if (testing) {
				R2HTML::HTML(R2HTML::as.title("Testing Regression Coefficients"), HR=2, file=stdout()) ;

				if (variables[3]=="1"){
					R2HTML::HTML("Testing is not available when the model contains only an intercept.", file=stdout(), na="")
				} else {
					if (ss_type %in% c(2, "2", "II", 3, "3", "III")) {
						if (chisq_mode %in% c("wald", "both")) {
							table_tmp	<- data.frame(Anova(fit_glm, test="Wald", type=ss_type)) ;
							colnames(table_tmp)	<- c("DF", "Chi-squared (Wald)", "P-value (Wald)") ;
							
							if (ss_type %in% c(3, "3", "III")) {
								test_print	<- table_tmp[-1,] ;
							} else {
								test_print	<- table_tmp ;
							}
						}
						if (chisq_mode %in% c("likelihood", "both")) {
							table_tmp	<- data.frame(Anova(fit_glm, test="LR", type=ss_type)) ;
							colnames(table_tmp)	<- c("Chi-squared (LR)", "DF", "P-value (LR)") ;
							table_tmp	<- table_tmp[c(2,1,3)] ;
							
							if (chisq_mode == "likelihood") { 
								test_print	<- table_tmp ;
							} else {
								test_print	<- cbind(test_print, table_tmp[-1]) ;
							}
						}
						test_print[-1]	<- format(round(test_print[-1], num_digit), nsmall=num_digit) ;

					} else if (ss_type %in% c(1, "1", "I")) {
						if (chisq_mode %in% c("wald", "both")) {
							nowald_flag	<- TRUE ;
						}
						if (chisq_mode %in% c("likelihood", "both")) {
							table_tmp	<- data.frame(anova(fit_glm, test="LRT"))[-1,-c(3,4)] ;
							colnames(table_tmp)	<- c("DF", "Chi-squared (LR)", "P-value (LR)") ;
							test_print	<- table_tmp ;
							test_print[-1]	<- format(round(test_print[-1], num_digit), nsmall=num_digit) ;
						}
					}
				}
				if (exists("test_print")) R2HTML::HTML(test_print, file=stdout(), na="", innerBorder = 1, align="left") ;
				if (nowald_flag) R2HTML::HTML("Wald test for type I SS is not available.", file=stdout(), na="") ;
				
			}
			if (covar_mat) {
				covar_print	<- vcov(fit_glm) ;
				rownames(covar_print)[1]	<- colnames(covar_print)[1]	<- "Intercept" ;
				
				R2HTML::HTML(R2HTML::as.title("Covariance Matrix for Estimates of Regression Coefficients"), HR=2, file=stdout()) ;
				R2HTML::HTML(format(round(covar_print, num_digit), nsmall=num_digit), file=stdout(), na="", innerBorder = 1, align="left") ;
			}
			
			if (corr_mat) {
				corr_print	<- summary(fit_glm, correlation=TRUE)$correlation ;
				rownames(corr_print)[1]	<- colnames(corr_print)[1]	<- "Intercept" ;
				
				R2HTML::HTML(R2HTML::as.title("Correlation Matrix for Estimates of Regression Coefficients"), HR=2, file=stdout()) ;
				R2HTML::HTML(format(round(corr_print, num_digit), nsmall=num_digit), file=stdout(), na="", innerBorder = 1, align="left") ;
			}

			getOutput <- FALSE
			if (fitted_value) {
				PoisReg_output	<- data.frame(fitted(fit_glm)) ;
				names(PoisReg_output)	<- "PoisReg_Fitted" ;
				getOutput <- TRUE
			}
			
			if (resid) {
				if ("original" %in% resid_mode) {
					if (exists("PoisReg_output")) {
						PoisReg_output	<- cbind(PoisReg_output, fit_glm$y - fitted(fit_glm)) ;
						names(PoisReg_output)[ncol(PoisReg_output)]	<- "PoisReg_ResidOriginal" ;
					} else {
						PoisReg_output	<- data.frame(fit_glm$y - fitted(fit_glm)) ;
						names(PoisReg_output)	<- "PoisReg_ResidOriginal" ;
					}
				getOutput <- TRUE
				}
				
				if ("pearson" %in% resid_mode) {
					if (exists("PoisReg_output")) {
						PoisReg_output	<- cbind(PoisReg_output, residuals(fit_glm, type="pearson")) ;
						names(PoisReg_output)[ncol(PoisReg_output)]	<- "PoisReg_ResidPearson" ;
					} else {
						PoisReg_output	<- data.frame(residuals(fit_glm, type="pearson")) ;
						names(PoisReg_output)	<- "PoisReg_ResidPearson" ;
					}
				getOutput <- TRUE
				}
				
				if ("pearson_std" %in% resid_mode) {
					if (exists("PoisReg_output")) {
						PoisReg_output	<- cbind(PoisReg_output, rstandard(fit_glm, type="pearson")) ;
						names(PoisReg_output)[ncol(PoisReg_output)]	<- "PoisReg_ResidPearsonStd" ;
					} else {
						PoisReg_output	<- data.frame(rstandard(fit_glm, type="pearson")) ;
						names(PoisReg_output)	<- "PoisReg_ResidPearsonStd" ;
					}
				getOutput <- TRUE
				}
				
				if ("deviance" %in% resid_mode) {
					if (exists("PoisReg_output")) {
						PoisReg_output	<- cbind(PoisReg_output, residuals(fit_glm, type="deviance")) ;
						names(PoisReg_output)[ncol(PoisReg_output)]	<- "PoisReg_ResidDev" ;
					} else {
						PoisReg_output	<- data.frame(residuals(fit_glm, type="deviance")) ;
						names(PoisReg_output)	<- "PoisReg_ResidDev" ;
					}
				getOutput <- TRUE
				}
				
				if ("deviance_std" %in% resid_mode) {
					if (exists("PoisReg_output")) {
						PoisReg_output	<- cbind(PoisReg_output, rstandard(fit_glm, type="deviance")) ;
						names(PoisReg_output)[ncol(PoisReg_output)]	<- "PoisReg_ResidDevStd" ;
					} else {
						PoisReg_output	<- data.frame(rstandard(fit_glm, type="deviance")) ;
						names(PoisReg_output)	<- "PoisReg_ResidDevStd" ;
					}
				getOutput <- TRUE
				}
				
				if ("likelihood" %in% resid_mode) {
					if (exists("PoisReg_output")) {
						PoisReg_output	<- cbind(PoisReg_output, rstudent(fit_glm)) ;
						names(PoisReg_output)[ncol(PoisReg_output)]	<- "PoisReg_ResidLikeli" ;
					} else {
						PoisReg_output	<- data.frame(rstudent(fit_glm)) ;
						names(PoisReg_output)	<- "PoisReg_ResidLikeli" ;
					}
				getOutput <- TRUE
				}
			}
			
			if (linear_pred) {
				if (exists("PoisReg_output")) {
					PoisReg_output	<- cbind(PoisReg_output, fit_glm$linear.predictors) ;
					names(PoisReg_output)[ncol(PoisReg_output)]	<- "PoisReg_LinPred" ;
				} else {
					PoisReg_output	<- data.frame(fit_glm$linear.predictors) ;
					names(PoisReg_output)	<- "PoisReg_LinPred" ;
				}
				getOutput <- TRUE
			}
			
			if (hat_value) {
				if (exists("PoisReg_output")) {
					PoisReg_output	<- cbind(PoisReg_output, hatvalues(fit_glm)) ;
					names(PoisReg_output)[ncol(PoisReg_output)]	<- "PoisReg_HatValue" ;
				} else {
					PoisReg_output	<- data.frame(hatvalues(fit_glm)) ;
					names(PoisReg_output)	<- "PoisReg_HatValue" ;
				}
				getOutput <- TRUE
			}
			
			if (cook_distance) {
				if (exists("PoisReg_output")) {
					PoisReg_output	<- cbind(PoisReg_output, cooks.distance(fit_glm)) ;
					names(PoisReg_output)[ncol(PoisReg_output)]	<- "PoisReg_CookDist" ;
				} else {
					PoisReg_output	<- data.frame(cooks.distance(fit_glm)) ;
					names(PoisReg_output)	<- "PoisReg_CookDist" ;
				}
				getOutput <- TRUE
			}
		}

		R2HTML::HTMLhr(file=stdout()) ;
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Poisson Regression",sep=""),file=stdout()) ;
		R2HTML::HTMLhr(file=stdout()) ;
	})

	if(getOutput) {
		return(list(html=html.output,Output=PoisReg_output))
	} else {
		return(html.output)
	}
}

REx_TSPT <- function(dataset,dataid,res_var,res_var1,res_var2,event,event1,event2,x1,x2,n1,n2,pool=FALSE,alternative,conf.level){
	#### 변수 설명 ####
	## dataset : 데이터셋이름
	## dataid : 데이터 그룹을 구분할 수 있는 2개의 카테고리로 구성된 변수. 필수로 입력되어야 하며 최대 1개 입력가능
	## res_var : 2개의 카테고리로 구성된 변수. 필수로 입력되어야 하며 최대 1개 입력가능. Samples에 해당하는 변수 이름
	## res_var1 : 2개의 카테고리로 구성된 변수. 필수로 입력되어야 하며 최대 1개 입력가능. Sample 1에 해당하는 변수 이름
	## res_var2 : 2개의 카테고리로 구성된 변수. 필수로 입력되어야 하며 최대 1개 입력가능. Sample 2에 해당하는 변수 이름
	## event : dataset의 2개의 값 중 사건으로 정의하는 값 지정, dataset이 있을 때만 활성화
	## event1 : dataset의 2개의 값 중 사건으로 정의하는 값 지정, dataset이 있을 때만 활성화
	## event2 : dataset의 2개의 값 중 사건으로 정의하는 값 지정, dataset이 있을 때만 활성화
	## x1 : Sample 1의 사건횟수, dataset이 들어오면 직접 입력을 받지 않음
	## x2 : Sample 2의 사건횟수, dataset이 들어오면 직접 입력을 받지 않음
	## n1 : Sample 1의 시행횟수, dataset이 들어오면 직접 입력을 받지 않음
	## n2 : Sample 2의 시행횟수, dataset이 들어오면 직접 입력을 받지 않음
	## pool : TRUE, FALSE. 합동 추정치 사용 여부
	## alternative : "two.sided", "less", "greater"
	## conf.level : confidence level, default:0.95

	###################
	#### Required packages : R2HTML, testit ####
	###################
	load.pkg(c("R2HTML", "testit"))

	global_op<-options()	#수정
	options("scipen"=0,"digits"=7)	#수정

	html.output <- capture.output({
			if(class(try(nrow(dataset),s=T))!='try-error'){
				if(class(try(nrow(dataid),s=T))=='try-error'){
					dataset <- dataset[,c(res_var1,res_var2),drop=F]
					n1<-length(which(!is.na(dataset[,res_var1])))
					n2<-length(which(!is.na(dataset[,res_var2])))
					x1<-length(which(dataset[,res_var1]==event1))
					x2<-length(which(dataset[,res_var2]==event2))
	
					data_desc<-c("Data type : Original data (Samples in different columns)",
					paste(c("Sample 1 : ","Sample 2 : "),colnames(dataset),sep=''),
					paste("Event from Sample 1 : ",event1,"-Success, ",unique(dataset[,res_var1])[which(!unique(dataset[,res_var1])%in% c(NA,event1))],"-Failure",sep=''),
					paste("Event from Sample 2 : ",event2,"-Success, ",unique(dataset[,res_var2])[which(!unique(dataset[,res_var2])%in% c(NA,event2))],"-Failure",sep=''))

				    }
				if(class(try(length(dataid),s=T))!='try-error'){
					dataset <- dataset[,c(dataid,res_var),drop=F]
					n1<-length(which(dataset[,dataid]==names(table(dataset[,dataid]))[1] & !dataset[,res_var]%in% c(NA,event))) 
					n2<-length(which(dataset[,dataid]==names(table(dataset[,dataid]))[2] & !dataset[,res_var]%in% c(NA,event))) 
					x1<-length(which(dataset[,dataid]==names(table(dataset[,dataid]))[1] & dataset[,res_var]==event))
					x2<-length(which(dataset[,dataid]==names(table(dataset[,dataid]))[2] & dataset[,res_var]==event))

					data_desc<-c("Data type : Original data (Samples in one column)",
					paste(c("Sample 1 : ","Sample 2 : "),dataid,"=",names(table(dataset[,dataid])),sep=''),
					paste("Event: ",event,"-Success, ",unique(dataset[,res_var])[which(!unique(dataset[,res_var])%in% c(NA,event))],"-Failure",sep='')) 

				    }

			}

			if(class(try(nrow(dataset),s=T))=='try-error'){
				data_desc<-c("Data type : Summarized data",
					paste("Number of trials from Sample 1: ",n1,sep=''),
					paste("Number of Successes from Sample 1: ",x1,sep=''),
					paste("Number of trials from Sample 2: ",n2,sep=''),
					paste("Number of Successes from Sample 2: ",x2,sep=''))
		}

 		 numerator = (x1/n1) - (x2/n2)
		 p1<-x1/n1
  		 p2<-x2/n2
   		 p.common = (x1+x2) / (n1+n2)
		 denominator2 <- sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2)  
  		if(pool){
      		   denominator <- sqrt(p.common * (1-p.common) * (1/n1 + 1/n2))
			   method_desc<-"Equality of variance : True"
            }
  		if(!pool){
    			   denominator <- sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2)
			   method_desc<-"Equality of variance : False"
    		}

		z.prop.ris = numerator / denominator
		
		if(alternative=="two.sided"){
							hypoth<-"Difference of proportions between Sample 1 and Sample 2 is not equal to 0."
							qt<-qnorm(1-(1-conf.level)/2)
							ci<-paste("( ",format(numerator-qt*denominator2)," , ",format(numerator+qt*denominator2)," )",sep='')  #수정
							pval<-2*pnorm(abs(z.prop.ris),lower.tail=(alternative=='less'))}
		if(alternative=="less"){
							hypoth<-"Proportion of Sample 1 is less than that of Sample 2."
							qt<-qnorm(conf.level)
							ci<-paste("( ",format(-1.00000000)," , ",format(numerator+qt*denominator2)," )",sep='')	#수정  
							pval<-pnorm(z.prop.ris,lower.tail=(alternative=='less'))}
		if(alternative=="greater"){
							hypoth<-"Proportion of Sample 1 is greater than that of Sample 2."
							qt<-qnorm(conf.level)
							ci<-paste("( ",format(numerator-qt*denominator2)," , ",format(1.00000000)," )",sep='')	#수정  
							pval<-pnorm(z.prop.ris,lower.tail=(alternative=='less'))}
		result_desc<-c(paste("alternative hypothesis : ",hypoth,sep=''),
				   paste("p-value : ",format(pval),sep=''),	#수정
				   paste("Sample 1 estimates : ",format(p1),sep=''),	#수정
				   paste("Sample 2 estimates : ",format(p2),sep=''),	#수정
				   paste(conf.level*100," percent confidence interval	: ",format(ci),sep=''))	#수정
		warn<-c()
		if(x1<5|x2<5|n1-x1<5|n2-x2<5){warn<-'Warning : z test may be incorrect.'}	#수정

		res_TSPT<-c(data_desc,method_desc,result_desc,warn)

		f_res_TSPT<-matrix(unlist(lapply(res_TSPT,FUN=strsplit,split=': | = ')),ncol=2,byrow=T)

		R2HTML::HTML(R2HTML::as.title("Two sample proportion test"),HR=1,file=stdout(),append=FALSE) 
		R2HTML::HTML(R2HTML::as.title("Analysis description"),HR=2,file=stdout())
		R2HTML::HTML(f_res_TSPT[1:which(f_res_TSPT[,1]=='Equality of variance '),],file=stdout(), innerBorder = 1,align="left")

		R2HTML::HTML(R2HTML::as.title("Result"),HR=2,file=stdout())
		R2HTML::HTML(f_res_TSPT[(which(f_res_TSPT[,1]=='Equality of variance ')+1):nrow(f_res_TSPT),],file=stdout(), innerBorder = 1,align="left")

		R2HTML::HTMLhr(file=stdout())
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Two sample proportion test",sep=""),file=stdout())
		R2HTML::HTMLhr(file=stdout())
	})
	options(global_op)	#수정
	return(html.output)
}

REx_CA <- function(dataset,Row,Column,dim=2,CT=TRUE,RP=TRUE,CP=TRUE,RPR=TRUE,CPR=TRUE,BPlot=TRUE,CPlot=TRUE,RPlot=TRUE,RLPlot=TRUE,CLPlot=TRUE,html_name="Correspondence_analysis",CAprint=FALSE)
{
  ####Arguments ####
  ## Row: Selected row from the dataset / GUI의 '행'과 대응
  ## Column: Selected column from the dataset /  GUI의 '열'과 대응
  ## dim: choose as few dimensions as needed to explain most of the variation.over 2 and less then number of categories of variables / GUI의 '차원'과 대응
  ## CT:the cross-tabulation of the row and column variables/ GUI의 '대응일치표'과 대응 
  ## RP: for each row point displays the mass, scores in dimension, inertia, contribution of the point to the inertia of the dimension./GUI의 '행포인트'과 대응 
  ## CP: for each column point displays the mass, scores in dimension, inertia, contribution of the point to the inertia of the dimension./GUI의 '열포인트'과 대응 
  ## RPR: the cell contents divided by their corresponding row total/GUI의 '행프로파일'과 대응 
  ## CPR : the cell contents divided by their corresponding column total/GUI의 '열프로파일'과 대응 
  ## BPlot: representing the rows and the columns of the matrix as points or vectors in a joint display /GUI의 'biplot'과 대응 
  ## RPlot: representing the rows of the matrix as points or vectors in a joint display /GUI의 'rowplot'과 대응 
  ## CPlot: representing the columns of the matrix as points or vectors in a joint display /GUI의 'columnplot'과 대응 
  ## RLPlot: plot of the transformation of the row category values into scores in dimension, with one plot per dimension. /GUI의 '행선도표'과 대응 
  ## CLPlot: plot of the transformation of the column category values into scores in dimension, with one plot per dimension. /GUI의 '열선도표'과 대응 
  ## html_name: name of .html output file / GUI의 '결과 파일 이름'에 대응
  ###################
  #### Required packages : R2HTML, FactoMineR, devtools, factoextra, MASS ####
  load.pkg(c("FactoMineR", "devtools", "factoextra", "MASS"))

  html.output <- capture.output({ 
    selected_row<-dataset[,Row]
    selected_column<-dataset[,Column]
    
    ## Warnings
    if(dim<2) {
      warn.msg1 <- "\a Error : Dimension should be larger than 1. Analysis has stopped."
    } 
    
    if(length(unique(selected_column)) < 3) {
      warn.msg2 <- "\a Error :Column dimension should be larger than 2. Analysis has stopped."
    }
    if(length(unique(selected_row))< 3) {
      warn.msg3 <- "\a Error :Row dimension should be larger than 2. Analysis has stopped."
    }
    
    
    
    if(exists('warn.msg1')) {
      R2HTML::HTML(warn.msg1,file=stdout())
    }
    if(exists('warn.msg2')) {
      R2HTML::HTML(warn.msg2,file=stdout())
    }
    if(exists('warn.msg3')) {
      R2HTML::HTML(warn.msg3,file=stdout())
    }  
    
    else{
      R2HTML::HTML(R2HTML::as.title("Correspondence  Analysis"),HR=2,file=paste0(html_name, ".html"),append=FALSE) # 분석이름
      
      data <- data.frame(selected_row,selected_column)
      data<-data[complete.cases(data),]
      newdata <- as.data.frame.matrix(table(data)) #원래 데이터(original data)에서 뽑은 두 범주형 변수로 categorical table 만들기
      nrow<-nrow(newdata)
      ncol<-ncol(newdata)
      
      
      if(exists('warn.msg2')) {
        R2HTML::HTML(warn.msg2,file=stdout())
      }else {
        # 분석내용 1
        R2HTML::HTML(R2HTML::as.title("Data Structure"),HR=3,file=paste0(html_name, ".html")) 
        DS <- matrix(c('Number of observations','Number of variables','Dimensions',paste(':',nrow(dataset)),paste(':',2),paste(':',dim)),nrow=3)
        R2HTML::HTML(DS,file=paste0(html_name, ".html"))
        
        R2HTML::HTML(R2HTML::as.title("Variable List"),HR=3,file=paste0(html_name, ".html")) 
        varlist <- data.frame(Variables=c(Row, Column))
        rownames(varlist) <- c('Row variable','Column variable')
        R2HTML::HTML(varlist,file=paste0(html_name, ".html"))
        
        # 분석내용 2
        R2HTML::HTML(R2HTML::as.title("Statistics"),HR=3,file=paste0(html_name, ".html")) # 분석내용 2
        
        ##통계량
        ###해법 차원을 제거함##
        res.ca <- CA(newdata  , graph = FALSE)
        row.sum <- apply(newdata , 1, sum)
        col.sum <- apply(newdata , 2, sum)
        n <- sum(newdata)
        average.rp <- col.sum/n 
        average.cp <- row.sum/n 
        row.profile <- newdata  /row.sum
        col.profile <- newdata /col.sum
        grand.total <- sum(newdata)
        row.mass <- row.sum/grand.total
        d2.row <- apply(row.profile, 1,function(row.p, av.p){sum(((row.p - av.p)^2)/av.p)}, average.rp)
        col.mass <- col.sum/grand.total
        d2.col <- apply(col.profile, 2,function(col.p, av.p){sum(((col.p - av.p)^2)/av.p)}, average.cp)
        chisq <- chisq.test(newdata)
        
        #1. 대응일치표
        if(CT){
          R2HTML::HTML(R2HTML::as.title("Corresponding analysis Table"),HR=4,file=paste0(html_name, ".html")) # 차원
          R2HTML::HTML(newdata,file=paste0(html_name, ".html"))
        }
        
        #2.행포인트 개요(mass. score. inerita, contribution)
        row.sum <- apply(newdata , 1, sum)
        grand.total <- sum(newdata )
        row.mass <- row.sum/grand.total #MASS	
        row.inertia <- row.mass * d2.row #INERTIA
        residuals <- chisq$residuals/sqrt(n)
        nb.axes<-dim
        res.svd <- svd(residuals, nu = nb.axes, nv = nb.axes)
        sv <- res.svd$d[1:nb.axes]
        u <-res.svd$u #SCORE
        cc <- t(apply(u, 1, '*', sv)) #
        row.coord <- apply(cc, 2, '/', sqrt(row.mass))
        eig <- sv^2
        variance <- eig*100/sum(eig)
        cumvar <- cumsum(variance)
        eig<- data.frame(eig = eig, variance = variance,
                         cumvariance = cumvar)
        cc <- apply(row.coord^2, 2, "*", row.mass)
        row.contrib <- t(apply(cc, 1, "/", eig[1:nb.axes,1])) *100 #CONTRIBUTION
        
        row <- cbind.data.frame(mass = row.mass, score=u, inertia = row.inertia,contribution=row.contrib)
        
        if(RP){
          R2HTML::HTML(R2HTML::as.title("Row point"),HR=4,file=paste0(html_name, ".html")) # 행 포인트
          R2HTML::HTML(row ,file=paste0(html_name, ".html"))
        } 
        
        #3.열포인트 개요
        col.sum <- apply(newdata , 2, sum)
        col.mass <- col.sum/grand.total
        col.inertia <- col.mass * d2.col
        v <- res.svd$v
        cc <- t(apply(v, 1, '*', sv)) # each row X sv
        col.coord <- apply(cc, 2, '/', sqrt(col.mass))
        cc <- apply(col.coord^2, 2, "*", col.mass)
        col.contrib <- t(apply(cc, 1, "/", eig[1:nb.axes,1])) *100
        
        col<- cbind.data.frame( mass = col.mass, score=v, inertia = col.inertia,contribution=col.contrib)
        
        if(CP){
          R2HTML::HTML(R2HTML::as.title("Column point"),HR=4,file=paste0(html_name, ".html")) # 열포인트
          R2HTML::HTML(col,file=paste0(html_name, ".html"))
        }
        
        
        ###4.대응일치표 순열 제거###
        #5.행프로파일 개요
        if(RPR){
          R2HTML::HTML(R2HTML::as.title("Row profile"),HR=4,file=paste0(html_name, ".html")) # 행프로파일
          R2HTML::HTML(row.profile,file=paste0(html_name, ".html"))
        }
        
        #6.열프로파일 개요
        if(CPR){
          R2HTML::HTML(R2HTML::as.title("Column profile"),HR=4,file=paste0(html_name, ".html")) # 열프로파일
          R2HTML::HTML(col.profile,file=paste0(html_name, ".html"))
        }
        
        #분석내용 3
        R2HTML::HTML(R2HTML::as.title("Graphs"),HR=3,file=paste0(html_name, ".html")) # 분석내용 3
        
        #그래프 저장해서 출력하는것
        ###산점도####
        R2HTML::HTML(R2HTML::as.title("Scatter plot"),HR=4,file=paste0(html_name, ".html")) # 산점도
        
        #1.bi-plot
        if (BPlot){
          R2HTML::HTML(R2HTML::as.title("BI-PLOT"),HR=5,file=paste0(html_name, ".html")) # BI-PLOT
          graph1 <- "BPplot.png"
          png(graph1,res = 80)
          plot(fviz_ca_biplot(res.ca))
          dev.off()
          R2HTML::HTMLInsertGraph(graph1,file=paste0(html_name, ".html"))
        }
        
        #2.행
        if (RPlot){
          R2HTML::HTML(R2HTML::as.title("Row"),HR=5,file=paste0(html_name, ".html")) # BI-PLOT
          graph1 <- "RPplot.png"
          png(graph1,res = 80)
          plot(fviz_ca_row(res.ca))
          dev.off()
          R2HTML::HTMLInsertGraph(graph1,file=paste0(html_name, ".html"))
        }
        
        #3.열
        if (CPlot){
          R2HTML::HTML(R2HTML::as.title("Column"),HR=5,file=paste0(html_name, ".html")) # BI-PLOT
          graph1 <- "CPplot.png"
          png(graph1,res = 80)
          plot(fviz_ca_col(res.ca))
          dev.off()
          R2HTML::HTMLInsertGraph(graph1,file=paste0(html_name, ".html"))
        }
        
        ##선도표
        R2HTML::HTML(R2HTML::as.title("Line graph"),HR=4,file=paste0(html_name, ".html")) # 선도표
        #1.행
        if (RLPlot){
          R2HTML::HTML(R2HTML::as.title("Row line graph"),HR=5,file=paste0(html_name, ".html")) # BI-PLOT
          for(i in seq(dim)){
            graph1=paste0("graph",i,".png")
            png(graph1, res = 80)
            plot(u[,i],type="l",ylab="Categories",xlab=Row,xaxt="n")
            axis(1,at=seq(1:nrow),labels=rownames(newdata))
            mtext(paste("Dimension", i))
            dev.off()
            
            R2HTML::HTMLInsertGraph(graph1,file=paste0(html_name, ".html"))
          }
        }
        
        #2.열
        if (CLPlot){
          R2HTML::HTML(R2HTML::as.title("Column line graph"),HR=5,file=paste0(html_name, ".html")) # BI-PLOT
          for(i in seq(dim)){
            graph1=paste0("graph",i,".png")
            png(graph1,res = 80)
            plot(v[,i],type="l",ylab="Categories",xlab=Row,xaxt="n")
            axis(1,at=seq(1:nrow),labels=rownames(newdata ))
            mtext(paste("Dimension", i))
            dev.off()
            
            R2HTML::HTMLInsertGraph(graph1,file=paste0(html_name, ".html"))
          }
        }
        R2HTML::HTMLhr(file=stdout())
        R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Correspondence Analysis",sep=""),file=stdout(),innerBorder = 1,align="left")
        R2HTML::HTMLhr(file=stdout())
      }
    }
  }
  
  )
  
  
  if(CAprint){
    return(list(html=html.output,Output=cm))
  } else {
    return(html.output)
  }
}

REx_KMsurv <- function(dataset, time, event, group="", confi=0.95, survival_plot=FALSE, a1_survival_plot=FALSE, CHazard_plot=FALSE, log_survival_plot=FALSE){
	#### 변수 설명 ####
	## time : 시간변수
	## event : 상태변수 (event=1, censoring=0)
	## group : 집단변수 (선택되지 않으면 "")
	## survival_plot : 출력옵션 탭의 '생존함수' 옵션
	## a1_survival_plot : 출력옵션 탭의 '1-생존함수' 옵션
	## CHazard_plot : 출력옵션 탭의 '누적위험함수' 옵션 
	## log_survival_plot : 출력옵션 탭의 '로그생존함수' 옵션
	###################
	#### Required packages : R2HTML, survival, KMsurv, markdown
	###################
	load.pkg(c("R2HTML", "KMsurv", "survival"))

	html.output <- capture.output({
		### Title
		R2HTML::HTML(R2HTML::as.title('Kaplan-Meier Analysis'),HR=1,file=stdout(),append=FALSE) 

		### Warnings
		if(class(dataset[,time])!='numeric') warn.msg1 <- 'Time variable should be numeric.'
		if(!all(unique(dataset[,event])%in%c(0,1))) warn.msg2 <- 'Invalid status values are observed (Expected value : 0 - censoring, 1 - event).'
		if(exists("warn.msg1")|exists("warn.msg2")) {
			R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
			if(exists("warn.msg1")) R2HTML::HTML(warn.msg1,file=stdout())
			if(exists("warn.msg2")) R2HTML::HTML(warn.msg2,file=stdout())
		} else {

			if(group=="") { 
				#### when group variable dosen't exist

				my.surv = Surv(dataset[,time], dataset[,event]) ;
				my.fit = survfit(formula = my.surv ~ 1, data=dataset, conf.int=confi) ;
				a <- summary(my.fit); 
				b <- as.data.frame(with(a,cbind(time,n.risk,n.event,surv,std.err,lower,upper)))
				colnames(b)[6:7] <- paste(confi*100,'% ',c('lower CI','upper CI'),sep="")

				### Description 
				R2HTML::HTML(R2HTML::as.title("Data structure"),HR=2,file=stdout(),append=TRUE) 
				R2HTML::HTML(matrix(c("Number of observations", paste(a$n.risk[1], "(Event:", paste0(sum(a$n.event), ","), "Censoring:", paste0(a$n.risk[1]-sum(a$n.event),")"))),ncol=2), HR=2,file=stdout(),append=TRUE,innerBorder = 1,align="left") 
				  
				R2HTML::HTML(R2HTML::as.title("Variable List"),HR=2,file=stdout(),append=TRUE) 
				R2HTML::HTML(matrix(c("Time variable", "Censoring variable", time, event),ncol=2),HR=2,file=stdout(),append=TRUE,innerBorder = 1,align="left") 

				### Default output  
				R2HTML::HTML(R2HTML::as.title("Kaplan-Meier estimator"),HR=2,file=stdout(),append=TRUE) 	
				R2HTML::HTML(round(b,4), file=stdout(), innerBorder = 1, align="left",row.names=F,digits=4,
				caption="<div style='text-align:left'> time = time point <br> n.risk = number of subjects at risk <br> n.event = number of events  <br> surv = survival probability <br> std.err = standard error of survival probability <br> lower CI = lower bound of confidence interval of survival probability <br> upper CI = upper bound of confidence interval of survival probability")
			
				### Option_plot : Survival function
				if(survival_plot) {
					R2HTML::HTML(R2HTML::as.title("Survival function"),HR=2,file=stdout(),append=TRUE)
					REx_ANA_PLOT()
					plot(my.fit, xlab="Time", ylab="Survival probability")
					REx_ANA_PLOT_OFF(paste(paste0(confi*100,"%"), "confidence interval bounded"))		 	
				}

				### Option_plot : a1-Survival function
				if(a1_survival_plot) {
					R2HTML::HTML(R2HTML::as.title("1-Survival function"),HR=2,file=stdout(),append=TRUE)
					REx_ANA_PLOT()
					plot(my.fit, xlab="Time", ylab="1-Survival probability", fun="event", ylim=c(0,1))
					REx_ANA_PLOT_OFF(paste(paste0(confi*100,"%"), "confidence interval bounded"))		 	
				}

				
				### CHazard_plot : Cumulative Hazard function
				if(CHazard_plot) {
					R2HTML::HTML(R2HTML::as.title("Cumulative Hazard function"),HR=2,file=stdout(),append=TRUE)
					REx_ANA_PLOT()
					H.hat <- -log(my.fit$surv); H.hat <- c(H.hat, H.hat[length(H.hat)])
					H.hat[!is.finite(H.hat)] <- NA
					h.sort.of <- my.fit$n.event / my.fit$n.risk
					H.tilde <- vector()
					for(i in 1:length(h.sort.of)) H.tilde[i] <- sum(h.sort.of[1:i])
					H.tilde <- c(H.tilde, H.tilde[length(H.tilde)])
					plot(my.fit, xlab="Time", ylab="Cumulative hazard", fun="cumhaz", ylim=range(c(na.omit(H.hat), na.omit(H.tilde))))
					REx_ANA_PLOT_OFF(paste(paste0(confi*100,"%"), "confidence interval bounded"))
				}

				### log_survival_plot : Log Survival function
				if(log_survival_plot) {
					R2HTML::HTML(R2HTML::as.title("Log Survival function"),HR=2,file=stdout(),append=TRUE)
					REx_ANA_PLOT()
					plot(my.fit, xlab="Time", ylab="Log survival probability", fun="log")
					REx_ANA_PLOT_OFF(paste(paste0(confi*100,"%"), "confidence interval bounded"))
				}
			} else {
				#### when group variable exists
				survobj = Surv(time = dataset[,time], event = dataset[,event])
				fit = survdiff(survobj ~ as.factor(dataset[,group]), rho=0) # log rank test
				p.val <- round(1 - pchisq(fit$chisq, length(fit$n) - 1),3)
				chisq <- round(fit$chisq,3)
				fit2 = survfit(formula(paste0('survobj ~ as.factor(',group,')')), data=dataset, conf.int=confi) ;
				fit3 = summary(fit2)	

				b <- as.data.frame(with(fit, cbind(n, obs, exp)))
				rownames(b) <- gsub(".+=","",rownames(b))
				colnames(b) <- c('n','observed events','expected events')
				grp <- gsub(".+=","",as.character(summary(fit2)$strata))
				c <- with(summary(fit2),data.frame(time, n.risk, n.event, surv, std.err, lower, upper))
				colnames(c)[6:7] <- paste(confi*100,'% ',c('lower CI','upper CI'),sep="")
				n <- length(unlist(fit[1]))

				### Description 
				R2HTML::HTML(R2HTML::as.title("Data structure"),HR=2,file=stdout(),append=TRUE) 
				R2HTML::HTML(matrix(c("Number of observations", paste(sum(unlist(fit[1])), "(Event:", paste0(sum(unlist(fit[2])), ","), "Censoring:", paste0(sum(unlist(fit[1]))-sum(unlist(fit[2])),")")),
					       "Number of groups", length(unlist(fit[1]))),ncol=2,byrow=T), HR=2,file=stdout(),append=TRUE,innerBorder = 1,align="left") 	
						
				R2HTML::HTML(R2HTML::as.title("Variable List"),HR=2,file=stdout(),append=TRUE) 
				R2HTML::HTML(matrix(c("Time variable", "Censoring variable", "Factor variable", time, event, group), ncol=2),HR=2,file=stdout(),append=TRUE,innerBorder = 1,align="left")


				### Default output  
				R2HTML::HTML(R2HTML::as.title("Overall comparison between groups"),HR=2,file=stdout(),append=TRUE) 	
				R2HTML::HTML(R2HTML::as.title("Pearson's chi-squared test"),HR=3,file=stdout(),append=TRUE) 	
				R2HTML::HTML(b,file=stdout(), innerBorder = 1, align="left",digits=4,
				caption=paste0('<div style="text-align:left"> Chisq=',chisq,' on ',nrow(b)-1,' degrees of fredom, <br>p-value=',p.val))

				R2HTML::HTML(R2HTML::as.title("Kaplan-Meier estimator"),HR=2,file=stdout(),append=TRUE) 
				for(g in unique(grp)){
					R2HTML::HTML(R2HTML::as.title(paste(group,g,sep=' - ')),HR=3,file=stdout(),append=TRUE)
					cc <- c[grp==g,]
					R2HTML::HTML(cc, file=stdout(), innerBorder = 1, align="left",digits=4,row.names=F,
					caption="<div style='text-align:left'>time = time point <br> n.risk = number of subjects at risk <br> n.event = number of events  <br> surv = survival probability <br> std.err = standard error of survival probability <br> lower CI = lower bound of confidence interval of survival probability <br> upper CI = upper bound of confidence interval of survival probability")
				}
				
				### Option_plot : Survival function
				if(survival_plot) {
					R2HTML::HTML(R2HTML::as.title("Survival function"),HR=2,file=stdout(),append=TRUE)
					REx_ANA_PLOT()
					plot(fit2, xlab="Time", ylab="Survival probability", col=rainbow(n))
					legend("topright", col=rainbow(n), lty=c(rep(1, length(rainbow(n)))), c(names(fit2$strata)) )
					REx_ANA_PLOT_OFF("")
				}

				### Option_plot : a1-Survival function
				if(survival_plot) {
					R2HTML::HTML(R2HTML::as.title("1-Survival function"),HR=2,file=stdout(),append=TRUE)
					REx_ANA_PLOT()
					plot(fit2, xlab="Time", ylab="1-Survival probability", fun="event", col=rainbow(n))
					legend("bottomright", col=rainbow(n), lty=c(rep(1, length(rainbow(n)))), c(names(fit2$strata)) )
					REx_ANA_PLOT_OFF("")	
				}

				### CHazard_plot : Cumulative Hazard function
				if(CHazard_plot) {
					R2HTML::HTML(R2HTML::as.title("Cumulative Hazard function"),HR=2,file=stdout(),append=TRUE)
					REx_ANA_PLOT()
					#H.hat <- -log(fit2$surv); H.hat <- c(H.hat, H.hat[length(H.hat)])
					#H.hat[!is.finite(H.hat)] <- NA
					#h.sort.of <- fit2$n.event / fit2$n.risk
					#H.tilde <- vector()
					#for(i in 1:length(h.sort.of)) H.tilde[i] <- sum(h.sort.of[1:i])
					#H.tilde <- c(H.tilde, H.tilde[length(H.tilde)])
					plot(fit2, xlab="Time", ylab="Cumulative hazard", fun="cumhaz", col=rainbow(n))
					legend("bottomright", col=rainbow(n), lty=c(rep(1, length(rainbow(n)))), c(names(fit2$strata)) )
					REx_ANA_PLOT_OFF("")
				}


				### log_survival_plot : Log Survival function

				if(log_survival_plot) {
					R2HTML::HTML(R2HTML::as.title("Log Survival function"),HR=2,file=stdout(),append=TRUE)
					REx_ANA_PLOT()
					plot(fit2, xlab="Time", ylab="Log survival probability", fun="log", col=rainbow(n))
					legend("topright", col=rainbow(n), lty=c(rep(1, length(rainbow(n)))), c(names(fit2$strata)) )
					REx_ANA_PLOT_OFF("")
				}
			}
		}
		### TIME 
		R2HTML::HTMLhr(file=stdout())
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Kaplan-Meier Analysis",sep=""),file=stdout())
		R2HTML::HTMLhr(file=stdout())		
	})
	return(html.output)
}

### 결측치 대체

REx_MIMP <- function(dataset,quan_var=NULL,unorder.qual_var=NULL,order.qual_var=NULL,level=NULL,predictor=NULL,method.quan=NULL,method.unorderqual=NULL,method.orderqual=NULL,method.any=NULL){
	#### 변수 설명 ####
	## quan_var : 변수선택 탭의 양적변수(2개 이상의 변수가 있는 경우 c('','')와 같이 입력 (ex. c('a','b','c')), 필수아님. 다른변수와 중복되어 사용될 수 없음)
	## unorder.qual_var : 변수선택 탭의 명목형 질적변수(2개 이상의 변수가 있는 경우 c('','')와 같이 입력 (ex. c('a','b','c')), 필수아님. 다른변수와 중복되어 사용될 수 없음)
	## order.qual_var : 변수선택 탭의 순서형 질적변수(최대 1개 입력가능, 필수아님. 다른변수와 중복되어 사용될 수 없음)
	## level : 변수선택 탭의 
	## 양적변수, 명목형 질적변수, 순서형 질적변수 중 하나는 반드시 변수를 입력받아야 함.
	## predictor : 변수선택 탭의 예측변수 (2개 이상의 변수가 있는 경우 c('','')와 같이 입력 (ex. c('a','b','c')), 필수아님. 다른변수와 중복되어 사용될 수 없음)
	## method.quan : 분석옵션 탭의 양적변수에서 선택된 방법. 2개 이상의 방법이 있는경우 c('','')와 같이 입력 (ex. c('a','b','c'))
	## method.unorderqual : 분석옵션 탭의 명목형 질적변수에서 선택된 방법. 2개 이상의 방법이 있는경우 c('','')와 같이 입력 (ex. c('a','b','c'))
	## method.orderqual : 분석옵션 탭의 명목형 양적변수에서 선택된 방법. 2개 이상의 방법이 있는경우 c('','')와 같이 입력 (ex. c('a','b','c'))
	## method.any : 분석옵션 탭의 모든 변수타입에서 선택된 방법. 2개 이상의 방법이 있는경우 c('','')와 같이 입력 (ex. c('a','b','c'))
	###################
	#### Required packages : mice, VIM, R2HTML, randomForest ####
	###################
	load.pkg(c("mice", "R2HTML", "VIM"))

	html.output <- capture.output({
		### HTML Output
		## Title
		HTML(as.title("Missing Data Imputation"),HR=1,file=stdout(),append=FALSE)

		## Warnings
		if(length(order.qual_var)!=0){
			exp.level <- as.character(unique(dataset[,order.qual_var]))
			exp.level <- exp.level[!is.na(exp.level)]
			if(!all(exp.level%in%level)) warn.msg1 <- paste0("Levels input by the user differ from the observed levels in the data. (input:",paste(level,collapse=','),", observed:",paste(exp.level,collapse=','),")")
		}
		
		if(exists("warn.msg1")){
			HTML(as.title("Warnings"),HR=2,file=stdout())
			HTML(warn.msg1,file=stdout())
		} else {
			## Data Structure
			HTML(as.title("Data structure"),HR=2,file=stdout())
			DS <- matrix(c('Number of observations',nrow(dataset),'Number of variables',length(c(quan_var,unorder.qual_var,order.qual_var,predictor))),ncol=2,byrow=T)
			HTML(DS,file=stdout(), innerBorder = 1,align="left")

			## Variable list
			HTML(as.title("Variable list"),HR=2,file=stdout())
			VL <- matrix(0,nrow=1,ncol=2)
			if(length(quan_var)!=0) VL <- rbind(VL,c('Numeric variable',paste(quan_var,collapse=', ')))
			if(length(unorder.qual_var)!=0) VL <- rbind(VL,c('Unordered categorical variable',paste(unorder.qual_var,collapse=', ')))
			if(length(order.qual_var)!=0) VL <- rbind(VL,c('Ordered categorical variable',paste(order.qual_var,collapse=', ')))
			if(length(predictor)!=0) VL <- rbind(VL,c('Predictor',paste(predictor,collapse=', ')))
			VL <- VL[-1,,drop=F]
			HTML(VL,file=stdout(), innerBorder = 1,align="left")

			## Analysis Description
			HTML(as.title("Analysis description"),HR=2,file=stdout())
			HTML(as.title("Selected Method"),HR=3,file=stdout())
			methods <- data.frame(abbr=c('pmm','norm','norm.nob','norm.boot','norm.predict','mean','2l.norm','2l.pan','2lonly.mean','2lonly.norm','2lonly.pmm','quadratic','logreg','polr','lda','cart','rf','ri','sample')
			,full=c('Predictive mean matching (pmm)','Bayesian linear regression (norm)','Linear regression ignoring model error (norm.nob)','Linear regression using bootstrap (norm.boot)','Linear regression, predicted values (norm.predict)','Unconditional mean imputation (mean)','Two-level normal imputation (2l.norm)','Two-level normal imputation using pan (2l.pan)','Imputation at level-2 of the class mean (2lonly.mean)',
			'Imputation at level-2 by Bayesian linear regression (2lonly.norm)','Imputation at level-2 by Predictive mean matching (2lonly.pmm)','Imputation of quadratic term (quadratic)',
			'Logistic regression (logreg)','Proportional odds model (polr)','Linear discriminant analysis (lda)','Classification and regression trees (cart)','Random forest imputation (rf)',
			'Random indicator method for nonignorable data (ri)','Random sample from the observed values (sample)'))
			AD <- matrix(0,nrow=1,ncol=2)
			if(length(method.quan)!=0) AD <- rbind(AD,c('Numeric variable',paste(methods$full[methods$abbr%in%method.quan],collapse=',<br>')))
			if(length(method.unorderqual)!=0) AD <- rbind(AD,c('Unordered categorical variable',paste(unique(methods$full[methods$abbr%in%method.unorderqual]),collapse=',<br>')))
			if(length(method.orderqual)!=0) AD <- rbind(AD,c('Ordered categorical variable',paste(unique(methods$full[methods$abbr%in%method.orderqual]),collapse=',<br>')))
			if(length(method.any)!=0) AD <- rbind(AD,c('Any type',paste(methods$full[methods$abbr%in%method.any],collapse=',<br>')))
			AD <- as.data.frame(AD[-1,,drop=F])
			colnames(AD) <- c('Variable type','Method')
			HTML(AD,file=stdout(), innerBorder = 1,align="left",row.names=F)

			## Inspecting the missing data
			HTML(as.title("Inspecting the missing data"),HR=2,file=stdout())
			newdat <- dataset[,c(quan_var,unorder.qual_var,order.qual_var),drop=F]
			if(ncol(newdat)==1){
				newdat1 <- newdat
				newdat1$obs <- 1:nrow(newdat)
				MP <- md.pattern(newdat1)
				MP <- MP[,colnames(MP)!='obs',drop=F]
				
			} else {
				MP <- md.pattern(newdat)
			}
			prop <- paste0(round(MP[nrow(MP),1:(ncol(MP)-1)]/nrow(dataset),3)*100,'%')
			MP1 <- MP[1:(nrow(MP)-1),1:(ncol(MP)-1),drop=F]
			MP1[MP1==1] <- 'Observed'; MP1[MP1==0] <- 'Missing'
			MP[1:(nrow(MP)-1),1:(ncol(MP)-1)] <- MP1
			MP[1:(nrow(MP)-1),ncol(MP)] <- gsub(' ','',rownames(MP)[1:(nrow(MP)-1)])
			colnames(MP)[ncol(MP)] <- '# of observations'
			rownames(MP) <- c(rep('',(nrow(MP)-1)),'# of missing')
			form <- ifelse(ncol(newdat)==1,MP[nrow(MP),ncol(MP)],paste0(paste(MP[nrow(MP),1:(ncol(MP)-1)],collapse='+'),'=',MP[nrow(MP),ncol(MP)]))
			MP[nrow(MP),1:(ncol(MP)-1)] <- paste0(MP[nrow(MP),1:(ncol(MP)-1)],' (',prop,')')
			HTML(MP,file=stdout(), innerBorder = 1,align="left",caption=paste0('<div style="text-align:left"> The total number of missing values is ',form,'.'))

			REx_ANA_PLOT(w=900,h=400)
			a <- capture.output(aggr(newdat,col=c('navyblue','yellow'),numbers=TRUE, sortVars=TRUE,labels=names(newdat), gap=3, ylab=c("Missing data","Pattern")))
			REx_ANA_PLOT_OFF('Yellow : Observed, Navy : Missing')
			
			## Imputation
			HTML(as.title("Imputation results"),HR=2,file=stdout())
			HTML("All results are stored in Excel sheet with following column name.",file=stdout())
			for(i in predictor) if(class(dataset[,i])=="character") dataset[,i] <- factor(dataset[,i])
			IR <- matrix(0,nrow=1,ncol=2)
			O <- data.frame(obs=1:nrow(dataset))
			
			# method : numeric variable
			if(length(method.quan)!=0){
				getIMP <- function(mt){
					mt1 <- mt
					mt <- rep('polyreg',ncol(newdat))
					pmm <- sapply(1:ncol(newdat),function(i) any(class(newdat[,i])=="numeric"))
					mt[pmm] <- "pmm"
					imp <- mice(newdat,pred=pred,meth=mt,m=1,maxit=10,print=F)
					completeData <- complete(imp)[,quan_var,drop=F]
					colnames(completeData) <- paste0(colnames(completeData),'_',mt1)
					return(completeData)
				}
				if(length(predictor)!=0){
					newdat <- dataset[,c(predictor,quan_var),drop=F]
					ini <- mice(newdat,maxit=0)
					pred <- ini$predictorMatrix
					pred[predictor,] <- pred[quan_var,quan_var] <- 0
				} else {
					newdat <- dataset[,quan_var,drop=F]
					newdat <- data.frame(obs=1:nrow(newdat),newdat)
					ini <- mice(newdat,maxit=0)
					pred <- ini$predictorMatrix
					pred[quan_var,] <- 0
				}
				resIMP <- lapply(method.quan,getIMP)
				resIMP.quan <- do.call(cbind,resIMP)
				O <- data.frame(O,resIMP.quan)
				IR <- rbind(IR,c('Numeric variable',paste(colnames(resIMP.quan),collapse=', ')))
			}

			# method : unordered categorical variable
			if(length(method.unorderqual)!=0){
				for(i in unorder.qual_var) dataset[,i] <- factor(dataset[,i])
				getIMP <- function(mt){
					mt1 <- mt
					mt <- rep(mt,ncol(newdat))
					if(mt1=='logreg') mt <- apply(newdat,2,function(i) ifelse(length(unique(i))==2,'logreg','polyreg'))
					pmm <- sapply(1:ncol(newdat),function(i) any(class(newdat[,i])=="numeric"))
					mt[pmm] <- "pmm"
					imp <- mice(newdat,pred=pred,meth=mt,m=1,maxit=10,print=F)
					completeData <- complete(imp)[,unorder.qual_var,drop=F]
					colnames(completeData) <- paste0(colnames(completeData),'_',mt1)
					return(completeData)
				}
				if(length(predictor)!=0){
					newdat <- dataset[,c(predictor,unorder.qual_var),drop=F]
					ini <- mice(newdat,maxit=0)
					pred <- ini$predictorMatrix
					pred[predictor,] <- pred[unorder.qual_var,unorder.qual_var] <- 0
				} else {
					newdat <- dataset[,unorder.qual_var,drop=F]
					newdat <- data.frame(obs=1:nrow(newdat),newdat)
					ini <- mice(newdat,maxit=0)
					pred <- ini$predictorMatrix
					pred[unorder.qual_var,] <- 0
				}
				resIMP <- lapply(method.unorderqual,getIMP)
				resIMP.unorderqual <- do.call(cbind,resIMP)
				O <- data.frame(O,resIMP.unorderqual)
				IR <- rbind(IR,c('Unordered categorical variable',paste(colnames(resIMP.unorderqual),collapse=', ')))
			}

			# method : ordered categorical variable
			if(length(method.orderqual)!=0){
				for(i in order.qual_var) dataset[,i] <- factor(dataset[,i],order=T,levels=level)
				getIMP <- function(mt){
					imp <- mice(newdat,pred=pred,m=1,maxit=10,print=F)
					completeData <- complete(imp)[,order.qual_var,drop=F]
					colnames(completeData) <- paste0(colnames(completeData),'_',mt)
					return(completeData)
				}
				if(length(predictor)!=0){
					newdat <- dataset[,c(predictor,order.qual_var),drop=F]
					ini <- mice(newdat,maxit=0)
					pred <- ini$predictorMatrix
					pred[predictor,] <- pred[order.qual_var,order.qual_var] <- 0
				} else {
					newdat <- dataset[,order.qual_var,drop=F]
					newdat <- data.frame(obs=1:nrow(newdat),newdat)
					ini <- mice(newdat,maxit=0)
					pred <- ini$predictorMatrix
					pred[order.qual_var,] <- 0
				}
				resIMP <- lapply(method.orderqual,getIMP)
				resIMP.orderqual <- do.call(cbind,resIMP)
				O <- data.frame(O,resIMP.orderqual)
				IR <- rbind(IR,c('Ordered categorical variable',paste(colnames(resIMP.orderqual),collapse=', ')))
			}

			# method : any variable
			if(length(method.any)!=0){
				getIMP <- function(mt){
					if(mt=="cart"|mt=="rf") {
						load.pkg("randomForest")

						imp <- mice(newdat,meth=mt,m=1,maxit=10,print=F)
					} else {
						imp <- mice(newdat,pred=pred,meth=mt,m=1,maxit=10,print=F)
					}
					completeData <- complete(imp)[,c(quan_var,unorder.qual_var,order.qual_var),drop=F]
					colnames(completeData) <- paste0(colnames(completeData),'_',mt)
					return(completeData)
				}
				if(length(predictor)!=0){
					newdat <- dataset[,c(predictor,quan_var,unorder.qual_var,order.qual_var),drop=F]
					ini <- mice(newdat,maxit=0)
					pred <- ini$predictorMatrix
					pred[predictor,] <- pred[quan_var,quan_var] <- 0
				} else {
					newdat <- dataset[,c(quan_var,unorder.qual_var,order.qual_var),drop=F]
					newdat <- data.frame(obs=1:nrow(newdat),newdat)
					ini <- mice(newdat,maxit=0)
					pred <- ini$predictorMatrix
					pred[quan_var,] <- 0
				}
				resIMP <- lapply(method.any,getIMP)
				resIMP.any <- do.call(cbind,resIMP)
				O <- data.frame(O,resIMP.any)
				IR <- rbind(IR,c('Any type',paste(colnames(resIMP.any),collapse=', ')))
			}

			O <- O[,-1]		
			IR <- as.data.frame(IR[-1,,drop=F])
			colnames(IR) <- c('Variable type','Name of imputed variable')
			HTML(IR,file=stdout(), innerBorder = 1,align="left",row.names=F)
		}
		# Analysis end time
		R2HTML::HTMLhr(file=stdout())
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Missing Data Imputation",sep=""),file=stdout())
		R2HTML::HTMLhr(file=stdout())
	})
	if(exists('O')){
		return(list(html=html.output,Output=O))
	} else {
		return(html.output)
	}
}

## 이변량 상관분석

REx_CORR <- function(dataset, var.x, var.y, Method = "pearson", Test.side = "two",exact=TRUE,continuity=FALSE) {
  #### 변수 설명
  ## dataset: dataset 이름
  ## var.x: 반응변수(Response variable). 필수로 한개 선택, 반응변수 중 첫번째 변수를 x로 대응함
  ## var.y: 또 다른 반응변수 (x와 correlation을 확인하고 싶은 변수). 필수로 한개 선택. 반응변수 중 두번째 변수를 y로 대응함
  ## Method : 통계량 탭에서 correlation 방법 중	
  ##	  "Pearson"을 선택하였을 경우, Method 입력값은 "pearson",
  ##	  "Kendall"을 선택하였을 경우, Method 입력값은 "kendall",
  ##	  "Spearman"을 선택하였을 경우, Method 입력값은 "spearman",
  ##	  선택하지 않았을 때의 디폴트 값은 "Pearson", 중복 선택 불가
  ## Test.side : 양측검정을 선택하였을 경우 'two'
  ##	       좌측 단측 검정을 선택하였을 경우 'less' 
  ##	       우측 단측 검정을 선택하였을 경우 'greater' 
  ##	       선택하지 않았을 경우 디폴트 값은 "two" 검정, 중복 선택 불가
  ## 변수가 숫자 변수가 아닐 경우에는 Warning 메시지만 나타남(All variables should be numeric.).
  ## exact : UI의 정확검정에 해당. 체크되면 TRUE, 그렇지 않으면 FALSE
  ## continutity :UI의 연속성보정에 해당. 체크되면 TRUE, 그렇지 않으면 FALSE
  #################################
  
  #################################
  #### Required packages : R2HTML
  #################################
  load.pkg("R2HTML")
  
  html.output <- capture.output({
    ## Analysis Title
    R2HTML::HTML(R2HTML::as.title("Correlation Analysis"), HR = 1, file=stdout(), append = FALSE) 
    
    ## Warning
    if (!(is.numeric(as.matrix(dataset[,c(var.x, var.y)])) || is.logical(as.matrix(dataset[,c(var.x, var.y)])))) 
      warn.msg1 <- "All variables should be numeric."
    
    if(exists('warn.msg1')) {
      R2HTML::HTML(R2HTML::as.title("Warnings"), HR = 2, file=stdout()) 
      R2HTML::HTML(warn.msg1, file = stdout())
    } else {
      
      corr.side.test <- function (x, y, method = c("pearson", "kendall", "spearman"), test.side = c("two", "less","greater"),exact,continuity,teststat) {
        method <- match.arg(method)
        test.side <- match.arg(test.side)
        
        if(method%in%c('kendall','spearman')){
          corr <- cor.test(x = x, y = y, alternative = test.side, method = method,continuity=continuity,exact=exact)
        }else{
          corr <- cor.test(x = x, y = y, alternative = test.side, method = method)
        }
        
        corr.out <- list(estimate = corr$estimate[[1]], statistic = corr$statistic[[1]], p.value = corr$p.value[[1]])
        if(teststat=="T-value") corr.out$df <- sum(apply(dataset[,c(var.x,var.y)],1,function(x) !any(is.na(x))))-2
        return(corr.out)
      }
      
      ## Data structure
      R2HTML::HTML(R2HTML::as.title("Data Structure"), HR=2, file=stdout())
      DS <- matrix(c('Number of observations', nrow(dataset), 'Number of variables', length(c(var.x, var.y))), ncol=2, byrow=T)
      R2HTML::HTML(DS, file = stdout(), align = "left", innerBorder = 1)
      
      ## Variable list
      R2HTML::HTML(R2HTML::as.title("Variable List"), HR = 2, file = stdout())
      VL <- matrix(c('Response variables', paste(var.x, var.y, sep=', ')), ncol =2, byrow=T)
      R2HTML::HTML(VL, file = stdout(), align = "left", innerBorder = 1)
      
      ## Analysis Description
      R2HTML::HTML(R2HTML::as.title("Analysis Description"), HR = 2, file = stdout())
      method <- switch(Method, pearson = "Pearson correlation", kendall = "Kendall rank correlation", spearman = "Spearman rank correlation")
      AD <- matrix(c('Correlation coefficient',method),ncol=2,byrow=T)
      teststat <- switch(Method, pearson = 'T-value', spearman = 'S-value', kendall = 'Tau-value')
      if(Method%in%c('kendall','spearman')) {
        # check exact
        if(!exact) {
          Dist <- ifelse(Method=="kendall",'standard normal approximation.','t approximation.')
          warn.msg3 <- paste0("<div style='text-align:left'>Since exact test was not selected, p-value was computed via asymptotic ",Dist)
        }
        # check tie
        if((any(duplicated(dataset[,var.x])) | any(duplicated(dataset[,var.y]))) & exact) {
          Dist <- ifelse(Method=="kendall",'standard normal approximation.','t approximation.')
          warn.msg2 <- paste0("<div style='text-align:left'>Since ties were observed in input variables, p-value was computed via asymptotic ",Dist,"<br>(Exact p-value can be computed when there are no ties.)")
          exact=exact
        }
        if(!exact) teststat <- ifelse(Method=="kendall",'Z-value','T-value')
        AD <- rbind(AD,matrix(c('Exact test', exact,'Continuity correction',continuity),ncol=2,byrow=T))
      }
      AD <- rbind(AD,c('Test statistics',teststat))
      # Alternative hypothesis
      H1 <- switch(Test.side,two = paste0('Two response variables are correlated.'),
                   less = paste0('Two response variables are negatively correlated.'),
                   greater = paste0('Two response variables are positively correlated.'))
      AD <- rbind(AD,c('Alternative hypothesis',H1))
      if(exists('warn.msg2')) {
        R2HTML::HTML(AD, file = stdout(), align = "left", innerBorder = 1, caption=warn.msg2)
      } else if(exists('warn.msg3')) {
        R2HTML::HTML(AD, file = stdout(), align = "left", innerBorder = 1, caption=warn.msg3)
      } else {
        R2HTML::HTML(AD, file = stdout(), align = "left", innerBorder = 1)
      }
      
      ## Anlaysis Result
      R2HTML::HTML(R2HTML::as.title("Analysis Result"), HR = 2, file = stdout())
      result <- corr.side.test(x = as.numeric(as.matrix(dataset[, var.x, drop=F])), 
                               y = as.numeric(as.matrix(dataset[, var.y, drop=F])),
                               method = Method, test.side = Test.side,exact=exact,continuity=continuity,teststat=teststat)
      
      nullToNA <- function(x) {
        x[sapply(x, is.null)] <- NA
        return(x)
      }
      
      result1 <- nullToNA(result)
      if(teststat=="T-value"){
        result2 <- as.data.frame(with(result1, matrix(c(method, round(estimate, 4), teststat, round(statistic,4), "Df", df, "P-value", round(p.value, 4)), ncol=2, byrow=T)))
      } else {
        result2 <- as.data.frame(with(result1, matrix(c(method, round(estimate, 4), teststat, round(statistic,4), "P-value", round(p.value, 4)), ncol=2, byrow=T)))
      }
      colnames(result2) <- c('Result','Value')
      R2HTML::HTML(result2, innerBorder = 1, file = stdout(), row.names = F, align = "left")
      
    }
    
    ## 분석종료시간
    R2HTML::HTMLhr(file=stdout())
    R2HTML::HTML(paste("Analysis is finished at ", Sys.time(),". REx : Correlation Analysis", sep=""), file=stdout())
    R2HTML::HTMLhr(file=stdout())
  })
  return(html.output)
}

## 편상관분석

REx_PCORR <- function(dataset, x, y, z, Method = "pearson", Test.side = "two") {
  #### 변수 설명
  ## dataset: dataset 이름
  ## x: 반응변수(Response variable). 필수로 한개 선택
  ## y: 또 다른 반응변수 (x와 correlation을 확인하고 싶은 변수). 필수로 한개 선택
  ## z : 제어변수(control variable) : x와 y의 관계를 파악하고자 할 때, 보정하고자 하는 변수, 필수 아님
  ## Method : 통계량 탭에서 correlation 방법 중 "Pearson"을 선택하였을 경우, Method 입력값은 "pearson",
  ##		"Kendall"을 선택하였을 경우, Method 입력값은 "kendall",
  ##		"Spearman"을 선택하였을 경우, Method 입력값은 "spearman",
  ##		선택하지 않았을 때의 디폴트 값은 "Pearson"
  ## Test.side : 유의성 검정 시에 양측 / 단측 검정 설정, 선택하지 않았을 경우 디폴트 값은 Two-sided 검정. 선택할 수 있는 값은 'two','one_left','one_right' 중 하나
  #################################
  #################################
  #### Required packages : R2HTML
  #################################
  load.pkg("R2HTML")
  
  html.output <- capture.output({
    ## Analysis Title
    R2HTML::HTML(R2HTML::as.title("Partial Correlation Analysis"), HR = 1, file=stdout(), append = FALSE) 
    
    ## Warning
    if (!(is.numeric(as.matrix(dataset[,c(x,y,z)])) || is.logical(as.matrix(dataset[,c(x,y,z)])))) 
      warn.msg1 <- "All variables should be numeric."
    
    if(exists('warn.msg1')) {
      R2HTML::HTML(R2HTML::as.title("Warnings"), HR = 2, file=stdout()) 
      R2HTML::HTML(warn.msg1, file = stdout())
    } else {
      ## Required user-defined function
      pcor.side <- function (x, method = c("pearson", "kendall", "spearman"), test.side = c("two", "one_left","one_right")) {
        method <- match.arg(method)
        test.side <- match.arg(test.side)
        if (is.data.frame(x)) 
          x <- as.matrix(x)
        n <- dim(x)[1]
        gp <- dim(x)[2] - 2
        cvx <- cov(x, method = method)
        if (det(cvx) < .Machine$double.eps) {
          warn.msg2 <- "The inverse of variance-covariance matrix is calculated using Moore-Penrose generalized matrix invers due to its determinant of zero."
          icvx <- ginv(cvx)
        }
        else icvx <- solve(cvx)
        pcor <- -cov2cor(icvx)
        diag(pcor) <- 1
        if (method == "kendall") {
          statistic <- pcor/sqrt(2 * (2 * (n - gp) + 5)/(9 * (n - gp) * (n - 1 - gp)))
          if (test.side == "two") {
            p.value <- 2 * pnorm(-abs(statistic))
          } else if (test.side == "one_left") {
            p.value <- pnorm(statistic)
          } else if (test.side == "one_right") {
            p.value <- pnorm(statistic,lower.tail=F)
          }
        } else {
          statistic <- pcor * sqrt((n - 2 - gp)/(1 - pcor^2))
          if (test.side == "two") {
            p.value <- 2 * pt(-abs(statistic), (n - 2 - gp))
          } else if (test.side == "one_left") {
            p.value <- pt(statistic, (n - 2 - gp))
          } else if (test.side == "one_right") {
            p.value <- pt(statistic, (n - 2 - gp),lower.tail=F)
          }
        }
        diag(statistic) <- 0
        diag(p.value) <- 0
        if(exists('warn.msg2')){
          list(estimate = pcor, p.value = p.value, statistic = statistic, n = n, gp = gp, method = method, warn.msg2=warn.msg2)
        } else {
          list(estimate = pcor, p.value = p.value, statistic = statistic, n = n, gp = gp, method = method, warn.msg2=NULL)
        }
      }
      
      pcor.side.test <- function (x, y, z, method = c("pearson", "kendall", "spearman"), test.side = c("two", "one_left","one_right")) {
        method <- match.arg(method)
        test.side <- match.arg(test.side)
        x <- c(x)
        y <- c(y)
        z <- as.data.frame(z)
        xyz <- data.frame(x, y, z)
        xyz <- xyz[complete.cases(xyz), ]
        pcor = pcor.side(xyz, method = method, test.side = test.side)
        list(estimate = pcor$est[1, 2], p.value = pcor$p.value[1, 2], statistic = pcor$statistic[1, 2], n = pcor$n, gp = pcor$gp, warn.msg2=pcor$warn.msg2)
      }
      
      ## Data structure
      R2HTML::HTML(R2HTML::as.title("Data Structure"), HR=2, file=stdout())
      DS <- matrix(c('Number of observations',nrow(dataset),'Number of variables',length(c(x,y,z))),ncol=2,byrow=T)
      R2HTML::HTML(DS, file = stdout(), align = "left", innerBorder = 1)
      
      ## Variable list
      R2HTML::HTML(R2HTML::as.title("Variable List"), HR = 2, file = stdout())
      VL <- matrix(c('Response variables',paste(x,y,sep=', '),'Control variables',paste(z,collapse=', ')),ncol=2,byrow=T)
      R2HTML::HTML(VL, file = stdout(), align = "left", innerBorder = 1)
      
      ## Analysis Description
      R2HTML::HTML(R2HTML::as.title("Analysis Description"), HR = 2, file = stdout())
      method <- switch(Method,pearson="Pearson r correlation",kendall="Kendall rank correlation",spearman="Spearman rank correlation")
      H1 <- switch(Test.side,two=paste0('With ',paste0(z,collapse=", "),' controlling, two response variables are correlated.'),one_left=paste0('With ',paste0(z,collapse=", "),' controlling, two response variables are negatively correlated.'),one_right=paste0('With ',paste0(z,collapse=", "),' controlling, two response variables are positively correlated.'))
      AD <- matrix(c('Correlation coefficient',method,'Alternative hypothesis',H1),ncol=2,byrow=T)
      R2HTML::HTML(AD, file = stdout(), align = "left", innerBorder = 1)
      
      ## Anlaysis Result
      R2HTML::HTML(R2HTML::as.title("Analysis Result"), HR = 2, file = stdout())
      result <- pcor.side.test(x = dataset[, x,drop=F], y = dataset[, y,drop=F], z = dataset[, z,drop=F], method = Method, test.side = Test.side)
      teststat <- switch(Method,pearson='T-value',spearman='T-value',kendall='Z-value')
      result1 <- with(result,matrix(c(method,round(estimate,4),teststat,round(statistic,4)),ncol=2,byrow=T))
      if(Method%in%c('pearson','spearman')) result1 <- rbind(result1,with(result,c('Degrees of freedom',n-2-gp)))
      result1 <- with(result,rbind(result1,c('P-value',round(p.value,4))))
      result2 <- as.data.frame(result1)
      colnames(result2) <- c('Result','Value')
      if(length(result$warn.msg2)==0){
        R2HTML::HTML(result2, innerBorder = 1, file = stdout(), row.names = F, align = "left")
      } else {
        cap <- paste("<div style='text-align:left'>",result$warn.msg2)
        R2HTML::HTML(result2, innerBorder = 1, file = stdout(), row.names = F, align = "left",caption=cap)
      }	
    }
    ## 분석종료시간
    R2HTML::HTMLhr(file=stdout())
    R2HTML::HTML(paste("Analysis is finished at ", Sys.time(),". REx : Partial Correlation Analysis", sep=""), file=stdout())
    R2HTML::HTMLhr(file=stdout())
  })
  return(html.output)
}

REx_PCA <- function(dataset,variables,scaling=TRUE, SP_PC=FALSE, PCprint=FALSE, npc=2){
  ####Arguments ####
  ## datset : 데이터셋 이름
  ## variables: 선택된 변수들. UI의 변수선택 탭의 선택변수에 해당
  ## scaling: UI의 분석옵션 탭의 
  ## SP-PC: Scree plot/ GUI의 스크리 도표에 해당
  ## PCprint: 
  ## npc: 엑셀창에 출력할 PC의 수
  ###########################################################################
  #### Required packages : R2HTML, psy, psych ###
  ###########################################################################
  load.pkg(c("psy", "psych", "R2HTML"))
  
  html.output <- capture.output({
    R2HTML::HTML(R2HTML::as.title("Principal Component Analysis"),HR=1,file=stdout(),innerBorder = 1,align="left",append=FALSE) 
    
    ## Warnings
    if(length(variables)==1) {
      warn.msg1 <- "\a Error : At least 2 variables should be selected. Analysis has stopped."
    } else {
      isnom <- !sapply(variables,function(i) is.numeric(dataset[,i]))
      if(sum(isnom)>0) {
        nn <- length(variables) - sum(isnom)
        msg <- 'PCA supports only numeric variables but categorical variables were selected. They were excluded from the analysis.'
        if(nn>1) {
          warn.msg2 <- paste0('\a Warning : ',msg)
          variables <- variables[!isnom]
        } else if(nn==1){
          warn.msg3 <- paste0('\a Error : ',msg,' After excluding them, there is a numeric variable. PCA can be conducted with at least 2 variables. Analysis has stopped.')
        } else {
          warn.msg4 <- paste0('\a Error : ',msg,' After excluding them, there is no variable. Analysis has stopped.')
        }
      }
    }
    
    if(exists('warn.msg1')|exists('warn.msg3')|exists('warn.msg4')) {
      R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout(),innerBorder = 1,align="left") 
      if(exists('warn.msg1')) {
        R2HTML::HTML(warn.msg1,file=stdout())
      } else {
        warn.msg <- ifelse(exists('warn.msg3'),warn.msg3,warn.msg4)
        R2HTML::HTML(warn.msg,file=stdout())
      }
    } else {
      
      # Data structure
      R2HTML::HTML(R2HTML::as.title("Data Structure"),HR=2,file=stdout(),innerBorder = 1,align="left") 
      data <- dataset[,variables,drop=F]
      data<-data[complete.cases(data),]
      DS <- matrix(c('Number of observations','Number of variables',nrow(dataset),length(variables)),nrow=2)
      if(exists('warn.msg2')) {
        R2HTML::HTML(DS,file=stdout(),innerBorder = 1,align="left",caption=paste0('<div style="text-align:left">',warn.msg2))
      } else {
        R2HTML::HTML(DS,file=stdout(),innerBorder = 1,align="left")
      }
      
      # Analysis description	
      R2HTML::HTML(R2HTML::as.title("Analysis Description"),HR=2,file=stdout(),innerBorder = 1,align="left")
      AD <- matrix(c('Variables',paste(variables,collapse=', '),'Scaling Variances',scaling),ncol=2,byrow=T)
      if(PCprint) AD <- rbind(AD,c('Number of printed PCs',npc))
      R2HTML::HTML(AD,file=stdout(),innerBorder = 1,align="left")
      
      # Fit PCA
      fit<-princomp(data, cor=scaling)
      R2HTML::HTML(R2HTML::as.title("PCA results"),HR=2,file=stdout(),innerBorder = 1,align="left")
      
      # Descriptive stat
      R2HTML::HTML(R2HTML::as.title("Descriptive Statistics"),HR=3,file=stdout(),innerBorder = 1,align="left")
      d <- as.data.frame(cbind(fit$center,sqrt(diag(cov(data)))))
      colnames(d)<-c("Mean","SD")
      R2HTML::HTML(d,file=stdout(),innerBorder = 1,align="left",digits=4)
      
      ## variable loadings
      R2HTML::HTML(R2HTML::as.title("Variable Loadings"),HR=3,file=stdout(),innerBorder = 1,align="left")
      FL <- loadings(fit)[1:(ncol(data)),]
      colnames(FL) <- paste0('PC',1:ncol(data))
      R2HTML::HTML(FL,file=stdout(),innerBorder = 1,align="left",digits=4)
      if(npc>1){
        for(i in 1:(npc-1)){
          for(j in (i+1):npc){
            R2HTML::HTML(R2HTML::as.title(paste0('PC',i,' & PC',j)),HR=4,file=stdout(),innerBorder = 1,align="left")
            xmax <- max(abs(FL[,i]))
            ymax <- max(abs(FL[,j]))
            REx_ANA_PLOT()
            plot(1,type="n",xlim=c(-xmax,xmax),ylim=c(-ymax,ymax),main="Variable Loadings",xlab=paste0('PC',i),ylab=paste0('PC',j))
            arrows(rep(0,length(variables)),rep(0,length(variables)),FL[,i],FL[,j],col=rep('red',length(variables)))
            text(x=FL[,i],y=FL[,j],labels=variables)
            REx_ANA_PLOT_OFF("")
          }
        }
      }
      
      #### Importance of Components
      R2HTML::HTML(R2HTML::as.title("Importance of Components"),HR=3,file=stdout(),innerBorder = 1,align="left") 
      IC <- capture.output(summary(fit))[-c(1:2)]
      IC.1 <- strsplit(IC,"[ ]+")
      IC.2 <- t(sapply(1:length(IC.1),function(i) as.numeric(IC.1[[i]][(length(IC.1[[i]])-ncol(data)+1):length(IC.1[[i]])])))
      colnames(IC.2) <- paste0('PC',1:ncol(data))
      rownames(IC.2) <- c('Standard Deviation','Proportion of Variance','Cumulative Proportion')
      R2HTML::HTML(IC.2,file=stdout(),innerBorder = 1,align="left",digits=4)
      
      #### Scree plot   
      if (SP_PC==TRUE){
        R2HTML::HTML(R2HTML::as.title("Scree plot"),HR=3,file=stdout(),innerBorder = 1,align="left") 
        REx_ANA_PLOT()
        plot(fit$sdev,type="b",ylab="Standard Deviation",xlab="Principal Components",xaxt="n",main="Scree Plot")
        axis(1,at=1:ncol(data),labels=paste0('PC',1:ncol(data)))
        REx_ANA_PLOT_OFF("")
      }
      
      #### The coordinates of the individuals (observations) on the principal components
      if(PCprint){
        cm <- as.matrix(fit$scores)[,1:npc,drop=F]
        colnames(cm) <- paste0('PC',1:npc)
        cm <- data.frame(cm)
      }
    }
    
    R2HTML::HTMLhr(file=stdout())
    R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Principal Component Analysis",sep=""),file=stdout(),innerBorder = 1,align="left")
    R2HTML::HTMLhr(file=stdout())
  })
  if(PCprint){
    return(list(html=html.output,Output=cm))
  } else {
    return(html.output)
  }
}

REx_cluster_hierachical<-function(datsett,conti_varname=NULL,dummy_varname=NULL,centers=2,maxcluster=20,
                      standardize=TRUE,dist_method='euclidean',linkage_method='complete',Silhouette=FALSE){
	#datsett : 데이터셋 이름
	#conti_varname=연속형 변수로 선택된 변수들 이름
	#		만약 선택된 변수가 두개 이상이면, c('A','B','C') 이런식으로.
	#		만약 선택된 변수가 하나라면, 'A' 이렇게표현
	#dummy_varname=범주형 변수로 선택된 변수들 이름
	#		만약 선택된 변수가 두개 이상이면, c('A','B','C') 이런식으로.
	#		만약 선택된 변수가 하나라면, 'A' 이렇게표현
	#dist_method :데이터간의 거리계산 방법. 'euclidean','manhattan','maximum','gower' 중 하나 선택가능/ Default : 'euclidean'
	#linkage_method :	그림선택기능 : h-clustering	알고리즘 'complete','ward.D','ward.D2','single','average','mcquitty','median','centroid' 중 하나 선택가능/ Default : 'complete',
	#standardize : 데이터 표준화 / Default : TRUE
	#centers : 군집수.
	#maxcluster : 최대군집수 (입력된 숫자값.) / Default : 10
	#Silhouette : 실루엣 
	
	#######################################################
	#### Required packages : R2HTML,cluster ,fpc,markdown
	#######################################################
	load.pkg(c("R2HTML", "cluster", "fpc"))
	
	datsett$EXR_numord<-1:nrow(datsett)
  
  html.output <- capture.output({
    ## Analysis title
    R2HTML::HTML(R2HTML::as.title('Hierarchical Cluster Analysis'),HR=1,file=stdout(),append=FALSE)
    
    ## Warnings
    # character conti_varname : analysis stopped
    if(!is.null(conti_varname)){
      is.chr <- !sapply(datsett[,conti_varname,drop=F],is.numeric)
      if(any(is.chr)) warn.msg1 <- paste0('\a Analysis has been stopped because non-numeric variable ',paste(conti_varname[is.chr],collapse=", "),' was included in quantitative variable.')
    }
    
    # numeric dummy_varname & method
    if(!is.null(dummy_varname)){
      is.num <- sapply(datsett[,dummy_varname,drop=F],is.numeric)
      for(i in 1:length(dummy_varname)) if(is.num[i]) datsett[,dummy_varname[i]] <- as.character(datsett[,dummy_varname[i]])
      if(any(is.num)) warn.msg2 <- paste0("<div style='text-align:left'> Warning message : Among qualitative variables you selected, the type of variable '",paste(dummy_varname[is.num],collapse=', '),"' was numeric. It was coerced into character.")
      
      if(!dist_method%in%'gower'){
        dist_method <- 'gower'
        warn.msg3 <- paste0("<div style='text-align:left'> \a Warning message : Only 'gower' can be supported for categorical variable. 'gower' will be used since categorical variable is selected.")
      }
    }
    
    # Maximum number of clusters
    temp.dat <- datsett[c(conti_varname,dummy_varname)]
    num.na <- sum(apply(temp.dat,1,function(i) any(is.na(i))))
    if(Silhouette & maxcluster>min(c(nrow(temp.dat)-num.na,20))) {
      warn.msg4 <- paste0("<div style='text-align:left'> \a Warning message : After removing observations with missing data, we had ",nrow(temp.dat)-num.na," observations. Maximum number of clusters for silhouette was replaced with ",nrow(temp.dat)-num.na,".")
      maxcluster <- min(c(nrow(temp.dat)-num.na,20))
    }

		#####################################################################################	
		#######################hierarchical clustering########################################
		######################################################################################
			
		## Data structure
		R2HTML::HTML(R2HTML::as.title('Data Structure'),HR=2,file=stdout())
    DS <- matrix(c("Number of observations","Number of variables",nrow(datsett),length(c(conti_varname,dummy_varname))),ncol=2)
    R2HTML::HTML(DS,file=stdout(), innerBorder = 1,align="left")
    
    ## Variable List
    R2HTML::HTML(R2HTML::as.title('Variable List'),HR=2,file=stdout())
    if(is.null(conti_varname)) {
      varlist <- matrix(c('Qualitative Variable',paste(dummy_varname,collapse=', ')),ncol=2)
    } else {
      varlist <- matrix(c('Quantitative Variable',paste(conti_varname,collapse=', ')),ncol=2)
      if(!is.null(dummy_varname)) varlist <- rbind(varlist,c('Qualitative Variable',paste(dummy_varname,collapse=', ')))
    }
    if(exists('warn.msg2')) {
      R2HTML::HTML(varlist,file=stdout(), innerBorder = 1,align="left",caption=warn.msg2)
    } else {
      R2HTML::HTML(varlist,file=stdout(), innerBorder = 1,align="left")
    }
    
    ## Descriptions of Analysis
    R2HTML::HTML(R2HTML::as.title('Descriptions of Analysis'),HR=2,file=stdout())
    DA <- matrix(c("Number of clusters","Method of distance","Method of linkage",centers,dist_method,linkage_method),ncol=2)
    if(Silhouette)	DA <- rbind(DA,matrix(c('Silhouette',Silhouette,'Maximum number of clusters for Silhouette',maxcluster),ncol=2,byrow=T))
    if(exists('warn.msg3')&exists('warn.msg4')) {
      R2HTML::HTML(DA,file=stdout(), innerBorder = 1,align="left",caption=paste(warn.msg3,warn.msg4))
    } else if(exists('warn.msg3')) {
      R2HTML::HTML(DA,file=stdout(), innerBorder = 1,align="left",caption=warn.msg3)
    } else if(exists('warn.msg4')) {
      R2HTML::HTML(DA,file=stdout(), innerBorder = 1,align="left",caption=warn.msg4)
    } else {
      R2HTML::HTML(DA,file=stdout(), innerBorder = 1,align="left")
    }

		############################################
		#############making data####################
		############################################

		if(!is.null(conti_varname)) {
      datt<-datsett[,c(conti_varname)]
      # Standardized
      if(standardize){
        datt<-data.frame(scale(datt))
      }
    }
    
    if(!is.null(dummy_varname)){
      datt2<-datsett[,dummy_varname,drop=F]
      for(i in 1:ncol(datt2)){
        if(class(datt2[,i])!='factor'){
          datt2[,i]<-as.factor(datt2[,i])
        } 
      }
      if(exists('datt')){
        datt <- cbind(datt,datt2)
      } else {
        datt <- datt2
      }
    }
    datt$EXR_numord<-1:nrow(datt)
			
		############################################		
		#########Calculate distance#################
		############################################
		if(length(dummy_varname)>=1){
      D<-daisy(na.omit(datt[,-ncol(datt)]),metric='gower',stand=FALSE) 
    }else{
      D<-dist(na.omit(datt[,-ncol(datt)]),method=dist_method)
    }
			
		#############################################
		#######Clustering by hclust##################
		##############################################		
		R2HTML::HTML(R2HTML::as.title("Results - Hierarchical Clustering"),HR=2,file=stdout()) # k-means 결과(군집간 할당된 샘플수, 제곱합 정보)
		R2HTML::HTML(R2HTML::as.title("Hierarchical Cluster Plot"),HR=3,file=stdout()) # k-means 결과(군집간 할당된 샘플수, 제곱합 정보)
		H.fit <- hclust(D, method=linkage_method)
    
    REx_ANA_PLOT(w=1000,h=1000)
    plot(H.fit)
    rect.hclust(H.fit, k=centers, border="red") 
    REx_ANA_PLOT_OFF("")

		R2HTML::HTML(R2HTML::as.title("Scatter Plot by Clustering Result"),HR=3,file=stdout()) # k-means 결과(군집간 할당된 샘플수, 제곱합 정보)
		groups <- cutree(H.fit, k=centers) 
    REx_ANA_PLOT(w=700,h=700)
    plot(datt[,-ncol(datt)],col=groups+1)
    REx_ANA_PLOT_OFF("")
			
		#############################################
		#######Silhouette Plot ######################
		##############################################			 
		if(Silhouette==TRUE){
      R2HTML::HTML(R2HTML::as.title("Silhouette"),HR=3,file=stdout()) 
      silhouette_hclust <- silhouette(cutree(H.fit, k = centers), D)
      
      REx_ANA_PLOT()
      plot(silhouette_hclust,main='Current Silhouette',border=NA,col=1:length(unique(silhouette_hclust[,1])))
      REx_ANA_PLOT_OFF("")
      
      #理쒖쟻?쓽 clustering number by silhouette
      ks<-2:maxcluster 
      set.seed(1234)
      ASW <- sapply(ks, FUN=function(k) {
        hh.fit<-hclust(D, method=linkage_method)
        out<-cutree(hh.fit,k=k)
        avg_sil<-fpc::cluster.stats(D,clustering=out)$avg.silwidth
        Dunn<-fpc::cluster.stats(D,clustering=out)$dunn
        return(c(avg_sil,Dunn))
      })
      
      R2HTML::HTML(R2HTML::as.title("Average Silhouette by cluster numbers"),HR=3,file=stdout()) 
      REx_ANA_PLOT()
      plot(ks, ASW[1,], type="l",main='Optimal cluster numbers',xlab='Cluster',ylab='Average Sihouette')
      maxout<-ks[which.max(ASW[1,])]
      abline(v=ks[which.max(ASW[1,])], col="red", lty=2)
      REx_ANA_PLOT_OFF("Red vertical line indicates optimal clustering numbers by average shilhouette.")
      
      R2HTML::HTML(R2HTML::as.title("Dunn Index by cluster numbers"),HR=3,file=stdout()) 
      REx_ANA_PLOT()
      plot(ks, ASW[2,], type="l",main='Optimal cluster numbers',xlab='Cluster',ylab='Dunn Index')
      abline(v=ks[which.max(ASW[2,])], col="red", lty=2)
      REx_ANA_PLOT_OFF("Red vertical line indicates optimal clustering numbers by Dunn Index.")
    }
			
		##################################################		
		########cluster 결과: 최종데이터에 붙이기#########
		##################################################
		datt<-cbind(datt,clustering=groups)
    result_cluster_dat<-merge(datsett,datt[,c('EXR_numord','clustering')],by='EXR_numord',all=T)
    result_cluster_dat$EXR_numord<-NULL
    O <- data.frame(Hcluster = result_cluster_dat$clustering)
    
    # Analysis end time
    R2HTML::HTMLhr(file=stdout())
    R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Hierarchical Cluster Analysis",sep=""),file=stdout())
    R2HTML::HTMLhr(file=stdout())
  })
  return(list(html=html.output,Output=O))
}

REx_oneTtest <- function(dataset, res_var ,alternative = "two.sided", mu=0,CI=F, conf.level = 0.95, digit=5)  {

	## Variables 
   	# res_var		: 변수선택 탭의 반응변수.
	# mu			: 분석옵션 탭의 모평균.
	# alternative		: 분석옵션 탭의 검정방법. 'two-sided' or 'less' or 'greater'
	# CI			: 분석옵션 탭의 신뢰구간 출력
	# conf.level		: 분석옵션 탭의 유의수준
	# digit 		: # of digits to be printed


	## Required packages : R2HTML
	load.pkg("R2HTML") 
	
	# Output  
	html.output <- capture.output({
		
		# HTML file title
		R2HTML::HTML(R2HTML::as.title("One Sample T-test"), HR = 1, file = stdout(), append = F)
				
		# Warnings
		if(!is.numeric(dataset[,res_var])) {
			R2HTML::HTML(R2HTML::as.title("Warnings"), HR = 2, file = stdout())
			warn.msg <- "Analysis has been stopped since non-numeric response variable was selected."
			R2HTML::HTML(warn.msg, file=stdout(), innerBorder = 1, align = "left") 
		} else {
			options(digits=digit)
			
			data <- dataset[,res_var,drop=F]
			name.var <- names(data)
			n.sample <- nrow(dataset)
			n.miss <- sum(!complete.cases(data))
			data <- data[complete.cases(data), ]
		 
			# Descriptive statistics
			mean.var<- mean(data)
			sd.var <- sd(data)
			
			# One sample t-test
			t.test <- t.test(x=data, alternative = alternative, mu=mu,conf.level = conf.level)
			t.test.t <- t.test$statistic
			t.test.df <- t.test$parameter
			t.test.p <- t.test$p.value
			t.test.mean <- t.test$estimate
			t.test.l.ci <- t.test$conf.int[1]
			t.test.u.ci <- t.test$conf.int[2]

			# Data structure
			R2HTML::HTML(R2HTML::as.title("Data Structure"), HR = 2, file = stdout())
			data_str <- matrix(c("Total number of observations", n.sample,"Number of observations with missing",n.miss), ncol=2,byrow=T)
			R2HTML::HTML(data_str, file=stdout(), row.names = FALSE, col.names = TRUE, innerBorder = 1, align = "left") ;

			# Analysis Description
			R2HTML::HTML(R2HTML::as.title("Analysis Description"), HR = 2, file = stdout())
			H1 <- capture.output(t.test)[grep("alternative",capture.output(t.test))]
			H1 <- gsub(".+: ","",H1)
			AD <- matrix(c('Response variable',res_var,'Alternative hypothesis',H1),ncol=2,byrow=T)
			if(CI) AD <- rbind(AD,c('Significance level for CI',conf.level))
			R2HTML::HTML(AD, file=stdout(), row.names = FALSE, col.names = TRUE, innerBorder = 1, align = "left") ;

			# Descriptive statistics
			R2HTML::HTML(R2HTML::as.title("Descriptive Statistics"), HR = 2, file = stdout())
			VL <- matrix(c(n.sample-n.miss,mean.var,sd.var), nrow=1, byrow=T)	
			colnames(VL) <- c('N','Mean','Std')
			rownames(VL) <- c(name.var)
			R2HTML::HTML(VL, file = stdout(), row.names = T,col.names = T, innerBorder = 1, align = "left")

			# One Sample T-Test
			R2HTML::HTML(R2HTML::as.title("One Sample T-Test Result"), HR = 2, file = stdout())
			result <- data.frame(t.test.mean,t.test.t, t.test.df,t.test.p)
			colnames(result) <- c('Estimate','T-value','DF','P-value')
			rownames(result) <- res_var
			if(CI) {
				result <- data.frame(result,t.test.l.ci,t.test.u.ci)
				colnames(result)[5:6] <- c(paste0('Lower bound of ',conf.level*100,'% CI'), paste0('Upper bound of ',conf.level*100,'% CI'))
			}
			R2HTML::HTML(result, file = stdout(), innerBorder = 1, align = "left")
		}

		#### End of output
		R2HTML::HTMLhr(file=stdout())
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : One Sample T-test",sep=""),file=stdout())
		R2HTML::HTMLhr(file=stdout())
	})
	
	return(html.output)
}

REx_twoTtest <- function(dataset, res_var, group_var, alternative = "two.sided", CI=F, conf.level = 0.95, levenes.test=T, var.equal=F, digit=4, mu=0)  {

	## Variables 
	# res_var		: 변수선택 탭의 반응변수.
	# group_var		: 변수선택 탭의 집단변수.
	# mu			: 분석옵션 탭의 모평균의 차이.
	# alternative		: 분석옵션 탭의 검정방법. 'two-sided' or 'less' or 'greater'
	# levenes.test		: 분석옵션 탭의 등분산 검정. a levene's test to test equalness of variance b/w groups
	# val.equal		: 분석옵션 탭의 등분산 가정. t-test with equal or unequal variance assumption
	# CI			: 분석옵션 탭의 신뢰구간 출력
	# conf.level		: 분석옵션 탭의 유의수준
	# digit 		: # of digits to be printed

	# Required packages : R2HTML, car
	load.pkg(c("R2HTML", "car"))
	options(digits=digit)

	html.output <- capture.output({
		# HTML file title
		R2HTML::HTML(R2HTML::as.title("Two Sample T-test"), HR = 1, file = stdout(), append = F)
		dataset <- dataset[complete.cases(dataset[,c(res_var,group_var),drop=F]),]
	
		# warnings
		if(!is.numeric(dataset[,res_var])) warn.msg1 <- "Analysis has been stopped since non-numeric response variable was selected."
		if(length(unique(dataset[,group_var]))!=2) warn.msg2 <- paste0("Dataset should consist of two groups. (Observed groups: ",paste(as.character(unique(dataset[,group_var])),collapse=", "),") Analysis has been stopped.")
		if(is.numeric(dataset[,group_var])) dataset[,group_var] <- factor(dataset[,group_var])
		if(exists('warn.msg1')|exists('warn.msg2')){
			R2HTML::HTML(R2HTML::as.title("Warnings"), HR = 2, file = stdout())
			if(exists('warn.msg1')) R2HTML::HTML(warn.msg1, file=stdout(), innerBorder = 1, align = "left") 
			if(exists('warn.msg2')) R2HTML::HTML(warn.msg2, file=stdout(), innerBorder = 1, align = "left") 
		} else {

			# Data structure
			group1 <- as.character(unique(dataset[,group_var])[1]); n.group1 <- sum(dataset[,group_var]==group1)
			group2 <- as.character(unique(dataset[,group_var])[2]); n.group2 <- sum(dataset[,group_var]==group2)

			data_str <- matrix(c("Number of observations", paste("Number of observations in group ",group1,sep=""),paste("Number of observations in group ",group2,sep=""), nrow(dataset), n.group1,n.group2), ncol=2)
			R2HTML::HTML(R2HTML::as.title("Data Structure"), HR = 2, file = stdout())
			R2HTML::HTML(data_str, file=stdout(), row.names = FALSE, col.names = TRUE, innerBorder = 1, align = "left") ;
			
			# Variable list
			R2HTML::HTML(R2HTML::as.title("Variable List"), HR = 2, file = stdout())
			VL <- matrix(c('Response variable','Group variable',res_var,group_var), ncol=2)
			R2HTML::HTML(VL, file = stdout(), row.names = FALSE, col.names = TRUE, innerBorder = 1, align = "left")
			
			# Analysis Description
			R2HTML::HTML(R2HTML::as.title("Analysis Description"), HR = 2, file = stdout())
			t.test.res <- t.test(formula(paste(res_var,'~',group_var)),data=dataset, alternative = alternative, var.equal=var.equal, conf.level = conf.level, mu=mu)
			H1 <- capture.output(t.test.res)[grep("alternative",capture.output(t.test.res))]
			H1 <- gsub(".+: ","",H1)
			AD <- matrix(c('Alternative hypothesis',H1,'Test for homogeneity variance',levenes.test,'Assumption of homogeneity variance',var.equal),ncol=2,byrow=T)
			if(CI) AD <- rbind(AD,c('Significance level for CI',conf.level))
			R2HTML::HTML(AD, file = stdout(), row.names = FALSE, col.names = TRUE, innerBorder = 1, align = "left")

			# Descriptive Statistics
			R2HTML::HTML(R2HTML::as.title("Descriptive Statistics"), HR = 2, file = stdout())
			DS <- t(sapply(c(group1,group2),function(i) return(c(length(dataset[dataset[,group_var]==i,res_var]),mean(dataset[dataset[,group_var]==i,res_var]),sd(dataset[dataset[,group_var]==i,res_var])))))
			DS <- as.data.frame(DS)
			colnames(DS) <- c('# of obs.','Mean','SD')
			rownames(DS) <- paste('Group',rownames(DS))
			R2HTML::HTML(DS, file = stdout(), row.names = T, col.names = TRUE, innerBorder = 1, align = "left",digits=6)

			# Levene's Test for Homogeneity of Variance between Groups
			if (levenes.test) {	
				R2HTML::HTML(R2HTML::as.title("Levene's Test for Homogeneity of Variance between Groups"), HR = 2, file = stdout())
				res.lev <- leveneTest(y=dataset[,res_var],group=dataset[,group_var])
				res.lev.F <- res.lev$F[1]
				res.lev.DF <- c(res.lev[1])[[1]]
				res.lev.P <- res.lev$"Pr(>F)"[1]
				LT <- data.frame(res.lev.F,res.lev.DF[1],res.lev.DF[2],res.lev.P)
				colnames(LT) <- c('F-value','DF1','DF2','P-value')
				rownames(LT) <- res_var
				R2HTML::HTML(LT, file = stdout(),col.names = TRUE, innerBorder = 1, align = "left")
			}

			# Test results
			R2HTML::HTML(R2HTML::as.title("Two Sample T-Test Result"), HR = 2, file =stdout())
			t.test.res.I <- as.data.frame(matrix(with(t.test.res,c(estimate[1],estimate[2],estimate[1]-estimate[2],statistic,parameter,p.value)),nrow=1))
			colnames(t.test.res.I) <- c(names(t.test.res$estimate),'Difference of the means','T-value','DF','P-value')
			rownames(t.test.res.I) <- res_var
			if(CI) {
				t.test.res.I <- cbind(t.test.res.I,t.test.res$conf.int[1],t.test.res$conf.int[2])
				data.frame(t.test.res.I,t.test.res$conf.int[1],t.test.res$conf.int[2])
				colnames(t.test.res.I)[7:8] <- 	c(paste0('Lower bound of ',conf.level*100,'% CI'),paste0('Upper bound of ',conf.level*100,'% CI'))
			}
			R2HTML::HTML(t.test.res.I, file = stdout(), col.names = TRUE, innerBorder = 1, align = "left",digits=5)
		}

		#### End of output
		R2HTML::HTMLhr(file=stdout())
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Two Sample T-test",sep=""),file=stdout())
		R2HTML::HTMLhr(file=stdout())
	})
	
	return(html.output)
}

REx_pairedTtest <- function(dataset, res_var1, res_var2, alternative = "two.sided", CI=F, conf.level = 0.95, levenes.test=T, var.equal=F, digit=4, mu=0)  {
	## Variables 
	# res_var1		: 변수선택 탭의 첫번째 반응변수.
	# res_var2		: 변수선택 탭의 두번째 반응변수.
	# mu			: 분석옵션 탭의 모평균의 차이.
	# alternative		: 분석옵션 탭의 검정방법. 'two-sided' or 'less' or 'greater'
	# CI			: 분석옵션 탭의 신뢰구간 출력
	# conf.level		: 분석옵션 탭의 유의수준
	# digit 		: # of digits to be printed

	# Required packages : R2HTML
	load.pkg("R2HTML")
	options(digits=digit)

	html.output <- capture.output({
		# HTML file title
		R2HTML::HTML(R2HTML::as.title("Paired Sample T-test"), HR = 1, file = stdout(), append = F)

		n.sample <- nrow(dataset)
		n.miss.var1 <- sum(!complete.cases(dataset[,res_var1,drop=F]))
		n.miss.var2 <- sum(!complete.cases(dataset[,res_var2,drop=F]))
		dataset <- dataset[complete.cases(dataset[,c(res_var1,res_var2),drop=F]),]
	
		# warnings
		if(!is.numeric(dataset[,res_var1])|!is.numeric(dataset[,res_var2])) {
			R2HTML::HTML(R2HTML::as.title("Warnings"), HR = 2, file = stdout())
			warn.msg1 <- "Analysis has been stopped since non-numeric response variable was selected."
			R2HTML::HTML(warn.msg1, file=stdout(), innerBorder = 1, align = "left") 
		} else {
			# Data structure
			#data_str <- matrix(c("Number of observations", nrow(dataset)), ncol=2)
			data_str <- matrix(c("Total number of paired observations", n.sample, paste0("Number of missing observations of ",res_var1),n.miss.var1,paste0("Number of missing observations of ",res_var2),n.miss.var2), nrow=3,byrow=T)

			R2HTML::HTML(R2HTML::as.title("Data Structure"), HR = 2, file = stdout())
			R2HTML::HTML(data_str, file=stdout(), row.names = FALSE, col.names = TRUE, innerBorder = 1, align = "left") ;
			
			# Analysis Description
			R2HTML::HTML(R2HTML::as.title("Analysis Description"), HR = 2, file = stdout())
			t.test.res <- t.test(x=dataset[,res_var1],y=dataset[,res_var2],paired=T,alternative = alternative,mu=mu)
			H1 <- capture.output(t.test.res)[grep("alternative",capture.output(t.test.res))]
			H1 <- gsub(".+: ","",H1)
			AD <- matrix(c('Response variables',paste(c(res_var1,res_var2),collapse=', '),'Alternative hypothesis',H1),ncol=2,byrow=T)
			if(CI) AD <- rbind(AD,c('Significance level for CI',conf.level))
			R2HTML::HTML(AD, file = stdout(), row.names = FALSE, col.names = TRUE, innerBorder = 1, align = "left")

			# Descriptive Statistics
			R2HTML::HTML(R2HTML::as.title("Descriptive Statistics"), HR = 2, file = stdout())
			R2HTML::HTML(R2HTML::as.title("Mean and Standard Deviation"), HR = 3, file = stdout())
			DS <- t(sapply(c(res_var1,res_var2),function(i) return(c(mean(dataset[,i]),sd(dataset[,i])))))
			DS <- as.data.frame(DS)
			colnames(DS) <- c('Mean','SD')
			rownames(DS) <- c(res_var1,res_var2)
			R2HTML::HTML(DS, file = stdout(), row.names = T, col.names = TRUE, innerBorder = 1, align = "left",digits=6)

			R2HTML::HTML(R2HTML::as.title("Sample Correlation Matrix"), HR = 3, file = stdout())
			Cor <- cor(dataset[,c(res_var1,res_var2)])
			R2HTML::HTML(Cor, file = stdout(), row.names = T, col.names = TRUE, innerBorder = 1, align = "left",digits=6)

			# Test results
			R2HTML::HTML(R2HTML::as.title("Paired Sample T-Test Result"), HR = 2, file =stdout())
			t.test.res.I <- as.data.frame(matrix(with(t.test.res,c(estimate,statistic,parameter,p.value)),nrow=1))
			colnames(t.test.res.I) <- c('Mean of the differences','T-value','DF','P-value')
			rownames(t.test.res.I) <- paste(res_var1,res_var2,sep=" & ")
			if(CI) {
				t.test.res.I <- cbind(t.test.res.I,t.test.res$conf.int[1],t.test.res$conf.int[2])
				data.frame(t.test.res.I,t.test.res$conf.int[1],t.test.res$conf.int[2])
				colnames(t.test.res.I)[5:6] <- 	c(paste0('Lower bound of ',conf.level*100,'% CI'),paste0('Upper bound of ',conf.level*100,'% CI'))
			}
			R2HTML::HTML(t.test.res.I, file = stdout(), col.names = TRUE, innerBorder = 1, align = "left",digits=5)
		}

		#### End of output
		R2HTML::HTMLhr(file=stdout())
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Paired Sample T-test",sep=""),file=stdout())
		R2HTML::HTMLhr(file=stdout())
	})
	
	return(html.output)
}

REx_FA <- function(dataset, variables, n.factors=2, Method = "ml", Scores="none", Rotation="varimax", D=TRUE, CORM=TRUE, SP=TRUE) {
  ####Arguments ####
  ## datset : 데이터셋 이름
  ## variables: 선택된 변수들 나열
  ## n.factors: The number of factors to be fitted./ GUI의 요인갯수에 해당 / default=1
  ## Method : "minres", "ols", "wls", "gls", "ml" / GUI의 방법에 해당 / default = "minres"
  ## Scores: "none", "regression", "Bartlett" / GUI의 점수에 해당 / default: "none"
  ## Rotation: "none", "varimax", "quatimax", "promax", "oblimin", "simplimax", "cluster"/ GUI의 회전에 해당/ default: "varimax"
  
  ####요인 분석##########################
  ## D: Descriptive statistics/ GUI의 일변량 기술통계에 해당
  ## CORM : Correlation Matrix, Uniqueness, Proportion of Variance / GUI 의 초기 해법에 해당
  ## SP: Scree plot/ GUI의 요인분석 스크리 도표에 해당
  ###########################################################################
  #### Required packages : R2HTML, psy, psych, GPArotation
  ###########################################################################
  
  #### Load Required Packages
  load.pkg(c("psy", "psych", "R2HTML", "GPArotation"))
  
  ## Catch Warning and Error messages
  myTryCatch <- function(expr) {
    myerror <- NULL
    myWarnings <- NULL
    
    wHandler <- function(w) {
      myWarnings <<- c(myWarnings, list(w))
      invokeRestart("muffleWarning")
    }
    
    val <- withCallingHandlers(tryCatch(expr, 
                                        error = function(e) {myerror <<- e; NULL}), 
                               warning = wHandler)
    
    return(list(warning = myWarnings, error = myerror))
  } 
  
  html.output <- capture.output({
    
    ##### Factor Analysis Start
    R2HTML::HTML(R2HTML::as.title("Factor Analysis"),HR=1,file=stdout(),innerBorder = 1,align="left",append=FALSE) 
    
    
    ##### Warnings
    ## Warning regarding the number of or types of selected variables
    if(length(variables) <= 2) {
      warn.msg1 <- "\a Error : At least 3 variables should be selected. Analysis has been stopped."
    } else {
      isnom <- !sapply(variables,function(i) is.numeric(dataset[,i]))
      if(sum(isnom)>0) {
        nn <- length(variables) - sum(isnom)
        msg <- 'FA supports only numeric variables but categorical variables were selected. They were excluded from the analysis.'
        if(nn > 2) {
          warn.msg2 <- paste0('\a Warning : ',msg)
          variables <- variables[!isnom]
        } else if (nn == 2) {
          warn.msg3 <- paste0('\a Error : ',msg,' After excluding them, there are two numeric variable. FA can be conducted with at least 3 variables. Analysis has been stopped.')
        } else if (nn == 1) {
          warn.msg4 <- paste0('\a Error : ',msg,' After excluding them, there is one numeric variable. Analysis has been stopped.')
        } else     if (n.factors<2) {
          warn.msg7 <-  paste0('\a Error : ',msg,' The number of factors should be more than one. Analysis has stopped.')
        }else {
          warn.msg5 <- paste0('\a Error : ',msg,' After excluding them, there is no variable. Analysis has been stopped.')
        }
      }
    }
    data <- dataset[,variables,drop=F]
    
    ## Warning regarding the number of factors
    if (length(variables) < n.factors) {
      warn.msg6 <- "\a Error : The number of factor variables are less than the number of total variables. Analysis has stopped."
    }

    ## Warning regarding the model
    fit.msg <- suppressMessages(myTryCatch(factanal(data, factors = n.factors, fm = Method, rotate = Rotation, scores = Scores, missing = T)))
    fit.warn <- sapply(fit.msg$warning, function(x) x$message)
    fit.error <- fit.msg$error$message
    
    if (exists('warn.msg1')|exists('warn.msg3')|exists('warn.msg4')|exists('warn.msg5')|exists('warn.msg6')|!is.null(fit.error)) {
      R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout(),innerBorder = 1,align="left") 
      
      if(exists('warn.msg1')) {
        R2HTML::HTML(warn.msg1,file=stdout())
      } else {
        if(exists('warn.msg3')) R2HTML::HTML(warn.msg3, file=stdout())
        if(exists('warn.msg4')) R2HTML::HTML(warn.msg4, file=stdout())
        if(exists('warn.msg5')) R2HTML::HTML(warn.msg5, file=stdout())
      }
      
      if(exists('warn.msg6')) {
        R2HTML::HTML(warn.msg6, file=stdout())
      }
      if(exists('warn.msg7')) {
        R2HTML::HTML(warn.msg7, file=stdout())
      }
      if(!is.null(fit.error)) {
        R2HTML::HTML(paste("\a Error :", fit.error), file=stdout())
      }
      
    } else {
      
      ##### Data structure
      R2HTML::HTML(R2HTML::as.title("Data Structure"),HR=2,file=stdout(),innerBorder = 1,align="left") 
      DS <- matrix(c('Number of observations','Number of variables',nrow(dataset),length(variables)),nrow=2)
      if(exists('warn.msg2')) {
        R2HTML::HTML(DS,file=stdout(),innerBorder = 1,align="left",caption=paste0('<div style="text-align:left">',warn.msg2))
      } else {
        R2HTML::HTML(DS,file=stdout(),innerBorder = 1,align="left")
      }
      
      
      ##### Analysis description	
      R2HTML::HTML(R2HTML::as.title("Analysis Description"),HR=2,file=stdout(),innerBorder = 1,align="left")
      
      methodName <- data.frame(abs = c("minres", "ols", "wls", "gls", "ml"),
                               full = c("Minimum residual", "Ordinary least square", "Weighted least square",
                                        "Generalized weighted least square", "Maximum likelihood"),
                               stringsAsFactors = F)
      
      AD <- matrix(c('Variables',paste(variables, collapse=', '),
                     'Number of factors', n.factors,
                     'Method', methodName[which(methodName$abs == Method), "full"],
                     'Scores', Scores,
                     'Rotation', Rotation),
                   ncol=2, byrow=T)
      R2HTML::HTML(AD,file=stdout(),innerBorder = 1,align="left")
      
      
      ##### Factor analysis
      data<-data[complete.cases(data),]
      fit <- suppressWarnings(factanal(data, factors = n.factors, fm = Method, rotate = Rotation, scores = Scores, missing = T))
      if(fit$dof == 0 ) {
        warn.msg10 <- "\a Error : The degrees of freedom for the model is 0. Analysis has been stopped."
      }
      if(exists('warn.msg10')) {
        R2HTML::HTML(warn.msg10, file=stdout())
      } else {
       O <- fit$scores
      colnames(O) <- paste0("Factor", 1:n.factors)
      
      R2HTML::HTML(R2HTML::as.title("Factor Analysis results"),HR=2,file=stdout(),innerBorder = 1,align="left")
      if (!is.null(fit.warn)) sapply(fit.warn, function(x) R2HTML::HTML(x, file=stdout()))
      
      
      ##### Descriptive stat
      R2HTML::HTML(R2HTML::as.title("Descriptive Statistics"),HR=3,file=stdout(),innerBorder = 1,align="left")
      
      if (D) {
        ## Mean and variance
        R2HTML::HTML(R2HTML::as.title("Univariate descriptives"),HR=4,file=stdout(),innerBorder = 1,align="left")
        d <- data.frame(Mean = colMeans(data), SD = apply(data, 2, sd))
        R2HTML::HTML(d,file=stdout(),innerBorder = 1,align="left",digits=4)
      }
      
      if (CORM) {
        ## Correlation Matrix
        R2HTML::HTML(R2HTML::as.title("Correlation Matrix"),HR=4,file=stdout(),innerBorder = 1,align="left")
        cor.mat <- cor(data, use = "pairwise.complete.obs")
        R2HTML::HTML(cor.mat,file=stdout(),innerBorder = 1,align="left",digits=4)
      }
      
      ## Uniqueness
      R2HTML::HTML(R2HTML::as.title("Uniqueness"),HR=4,file=stdout(),innerBorder = 1,align="left")
      uniq <- data.frame(Uniqueness = fit$uniquenesses)
      R2HTML::HTML(uniq,file=stdout(),innerBorder = 1,align="left",digits=4)
      
      
      ##### Factor analysis results
      ## variable loadings
      R2HTML::HTML(R2HTML::as.title("Variable Loadings"),HR=3,file=stdout(),innerBorder = 1,align="left")
      FL <- loadings(fit)[]
      colnames(FL) <- paste0("Factor", 1:n.factors)
      R2HTML::HTML(FL,file=stdout(),innerBorder = 1,align="left",digits=4)
      
      ## Loading plots
      if(n.factors > 1){
        for(i in 1:(n.factors-1)){
          for(j in (i+1):n.factors){
            R2HTML::HTML(R2HTML::as.title(paste0('Factor',i,' & Factor',j)),HR=4,file=stdout(),innerBorder = 1,align="left")
            xmax <- max(abs(FL[,i]))
            ymax <- max(abs(FL[,j]))
            REx_ANA_PLOT()
            plot(1,type="n",xlim=c(-xmax,xmax),ylim=c(-ymax,ymax),main="Variable Loadings",xlab=paste0('Factor',i),ylab=paste0('Factor',j))
            arrows(rep(0,length(variables)),rep(0,length(variables)),FL[,i],FL[,j],col=rep('red',length(variables)))
            text(x=FL[,i],y=FL[,j],labels=variables)
            REx_ANA_PLOT_OFF("")
          }
        }
      }
      
      ## Summary of loadings
      R2HTML::HTML(R2HTML::as.title("Summary of Loadings"),HR=3,file=stdout(),innerBorder = 1,align="left")
      
      extract.value <- function(vec) {
        trim <- sapply(vec, function(x) gsub("^\\s+|\\s+$", "", x))
        val <- trim[trim != ""]
        
        return(val)
      }
      
      cap.load <- capture.output(loadings(fit))
      
      ss.pos <- grep("SS loadings", cap.load)
      ss.load.raw <- strsplit(cap.load[ss.pos], "SS loadings")
      
      prop.pos <- grep("Proportion Var", cap.load)
      prop.var.raw <- strsplit(cap.load[prop.pos], "Proportion Var")
      
      if(n.factors == 1) {
        ss.load <- as.numeric(extract.value(ss.load.raw))
        prop.var <- as.numeric(extract.value(prop.var.raw))
        cum.var <- prop.var  
      } else if (n.factors > 1) {
        ss.load <- as.numeric(strsplit(extract.value(ss.load.raw), split = " ")[[1]])
        prop.var <- as.numeric(strsplit(extract.value(prop.var.raw), split = " ")[[1]])
        cum.pos <- grep("Cumulative Var", cap.load)
        cum.var.raw <- strsplit(cap.load[cum.pos], "Cumulative Var")
        cum.var <- as.numeric(strsplit(extract.value(cum.var.raw), split = " ")[[1]])
      }
      
      summary.fit <- rbind(ss.pos, prop.var, cum.var)
      rownames(summary.fit) <- c("Sums of Squared Loadings", "Proportion of Variance", "Cumulative Proportion of Variance")
      colnames(summary.fit) <- paste0("Factor", 1:n.factors)
      
      R2HTML::HTML(summary.fit, file=stdout(), innerBorder = 1, row.names = T, align="left",digits=4)
      
      ## Scree plot
      if (SP){
        R2HTML::HTML(R2HTML::as.title("Scree plot"),HR=2,file=stdout(),innerBorder = 1,align="left")
        REx_ANA_PLOT()
        plot(fit$e.values, type = "b", pch = 16, bty = "o", main = "Scree plot", xlab = "Dimension", ylab = "Eigenvalue")
        REx_ANA_PLOT_OFF("")
      }
      
      R2HTML::HTMLhr(file=stdout())
      R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Factor Analysis",sep=""),file=stdout(),innerBorder = 1,align="left")
      R2HTML::HTMLhr(file=stdout())
    }
  }
    }
  )
  
  if(exists('O')){
    return(list(html=html.output,Output=O))
  } else {
    return(html.output)
  }
}

REx_KMC <-function(datsett,conti_varname=NULL,centers,nstart=1L,iter.max=10L,
                             algorithm="Hartigan-Wong",maxcluster=10L,standardize=TRUE,Silhouette=TRUE){
  #datsett : 데이터셋 이름
  #conti_varname=연속형 변수로 선택된 변수들 이름
  #               만약 선택된 변수가 두개 이상이면, c('A','B','C') 이런식으로.
  #               만약 선택된 변수가 하나라면, 'A' 이렇게표현
  #centers : 그림선택기능 - 군집수 (입력된 숫자값.)
  #nstart : 그림선택기능 - 초기값 설정횟수 (입력된 숫자값.) / Default : 1L
  #iter.max : 그림선택기능 - 알고리즘 반복횟수 (입력된 숫자값.) / Default : 10L
  #algorithm : 'Hartigan-Wong', 'Lloyd','Forgy','MacQueen' / default : 'Hartigan-Wong'
  #standardize : 데이터 표준화 / Default : TRUE
  
  ## OPTION부분
  #maxcluster : 그림선택기능-옵션-최대군집수 (입력된 숫자값.) / Default : 10
  #Silhouette : 그림선택기능-옵션-실루엣 / Default : TRUE
  
  #######################################################
  #### Required packages : R2HTML,cluster ,fpc,markdown
  #######################################################
  load.pkg(c("R2HTML", "cluster", "fpc"))
  
  html.output <- capture.output({
    ############################################
    ############# Making data ##################
    ############################################
    
    ## Analysis title
    R2HTML::HTML(R2HTML::as.title('K-means Cluster Analysis'),HR=1,file=stdout(),append=FALSE)
    
    if(is.null(conti_varname)) warn.msg0 <- '\a Any variables are not selected. Analysis has been stopped.'
    if(centers <= 1) warn.msg0.1 <- '\a  The number of centers should be at least 2. Analysis has been stopped.'
    
    if (exists("warn.msg0") | exists("warn.msg0.1")) {
      R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout(),innerBorder = 1,align="left") 
      if(exists("warn.msg0")) R2HTML::HTML(warn.msg0, file=stdout())
      if(exists("warn.msg0.1")) R2HTML::HTML(warn.msg0.1, file=stdout())
    }  else {
      datset<-datsett[,conti_varname]
      if(standardize) datset <- data.frame(scale(datset))
      dat <- na.exclude(datset)
      
      datsett$EXR_numord <- 1:nrow(datsett) #분석하다 데이터 번호가 바뀌는 것을 방지하기 위해 열을 미리 추가해둠
      datset$EXR_numord <- 1:nrow(datset)
      dat$EXR_numord <- 1:nrow(dat)
      
      ####################################################################################  
      ######################## Warning Messages ##########################################
      ####################################################################################
      
      ## Warnings
      # character conti_varname : analysis stopped
      if(!is.null(conti_varname)){
        is.chr <- !sapply(datsett[,conti_varname,drop=F],is.numeric)
        if(any(is.chr)) warn.msg1 <- paste0("\a Analysis has been stopped because non-numeric variable '",paste(conti_varname[is.chr],collapse=", "),"' was included in quantitative variable.")
      }
      
      # Maximum number of clusters
      temp.dat <- datsett[conti_varname]
      num.na <- sum(apply(temp.dat,1,function(i) any(is.na(i))))
      if(Silhouette & maxcluster>min(c(nrow(temp.dat)-num.na,20))) {
        warn.msg2 <- paste0("<div style='text-align:left'> \a Warning message : After removing observations with missing data, we had ",nrow(temp.dat)-num.na," observations. Maximum number of clusters for silhouette was replaced with ",nrow(temp.dat)-num.na,".")
        maxcluster <- min(c(nrow(temp.dat)-num.na,20))
      }
      rm("temp.dat")
      
      # NA elements
      if (nrow(dat) < nrow(datset)) {
        exclude <- sum(!datset$EXR_numord %in% dat$EXR_numord)
        if (exclude == 1) {
          warn.msg3 <- c("<div style='text-align:left'> \a Warning messge : one row has NA and eliminated for analysis.")
        } else warn.msg3 <- paste0("<div style='text-align:left'> \a Warning messge : ", exclude, " rows have NA and eliminated for analysis.")  
      }
      if(exists("warn.msg1")) {
        R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout(),innerBorder = 1,align="left") 
        R2HTML::HTML(warn.msg1, file=stdout())
      } else {
        
        ############################################
        ############# Description ##################
        ############################################
        
        ## Data structure
        R2HTML::HTML(R2HTML::as.title('Data Structure'),HR=2,file=stdout())
        DS <- matrix(c("Number of observations",nrow(dat),
                       "Number of variables",length(conti_varname),
                       "Variables", paste(conti_varname,collapse=", ")), byrow = T, ncol=2)
        if(exists("warn.msg3")) {
          R2HTML::HTML(DS, file=stdout(), innerBorder = 1, align="left", caption = warn.msg3)
        } else R2HTML::HTML(DS, file=stdout(), innerBorder = 1, align="left")
        
        ## Descriptions of Analysis
        R2HTML::HTML(R2HTML::as.title('Descriptions of Analysis'),HR=2,file=stdout())
        DA <- matrix(c("Number of clusters", centers,
                       "Algorithm", algorithm,
                       "Maximum number of interation", iter.max,
                       "Number of random sets", nstart), byrow = T, ncol=2)
        if(Silhouette)	DA <- rbind(DA,matrix(c('Silhouette',Silhouette,
                                               'Maximum number of clusters for Silhouette',maxcluster),ncol=2,byrow=T))
        if(exists("warn.msg2")) {
          R2HTML::HTML(DA,file=stdout(), innerBorder = 1, align="left", caption = warn.msg2)
        } else {
          R2HTML::HTML(DA,file=stdout(), innerBorder = 1, align="left")
        }
        
        
        #####################################################
        ############# Clustering by hclust ##################
        #####################################################
        
        km <- kmeans(x=dat[,-ncol(dat)], centers=centers, nstart=nstart, algorithm=algorithm, iter.max=iter.max) #K-MEANS 수행
        
        R2HTML::HTML(R2HTML::as.title("Results - K-means Clustering"),HR=2,file=stdout()) # k-means 결과(군집간 할당된 샘플수, 제곱합 정보)
        res <- matrix(c(paste0('Sizes of ',centers,' clusters'), paste(km$size,collapse=', '),                         # 각 군집단 할당된 샘플수
                        'Total sum of squares', round(km$totss,3),                                                     # Total sum of squares
                        'Within Cluster SS for each cluster', paste(round(km$withinss,3),collapse=', '),               # Within Cluster sum of squares
                        'Between Cluster SS', round(km$betweenss,3),                                                   # Between Cluster sum of squares
                        'Explanation Ratio (Between SS / Total SS)', paste0(round(km$betweenss/km$totss,5)*100,'%')),  # Explanation Ratio
                      byrow = T, ncol = 2)
        R2HTML::HTML(res, file=stdout(), innerBorder = 1, align="left")
        
        
        #####################################################
        ############# K-means cluster plot ##################
        #####################################################
        R2HTML::HTML(R2HTML::as.title("K-means Cluster Plot"), HR=3, file=stdout()) # k-means 결과(군집간 할당된 샘플수, 제곱합 정보)
        REx_ANA_PLOT()
        plot(dat[,-ncol(dat)],col=km$cluster+1,pch=km$cluster)
        REx_ANA_PLOT_OFF("")
        
        
        ###########################################
        ############# Silhouette ##################
        ###########################################
        if(Silhouette) {
          d<-dist(dat[,-ncol(dat)])
          R2HTML::HTML(R2HTML::as.title("Silhouette"),HR=3,file=stdout())  
          
          #현재에서의 clustering silhouette plot
          REx_ANA_PLOT()
          plot(silhouette(km$cluster, d),main='Current Silhouette',border=NA,col=1:length(unique(km$cluster))) 
          REx_ANA_PLOT_OFF("")
          
          #최적의 clustering number by silhouette
          ks<-2:maxcluster 
          set.seed(1234)
          ASW <- sapply(ks, FUN=function(k) {
            avg_sil<-fpc::cluster.stats(d,kmeans(dat[,-ncol(dat)], centers=k, nstart=25)$cluster)$avg.silwidth
            Dunn<-fpc::cluster.stats(d,kmeans(dat[,-ncol(dat)], centers=k, nstart=25)$cluster)$dunn
            return(c(avg_sil,Dunn))
          })
          
          R2HTML::HTML(R2HTML::as.title("Average Silhouette by cluster numbers"),HR=3,file=stdout()) 
          REx_ANA_PLOT()
          plot(ks, ASW[1,], type="l",main='Optimal cluster numbers',xlab='Cluster',ylab='Average Sihouette')
          abline(v=ks[which.max(ASW[1,])], col="red", lty=2)
          REx_ANA_PLOT_OFF("Red vertical line indicates optimal clustering numbers by average shilhouette.")
          
          ### 최적의 clustering number by Dunn Index
          R2HTML::HTML(R2HTML::as.title("Dunn Index by cluster numbers"),HR=3,file=stdout()) 
          REx_ANA_PLOT()
          plot(ks, ASW[2,], type="l",main='Optimal cluster numbers',xlab='Cluster',ylab='Dunn Index')
          abline(v=ks[which.max(ASW[2,])], col="red", lty=2)
          REx_ANA_PLOT_OFF("Red vertical line indicates optimal clustering numbers by Dunn Index.")
        }
        
        #cluster 결과: 최종데이터에 붙이기.
        datt<-cbind(dat, clustering=km$cluster)
        result_cluster_dat<-merge(datsett,datt[,c('EXR_numord','clustering')],by='EXR_numord',all=T)
        result_cluster_dat$EXR_numord <- NULL
        O <- data.frame(KMcluster = result_cluster_dat$clustering)
      }
    }
    
    # Analysis end time
    R2HTML::HTMLhr(file=stdout())
    R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : K-means Cluster Analysis",sep=""),file=stdout())
    R2HTML::HTMLhr(file=stdout())
  })
  
  if(exists('O')){
    return(list(html=html.output,Output=O))
  } else {
    return(html.output)
  }
}

REx_PAM <-function(datsett,conti_varname=NULL,dummy_varname=NULL,k=2,maxcluster=10L,diss=c('original_data','distance'),metric='euclidean',
                          standardize=FALSE,Silhouette=FALSE){
	#datsett : 데이터셋 이름
	#conti_varname=연속형 변수로 선택된 변수들 이름
	#		만약 선택된 변수가 두개 이상이면, c('A','B','C') 이런식으로.
	#		만약 선택된 변수가 하나라면, 'A' 이렇게표현
	#dummy_varname=범주형 변수로 선택된 변수들 이름
	#		만약 선택된 변수가 두개 이상이면, c('A','B','C') 이런식으로.
	#		만약 선택된 변수가 하나라면, 'A' 이렇게표현
	#k : 그림선택기능 - 군집수 (입력된 숫자값 (자연수만 설정가능))
	#maxcluster : 그림선택기능-옵션-최대군집수 (입력된 숫자값.) / Default : 10
	#diss : 계산 데이터 - 'original_data' 또는 'distance' / 만약 범주형변수 (dummy_varname)으로 할당된 값이 있으면, 무조건 'distance'로만 체크될수있게.
	#metric : 데이터간의 거리계산방법 : 'euclidean','manhattan','maximum','gower'중 하나 선택가능. // 만약 범주형 변수로 할당된 값이 있으면, 무조건 'gower'가 활성화됨.
	#standard : 데이터 표준화 
	## OPTION부분
	#Silhouette : 그림선택기능-옵션-실루엣 

	#######################################################
	#### Required packages : cluster
	#######################################################
	load.pkg(c("cluster", "R2HTML"))
		
	html.output <- capture.output({
    ## Analysis title
    R2HTML::HTML(R2HTML::as.title('Partitioning Around Medoids Clustering'),HR=1,file=stdout(),append=FALSE)
    
    ## Warnings
    # character conti_varname
    if(!is.null(conti_varname)){
      is.chr <- !sapply(datsett[,conti_varname,drop=F],is.numeric)
      if(any(is.chr)) warn.msg1 <- paste0("\a Warning : The type of variable '",paste(conti_varname[is.chr],collapse=', '),"' is not numeric but selected as the quantitative variable. It was excluded from the analysis.")
      conti_varname <- conti_varname[!is.chr]
      if(length(conti_varname)==0) conti_varname <- NULL
    }
    
    # numeric dummy_varname & method
    if(!is.null(dummy_varname)){
      is.num <- sapply(datsett[,dummy_varname,drop=F],is.numeric)
      for(i in 1:length(dummy_varname)) if(is.num[i]) datsett[,dummy_varname[i]] <- as.character(datsett[,dummy_varname[i]])
      if(any(is.num)) {
        waswere1 <- ifelse(sum(is.num)>1,'were','was')
        waswere2 <- ifelse(sum(is.num)>1,'They were','It was')
        warn.msg2 <- paste0("<div style='text-align:left'> \a Warning : Among qualitative variables you selected, the type of variable '",paste(dummy_varname[is.num],collapse=', '),"' ",waswere1," numeric. ",waswere2," coerced into character.")
      }
      if(!metric%in%'gower' & diss=='distance'){
        metric <- 'gower'
        warn.msg3 <- paste0("<div style='text-align:left'> \a Warning : Only 'gower' can be used for categorical variable. 'gower' will be used since categorical variable is selected.")
      }
    }
    
    if(is.null(c(conti_varname,dummy_varname))) {
      warn.msg4 <- '\a Error : At least 1 independent variable should be selected. Analysis has been stopped.'
    } else {
      # Maximum number of clusters
      temp.dat <- datsett[c(conti_varname,dummy_varname)]
      num.na <- sum(apply(temp.dat,1,function(i) any(is.na(i))))
      if(Silhouette & maxcluster>min(c(nrow(temp.dat)-num.na,20))) {
        warn.msg5 <- paste0("<div style='text-align:left'> \a Warning : After removing observations with missing data, we had ",nrow(temp.dat)-num.na," observations. Maximum number of clusters for silhouette was replaced with ",nrow(temp.dat)-num.na,".")
        maxcluster <- min(c(nrow(temp.dat)-num.na,20))
      }
    }
    if(k>min(c(nrow(datsett),20))) warn.msg6 <- '\a Error : Too large number of cluster. It should be less than min(# of observations, 20). Analysis has been stopped.'
    
    if(exists('warn.msg4')|exists('warn.msg6')){
      if(exists('warn.msg4')){
        if(exists('warn.msg1')) R2HTML::HTML(warn.msg1,file=stdout(), innerBorder = 1,align="left")
        R2HTML::HTML(warn.msg4,file=stdout(), innerBorder = 1,align="left")
      }
      if(exists('warn.msg6')) R2HTML::HTML(warn.msg6,file=stdout(), innerBorder = 1,align="left")
    } else {
      ## Data structure
      R2HTML::HTML(R2HTML::as.title('Data Structure'),HR=2,file=stdout())
      DS <- matrix(c("Number of observations","Number of variables",nrow(datsett),length(c(conti_varname,dummy_varname))),ncol=2)
      R2HTML::HTML(DS,file=stdout(), innerBorder = 1,align="left")
      
      ## Variable List
      R2HTML::HTML(R2HTML::as.title('Variable List'),HR=2,file=stdout())
      if(is.null(conti_varname)) {
        varlist <- matrix(c('Qualitative Variable',paste(dummy_varname,collapse=', ')),ncol=2)
      } else {
        varlist <- matrix(c('Quantitative Variable',paste(conti_varname,collapse=', ')),ncol=2)
        if(!is.null(dummy_varname)) varlist <- rbind(varlist,c('Qualitative Variable',paste(dummy_varname,collapse=', ')))
      }
      R2HTML::HTML(varlist,file=stdout(), innerBorder = 1,align="left")
      if(exists('warn.msg1')) R2HTML::HTML(warn.msg1,file=stdout(), innerBorder = 1,align="left")
      if(exists('warn.msg2')) R2HTML::HTML(warn.msg2,file=stdout(), innerBorder = 1,align="left")
      
      ## Descriptions of Analysis
      R2HTML::HTML(R2HTML::as.title('Descriptions of Analysis'),HR=2,file=stdout())
      DA<- matrix(c("Number of clusters",k,"Data type",diss),ncol=2,byrow=T)
      if(diss=='distance') DA <- rbind(DA,c('Distance measure',metric))
      if(Silhouette)	DA <- rbind(DA,matrix(c('Silhouette',Silhouette,'Maximum number of clusters for Silhouette',maxcluster),ncol=2,byrow=T))
      R2HTML::HTML(DA,file=stdout(), innerBorder = 1,align="left")
      if(exists('warn.msg3')) R2HTML::HTML(warn.msg3,file=stdout(), innerBorder = 1,align="left")
      if(exists('warn.msg5')) R2HTML::HTML(warn.msg5,file=stdout(), innerBorder = 1,align="left")
      
      
      ############################################
      #############making data####################
      ############################################
      datsett$EXR_numord<-1:nrow(datsett)
      datt<-datsett[,conti_varname]
      
      if(standardize==TRUE) datt<-data.frame(scale(datt))			
      
      if(length(dummy_varname)>0){
        if(length(dummy_varname)==1){
          dat2<-data.frame(datsett[,dummy_varname])
          colnames(dat2)<-dummy_varname
        }else{
          dat2<-datsett[,dummy_varname]
        }
        for(i in 1:ncol(dat2)){
          if(class(dat2[,i])!='factor'){
            dat2[,i]<-as.factor(dat2[,i])
          } 
        }
        datt<-cbind(datt,dat2)
      }
      
      datt$EXR_numord<-1:nrow(datt)
      datt2<-na.omit(datt) #오류 방지 위해 미리 NA 제거
			
			
			############################################		
			######### Data input 만들기#################
			############################################
			
			if(diss=='original_data'){
        x<-datt2; diss=FALSE
      }
      if(diss=='distance'){
        if(length(dummy_varname)>=1|metric=='gower'){
          D<-daisy(datt2[,-ncol(datt2)],metric='gower',stand=FALSE) 
        }else {
          D<-dist(datt2[,-ncol(datt2)],method=metric)
        }
        x<-D ; diss=TRUE
      }
			
			########################################################################
			###### Recommended number of clustering by Silhouette width ############
			#########################################################################
			R2HTML::HTML(R2HTML::as.title("Results - PAM Clustering"),HR=2,file=stdout()) 

			##############################################################
			################## PAM clustering ############################
			##############################################################
			pam_fitting <- pam(x,diss = diss,k = k)
			
			########################################
			######result############################
			########################################
			R2HTML::HTML(R2HTML::as.title("Results - Medoids"),HR=3,file=stdout()) 
      if(diss) res <- matrix(c('medoids', paste0(paste0(pam_fitting$medoids,'th',collapse=','),' observations are centers of each clustering')),nrow=1)
      if(!diss) res <- pam_fitting$medoids
      R2HTML::HTML(res, file=stdout(), innerBorder = 1, align="left")
		 
			
			#############################################################
			### PLOT은 OPTION과 상관없이 군집분석에서는 항상 출력되도록###
			##############################################################
			R2HTML::HTML(R2HTML::as.title("PAM Clustering Plot 1"),HR=3,file=stdout())
      REx_ANA_PLOT()
      plot(datt2,col=pam_fitting$cluster+1)
      REx_ANA_PLOT_OFF("")
      inform_plot1<-'The plot 1 is plotted based on original data.'
      R2HTML::HTML(inform_plot1,file=stdout(), innerBorder = 1,align="left")
      
      R2HTML::HTML(R2HTML::as.title("PAM Clustering Plot 2"),HR=3,file=stdout())
      REx_ANA_PLOT()
      clusplot(pam_fitting,main=paste0('PAM clustering with ',k,' groups'))
      REx_ANA_PLOT_OFF("")
      inform_plot2<-'The plot 2 is plotted based on distance.'
      R2HTML::HTML(inform_plot2,file=stdout(), innerBorder = 1,align="left")
				

			#############################################
			#######Silhouette Plot ######################
			##############################################			 
			if(Silhouette==TRUE){
        R2HTML::HTML(R2HTML::as.title("Silhouette"),HR=3,file=stdout()) 
        REx_ANA_PLOT()
        plot(silhouette(pam_fitting),border=NA,col=1:length(unique(pam_fitting$clustering)))		 
        REx_ANA_PLOT_OFF("")
        
      }
      
      R2HTML::HTML(R2HTML::as.title("Silhouette Width by Number of Clusters"),HR=3,file=stdout()) 
      sil_width <- c(NA)
      for(i in 2:maxcluster){
        pam_fit <- pam(x,diss = diss,k = i)
        sil_width[i] <- pam_fit$silinfo$avg.width
      }
      kk<-which(sil_width==max(sil_width,na.rm=T))
      
      REx_ANA_PLOT()
      # Plot sihouette width (higher is better)
      plot(1:maxcluster, sil_width, xlab = "Number of clusters", ylab = "Silhouette Width")
      lines(1:maxcluster, sil_width)
      abline(h=k,lty=2,col="red")
      REx_ANA_PLOT_OFF("")
      inform_mention1<-paste0('Since a higher silhouette width is better, the recommended number of clustering groups is ',kk,'. <br> (Our analyis was conducted under number of ',k,' clustering groups.)')
      R2HTML::HTML(inform_mention1,file=stdout(), innerBorder = 1,align="left")
			
			##################################################		
			########cluster 결과: 최종데이터에 붙이기#########
			##################################################
			fdat<-cbind(datt2,clustering=pam_fitting$cluster)
      result_cluster_dat<-merge(datsett,fdat[,c('EXR_numord','clustering')],by='EXR_numord',all=T)
      result_cluster_dat$EXR_numord<-NULL
      O <- data.frame(PAMcluster = result_cluster_dat$clustering)
    }
    
    # Analysis end time
    R2HTML::HTMLhr(file=stdout())
    R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Partitioning Around Medoids Clustering",sep=""),file=stdout())
    R2HTML::HTMLhr(file=stdout())
  })
  
  if(exists('O')){
    return(list(html=html.output,Output=O))
  } else {
    return(html.output)
  }
}

REx_DBSCAN <-function(datsett, conti_varname=NULL, standardize=TRUE, eps=0.15, minPts=5, method = "dist", Silhouette=TRUE) {
	# datsett : 데이터셋 이름
	# conti_varname : 연속형 변수로 선택된 변수들 이름
	#		  만약 선택된 변수가 두개 이상이면, c('A','B','C') 이런식으로.
	#		  만약 선택된 변수가 하나라면, 'A' 이렇게표현
	# standardize : 데이터 표준화 / Default : TRUE
	# eps : 군집 설정거리 : Default : 0.15 (설정가능 범위는 0 초과 1 미만으로)
	# minPts : 군집당 최소 샘플 : Default : 5 ( 설정가능 범위는 1 이상부터 총 데이터 수 (행 개수) 미만으로, 자연수만 선택가능)
	# method : nearest search method. "dist", "kdtree", "linear" 중 선택 / Default : "dist"
	# Silhouette : 그림선택기능-옵션-실루엣 / Default : TRUE
	
	#######################################################
	#### Required packages : dbscan, cluster
	#######################################################
	load.pkg(c("dbscan", "cluster"))
		
	findColorName <- function(colornum) {
    if (colornum == 1) {
		colornm <- "Noise : black" 
    } else colornm <-	paste(paste("Cluster", colornum-1), ":", gsub("[0-9]+", "", palette()[colornum %% 7]))
		return(colornm)
	}
	
	html.output <- capture.output ({
		
		############################################
		############# Making data ##################
		############################################
		
		## Analysis title
		R2HTML::HTML(R2HTML::as.title('DBSCAN Clustering'),HR=1,file=stdout(),append=FALSE)
    
		## Warnings
		# character conti_varname
		if(!is.null(conti_varname)){
		  is.chr <- !sapply(datsett[,conti_varname,drop=F],is.numeric)
		  if(any(is.chr)) warn.msg1 <- paste0("\a Warning : The type of variable '",paste(conti_varname[is.chr],collapse=', '),"' is not numeric but selected as the quantitative variable. It was excluded from the analysis.")
		  conti_varname <- conti_varname[!is.chr]
		  if(length(conti_varname)==0) conti_varname <- NULL
		}
    
		if(is.null(conti_varname)) warn.msg4 <- '\a Error : At least 1 independent variable should be selected. Analysis has been stopped.'
    
		if(exists('warn.msg4')){
		  if(exists('warn.msg1')) R2HTML::HTML(warn.msg1,file=stdout(), innerBorder = 1,align="left")
		  R2HTML::HTML(warn.msg4,file=stdout(), innerBorder = 1,align="left")
		} else {
      
		  datset<-datsett[,conti_varname]
		  if(standardize) datset <- data.frame(scale(datset))
		  dat <- na.exclude(datset)
		  usedat <- dat
			
			datsett$EXR_numord <- 1:nrow(datsett) #분석하다 데이터 번호가 바뀌는 것을 방지하기 위해 열을 미리 추가해둠
			datset$EXR_numord <- rownames(datset)
			dat$EXR_numord <- rownames(dat)
			
			if (minPts > nrow(usedat)) warn.msg2 <- '\a The number of minimum points in the epsilon region should be less than the number of observation.'
			  if(exists("warn.msg2")) {
				R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout(),innerBorder = 1,align="left") 
				R2HTML::HTML(warn.msg2, file=stdout())
			  } else {
				## Warnings
				# NA elements
				if (nrow(dat) < nrow(datset)) {
				  exclude <- sum(!datset$EXR_numord %in% dat$EXR_numord)
				  if (exclude == 1) {
					warn.msg3 <- c("<div style='text-align:left'> \a Warning messge : one row has NA and eliminated for analysis.")
				  } else warn.msg3 <- paste0("<div style='text-align:left'> \a Warning messge : ", exclude, " rows have NA and eliminated for analysis.")	
				}
				
				############################################
				############# Description ##################
				############################################
					
				## Data structure
				R2HTML::HTML(R2HTML::as.title('Data Structure'),HR=2,file=stdout())
				DS <- matrix(c("Number of observations",nrow(usedat),
							   "Number of variables",length(conti_varname),
							   "Variables", paste(conti_varname,collapse=", ")), byrow = T, ncol=2)
        
				R2HTML::HTML(DS,file=stdout(), innerBorder = 1, align="left")
				if(exists('warn.msg3')) R2HTML::HTML(warn.msg3,file=stdout(), innerBorder = 1,align="left")
				if(exists('warn.msg1')) R2HTML::HTML(warn.msg1,file=stdout(), innerBorder = 1,align="left")
        
				## Descriptions of Analysis
				R2HTML::HTML(R2HTML::as.title('Descriptions of Analysis'),HR=2,file=stdout())
				method.table <- data.frame(abs = c("dist", "kdtree", "linear"),
										   full = c("Euclidean Distances", "KD-tree", "Linear"),
										   stringsAsFactors = F)
				DA <- matrix(c("Nearest neighbor search", method.table[method.table$abs == method, "full"],
							   "Epsilon", eps,
							   "Minimum points", minPts), byrow = T, ncol=2)
				if(Silhouette)	DA <- rbind(DA,matrix(c('Silhouette',Silhouette),ncol=2,byrow=T))
				R2HTML::HTML(DA,file=stdout(), innerBorder = 1, align="left")
					
				#######################################################
				###### Inform of K-Nearest Neighbor Distance ##########
				#######################################################
					
				R2HTML::HTML(R2HTML::as.title("Results - DBSCAN Clustering"),HR=2,file=stdout()) 
        
				R2HTML::HTML(R2HTML::as.title("DBSCAN Cluster k-Nearest Neighbor plot"),HR=3,file=stdout()) 
				inform_mention <- paste('\a This plot is plotted under', minPts, 'minimum points setting.',
										'The recommended epsilon value is point that sharp change occurs along the curve.',
										'Our analyis will be conducted under epsilon =', eps, 'setting.')
				REx_ANA_PLOT()
				kNNdistplot(usedat, k=minPts)
				REx_ANA_PLOT_OFF(caption = inform_mention)
				
				#####################################################
				############# Clustering by dbscan ##################
				#####################################################
					
				db <- dbscan::dbscan(x = usedat, eps = eps, minPts = minPts, search = method)
					
				## The number of observations of each cluster
				R2HTML::HTML(R2HTML::as.title("The number of observations of each cluster"), HR=3, file=stdout()) # k-means 결과(군집간 할당된 샘플수, 제곱합 정보)
				cl <- table(db$cluster)
				cl.name <- sapply(names(cl), function(colnm) ifelse (colnm == "0", "Noise", paste0("Cluster", colnm)))
				cl.df <- data.frame(NumberofObs = data.frame(cl)$Freq, row.names = cl.name)
				R2HTML::HTML(cl.df, file=stdout(), innerBorder = 1, align="left", row.names = T)
					
				## Cluster plot
				R2HTML::HTML(R2HTML::as.title("DBSCAN clustering plot"),HR=3,file=stdout()) # k-means 결과(군집간 할당된 샘플수, 제곱합 정보)
				cl1 <- sort(unique(db$cluster+1))
				mylegend <- paste("\a", paste(sapply(cl1, findColorName), collapse = ", "))
        
				REx_ANA_PLOT()
				if(ncol(usedat)==1) {
				  plot(usedat[,1],db$cluster,col=db$cluster+1, yaxt="n",xlab = conti_varname,ylab = 'cluster')
				  axis(side=2, at=unique(db$cluster))
          
				} else {
				  plot(usedat, col=db$cluster+1L)  
				}
				REx_ANA_PLOT_OFF(mylegend)

				#############################################
				####### Silhouette Plot #####################
				#############################################
					
				if(Silhouette==TRUE) {
				  R2HTML::HTML(R2HTML::as.title("Silhouette plot of current clusters"),HR=3,file=stdout()) 
				  silhouette_dbscan <- silhouette(db$cluster, dist(usedat))
          
				  if (any(is.na(silhouette_dbscan))) {
					warn.msg4 <- '\a There are many noises. Silhouette plot is not supported'
					R2HTML::HTML(warn.msg4, file=stdout())
				  } else {
					REx_ANA_PLOT()
					plot(silhouette_dbscan, main = "",border=NA,col=1:length(unique(db$cluster)))
					REx_ANA_PLOT_OFF("")
				  }
				}
					
				##################################################		
				########cluster 결과: 최종데이터에 붙이기#########
				##################################################
				fdat <- cbind(dat, clustering=db$cluster)
				result_cluster_dat <- merge(datsett, fdat[,c('EXR_numord','clustering')], by='EXR_numord', all=T)
				result_cluster_dat$EXR_numord <- NULL
				O <- data.frame(cluster = result_cluster_dat$clustering)
			}
		}
			
		# Analysis end time
		R2HTML::HTMLhr(file=stdout())
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : DBSCAN Clustering",sep=""),file=stdout())
		R2HTML::HTMLhr(file=stdout()) 
	})

	if(exists('O')){
		return(list(html=html.output,Output=O))
	} else {
		return(html.output)
	}
}

REx_ETS <- function(dataset,ts_var,method,lambda=NULL,freq=1,timeplot=FALSE,pred=FALSE,pred.plot=FALSE,conf=NULL,h=10,diag=FALSE,decompplot=FALSE,alpha=NULL,beta=NULL,gamma=NULL,phi=NULL,fitted=FALSE,resid=FALSE){
  
  ######
  ## dataset : dataset 이름
  ## ts_var : 변수선택 탭의 시계열 변수: 결측값이 없는 numeric 변수 하나필수 선택.
  ## method : 분석옵션 탭의 모형선택.
  ## alpha : 분석옵션탭의 Level. 
  ## beta : 분석옵션탭의 Trend 
  ## gamma : 분석옵션탭의 Season
  ## phi : 분석옵션탭의 Damping
  ## freq : 분석옵션탭의 계절길이
  ## lambda : 분석옵션탭의 lambda
  ## timeplot : 출력옵션 탭의 시계열 도표
  ## decompplot : 출력옵션 탭의 분해 도표
  ## diag : 출력옵션 탭의 잔차진단 도표
  ## pred : 출력옵션 탭의 예측값
  ## pred.plot : 출력옵션 탭의 예측도표
  ## conf : 출력옵션 탭의 신뢰구간
  ## h : 출력옵션 탭의 예측시차
  ## fitted : 출력옵션 탭의 적합값
  ## resid : 출력옵션 탭의 잔차
  ######
  ## required packages : forecast, R2HTML
  
  load.pkg(c("forecast", "R2HTML"))
  
  html.output <- capture.output({
    ### HTML Output
    ## Title
    HTML(as.title('Exponential Smoothing Method'),HR=1,file=stdout(),append=FALSE)
    
    ## Warnings
    if (!(is.numeric(dataset[,ts_var]))) warn.msg1 <- "\a Error : Time series variable should be numeric. Analysis has been stopped."
    if (sum(is.na(dataset[,ts_var]))!=0 ) warn.msg2 <- "\a Error : Missing values encountered in Time series variable. Analysis has been stopped."
    data.positive <- (min(dataset[,ts_var]) > 0)
    if(!data.positive & method%in%c("Exponential trend method","Multiplicative damped trend method","Multiplicative Holt-Winters method","Holt-Winters damped method"))
      warn.msg3 <- "\a Error : Inappropriate model for data with negative or zero values. Analysis has been stopped."
    
    # Parameters
    warn.para <- NULL
    if(length(alpha)!=0){if(alpha<0.0001 | alpha>0.9999) warn.para <- paste(warn.para,"'Level' should be in [0.0001,0.9999].")}
    if(length(beta)!=0){if(beta<0.0001 | beta>0.9999) warn.para <- paste(warn.para,"'Trend' should be in [0.0001,0.9999].")}
    if(length(gamma)!=0){if(gamma<0.0001 | gamma>0.9999) warn.para <- paste(warn.para,"'Season' should be in [0.0001,0.9999].")}
    if(length(phi)!=0){if(phi<0.8 | phi>0.98) warn.para <- paste(warn.para,"'Damping' should be in [0.8,0.98].")}
    if(!is.null(freq) & freq>24 ) warn.para <- paste(warn.para,"'Seasonal period' should be an integer ranged from 1 to 24.")
    if(length(alpha)!=0 & length(beta)!=0 ){ if(alpha<beta ) warn.para <- paste(warn.para,"'Trend' should be larger than 'Level'.")}
    if(length(alpha)!=0 & length(gamma)!=0 ){ if(gamma>1-alpha ) warn.para <- paste(warn.para,"'Season' should be less than 1-'Level'.")}
    if(length(lambda)!=0 ){ if(lambda < -1|lambda > 2 ) warn.para <- paste(warn.para,"'lambda' should be in [-1,2].")}
    if(!is.null(warn.para)) warn.msg4 <- paste0('\a Error : Parameters out of range. Analysis has been stopped. (',warn.para,')')
    
    if(exists("warn.msg1")|exists("warn.msg2")|exists("warn.msg3")|exists("warn.msg4")) {
      R2HTML::HTML(R2HTML::as.title('Warnings'),HR=2,file=stdout())
      if(exists("warn.msg1")) R2HTML::HTML(warn.msg1,file=stdout())
      if(exists("warn.msg2")) R2HTML::HTML(warn.msg2,file=stdout())
      if(exists("warn.msg3")) R2HTML::HTML(warn.msg3,file=stdout())
      if(exists("warn.msg4")) R2HTML::HTML(warn.msg4,file=stdout())
    } else {
      ## Data Structure
      R2HTML::HTML(as.title('Data Structure'),HR=2,file=stdout())
      DS <- matrix(c('Number of observations',dim(dataset)[1],'Seasonal period',freq),ncol=2,byrow=T)
      HTML(DS,file=stdout(), innerBorder = 1,align="left")
      
      ## Analysis Description
      HTML(as.title('Analysis Description'),HR=2,file=stdout())
      AD <- matrix(c("Method",method,"Variable name",ts_var),ncol=2,byrow=T)
      if(!is.null(lambda)) AD <- rbind(AD,matrix(c('Box-Cox transformation',ifelse(lambda==0,paste0('f(y)=log(y)'),paste0('f(y)=(y^(',lambda,')-1)/(',lambda,')'))),ncol=2,byrow=T))
      R2HTML::HTML(AD,HR=2,file=stdout(),append=TRUE,innerBorder = 1,align="left") 
      
      ## freq
      y <- ts(dataset[,ts_var],freq=freq)
      
      ##Time series Analysis 
      if(method=="Simple exponential smoothing") model<-ets(y,"ANN",damped=F,alpha=alpha,beta=beta,gamma=gamma,phi=phi,lambda=lambda)
      if(method=="Holt's linear method") model<-ets(y,"AAN",damped=F,alpha=alpha,beta=beta,gamma=gamma,phi=phi,lambda=lambda)
      if(method=="Exponential trend method") model<-ets(y,"MMN",damped=F,alpha=alpha,beta=beta,gamma=gamma,phi=phi,lambda=lambda)
      if(method=="Additive damped trend method") model<-ets(y,"AAN",damped=T,alpha=alpha,beta=beta,gamma=gamma,phi=phi,lambda=lambda)
      if(method=="Multiplicative damped trend method") model<-ets(y,"MMN",damped=T,alpha=alpha,beta=beta,gamma=gamma,phi=phi,lambda=lambda)
      if(method=="Additive Holt-Winters method") model<-ets(y,"AAA",damped=F,alpha=alpha,beta=beta,gamma=gamma,phi=phi,lambda=lambda)
      if(method=="Multiplicative Holt-Winters method") model<-ets(y,"MAM",damped=F,alpha=alpha,beta=beta,gamma=gamma,phi=phi,lambda=lambda)
      if(method=="Holt-Winters damped method") model<-ets(y,"MAM",damped=T,alpha=alpha,beta=beta,gamma=gamma,phi=phi,lambda=lambda)
      
      ### Default output	
      
      HTML(R2HTML::as.title('Estimates of Smoothing Parameters'),HR=3,file=stdout(),append=TRUE) 	
      
      par<-data.frame(t(model$par))
      VL <- matrix(0,nrow=1,ncol=2)
      if(length(par$alpha)!=0) VL <- rbind(VL,c('Level Smooothing Weight (alpha)',round(par$alpha, 4)))
      if(length(par$beta)!=0) VL <- rbind(VL,c('Trend Smoothing Weight (beta)',round(par$beta, 4)))
      if(length(par$gamma)!=0) VL <- rbind(VL,c('Seasonal Smoothing Weight (gamma)',round(par$gamma,4)))
      if(length(par$phi)!=0) VL <- rbind(VL,c('Damped Smooothing Weight (phi)',round(par$phi,4)))
      VL <- VL[-1,,drop=F]
      colnames(VL) <- c('Parameter','Estimate')
      HTML(data.frame(VL),file=stdout(), innerBorder = 1,align="left",row.names=F)
      
      R2HTML::HTML(R2HTML::as.title('Goodness of Fit Measures'),HR=3,file=stdout(),append=TRUE)
      fit1<-matrix(round(c(model$aic,model$aicc,model$bic),4),ncol=3)
      colnames(fit1)<-c("AIC","AICc","BIC"); rownames(fit1) <- 'Criterion'
      HTML(fit1,file=stdout(), innerBorder = 1,align="left")
      
      R2HTML::HTML(R2HTML::as.title('Training Set Accuracy Measures'),HR=3,file=stdout(),append=TRUE)
      fit2<-round(accuracy(model),4)
      rownames(fit2)<-'Accuracy'
      HTML(fit2,file=stdout(), innerBorder = 1,align="left")
      
      ### Option_plot : Time-series Plot
      if(timeplot) {
        R2HTML::HTML(R2HTML::as.title('Time Series Plot'),HR=3,file=stdout(),append=TRUE)
        REx_ANA_PLOT()
        plot(y,ylab=ts_var)
        REx_ANA_PLOT_OFF("")
      }
      
      ### Option_plot : Decomposition Plot
      if(decompplot) {
        plot2<-function (x, ...) {
          if (!is.null(x$lambda)) 
            y <- BoxCox(x$x, x$lambda)
          else y <- x$x
          if (x$components[3] == "N" & x$components[2] == "N") {
            plot(cbind(observed = y, level = x$states[, 1]), main = paste("Decomposition by",x$method, "method"), ...)
          } 
          else if (x$components[3] == "N") {
            plot(cbind(observed = y, level = x$states[, 1], trend = x$states[,"b"]), main = paste("Decomposition by", x$method,"method"), ...)
          }
          else if (x$components[2] == "N") {
            plot(cbind(observed = y, level = x$states[, 1], season = x$states[,"s1"]), main = paste("Decomposition by", x$method, "method"), ...)
          }
          else {
            plot(cbind(observed = y, level = x$states[, 1], trend = x$states[,"b"], season = x$states[, "s1"]), main = paste("Decomposition by", x$method, "method"), ...)
          }
        }
        R2HTML::HTML(R2HTML::as.title('Decomposition Plot'),HR=3,file=stdout(),append=TRUE)
        REx_ANA_PLOT()
        plot2(model)
        REx_ANA_PLOT_OFF("")
      }
      
      ### Option : Prediction
      if(pred){
        step<-1:h
        fp<-forecast(model,h=h,level=conf)
        f.mean<-fp$mean;f.lo<-fp$lower;f.up<-fp$upper
        f.table<-matrix(round(c(f.mean,f.lo,f.up),4),ncol=3)
        f.table<-data.frame(step,f.table)
        colnames(f.table)<-c("Step","Point Forecast",paste("Lo",conf*100),paste("Hi",conf*100))
        rownames(f.table)<-NULL
        R2HTML::HTML(R2HTML::as.title('Forecast Intervals'),HR=3,file=stdout(),append=TRUE)
        HTML(f.table,file=stdout(), innerBorder = 1,align="left",row.names=F)
      }
      
      ### Option_plot : Forecast Plot
      if(pred.plot)
      {
        R2HTML::HTML(R2HTML::as.title('Forecast Plot'),HR=3,file=stdout(),append=TRUE)
        REx_ANA_PLOT()
        plot(fp,ylab=ts_var,xlab="Time")
        REx_ANA_PLOT_OFF("")
      }
      
      ### Option_plot : Diagnostic Plot
      if(diag)
      {
        R2HTML::HTML(R2HTML::as.title('Diagnostic Plots'),HR=3,file=stdout(),append=TRUE)
        REx_ANA_PLOT()
        tsdiag(model)
        REx_ANA_PLOT_OFF("")
      }
      
      if(fitted | resid)
      {
        O<-data.frame(obs=1:length(y))
        if(fitted) O<-data.frame(O,Fitted_ets=as.matrix(model$fitted))
        if(resid) O<-data.frame(O,Resid_ets=as.matrix(model$residuals))
        
      }
    }
    
    ### TIME 
    R2HTML::HTMLhr(file=stdout())
    R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Exponential Smoothing Method",sep=""),file=stdout())
    R2HTML::HTMLhr(file=stdout())
  })
  
  if(exists('O')){
    return(list(html=html.output,Output=O))
  } else {
    return(html.output)
  }
}

REx_MANOVA <- function(dataset,res_var,qual_var=NULL,vars=NULL,noint=FALSE,ss="I",test_mode="Pillai"){

	#### 변수 설명 ####
	## dataset : 데이터셋이름
	## res_var : 종속변수, 필수로 입력되어야 하며 2개 이상 입력가능
	## qual_var : 질적설명변수(independent variables) (반드시 한개 이상 입력. 2개 이상의 변수가 있는 경우 c(,,)로 구분)
	## vars : 최종모형에 선택된 변수들 (2개 이상의 변수가 있는 경우 c(,,)로 구분,교호작용의 경우 :를 이용하여 입력, 필수사항)	#수정
	## noint : 모형에 절편한 포함 여부. TRUE or FALSE.
	## ss: 제곱합 유형 지정. I,II,III 유형 지정가능
	## test_mode : c("Pillai", "Hotelling-Lawley", "Wilks", "Roy") 의 값을 가질 수 있음.
	## ref : http://users.stat.umn.edu/~helwig/notes/mvlr-Notes.pdf, http://ibgwww.colorado.edu/~carey/p7291dir/handouts/manova1.pdf
	###################
	#### Required packages : R2HTML, car, MASS ####
	###################
	load.pkg(c("R2HTML", "car", "MASS"))

	global_op<-options()	#수정
	options(contrasts=c("contr.sum", "contr.poly"))	#수정

	html.output <- capture.output({
		# Title
		R2HTML::HTML(R2HTML::as.title("Multivariate Analysis of Variance"),HR=1,file=stdout(),append=FALSE)

		## Warnings
		if(any(sapply(res_var,function(i) !is.numeric(dataset[,i])))) warn.msg1 <- '\a Error : Response variable should be numeric. Analysis has been stopped.'
		if(is.null(qual_var)) warn.msg4 <- '\a Error : At least 1 group variable should be selected. Analysis has been stopped.'	#수정
		if(!is.null(qual_var)) {
			is.num <- sapply(qual_var,function(i) is.numeric(dataset[,i]))
			if(any(is.num)) warn.msg5 <- paste0("\a Warning : The type of variable '",paste(qual_var[is.num],collapse=', '),"' is numeric but selected as the qualitative variable. It was coreced into character.")
		}

		if(exists('warn.msg1')|exists('warn.msg4')){
			R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
			if(exists('warn.msg1')) R2HTML::HTML(warn.msg1,file=stdout())
			if(exists('warn.msg4')) R2HTML::HTML(warn.msg4,file=stdout())
		
		} else {

			## Model formula
			y<-noquote(paste("cbind(",paste(res_var,collapse=','),")"))
			form.0 <- paste0("cbind(",paste(res_var,collapse=','),") ~ ",paste(vars,collapse=' + '))	#추가수정
			form <- as.formula(paste(form.0,ifelse(noint,'-1',''),sep=''))	#수정
			form.1 <- gsub('cbind','',form.0)
			newdat <- dataset[,c(res_var,qual_var),drop=F]
			newdat <- newdat[complete.cases(newdat),,drop=F]
			if(!is.null(qual_var)) for(i in qual_var) newdat[,i] <- factor(newdat[,i])

			# Data Structure
			R2HTML::HTML(R2HTML::as.title("Data Structure"),HR=2,file=stdout())
			DS <- matrix(c('Number of observations','Number of response variables','Number of group variables',nrow(newdat),length(res_var),length(qual_var)),ncol=2)
			R2HTML::HTML(DS,file=stdout(), innerBorder = 1,align="left")
			if(exists('warn.msg3')) R2HTML::HTML(warn.msg3,file=stdout())

			# Varibale list
			R2HTML::HTML(R2HTML::as.title("Variable List"),HR=2,file=stdout())
			varlist <- matrix(c('Response variables',paste(res_var,collapse=', '),'Group variable',paste(qual_var,collapse=', ')),ncol=2,byrow=T)
			R2HTML::HTML(varlist,file=stdout(), innerBorder = 1,align="left")
			if(exists('warn.msg5')) R2HTML::HTML(warn.msg5,file=stdout())

			# Analysis Description
			R2HTML::HTML(R2HTML::as.title("Analysis Description"),HR=2,file=stdout())
			AD <- matrix(c('Fitted model',ifelse(noint,paste0(gsub(" - 1","",form.1),' (Intercept is not included)'),form.1),'Test statistics',paste(test_mode,collapse=', '),'Sums of Squares',ss),ncol=2,byrow=T)
			R2HTML::HTML(AD,file=stdout(), innerBorder = 1,align="left")

			# Scatter plot matrix
			R2HTML::HTML(R2HTML::as.title("Scatter Plots"),HR=2,file=stdout()) 
			n <- ifelse(length(res_var)<4,4,length(res_var))
			for(i in qual_var){
				REx_ANA_PLOT(w=130*n,h=130*n)
				form.2 <- as.formula(paste0('~',paste(res_var,collapse='+'),'|',i))
				check.ellipse <- function(j) apply(newdat[newdat[,i]==j,res_var],2,function(ii) length(unique(ii))==1)
				if(any(sapply(levels(newdat[,i]),check.ellipse))) {					
					scatterplotMatrix(form.2,data=newdat, smooth=FALSE, reg.line=FALSE, by.groups=TRUE, diagonal="none",main=i)
					cap <- paste0("<div style='text-align:left'> \a Due to constant values in some groups in variable '",i,"' there is no way to draw an ellipse.")
				} else {
					scatterplotMatrix(form.2,data=newdat, smooth=FALSE, reg.line=FALSE, ellipse=TRUE, by.groups=TRUE, diagonal="none",main=i)
					cap <- ""
				}
				REx_ANA_PLOT_OFF(cap)
			}
	
			## Model Fitting
			R2HTML::HTML(R2HTML::as.title("Multivariate ANOVA Results"),HR=2,file=stdout())
			
			# Coefficient estimates
			R2HTML::HTML(R2HTML::as.title("Coefficient Estimates"),HR=3,file=stdout())
			mlm.fit <- lm(form, data=newdat)
			for(i in 1:length(res_var)){
				R2HTML::HTML(R2HTML::as.title(paste0(res_var[i])),HR=4,file=stdout())
				CE <- as.data.frame(summary(mlm.fit)[[i]]$coef)
				CE[,ncol(CE)] <- format(CE[,ncol(CE)],scientific=T,digits=4)
				R2HTML::HTML(CE,file=stdout(), innerBorder = 1,align="left",digits=4)
			}
				
			# ANOVA Table
			R2HTML::HTML(R2HTML::as.title(paste0("Multivariate ANOVA Table with Type ",ss," SS")),HR=3,file=stdout())
			if(ss=='I'){
				res<-try(manova(form, data = newdat),s=T)
				if(class(res)[1]!='try-error'){
					ANOVA <- summary(res,intercept=!noint)$stat[,-2]	#수정
					for(i in test_mode) ANOVA <- cbind(summary(res,intercept=!noint, test=i)$stat[,2,drop=F],ANOVA)	#수정
					if(nrow(ANOVA)==2){ANOVA <- as.data.frame(t(ANOVA[-nrow(ANOVA),]));rownames(ANOVA)<-qual_var}	#수정
					if(nrow(ANOVA)>2){ANOVA <- as.data.frame(ANOVA[-nrow(ANOVA),])}		#수정
					ANOVA[,ncol(ANOVA)] <- format(ANOVA[,ncol(ANOVA)],scientific=T,digits=4)
					R2HTML::HTML(ANOVA,file=stdout(), innerBorder = 1,align="left",digits=4)
				} else {
					warn.msg6 <- "\a Error : Fail to fit the multivariate Analysis of Variance."
					R2HTML::HTML(warn.msg6,file=stdout())
				}
			}

			if(ss %in% c('II','III')){
				fit<-try(lm(form,data=newdat),silent=T)
				if(ss=='II') res<-try(car::Manova(fit,type='II'))
				if(ss=='III'){res<-try(car::Manova(fit,type='III'))}
				if(class(res)[1]!='try-error'){
					outtests <- car:::print.Anova.mlm
					body(outtests)[[16]] <- quote(invisible(tests))
					body(outtests)[[15]] <- NULL
					ANOVA <- outtests(car::Manova(fit, type=ss, test.statistic="Pillai"))[,-2]	
					for(i in test_mode) {
						ANOVA <- cbind(outtests(car::Manova(fit, type=ss, test.statistic=i))[,2,drop=F],ANOVA)
						colnames(ANOVA)[1] <- i
					}
					ANOVA[,ncol(ANOVA)] <- format(ANOVA[,ncol(ANOVA)],scientific=T,digits=4)
					R2HTML::HTML(ANOVA,file=stdout(), innerBorder = 1,align="left",digits=4)
				} else {
					warn.msg7 <- "\a Error : Fail to fit the multivariate Analysis of Variance."
					R2HTML::HTML(warn.msg7,file=stdout())
				}
			}
		}

		R2HTML::HTMLhr(file=stdout())
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Multivariate Analysis of Variance",sep=""),stdout())
		R2HTML::HTMLhr(file=stdout())
	})

	options(global_op)	#수정
	return(html.output)
}

REx_MultiVarReg <- function(dataset,res_var,quan_var=NULL,qual_var=NULL,vars=NULL,noint=FALSE,ss="I",test_mode="Pillai"){

	#### 변수 설명 ####
	## dataset : 데이터셋이름
	## res_var : 종속변수, 필수로 입력되어야 하며 2개 이상 입력가능
	## qual_var : 질적설명변수(independent variables) (2개 이상의 변수가 있는 경우 c(,,)로 구분, 필수사항아님)
	## quan_var : 양적설명변수(independent variables) (2개 이상의 변수가 있는 경우 c(,,)로 구분, 필수사항아님)
	## qual_var와 quan_var 중 적어도 하나는 변수를 입력받아야 함.
	## inter_var : 교호작용 (2개 이상의 변수가 있는 경우 c(,,)로 구분, 필수사항아님)
	## noint : 모형에 절편한 포함 여부. TRUE or FALSE.
	## ss: 제곱합 유형 지정. I,II,III 유형 지정가능
	## test_mode : c("Pillai", "Hotelling-Lawley", "Wilks", "Roy") 의 값을 가질 수 있음.
	## ref : http://users.stat.umn.edu/~helwig/notes/mvlr-Notes.pdf, http://ibgwww.colorado.edu/~carey/p7291dir/handouts/manova1.pdf
	###################
	#### Required packages : R2HTML, car, MASS ####
	###################
	load.pkg(c("R2HTML", "car", "MASS", "Matrix"))

	html.output <- capture.output({
		# Title
		R2HTML::HTML(R2HTML::as.title("Multivariate Linear Regression"),HR=1,file=stdout(),append=FALSE)
		
		inter_var	<- vars[sapply(vars, function(s) grepl(":", s))]

		## Warnings
		if (length(res_var) <= 1) warn.msg0 <- '\a Error : Please select more than 1 response variables. Analysis has been stopped.' ;
		
		if(any(sapply(res_var,function(i) !is.numeric(dataset[,i])))) warn.msg1 <- '\a Error : Response variable should be numeric. Analysis has been stopped.'
		if(!is.null(quan_var)) {
			is.nom <- sapply(quan_var,function(i) !is.numeric(dataset[,i]))
			if(any(is.nom)){
				warn.msg3 <- paste0("\a Warning : The type of variable '",paste(quan_var[is.nom],collapse=', '),"' is not numeric but selected as the quantitative variable. It was excluded from the analysis.")
				# interaction
				ec_inter <- unique(do.call(c,lapply(quan_var[is.nom],function(i) grep(i,inter_var))))
				if(length(ec_inter)!=0) inter_var <- inter_var[-ec_inter]

				# main effect
				quan_var <- quan_var[!is.nom]
			}
		}
		if(is.null(quan_var)&is.null(qual_var)) warn.msg4 <- '\a Error : At least 1 independent variable should be selected. Analysis has been stopped.'
		if(!is.null(qual_var)) {
			is.num <- sapply(qual_var,function(i) is.numeric(dataset[,i]))
			if(any(is.num)) warn.msg5 <- paste0("\a Warning : The type of variable '",paste(qual_var[is.num],collapse=', '),"' is numeric but selected as the qualitative variable. It was coreced into character.")
		}
		
		newdat <- dataset[,c(res_var,quan_var,qual_var),drop=F]
		newdat <- newdat[complete.cases(newdat),,drop=F]
		newdat$inter	<- 1 ;
		if (rankMatrix(newdat) < ncol(newdat)) warn.msg6 <- "\a Error : Linear dependency between columns of the design matrix (including the intercept) detected. Please check the values of indepedent/response variables."

		if(exists('warn.msg0')|exists('warn.msg1')|exists('warn.msg2')|exists('warn.msg4')|exists('warn.msg6')){
			R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
			if(exists('warn.msg0')) R2HTML::HTML(warn.msg0,file=stdout())
			if(exists('warn.msg1')) R2HTML::HTML(warn.msg1,file=stdout())
			if(exists('warn.msg2')) R2HTML::HTML(warn.msg2,file=stdout())
			if(exists('warn.msg4')) {
				if(exists('warn.msg3')) R2HTML::HTML(warn.msg3,file=stdout())
				R2HTML::HTML(warn.msg4,file=stdout())
			}
			if(exists('warn.msg6')) R2HTML::HTML(warn.msg6,file=stdout())
		
		} else {

			## Model formula
			y <- noquote(paste("cbind(",paste(res_var,collapse=','),")"))
			form.0 <- ifelse(is.null(c(quan_var,qual_var)),paste0("cbind(",paste(res_var,collapse=','),") ~ 1"),paste0("cbind(",paste(res_var,collapse=','),") ~ ",paste(c(vars),collapse=' + ')))
			form <- as.formula(form.0) ;
			form.1 <- gsub('cbind','',form.0)
			newdat <- dataset[,c(res_var,quan_var,qual_var),drop=F]
			newdat <- newdat[complete.cases(newdat),,drop=F]
			if(!is.null(qual_var)) for(i in qual_var) newdat[,i] <- factor(newdat[,i])

			# Data Structure
			R2HTML::HTML(R2HTML::as.title("Data Structure"),HR=2,file=stdout())
			DS <- matrix(c('Number of observations','Number of variables',nrow(newdat),length(c(res_var,quan_var,qual_var))),ncol=2)
			R2HTML::HTML(DS,file=stdout(), innerBorder = 1,align="left")
			if(exists('warn.msg3')) R2HTML::HTML(warn.msg3,file=stdout())

			# Varibale list
			R2HTML::HTML(R2HTML::as.title("Variable List"),HR=2,file=stdout())
			varlist <- matrix(c('Response variables',paste(res_var,collapse=', ')),ncol=2)
			if(!is.null(quan_var)) varlist <- rbind(varlist,c('Independent variable (Quantitative)',paste(quan_var,collapse=', ')))
			if(!is.null(qual_var)) varlist <- rbind(varlist,c('Independent variable (Qualitative)',paste(qual_var,collapse=', ')))
			R2HTML::HTML(varlist,file=stdout(), innerBorder = 1,align="left")
			if(exists('warn.msg3')) R2HTML::HTML(warn.msg3,file=stdout())
			if(exists('warn.msg5')) R2HTML::HTML(warn.msg5,file=stdout())

			# Analysis Description
			R2HTML::HTML(R2HTML::as.title("Analysis Description"),HR=2,file=stdout())
			AD <- matrix(c('Fitted model',ifelse(noint,paste0(gsub(" - 1","",form.1),' (Intercept is not included)'),form.1),'Test statistics',paste(test_mode,collapse=', '),'Sums of Squares',ss),ncol=2,byrow=T)
			R2HTML::HTML(AD,file=stdout(), innerBorder = 1,align="left")

			# Scatter plot matrix
			R2HTML::HTML(R2HTML::as.title("Scatter Plots"),HR=2,file=stdout()) 
			n <- length(c(res_var,quan_var))
			if(!is.null(qual_var)){
				for(i in qual_var){
					REx_ANA_PLOT(w=125*n,h=125*n)
					form.2 <- as.formula(paste0('~',paste(c(res_var,quan_var),collapse='+'),'|',i))
					check.ellipse <- function(j) apply(newdat[newdat[,i]==j,c(res_var,quan_var)],2,function(ii) length(unique(ii))==1)
					if(any(sapply(levels(newdat[,i]),check.ellipse))) {					
						scatterplotMatrix(form.2,data=newdat, smooth=FALSE, reg.line=FALSE, by.groups=TRUE, diagonal="none",main=i)
						cap <- paste0("<div style='text-align:left'> \a Due to constant values in some groups in variable '",i,"' there is no way to draw an ellipse.")
					} else {
						scatterplotMatrix(form.2,data=newdat, smooth=FALSE, reg.line=FALSE, ellipse=TRUE, by.groups=TRUE, diagonal="none",main=i)
						cap <- ""
					}
					REx_ANA_PLOT_OFF(cap)
				}
			} else {
				REx_ANA_PLOT(w=130*n,h=130*n)
				form.2 <- as.formula(paste0('~',paste(c(res_var,quan_var),collapse='+')))
				scatterplotMatrix(form.2,data=newdat, smooth=FALSE, reg.line=FALSE, by.groups=TRUE, diagonal="none")
				REx_ANA_PLOT_OFF("")			
			} 
			
			## Model Fitting
			R2HTML::HTML(R2HTML::as.title("Multivariate Linear Regression Results"),HR=2,file=stdout())
			
			# Coefficient estimates
			R2HTML::HTML(R2HTML::as.title("Coefficient Estimates"),HR=3,file=stdout())
			mlm.fit <- lm(form, data=newdat)
			for(i in 1:length(res_var)){
				R2HTML::HTML(R2HTML::as.title(paste0(res_var[i])),HR=4,file=stdout())
				CE <- as.data.frame(summary(mlm.fit)[[i]]$coef)
				CE[,ncol(CE)] <- format(CE[,ncol(CE)],scientific=T,digits=4)
				R2HTML::HTML(CE,file=stdout(), innerBorder = 1,align="left",digits=4)
			}
				
			# ANOVA Table
			R2HTML::HTML(R2HTML::as.title(paste0("Multivariate ANOVA Table with Type ",ss," SS")),HR=3,file=stdout())
			if(ss=='I'){
				res<-try(manova(form, data = newdat),s=T)
				if(class(res)[1]!='try-error'){
					ANOVA <- summary(res)$stat[,-2]
					for(i in test_mode) ANOVA <- cbind(summary(res,test=i)$stat[,2,drop=F],ANOVA)
					ANOVA <- as.data.frame(ANOVA[-nrow(ANOVA),])
					ANOVA[,ncol(ANOVA)] <- format(ANOVA[,ncol(ANOVA)],scientific=T,digits=4)
					R2HTML::HTML(ANOVA,file=stdout(), innerBorder = 1,align="left",digits=4)
				} else {
					warn.msg6 <- "\a Error : Fail to fit the multivariate linear regression."
					R2HTML::HTML(warn.msg6,file=stdout())
				}
			}

			if(ss %in% c('II','III')){
				fit<-try(lm(form,data=newdat),silent=T)
				if(ss=='II') res<-try(car::Manova(fit,type='II'))
				if(ss=='III'){ options(contrasts=c("contr.sum", "contr.poly"))
							res<-try(car::Manova(fit,type='III'))}
				if(class(res)[1]!='try-error'){
					outtests <- car:::print.Anova.mlm
					body(outtests)[[16]] <- quote(invisible(tests))
					body(outtests)[[15]] <- NULL
					ANOVA <- outtests(car::Manova(fit, test.statistic="Pillai"))[,-2]
					for(i in test_mode) {
						ANOVA <- cbind(outtests(car::Manova(fit, test.statistic=i))[,2,drop=F],ANOVA)
						colnames(ANOVA)[1] <- i
					}
					ANOVA[,ncol(ANOVA)] <- format(ANOVA[,ncol(ANOVA)],scientific=T,digits=4)
					R2HTML::HTML(ANOVA,file=stdout(), innerBorder = 1,align="left",digits=4)
				} else {
					warn.msg7 <- "\a Error : Fail to fit the multivariate linear regression."
					R2HTML::HTML(warn.msg7,file=stdout())
				}
			}
		}

		R2HTML::HTMLhr(file=stdout())
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Multivariate Linear Regression",sep=""),stdout())
		R2HTML::HTMLhr(file=stdout())
	})
	return(html.output)
}

### ANOVA

REx_ANOVA	<- function(dataset, res_var, qual_var=NULL, vars=NULL, ss_type=c(2, "2", "II", 1, "1", "I", 3, "3", "III"), 
					posthoc=FALSE, posthoc_var=NULL, posthoc_method=c("Tukey", "LSD", "Scheffe"), p_adj=p.adjust.methods, posthoc_level=0.95,
					fitted_value=FALSE, resid=FALSE, resid_mode=c("original", "standard", "student"), cook_distance=FALSE, num_digit=3) {
	#### Arguments ####
	## res_var : GUI의 '반응변수', 필수로 1개 필요
	## main_var : 주효과 항(2개 이상의 변수가 있는 경우 +로 구분, 최소 1개 필수. 아무변수도 선택되지 않으면 1로 코딩, 양적변수와 중복되어 사용될 수 없음)
	## inter_var : 교호작용 항(2개 이상의 변수가 있는 경우 +로 구분, 필수아님. 아무변수도 선택되지 않으면 코딩하지 않음)
	## no_intercept: exclude the intercept term / GUI의 '상수항 포함하지 않음'과 대응
	# 굳이 넣어야 하나? R에서 분리가 안 되는 모양인데...
	## ss_type: type of F statistics / GUI의 '통계량 유형'과 대응
	## posthoc: post-hoc analysis / GUI의 '사후 분석 수행'과 대응
	## posthoc_var: 사후분석을 수행할 변수 / GUI의 '사후 분석 변수'과 대응
	## posthoc_method: default "Tukey" / GUI의 '분석 방법'과 대응
	## p_adj: default "none" / GUI의 '유의 확률 보정'에 대응
	## posthoc_level:  / GUI의 '유의수준'에 대응
	## fiited_value: print/save fitted values if TRUE / GUI의 '적합값'에 대응
	## resid: print/save residuals if TRUE / GUI의 '잔차'에 대응
	## resid_mode: modes of residuals, can designate multiple ("original", "standard", "student")
	## 	// GUI의 '기본 잔차', '표준화 잔차', '스튜던트화 잔차'에 대응, resid=TRUE이면 반드시 하나 이상 선택되어야 함.
	## cook_distance: Cook'd distance / GUI의 'Cook 거리'에 대응
	## num_digit: number of output digits
	###################

	###################
	#### Required packages : R2HTML, car, AICcmodavg, agricolae, markdown ####
	###################

	###################
	load.pkg(c("R2HTML", "car", "AICcmodavg", "agricolae", "Matrix"))
	
	ss_type	<- as.character(ss_type) ;
			
	getOutput <- FALSE
	
	html.output	<- capture.output({
	
		inter_var	<- vars[sapply(vars, function(s) grepl(":", s))]

		error_msg	<- c() ;
		warn_msg	<- c() ;

		ss_type		<- match.arg(ss_type) ;

		# Title
		R2HTML::HTML(R2HTML::as.title("Analysis of Variance (ANOVA)"),HR=1,file=stdout(),append=FALSE) ;
		
		## Warnings
		# Response variable type
		if(!is.numeric(dataset[,res_var])) {
			error_msg	<- c(error_msg, '\a Error : Response variable should be numeric. Analysis has been stopped.')
		}
		
		# explanatory variable type
		if(is.null(qual_var)) {
			error_msg	<- c(error_msg, '\a Error : At least 1 independent variable should be selected. Analysis has been stopped.') 
		} else {
			for(ql in qual_var) {
				if (length(intersect(c("character", "factor"), class(dataset[[ql]]))) < 1) {
					warn_msg	<- c(warn_msg, paste0("\a Warning : The variable '", ql, "' seems not to be categorical: the variable is coerced to the factor, and the fit may be unstable bacause of too many parameters.")) ;
				}
				dataset[,ql] <- as.factor(dataset[,ql]) ;
			}
		}
		
		if (posthoc) { 
			 if (!grepl(posthoc_var, qual_var)) error_msg	<- c(error_msg, '\a Error : The variable requested for post-hoc analysis is not in the main analysis variable.') 
		
		} ;
		
		formula <- as.formula(paste0(res_var, ' ~ ', paste(c(vars),collapse=' + '))) ;
		if (missing(dataset)) dataset <- environment(formula) ;

		newdat <- dataset[,c(res_var,qual_var),drop=F]
		newdat <- newdat[complete.cases(newdat),,drop=F]
		newdat_2 <- model.matrix(formula, data=newdat)
		newdat_2 <- cbind(newdat_2, newdat[,res_var]) ;
		if (rankMatrix(newdat_2) < ncol(newdat_2)) error_msg 	<- c(error_msg, "\a Error : Linear dependency between columns of the design matrix (including the intercept) detected. Please check the values of indepedent/response variables.")
		
		if(length(error_msg) > 0) {
			R2HTML::HTML(R2HTML::as.title("Error messages"),HR=2,file=stdout())
			for (msg in error_msg) {
				R2HTML::HTML(msg, file=stdout()) ;
			}
		} else {
		
			if (length(warn_msg) > 0) {
				R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
				for (msg in warn_msg) {
					R2HTML::HTML(msg, file=stdout()) ;
				}
			}
		
			fit_aov	<- aov(formula, data=dataset) ;
		
			variables	<- as.character(formula) ;

			# Data structure
			var_list	<- unlist(strsplit(variables[3], c(" \\* | \\+ |:"))) ;
			var_list	<- unique(var_list) ;
			if(is.null(c(qual_var))){
				summ_var	<- 1 ;
				var_lists	<- variables[2] ;
			} else {
				var_list	<- var_list[var_list != "1"] ;
				var_lists	<- c(variables[2], var_list) ;
				summ_var	<- length(var_lists) ;
			}
			data_str	<- matrix(c("Number of observations", "Number of variables", nobs(fit_aov), summ_var), ncol=2) ;
			
			R2HTML::HTML(R2HTML::as.title("Data Structure"), HR=2, file=stdout()) ;
			R2HTML::HTML(data_str, file=stdout(), row.names=FALSE, col.names=TRUE, innerBorder = 1, align="left") ;

			# Model formula
			fit_print	<- matrix(c(variables[2], 
						paste0(variables[2], " = ", variables[3], " + error")
						)) ;
			fit_print	<- cbind(c("Response variable", "Model"), fit_print) ;
			colnames(fit_print)	<- c("Category", "Value") ;
			
			R2HTML::HTML(R2HTML::as.title("Model Information"), HR=2, file=stdout()) ;
			R2HTML::HTML(fit_print, file=stdout(), row.names=FALSE, col.names=TRUE, innerBorder = 1, align="left") ;
			if (exists("warning_msg1")) R2HTML::HTML(warning_msg1, file=stdout()) ;

			# Data structure
			summ_print <- summary(dataset[,var_lists,drop=F]) ;
			summ_print <- data.frame(unclass(summ_print), stringsAsFactors=FALSE) ;
			colnames(summ_print)[1] <- variables[2] ; colnames(summ_print)[-1] <- var_list ;
			colnames(summ_print) <- gsub(" ", "", colnames(summ_print), fixed=TRUE) ;

			R2HTML::HTML(R2HTML::as.title("Summary Statistics"), HR=2, file=stdout()) ;
			R2HTML::HTML(summ_print, file=stdout(),row.names=FALSE,na="", innerBorder = 1, align="left") ;
		
			# Model fitness measures
			fitness_print	<- data.frame(numeric(0), numeric(0)) ;

			fitness_print	<- rbind(fitness_print, c(deviance(fit_aov), fit_aov$df.residual)) ;

			fitness_print	<- rbind(fitness_print, c(as.numeric(logLik(fit_aov)), attr(logLik(fit_aov), "df"))) ;
			
			fitness_print	<- rbind(fitness_print, c(AIC(fit_aov), NA)) ;
			fitness_print	<- rbind(fitness_print, c(AICc(fit_aov), NA)) ;
			fitness_print	<- rbind(fitness_print, c(BIC(fit_aov), NA)) ;
			names(fitness_print)	<- c("Value", "DF") ;
			row.names(fitness_print)	<- c("Deviance", "log-likelihood", "AIC", "AICc", "BIC") ;
			fitness_print[1]	<- format(round(fitness_print[1], num_digit), nsmall=num_digit) ;

			R2HTML::HTML(R2HTML::as.title("Model Fitness Measurements"), HR=2, file=stdout()) ;
			R2HTML::HTML(fitness_print, file=stdout(), na="", innerBorder = 1, align="left") ;
		
			# ANOVA table
			R2HTML::HTML(R2HTML::as.title("ANOVA Table"), HR=2, file=stdout()) ;
			if (ss_type %in% c(1, "1", "I")) {
				table_tmp	<- data.frame(anova(fit_aov)) ;
				colnames(table_tmp)	<- c("DF", "Sum of squares", "Mean squares", "F", "P-value") ;
				
				test_print	<- table_tmp ;
			} else {
				table_tmp	<- data.frame(Anova(fit_aov, type=ss_type)) ;
				table_tmp$Mean.Sq	<- table_tmp$Sum.Sq / table_tmp$Df ;
				table_tmp	<- table_tmp[c(2,1,5,3,4)] ;
				colnames(table_tmp)	<- c("DF", "Sum of squares", "Mean squares", "F", "P-value") ;
					
				if (ss_type %in% c(3, "3", "III")) {
					test_print	<- table_tmp[-1,] ;
				} else {
					test_print	<- table_tmp ;
				}
			}
			
			test_print[-1]	<- format(round(test_print[-1], num_digit), nsmall=num_digit) ;
			R2HTML::HTML(test_print, file=stdout(), na="", innerBorder = 1, align="left") ;
			
			# Post-hoc analysis
			if (posthoc) {
				
				R2HTML::HTML(R2HTML::as.title(paste0("Post-hoc Analysis for variable ", posthoc_var)), HR=2, file=stdout()) ;
			
				# Tukey's
				if (posthoc_method == "Tukey") {
				
					post_table	<- data.frame(TukeyHSD(fit_aov, posthoc_var, conf.level=posthoc_level)[[1]]) ;
					names(post_table)	<- c("Difference", paste0(100*posthoc_level, "% CI ", c("lower", "upper")), "P-value (adjusted)") ;
					
					post_table	<- format(round(post_table, num_digit), nsmall=num_digit) ;
					if (p_adj != "none") warning_msg2	<- "\a Warning : The p-value adjusting method will be ignored." ;
					
				# Scheffe's
				} else if (posthoc_method == "Scheffe") {
				
					post_tmp	<- scheffe.test(fit_aov, posthoc_var, alpha=1-posthoc_level) ;
					post_groups	<- post_tmp$groups ;
					n_allcomb	<- choose(nrow(post_groups), 2) ;
					tab_comb	<- combn(rownames(post_groups), 2) ;
					post_table	<- data.frame(matrix(numeric(n_allcomb*2), ncol=2)) ;
					names_row	<- character(n_allcomb) ;
					
					for (i in seq(n_allcomb)) {
						lev1	<- as.character(tab_comb[2,i]) ;
						lev2	<- as.character(tab_comb[1,i]) ;
						
						names_row[i]	<- paste0(lev1, "-", lev2) ;
						post_table[i,1]	<- post_groups[rownames(post_groups) == lev1, res_var] - post_groups[rownames(post_groups) == lev2, res_var] ;
						if (post_groups$groups[rownames(post_groups) == lev1] == post_groups$groups[rownames(post_groups) == lev2] |
							grepl(post_groups$groups[rownames(post_groups) == lev1], post_groups$groups[rownames(post_groups) == lev2]) |
							grepl(post_groups$groups[rownames(post_groups) == lev2], post_groups$groups[rownames(post_groups) == lev1])) {
							post_table[i,2]	<- "No" ;
						} else {
							post_table[i,2]	<- "Yes" ;
						}
					}
					
					names(post_table)	<- c("Difference", "Significance") ;
					row.names(post_table)	<- names_row ;
					post_table[1]	<- format(round(post_table[1], num_digit), nsmall=num_digit) ;
					
					diff_crit	<- format(round(post_tmp$statistics$CriticalDifference, num_digit), nsmall=num_digit) ;
					notice_msg	<- paste0("\a Notice : groups with the difference larger than ", diff_crit, " is considered to be significantly different.") ;
					if (p_adj != "none") warning_msg2	<- "\a Warning : The p-value adjusting method will be ignored." ;
					
				# Fisher's LSD
				} else if (posthoc_method == "LSD") {
				
					post_tmp	<- LSD.test(fit_aov, posthoc_var, p.adj=p_adj, alpha=1-posthoc_level) ;
					post_groups	<- post_tmp$groups ;
					n_allcomb	<- choose(nrow(post_groups), 2) ;
					tab_comb	<- combn(rownames(post_groups), 2) ;
					post_table	<- data.frame(matrix(numeric(n_allcomb*3), ncol=3)) ;
					names_row	<- character(n_allcomb) ;
					
					pval_tmp	<- pairwise.t.test(dataset[[res_var]], dataset[[posthoc_var]], p.adj=p_adj) ;
					
					for (i in seq(n_allcomb)) {
						lev1	<- as.character(tab_comb[2,i]) ;
						lev2	<- as.character(tab_comb[1,i]) ;
						
						names_row[i]	<- paste0(lev1, "-", lev2) ;
						post_table[i,1]	<- post_groups[rownames(post_groups) == lev1, res_var] - post_groups[rownames(post_groups) == lev2, res_var] ;
						post_table[i,2]	<- ifelse(lev1 %in% rownames(pval_tmp$p.value) & lev2 %in% colnames(pval_tmp$p.value),
												ifelse(is.na(pval_tmp$p.value[lev1, lev2]), pval_tmp$p.value[lev2, lev1], pval_tmp$p.value[lev1, lev2]),
												pval_tmp$p.value[lev2, lev1]) ;
						if (post_groups$groups[rownames(post_groups) == lev1] == post_groups$groups[rownames(post_groups) == lev2] |
							grepl(post_groups$groups[rownames(post_groups) == lev1], post_groups$groups[rownames(post_groups) == lev2]) |
							grepl(post_groups$groups[rownames(post_groups) == lev2], post_groups$groups[rownames(post_groups) == lev1])) {
							post_table[i,3]	<- "No" ;
						} else {
							post_table[i,3]	<- "Yes" ;
						}
					}
					
					names(post_table)	<- c("Difference", ifelse(p_adj == "none", "P-value", "P-value (adjusted)"), "Significance") ;
					row.names(post_table)	<- names_row ;
					post_table[1:2]	<- format(round(post_table[1:2], num_digit), nsmall=num_digit) ;
					
					if ("LSD" %in% names(post_tmp$statistics)) {
						diff_crit	<- format(round(post_tmp$statistics$LSD, num_digit), nsmall=num_digit) ;
						notice_msg	<- paste0("\a Notice : groups with the difference larger than ", post_tmp$statistics$LSD, " is considered to be significantly different.") ;
					}
					if (p_adj != "none") {
						p_adj_print	<- switch(p_adj, 
												holm = "Holm's",
												hochberg = "Hochberg's",
												hommel = "Hommel's",
												bonferroni = "Bonferroni's",
												BH = "Benjamini & Hochberg's", 
												BY = "Benjamini & Yekutieli's", 
												fdr = "Benjamini & Hochberg's") ;
						notice_msg2	<- paste0("\a The significance criterion in adjusted by ", p_adj_print, " method.") ;					
					}
					
				# } else {
				
					# post_tmp	<- pairwise.t.test(dataset[[res_var]], dataset[[posthoc_var]], p.adj=p_adj) ;
					# n_table		<- table(dataset[[posthoc_var]]) ;
					# mean_table	<- aggregate(as.formula(paste0(res_var, "~", posthoc_var)), data=dataset, FUN="mean") ;
					
					# aov_posthoc	<- data.frame(anova(aov(as.formula(paste0(res_var, "~", posthoc_var)), data=dataset))) ;
					
					# mse		<- aov_posthoc["Residuals", "Mean.Sq"] ;
					# df_mse	<- aov_posthoc["Residuals", "Df"] ;
					
					# tt	<- qt((1-posthoc_level)/2, df_mse, lower.tail=FALSE) ;
					
					# post_levels	<- fit_aov$xlevels[[posthoc_var]] ;
					
					# n_allcomb	<- choose(length(post_levels), 2) ;
					# tab_comb	<- combn(post_levels, 2) ;
					# post_table	<- data.frame(matrix(numeric(n_allcomb*4), ncol=4)) ;
					# names_row	<- character(n_allcomb) ;
					
					# for (i in seq(n_allcomb)) {
						# lev1	<- as.character(tab_comb[2,i]) ;
						# lev2	<- as.character(tab_comb[1,i]) ;
						
						# n1	<- n_table[[lev1]] ; n2	<- n_table[[lev2]] ;
						# range_CI	<- tt * sqrt((mse / n1) + (mse / n2)) ;
						
						# names_row[i]	<- paste0(lev1, "-", lev2) ;
						# post_table[i,1]	<- mean_table[[res_var]][mean_table[[posthoc_var]] == lev1] - mean_table[[res_var]][mean_table[[posthoc_var]] == lev2] ;
						# post_table[i,2]	<- post_table[i,1] - range_CI ; post_table[i,3]	<- post_table[i,1] + range_CI ;
						# post_table[i,4]	<- post_tmp$p.value[lev1, lev2] ;
					# }
					# names(post_table)	<- c("Difference", paste0(100*posthoc_level, "% CI ", c("lower", "upper")), ifelse(p_adj == "none", "P-value", "P-value (adjusted)")) ;
				
					# post_table	<- format(round(post_table, num_digit), nsmall=num_digit) ;
					# if (p_adj != "none") {
						# p_adj_print	<- switch(p_adj, 
												# holm = "Holm's",
												# hochberg = "Hochberg's",
												# hommel = "Hommel's",
												# bonferroni = "Bonferroni's",
												# BH = "Benjamini & Hochberg's", 
												# BY = "Benjamini & Yekutieli's", 
												# fdr = "Benjamini & Hochberg's") ;
						# notice_msg2	<- paste0("\a The significance criterion in adjusted by ", p_adj_print, " method.") ;					
					# }
					
				}
				
				post_table	<- cbind(rownames(post_table), post_table) ;
				names(post_table)[1]	<- "Comparison" ;
				
				R2HTML::HTML(post_table, file=stdout(), row.names=FALSE, col.names=TRUE, innerBorder = 1, align="left") ;
				if (exists("warning_msg2")) R2HTML::HTML(warning_msg2, file=stdout()) ;
				if (exists("notice_msg")) R2HTML::HTML(notice_msg, file=stdout()) ;
				if (exists("notice_msg2")) R2HTML::HTML(notice_msg2, file=stdout()) ;
				
				REx_ANA_PLOT()
				form_tmp <- as.formula(paste0(res_var, "~", posthoc_var)) ;
				boxplot(form_tmp, data=dataset, xlab=posthoc_var, ylab=res_var) ;
				stripchart(form_tmp, data=dataset, vertical=TRUE, method = "jitter", add=TRUE, pch=20, col='blue')
				REx_ANA_PLOT_OFF("")
			}

			if (fitted_value) {
				ANOVA_output	<- data.frame(fitted(fit_aov)) ;
				names(ANOVA_output)	<- "ANOVA_Fitted" ;
				getOutput <- TRUE ;
			}
			
			if (resid) {
				if ("original" %in% resid_mode) {
					if (exists("ANOVA_output")) {
						ANOVA_output	<- cbind(ANOVA_output, residuals(fit_aov)) ;
						names(ANOVA_output)[ncol(ANOVA_output)]	<- "ANOVA_ResidOriginal" ;
					} else {
						ANOVA_output	<- data.frame(residuals(fit_aov)) ;
						names(ANOVA_output)	<- "ANOVA_ResidOriginal" ;
					}
				getOutput <- TRUE
				}
				
				if ("standard" %in% resid_mode) {
					if (exists("ANOVA_output")) {
						ANOVA_output	<- cbind(ANOVA_output, rstandard(fit_aov)) ;
						names(ANOVA_output)[ncol(ANOVA_output)]	<- "ANOVA_ResidStandardized" ;
					} else {
						ANOVA_output	<- data.frame(rstandard(fit_aov)) ;
						names(ANOVA_output)	<- "ANOVA_ResidStandardized" ;
					}
				getOutput <- TRUE
				}
				
				if ("student" %in% resid_mode) {
					if (exists("ANOVA_output")) {
						ANOVA_output	<- cbind(ANOVA_output, rstudent(fit_aov)) ;
						names(ANOVA_output)[ncol(ANOVA_output)]	<- "ANOVA_ResidStudentized" ;
					} else {
						ANOVA_output	<- data.frame(rstudent(fit_aov)) ;
						names(ANOVA_output)	<- "ANOVA_ResidStudentized" ;
					}
				getOutput <- TRUE
				}
			}
			
			if (cook_distance) {
				if (exists("ANOVA_output")) {
					ANOVA_output	<- cbind(ANOVA_output, cooks.distance(fit_aov)) ;
					names(ANOVA_output)[ncol(ANOVA_output)]	<- "ANOVA_CookDist" ;
				} else {
					ANOVA_output	<- data.frame(cooks.distance(fit_aov)) ;
					names(ANOVA_output)	<- "ANOVA_CookDist" ;
				}
				getOutput <- TRUE
			}
			
		}
		
		# Analysis end time
		R2HTML::HTMLhr(file=stdout())
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Analysis of Variance (ANOVA).",sep=""),file=stdout())
		R2HTML::HTMLhr(file=stdout())
	
	})
	
	if(getOutput) {
		return(list(html=html.output, Output=ANOVA_output))
	} else {
		return(html.output)
	}
}

REx_Coxph <- function(dataset, time1, event, cov_set1=NULL, cov_set2=NULL, tdp_time=NULL, tdp_cov=NULL, tdp_var=NULL, confi=0.95, Mselec=FALSE, Mselel_option= "both", confi2=0.95, survival_plot=FALSE, a1_survival_plot=FALSE, CHazard_plot=FALSE, log_survival_plot=FALSE) {
	#### 변수 설명 ####
	## time1 : 시간변수
	## confi : confidence interval (디폴트는 0.95)
	## event : 상태변수 (event=1, censoring=0) 
	## cov_set1 : 양적변수 (2개 이상의 변수가 있는 경우 c(,,)로 입력, 필수아님. 질적변수와 중복되어 사용될 수 없음)
	## cov_set2 : 질적변수 (2개 이상의 변수가 있는 경우 c(,,)로 입력, 필수아님. 양적변수와 중복되어 사용될 수 없음) 
	## tdp_cov : 시간종속변수의 기준변수, 2개 이상의 변수가 있는 경우 c(,,)로 입력, 필수아님
	## tdp_time : 시간종속변수의 시간변수, 2개 이상의 변수가 있는 경우 c(,,)로 입력, 필수아님
	## tdp_var : 시간종속변수 새로운 변수이름, 2개 이상의 변수가 있는 경우 c(,,)로 입력, 필수아님
	## Mselec : 출력옵션 탭의 '모델선택' 옵션 
	## Mselel_option : 모델선택에서의 옵션 (디폴트는 both이고 backward나 forward 선택가능)
	## confi2 : 모델선택에서의 confidence interval (디폴트는 0.95)
	## cov_set1, cov_set2, tdp_var 중 적어도 하나는 변수를 입력받아야 함.
	## survival_plot : 출력옵션 탭의 '생존함수' 옵션
	## a1_survival_plot : 출력옵션 탭의 '1-생존함수' 옵션
	## CHazard_plot : 출력옵션 탭의 '누적위험함수' 옵션 
	## log_survival_plot : 출력옵션 탭의 '로그생존함수' 옵션
	###################
	#### Required packages : R2HTML, survival, KMsurv, markdown, MASS
	###################
	load.pkg(c("survival", "KMsurv", "MASS"))

	html.output <- capture.output({
		### Title
		R2HTML::HTML(R2HTML::as.title('Cox Proportional-Hazards Regression'),HR=1,file=stdout(),append=FALSE) 

		### Warnings
		# string time variable : analysis stopped
		Time <- time1
		if(length(tdp_time)!=0) Time <- c(Time,tdp_time)
		is.chr <- !apply(dataset[,Time,drop=F],2,is.numeric)
		if(any(is.chr)) warn.msg1 <- paste0('\a Analysis stopped due to the type of time variable ',paste(Time[is.chr],collapse=', '),'. (Time variable should be numeric.)')

		# invalid status values : analysis stopped
		Status <- event
		if(length(tdp_cov)!=0) Status <- c(Status,tdp_cov)
		is.invalid <- sapply(1:length(Status), function(i) !all(unique(dataset[,Status[i]])%in%c(0,1)))
		if(any(is.invalid)) warn.msg2 <- paste0('\a Analysis stopped due to the invalid status in the event variable ',paste(Status[is.invalid],collapse=', '),'. (Expected value : 0 - censoring, 1 - event)')

		# character cov_set1 : analysis stopped
		if(length(cov_set1)!=0){
			is.chr <- !apply(dataset[,cov_set1,drop=F],2,is.numeric)
			if(any(is.chr)) warn.msg3 <- paste0('\a Analysis stopped because non-numeric variable ',paste(cov_set1[is.chr],collapse=", "),' was included in Numeric variable.')
		}

		# character cov_set2
		if(length(cov_set2)!=0){
			is.num <- apply(dataset[,cov_set2,drop=F],2,is.numeric)
			for(i in 1:length(cov_set2)) if(is.num[i]) dataset[,cov_set2[i]] <- as.character(dataset[,cov_set2[i]])
			if(any(is.num)) warn.msg4 <- paste0("<div style='text-align:left'> Warning message : Among Categorical variables you selected, the type of variable ",paste(cov_set2[is.num],collapse=', ')," was numeric. It was coerced into character.")
		}

		### Warning & analysis stopped
		if(exists("warn.msg1")|exists("warn.msg2")|exists("warn.msg3")) {
			R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
			if(exists("warn.msg1")) R2HTML::HTML(warn.msg1,file=stdout())
			if(exists("warn.msg2")) R2HTML::HTML(warn.msg2,file=stdout())
			if(exists("warn.msg3")) R2HTML::HTML(warn.msg3,file=stdout())
		} else { 
			if(length(tdp_var)!=0){
				# make time-dependent covariates
				cut.point <- lapply(1:length(tdp_var),function(i) unique(with(dataset, dataset[dataset[,tdp_cov[i]]==1,tdp_time[i]])))
				cut.points <- unique(do.call(c,cut.point))
				# cut.points <- unique(with(dataset, dataset[dataset[,tdp_cov]==1,tdp_time]))
				survobj <- with(dataset, Surv(dataset[,time1], dataset[,event], type="right"))
				
				splited <- survSplit(survobj~., zero=0, id="ID", cut=cut.points, data=dataset)

				N <- dim(splited)[1]
				splited$episode_time = 0
				for( i in 1:N){
					 splited$episode_time[i] =splited$survobj[i][[1]] 
				}
				splited.ordered = splited[order(splited$ID),]

				tdp <- lapply(1:length(tdp_var),function(i) with(splited.ordered, as.numeric((episode_time>=splited.ordered[,tdp_time[i]])&(splited.ordered[,tdp_cov[i]]==1))))
				tdp1 <- do.call(cbind,tdp)
				colnames(tdp1) <- tdp_var
				dat <- cbind(splited.ordered,tdp1)
			} else {
				dat <- na.omit(dataset)
			}

			fm <- formula(paste0("Surv(dat[,time1], dat[,event]) ~ ",paste(c(cov_set1,cov_set2,tdp_var),collapse='+')))
			fit <- coxph(fm, data=dat)
			
			rn <- cov_set1
			if(length(cov_set2)!=0) for(i in cov_set2) rn <- c(rn,paste(i,levels(factor(dataset[,i]))[-1],sep=' - '))
			rn <- c(rn,tdp_var)

			a <- with(summary(fit, conf.int=confi), cbind(coefficients)); rownames(a) <- rn
			b <- with(summary(fit, conf.int=confi), cbind(conf.int)); rownames(b) <- rn
			con <- round(unname(summary(fit)$concordance[1]),3)
			con_se <- round(unname(summary(fit)$concordance[2]),3)
			rsq <- round(unname(summary(fit)$rsq[1]),3)
			rsq_max <- round(unname(summary(fit)$rsq[2]),3)
			logtest <- round(unname(summary(fit)$logtest[1]),3); logtest_df <- round(unname(summary(fit)$logtest[2]),3); logtest_pval <- round(unname(summary(fit)$logtest[3]),3)
			waldtest <- round(unname(summary(fit)$waldtest[1]),3); waldtest_df <- round(unname(summary(fit)$waldtest[2]),3); waldtest_pval <- round(unname(summary(fit)$waldtest[3]),3)
			sctest <- round(unname(summary(fit)$sctest[1]),3); sctest_df <- round(unname(summary(fit)$sctest[2]),3); sctest_pval <- round(unname(summary(fit)$sctest[3]),3)
		
			### Description 
			R2HTML::HTML(R2HTML::as.title("Data structure"),HR=2,file=stdout(),append=TRUE) 
			R2HTML::HTML(matrix(c("Number of observations", "Number of event", fit$n, fit$nevent), ncol=2), file=stdout(),innerBorder = 1,align="left")

			### Variable list
			R2HTML::HTML(R2HTML::as.title("Variable List"),HR=2,file=stdout(),append=TRUE) 
			VL <- matrix(c("Time variable", time1, "Status variable",event),ncol=2,byrow=T)
			if(length(cov_set1)!=0) VL <- rbind(VL,c('Numeric variable', paste(cov_set1, collapse=", ")))
			if(length(cov_set2)!=0) VL <- rbind(VL,c('Categorical variable', paste(cov_set2, collapse=", ")))
			if(length(tdp_var)!=0) VL <- rbind(VL,c('Time-dependent variable', paste(tdp_var, collapse=", ")))
			VL <- data.frame(VL); colnames(VL) <- c('Variable type','Variables')
			R2HTML::HTML(VL,file=stdout(),innerBorder = 1,align="left",row.names=F) 

			if(exists("warn.msg4")) R2HTML::HTML(warn.msg4, file=stdout())
			  
			### Default output  
			R2HTML::HTML(R2HTML::as.title("Cox Proportional-Hazards Regression"),HR=2,file=stdout(),append=TRUE) 	
			R2HTML::HTML(a, file=stdout(), innerBorder = 1, align="left",row.names=F,digits=4,
			caption=paste('<div style="text-align:left"> Concordance=',con,' (se=',con_se,') <br> Rsquare=',rsq,' (max possible=',rsq_max,') <br> Likelihood ratio test=',logtest, 'on', logtest_df,' degree of freedom,','p-value=',logtest_pval,
				'<br> Wald test=',logtest, 'on', waldtest_df,' degree of freedom,','p-value=',waldtest_pval,
				'<br> Score (logrank) test=',sctest, 'on', sctest_df,' degree of freedom,','p-value=',sctest_pval))
				#'<br> <br> coef=coefficients <br> exp(coef)=exponential of coefficients <br> se(coef)=standard error of coefficients <br> z=threshold <br>  Pr(>|z|)=p-value'))
		
			R2HTML::HTML(R2HTML::as.title("Confidence Interval of exp(coef)"),HR=2,file=stdout(),append=TRUE) 	
			R2HTML::HTML(b, file=stdout(), innerBorder = 1, align="left",row.names=F,digits=4,
			caption="<div style='text-align:left'> exp(coef)=exponential of coefficients <br> exp(-coef)=exponential of - coefficients <br> lower=lower confidence interval of exp(coef) <br> upper=upper confidence interval of exp(coef)")
		
			### Model selection 
			if(Mselec) {
				mo <- stepAIC(fit, direction=Mselel_option, conf.int=confi2)
				if(class(mo)[1]=="coxph"){
					m <- with(summary(mo, conf.int=confi2), cbind(coefficients))
					m2 <- with(summary(mo, conf.int=confi2), cbind(conf.int))
					con <- round(unname(summary(mo)$concordance[1]),3)
					con_se <- round(unname(summary(mo)$concordance[2]),3)
					rsq <- round(unname(summary(mo)$rsq[1]),3)
					rsq_max <- round(unname(summary(mo)$rsq[2]),3)
					logtest <- round(unname(summary(mo)$logtest[1]),3); logtest_df <- round(unname(summary(mo)$logtest[2]),3); logtest_pval <- round(unname(summary(mo)$logtest[3]),3)
					waldtest <- round(unname(summary(mo)$waldtest[1]),3); waldtest_df <- round(unname(summary(mo)$waldtest[2]),3); waldtest_pval <- round(unname(summary(mo)$waldtest[3]),3)
					sctest <- round(unname(summary(mo)$sctest[1]),3); sctest_df <- round(unname(summary(mo)$sctest[2]),3); sctest_pval <- round(unname(summary(mo)$sctest[3]),3)
					R2HTML::HTML(R2HTML::as.title("Model Selection"),HR=2,file=stdout(),append=TRUE) 
					R2HTML::HTML(R2HTML::as.title("(1) Cox Proportional-Hazards Regression"),HR=3,file=stdout(),append=TRUE) 	
					R2HTML::HTML(m, file=stdout(), innerBorder = 1, align="left",row.names=F,digits=4,
					caption=paste('<div style="text-align:left"> Concordance=',con,' (se=',con_se,') <br> Rsquare=',rsq,' (max possible=',rsq_max,') <br> Likelihood ratio test=',logtest, 'on', logtest_df,' degree of freedom,','p-value=',logtest_pval,
					'<br> Wald test=',logtest, 'on', waldtest_df,' degree of freedom,','p-value=',waldtest_pval,
					'<br> Score (logrank) test=',sctest, 'on', sctest_df,' degree of freedom,','p-value=',sctest_pval,
					'<br> <br> coef=coefficients <br> exp(coef)=exponential of coefficients <br> se(coef)=standard error of coefficients <br> z=threshold <br>  Pr(>|z|)=p-value'))

					R2HTML::HTML(R2HTML::as.title("(2) Confidence Interval of exp(coef)"),HR=3,file=stdout(),append=TRUE) 	
					R2HTML::HTML(m2, file=stdout(), innerBorder = 1, align="left",row.names=F,digits=4,
					caption="<div style='text-align:left'> exp(coef)=exponential of coefficients <br> exp(-coef)=exponential of - coefficients <b> lower=lower confidence interval of exp(coef) <br> upper=upper confidence interval of exp(coef)")
				} else { 	
					R2HTML::HTML(R2HTML::as.title("Model Selection"),HR=2,file=stdout(),append=TRUE) 
					R2HTML::HTML(R2HTML::as.title("'No variables were selected'"),HR=4,file=stdout(),append=TRUE) 

				}
				}


			### survival_plot : Survival function
			if(survival_plot) {
				R2HTML::HTML(R2HTML::as.title("Survival function"),HR=2,file=stdout(),append=TRUE)
				REx_ANA_PLOT()
				plot(survfit(fit, data=dataset), xlab="Time", ylab="Survival probability")
				REx_ANA_PLOT_OFF(paste(paste0(confi*100,"%"), "confidence interval bounded"))		 	
			}

			### a1_survival_plot : a1-Survival function
			if(a1_survival_plot) {
				R2HTML::HTML(R2HTML::as.title("1-Survival function"),HR=2,file=stdout(),append=TRUE)
				REx_ANA_PLOT()
				plot(survfit(fit, data=dataset), xlab="Time", ylab="1-Survival probability", fun="event", ylim=c(0,1))
				REx_ANA_PLOT_OFF(paste(paste0(confi*100,"%"), "confidence interval bounded"))		 	
			}

			### CHazard_plot : Cumulative Hazard function
			if(CHazard_plot) {
				R2HTML::HTML(R2HTML::as.title("Cumulative Hazard function"),HR=2,file=stdout(),append=TRUE)
				REx_ANA_PLOT()
				H.hat <- -log(survfit(fit, data=dataset)$surv); H.hat <- c(H.hat, H.hat[length(H.hat)])
				H.hat[!is.finite(H.hat)] <- NA
				h.sort.of <- survfit(fit, data=dataset)$n.event / survfit(fit, data=dataset)$n.risk
				H.tilde <- vector()
				for(i in 1:length(h.sort.of)) H.tilde[i] <- sum(h.sort.of[1:i])
				H.tilde <- c(H.tilde, H.tilde[length(H.tilde)])
				plot(survfit(fit, data=dataset), xlab="Time", ylab="Cumulative hazard", fun="cumhaz", ylim=range(c(na.omit(H.hat), na.omit(H.tilde))))
				REx_ANA_PLOT_OFF(paste(paste0(confi*100,"%"), "confidence interval bounded"))
			}
				
			### log_survival_plot : Log Survival function
			if(log_survival_plot) {
				R2HTML::HTML(R2HTML::as.title("Log Survival function"),HR=2,file=stdout(),append=TRUE)
				REx_ANA_PLOT()
				plot(survfit(fit, data=dataset), xlab="Time", ylab="Log survival probability", fun="log")
				REx_ANA_PLOT_OFF(paste(paste0(confi*100,"%"), "confidence interval bounded"))
			}
		}
		### TIME 
		R2HTML::HTMLhr(file=stdout())
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Cox Proportional-Hazards Regression",sep=""),stdout())
		R2HTML::HTMLhr(file=stdout())		
	})
	return(html.output)
}

REx_KNN=function(dataset, y, quan_x=NULL, qual_x=NULL, train_percent=0.7, anal_method="classification", k=5, tuning=FALSE, resample_method="cv",resample_repeat=3, new_predict=FALSE, new_quan_x = NULL, new_qual_x = NULL){

	###################
	#### Arguments ####
	###################
	##Require four inputs: 

	# dataset: the data which will be used in analysis.(데이터셋 이름) 
	# y: name of target variable (or true classification).(종속변수) 
	# quan_x: names of quantitative attribute.(설명변수: 질적변수) 
	# qual_x: names of qualitative attribute.(설명변수: 양적변수)
	# train_percent: percentage of training data to divide, default=0.7 (range of 0~1).(전체 데이터중 몇%를 훈련데이터로 사용할지) 
	# k: number of neighbours considered. (이웃 수 선택) 
	# anal_method: the method of model fitting. "classification" or "regression". (모형 적합 방법선택)   
	# tuning: logical, if TRUE, tuning the model by tuning the parameters as k.(모텔튜닝 여부)    
	# resample_method: resampling method with values like "boot","cv","LOOCV" etc. c("boot","cv","LOOCV","repeatedcv").(resampling 방법: "boot","cv","LOOV","repeatedcv"로 나눠짐)
	# resample_repeat: the paramter contians the complete sets of folds to compute for repeaded cross-validation (only for repeatedcv) 
	# new_predict: logical, if TRUE, predict new dataset with fitted model, default is FALSE. (새로운 데이터로 예측 여부) 
	# new_quan_x : (새로 입력받는 값 중 quan_x에 대응되는 새로운 변수)
	# new_qual_x : (새로 입력받는 값 중 qual_x에 대응되는 새로운 변수)


	###################
	#### Packages #####
	###################
	load.pkg(c("caret", "e1071"))

	###############################
	#### Functions in Analysis ####
	###############################
	
	html.output <- capture.output({

		## Generage HTML File
		R2HTML::HTML(R2HTML::as.title("K-Nearest Neighbour Classification"),HR=1,file=stdout(),append=FALSE) 

		## Warnings
		if(!is.null(quan_x)) {
			is.nom <- sapply(quan_x,function(i) !is.numeric(dataset[,i]))
			if(any(is.nom)){
				warn.msg1 <- paste0("\a Warning : The type of variable '",paste(quan_x[is.nom],collapse=', '),"' is not numeric but selected as the quantitative variable. It was excluded from the analysis.")
				quan_x <- quan_x[!is.nom]
				new_quan_x <- new_quan_x[!is.nom]
				if(length(quan_x)==0) quan_x <- new_quan_x <- NULL
			}
		}
		if(is.null(quan_x)&is.null(qual_x)) warn.msg2 <- '\a Error : At least 1 independent variable should be selected. Analysis has been stopped.'
		if(!is.null(qual_x)) {
			is.num <- sapply(qual_x,function(i) is.numeric(dataset[,i]))
			if(any(is.num)){
				warn.msg3 <- paste0("\a Warning : The type of variable '",paste(qual_x[is.num],collapse=', '),"' is numeric but selected as the qualitative variable. It was coreced into character.")
				for(i in qual_x[is.num]) dataset[,i] <- factor(dataset[,i])
				for(i in new_qual_x[is.num]) dataset[,i] <- factor(dataset[,i])
			}
		}

		if(exists('warn.msg2')){
			R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
			if(exists('warn.msg1')) R2HTML::HTML(warn.msg1,file=stdout())
			R2HTML::HTML(warn.msg2,file=stdout())
		} else {
			D <- dataset
			x <- c(quan_x,qual_x)
			dataset <- D[,c(y,x),drop=F]
			dataset <- dataset[complete.cases(dataset),,drop=F]

			if(new_predict){
				new_x <- c(new_quan_x,new_qual_x)
				new_dataset <- D[,c(new_x),drop=F]
				colnames(new_dataset) <- x
			}

			## Generate model formula
			x.vars=paste(x,collapse="+")
			model=as.formula(paste(y,x.vars,sep="~"))

			## Data Slicing (Training & Testing data)
			if(anal_method=="classification"){ 
				dataset[,y] <- factor(dataset[,y])
			}else{ 
				dataset[,y] <- as.numeric(dataset[,y])
			}
			y.val=dataset[,y]

			Wmes="Some classes have a single record"
			tryCatch(
				{idx=createDataPartition(y=y.val,p=train_percent,list=F)},		
				warning=function(cond){ 
					R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
					mes=capture.output(cond)
					warn.msg5="\a WARNING : Some classes of response have a single record. Please check the response variable again." 
					if (grepl(Wmes,mes)) R2HTML::HTML(warn.msg5,file=stdout())		
				}
			)
			
			if(exists("idx")){	
				train=dataset[idx,]
				test=dataset[-idx,]
			
				# Data Structure
				R2HTML::HTML(R2HTML::as.title("Data Structure"),HR=2,file=stdout())
				DS <- matrix(c('Number of observations' ,nrow(D)),ncol=2)
				DS <- rbind(DS,c('Number of observations for model fitting',nrow(dataset)))
				if(new_predict) DS <- rbind(DS,c('Number of observations for prediction',nrow(new_dataset[complete.cases(new_dataset),])))
				DS <- rbind(DS,c('Number of variables',length(x)+1))
				R2HTML::HTML(DS,file=stdout(), innerBorder = 1,align="left")

				# Varibale list
				R2HTML::HTML(R2HTML::as.title("Variable List"),HR=2,file=stdout())
				varlist <- matrix(c('Response variable',y,'Attributes for model fitting',paste(x,collapse=', ')),ncol=2,byrow=T)
				if(new_predict) varlist <- rbind(varlist,c('Attributes for prediction',paste(new_x,collapse=', ')))
				R2HTML::HTML(varlist,file=stdout(), innerBorder = 1,align="left")
				if(exists('warn.msg1')) R2HTML::HTML(warn.msg1,file=stdout())
				if(exists('warn.msg3')) R2HTML::HTML(warn.msg3,file=stdout())

				# Analysis Description
				R2HTML::HTML(R2HTML::as.title("Analysis Description"),HR=2,file=stdout())
				if(anal_method=="classification"){ 
					anl_methd="classification"
				}else{ 
					anl_methd="regression"
				}
				AD <- matrix(c('Fitted model',paste(y,x.vars,sep="~"),'% of data that goes to training',paste0(train_percent*100,'%'),'Analysis method for KNN',anl_methd),ncol=2,byrow=T)
				if(tuning){ 	
					if(resample_method=="boot") methd="Bootstrap"
					if(resample_method=="cv") methd="10-fold cross-validation"
					if(resample_method=="LOOCV") methd="Leave one out cross validation"
					if(resample_method=="repeatedcv") methd="repeated 10-fold cross-validation"
					AD <- rbind(AD,c('Resampling method for KNN',methd),ncol=2,byrow=T)
					if(resample_method=="repeatedcv") AD <- rbind(AD,c('The number of repeated',resample_repeat))
				}
				R2HTML::HTML(AD,file=stdout(), innerBorder = 1,align="left")

				## Training control  
				trctrl= trainControl(method=resample_method,repeats=resample_repeat)

				## Fit model 
			
				if(tuning){ 
				knn_fit=suppressWarnings(train(model, data=train, method="knn", trControl=trctrl, preProcess=c("center","scale")))
				msg <- paste0('\a Accuracy was used to select the optimal model using the largest accuracy. <br> The final value used for the model was K = ',knn_fit$bestTune)
				R2HTML::HTML(R2HTML::as.title(paste("Model Fit with Training Datasets (",methd,")",sep="")),HR=2,file=stdout()) 
				R2HTML::HTML(R2HTML::as.title("Resampling Results Across Tuning Parameters"),HR=3,file=stdout()) 
				R2HTML::HTML(knn_fit$results,file=stdout(),innerBorder = 1,align="left",row.names=F,digits=4)
				R2HTML::HTML(msg,file=stdout())
				}else{ 
					if(anal_method=="classification"){ 
						knn_fit=knn3(model,data=train,k=k)
					}else if(anal_method=="regression"){
						knn_fit=knnreg(model,data=train,k=k)
					}
				}

				## Plot the Accuracy VS K-value 
				R2HTML::HTML(R2HTML::as.title("Plot of Accuracy vs K-value"),HR=3,file=stdout())
				REx_ANA_PLOT()
				#plot(knn_fit)
				if(tuning){
					plot(Accuracy~k, data=knn_fit$result, xlab="#Neighbors",ylab=paste0("Accuracy (",methd,")"),type="o",col="blue",main="Plot of Accuracy vs K-Neighbors")
					grid()
					}
		

				REx_ANA_PLOT_OFF("")

				## Predict with test dataset
				if(nrow(train)<nrow(dataset)){ 
					R2HTML::HTML(R2HTML::as.title("Predicted Result with Test Datasets"),HR=2,file=stdout()) 
		
					if(anal_method=="classification"){
						test_pred=predict(knn_fit, newdata=test,type="class")
					}else{
						test_pred=predict(knn_fit, newdata=test)
					}

					test_confs=confusionMatrix(test_pred,test[,y])
					
					R2HTML::HTML(R2HTML::as.title("Confusion Matrix"),HR=3,file=stdout()) 
					cft <- test_confs$table
					colnames(cft) <- rownames(cft) <- paste0('Class: ',colnames(cft))
					R2HTML::HTML(cft,file=stdout(),innerBorder = 1,align="left",caption='<div style="text-align:left"> Row : Predicted class <br> Column : True class')

					R2HTML::HTML(R2HTML::as.title("Overall Statistics"),HR=3,file=stdout()) 
					ov <- t(as.matrix(test_confs$overall))
					ov <- ov[,c('Accuracy','AccuracyLower','AccuracyUpper','AccuracyNull','AccuracyPValue','Kappa','McnemarPValue'),drop=F]
					colnames(ov) <- c("ACC","Lower bound of 95% CI for Acc","Upper bound of 95% CI for Acc","NIR","P-value [Acc > NIR]","Kappa","P-value for Mcnemar's test")
					R2HTML::HTML(round(ov,4),file=stdout(),innerBorder = 1,align="left",digits=4,row.names=F)
					
					R2HTML::HTML(R2HTML::as.title("Statistics by Class"),HR=3,file=stdout()) 
					R2HTML::HTML(test_confs$byClass,file=stdout(),innerBorder = 1,align="left",digits=4)
				}

				## predict with new dataset 
				if(new_predict){ 	
					tryCatch(
						{
						if(anal_method=="classification"){
						 new_dataset[complete.cases(new_dataset),"KNN_newpred"]=predict(knn_fit, newdata=new_dataset[complete.cases(new_dataset),],type="class")
						}else{
						 new_dataset[complete.cases(new_dataset),"KNN_newpred"]=predict(knn_fit, newdata=new_dataset[complete.cases(new_dataset),])
						}
						 out_dataset=new_dataset[rowSums(is.na(new_dataset))!=ncol(new_dataset),]
						 Ot=data.frame(KNN_newpred=out_dataset[,"KNN_newpred"])},		
						 error=function(cond){ 
							R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
							warn.msg4="\a Error : The variables in prediction should be similar with variables in the model fitting. Please check the variables again."	
							R2HTML::HTML(warn.msg4,file=stdout())
						 },
						 warning=function(cond){ 
							R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
							warn.msg4="\a Error : The variables in prediction should be similar with variables in the model fitting. Please check the variables again."	
							R2HTML::HTML(warn.msg4,file=stdout())
						}
					)
				}
			}else{
				R2HTML::HTML("Please check the input variables again.",file=stdout())
			}
		}
		## Time Track for Analysis 
		R2HTML::HTMLhr(file=stdout()) ;
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : K-Nearest Neighbour Classification",sep=""),file=stdout()) ;
		R2HTML::HTMLhr(file=stdout()) ;
	})

	if(exists('Ot')) {
		return(list(html=html.output,Output=Ot))
	} else {
		return(html.output)
	}
}

REx_Tree = function(dataset, y, quan_x=NULL, qual_x=NULL, train_percent=0.7, method=c("class") , plot_type=c("basic"),cv_plot=FALSE,Prune=FALSE, new_predict=FALSE, new_quan_x = NULL, new_qual_x = NULL) {

	###################
	#### Arguments ####
	###################
	##Require four inputs:
	
	# dataset:  the data which will be used in analysis. (데이터셋 이름)  
	# y: response variable. (종속변수) 
	# quan_x: names of quantitative attribute. (설명변수: 질적변수) 
	# qual_x: names of qualitative attribute. (설명변수: 양적변수) 
	# method: classification, regression or conditional inference (Three options: "class" or "reg"(anova) or "condinf (conditional inference)") (분석방법: classification, regression, conditional inference)  	
	# train_percent: percentage of training data to divide, default=0.7 (range of 0~1).(전체데이터중 몇%를 훈련데이터로 사용할지)
	# Prune : if true, prunes the tree according cross validation results, Purning not suitable for "condinf" condtion. only for "class" and "reg". (Pruning 여부)   
	# plot_type : returns one of two kinds of tree plot, one is "basic" tree plot, the other one is better("extend") tree plot. (tree 그림 종류선택)   
	# cv_plot : if true, returns the table of k-fold cross validation result that is used to find the best tree. only for "class" and "reg". (cv_plot 출력여부: classification 하고 regression 에만 적용)    
	# new_predict: logical, if TRUE, predict new dataset with fitted model, default is FALSE. (새로운 데이터에 fitted된 모형 적용 여부) 
	# new_quan_x : (새로 입력받는 값 중 quan_x에 대응되는 새로운 변수)
	# new_qual_x : (새로 입력받는 값 중 qual_x에 대응되는 새로운 변수)
	
	###################
	#### Packages #####
	###################
	## Required packages : R2HTML, tree, rpart, caret, partykit ##
	load.pkg(c("R2HTML", "tree", "rpart", "caret", "partykit", "markdown", "party"))
	
	###########################
	#### Functions of Tree ####
	###########################
	
		## tree fit 
		Tree.fit= function(model,method,data){
			method=ifelse(method=="reg","anova",method)

			if(method=="class"|method=="anova"){	
				rpart(model, method=method, data=data)  # fit the model
			}else if (method=="conditional inference"){
				ctree(model, data=na.omit(data))
			}else{ 
				cat('Please choose correct tree method:\n"class" or "anova" or "conditional inference"\n')
			}			
		}
		
				
		## pruning tree
		Tree.pfit= function(tree.fit,method){ 
			fit=tree.fit
			method=ifelse(method=="reg","anova",method)
			if(method=="class"|method=="anova"){
				pfit=prune(fit, cp=fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"])			
			}else{ 
				R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
				warn.msg10="\a Error : Please choose correct tree method: Pruning is suitable for 'class' or 'reg'."
				R2HTML::HTML(warn.msg10,file=stdout())
			}
		}


		## Result table : return complexity paramter table
		Tree.result= function(tree.fit,method){
			fit=tree.fit
			method=ifelse(method=="reg","anova",method)
			if(method=="class"){
				return(fit$cptable)   # display result
			}else if(method=="anova"){
				return(fit$cptable)   # display result
			}else{
				R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
				warn.msg11="\a Error : Please choose correct tree method: complexity parameter table is suitable for 'class' or 'reg'."
				R2HTML::HTML(warn.msg11,file=stdout())
			}
		}


		## Cross Validation Plot 
		CV.plot= function(tree.fit,method){ 
			fit=tree.fit
			method=ifelse(method=="reg","anova",method)
			if(method=="class"){
				plotcp(fit)    # visualize cross-validation results
			}else if(method=="anova"){
				par(mfrow=c(1,2))
				rsq.rpart(fit)
			}else{
				cat('Please choose correct tree method:\nPruning suitable for "class" or "anova"\n')
			}
		}
			

		###Tree Plot function 
		## plot function for rpart result 
		plot.rpart.tree <- function(fit.model, fit.title="Decision Tree", font.size = 0.8) {
		    # plot decision tree
			title=as.character(fit.title)
		    plot(fit.model,
			 uniform   = T,    # if 'TRUE', uniform vertical spacing of the nodes is used
			 branch    = 1,    # controls the shape of the branches from parent to child node
			 compress  = F,    # if 'FALSE', the leaf nodes will be at the horizontal plot
			 nspace    = 0.1,
			 margin    = 0.1, # an extra fraction of white space to leave around the borders
			 minbranch = 0.3,  # set the minimum length for a branch
			 main      = title )

		    # Add text
		    text(x      = fit.model,   #
			 splits = T,           # If tree are labeled with the criterion for the split
			 all    = T,           # If 'TRUE', all nodes are labeled, otherwise just terminal nodes
			 use.n  = T,           # Use numbers to annotate
			 cex    = font.size)   # Font size
		}


		## plot function for party result
		plot.party.tree <- function(fit.model,fit.title="Decision Tree") { 
			
			title=as.character(fit.title)
			if (any(class(fit.model)=="rpart")){ 
				party.tree  = as.party(fit.model)  # change to party like plot
			} else { 
				party.tree  = fit.model
			}

			plot(party.tree,main   = title)               # set main title 
		}


		## Combination of rpart and party plot function 
		Tree.plot=function(tree.fit,method,fit.title,plot.type="basic",font.size=0.8){ 
			method=ifelse(method=="reg","anova",method)
			if(plot.type=="basic"){
				if(method=="class"|method=="anova"){
					plot.rpart.tree(tree.fit,fit.title,font.size)
				}else{ 
					cat('Please choose correct tree method:\n"basic" plot type suitable for "class" or "anova"\n')
				}
			}else if(plot.type=="extend"){ 
				plot.party.tree(tree.fit,fit.title)
			}
		}


	html.output <- capture.output({ 

		## Generage HTML File
		R2HTML::HTML(R2HTML::as.title("Decision Tree Analysis"),HR=1,file=stdout(),append=FALSE) 

		## Warnings
		if(!is.null(quan_x)) {
			is.nom <- sapply(quan_x,function(i) !is.numeric(dataset[,i]))
			if(any(is.nom)){
				warn.msg1 <- paste0("\a Warning : The type of variable '",paste(quan_x[is.nom],collapse=', '),"' is not numeric but selected as the quantitative variable. It was excluded from the analysis.")
				quan_x <- quan_x[!is.nom]
				new_quan_x <- new_quan_x[!is.nom]
				if(length(quan_x)==0) quan_x <- new_quan_x <- NULL
			}
		}
		if(is.null(quan_x)&is.null(qual_x)) warn.msg2 <- '\a Error : At least 1 independent variable should be selected. Analysis has been stopped.'
		if(!is.null(qual_x)) {
			is.num <- sapply(qual_x,function(i) is.numeric(dataset[,i]))
			if(any(is.num)){
				warn.msg3 <- paste0("\a Warning : The type of variable '",paste(qual_x[is.num],collapse=', '),"' is numeric but selected as the qualitative variable. It was coreced into character.")
				for(i in qual_x[is.num]) dataset[,i] <- factor(dataset[,i])
				for(i in new_qual_x[is.num]) dataset[,i] <- factor(dataset[,i])
			}
		}
		if(Prune & method=="condinf") warn.msg4 <- "\a Error : Please choose correct tree method. Pruning is not available for 'conditional inference'. Analysis has been stopped."
		if(method=="reg" & !is.numeric(dataset[,y])) warn.msg5 <- "\a Error : Inappropriate type for response variable. (Regression tree is supported for only numeric variable.) Analysis has been stopped."

		if(exists('warn.msg2')|exists('warn.msg4')|exists('warn.msg5')){
			R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
			if(exists('warn.msg2')) {
				if(exists('warn.msg1')) R2HTML::HTML(warn.msg1,file=stdout())
				R2HTML::HTML(warn.msg2,file=stdout())
			}
			if(exists('warn.msg4')) R2HTML::HTML(warn.msg4,file=stdout())
			if(exists('warn.msg5')) R2HTML::HTML(warn.msg5,file=stdout())
		} else {
			D <- dataset
			x <- c(quan_x,qual_x)
			dataset <- D[,c(y,x),drop=F]
			dataset <- dataset[complete.cases(dataset),,drop=F]

			if(new_predict){
				new_x <- c(new_quan_x,new_qual_x)
				new_dataset <- D[,c(new_x),drop=F]
				colnames(new_dataset) <- x
			}


			## Generate model formula
			x.vars=paste(x,collapse="+")
			model=as.formula(paste(y,x.vars,sep="~"))

			## Data Slicing (Training & Testing data)
			dataset[,y] <- factor(dataset[,y])
			y.val=dataset[,y]
			idx=createDataPartition(y=y.val,p=train_percent,list=F) 
			train=dataset[idx,]
			test=dataset[-idx,]
			
			## Data Structure
			R2HTML::HTML(R2HTML::as.title("Data Structure"),HR=2,file=stdout())
			DS <- matrix(c('Number of observations' ,nrow(D)),ncol=2)
			DS <- rbind(DS,c('Number of observations for model fitting',nrow(dataset)))
			if(new_predict) DS <- rbind(DS,c('Number of observations for prediction',nrow(new_dataset[complete.cases(new_dataset),])))
			DS <- rbind(DS,c('Number of variables',length(x)+1))	
			R2HTML::HTML(DS,file=stdout(), innerBorder = 1,align="left",caption="<div style='text-align:left'> Observations with missing were discarded.")

			## Varibale list
			R2HTML::HTML(R2HTML::as.title("Variable List"),HR=2,file=stdout())
			varlist <- matrix(c('Response variable',y,'Attributes for model fitting',paste(x,collapse=', ')),ncol=2,byrow=T)
			if(new_predict) varlist <- rbind(varlist,c('Attributes for prediction',paste(new_x,collapse=', ')))
			R2HTML::HTML(varlist,file=stdout(), innerBorder = 1,align="left")
			if(exists('warn.msg1')) R2HTML::HTML(warn.msg1,file=stdout())
			if(exists('warn.msg3')) R2HTML::HTML(warn.msg3,file=stdout())

			## Analysis Description
			R2HTML::HTML(R2HTML::as.title("Analysis Description"),HR=2,file=stdout())
			method_print <- switch(method, class="Classification", reg="Regression", condinf="Conditional Inference")
			AD <- matrix(c('Fitted model',paste(y,x.vars,sep="~"),'% of data that goes to training',paste0(train_percent*100,'%'),'Method',method_print, 'Tree pruning',Prune),ncol=2,byrow=T)
			R2HTML::HTML(AD,file=stdout(), innerBorder = 1,align="left")

			## Fit model with trainning dataset 
			method=ifelse(method=="reg","anova",method)
			if(method=="class"|method=="anova"){	
				fit=rpart(model, method=method, data=train)  # fit the model
			}else if (method=="condinf"){
				fit=ctree(model, data=na.omit(train))
			}else{ 
				cat('Please choose correct tree method:\n"class" or "anova" or "conditional inference"\n')
			} 
			R2HTML::HTML(R2HTML::as.title(paste("Model Fit with Training Datasets (",method_print,")",sep="")),HR=2,file=stdout()) 
			REx_ANA_PLOT()
			par(mar=c(0,0,5,0))
			Emes1=c("'x' is a list, but does not have components 'x' and 'y'")
			Emes2=c("fit is not a tree, just a root")
			tryCatch(
				{Tree.plot(fit,method,"Decision Tree before Pruning",plot_type,0.8)},		
				error=function(cond){ 
					R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
					mes=capture.output(cond)
					if(grepl(Emes1,mes)|grepl(Emes2,mes)){
						warn.msg6="\a Error : The fitted model is not a tree, just a root.Please select variables properly."	
						R2HTML::HTML(warn.msg6,file=stdout())
					}else{
						warn.msg7="\a Error : Please select variables properly."
						R2HTML::HTML(warn.msg7,file=stdout())
					}
					
				},
				warning=function(cond){ 
						R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
						warn.msg8="\a WARNING : Please select variables properly."
						R2HTML::HTML(warn.msg8,file=stdout())
					
				}
			)
			REx_ANA_PLOT_OFF("")

			## Complexity Parameter Table	
			if(method=="class"|method=="anova"){
				res=as.data.frame(Tree.result(fit,method))
				R2HTML::HTML(R2HTML::as.title("Complexity Parameter Table of Fitted Tree"),HR=3,file=stdout()) 
				R2HTML::HTML(res,file=stdout(),innerBorder = 1,align="left",caption='<div style="text-align:left"> CP : Complex parameter <br> nsplit : Number of split <br> rel_error : relative error <br> xerror : cross-validated error <br> xstd : Standard deviation of cross-validation')
			}

			## Cross Validation Plot
			if(method=="class"|method=="anova"){
				if (cv_plot) {
					R2HTML::HTML(R2HTML::as.title("Cross Validation Plot of Fitted Tree"),HR=3,file=stdout()) # Result Table
					REx_ANA_PLOT()
					capture.output(CV.plot(fit,method))
					REx_ANA_PLOT_OFF("")
				}
			}

			## Tree Pruning
			if(method=="class"|method=="anova"){
				if (Prune) { 
					R2HTML::HTML(R2HTML::as.title("Analysis Result after Pruning"),HR=2,file=stdout())

					## Tree Plot of Fitted Model
					pfit_tree=Tree.pfit(fit,method)
					R2HTML::HTML(R2HTML::as.title("Tree Plot of Pruned Model"),HR=3,file=stdout())
					REx_ANA_PLOT()
					par(mar=c(0,0,5,0))
					Tree.plot(pfit_tree,method,"Decision Tree after Pruning",plot_type,0.8)
					REx_ANA_PLOT_OFF("")

					## Complexity Parameter Table	
					pres=as.data.frame(Tree.result(pfit_tree,method))
					R2HTML::HTML(R2HTML::as.title("Complexity Parameter Table of Purned Tree"),HR=3,file=stdout()) # Complexity Parameter Table
					R2HTML::HTML(pres,file=stdout(),innerBorder = 1,align="left",row.names=F)
					
					## Cross Validation Plot
					R2HTML::HTML(R2HTML::as.title("Cross Validation Plot of Purned Tree"),HR=3,file=stdout()) 
					REx_ANA_PLOT()
					capture.output(CV.plot(pfit_tree,method))
					REx_ANA_PLOT_OFF("")
					
				}
			}else{ 
				if(Prune){ 
					cat('Please choose correct tree method:\nPruning is not available for "conditional inference" \n')
				}
			}
				
			## Predict with test dataset
			if(method=="class"){
				type="class"
			}else{
				type=NULL
			}

			if(nrow(train)<nrow(dataset)){ 
				R2HTML::HTML(R2HTML::as.title("Predicted Result with Test Datasets"),HR=2,file=stdout()) 

				if(Prune==FALSE) test_pred=predict(fit,newdata=test,type=type)
				
				if(method=="class"|method=="anova"){
					if(Prune) test_pred=predict(pfit_tree, newdata=test, type=type)	
				}else{
					if(Prune) cat('Please choose correct tree method:\nPruning is not available for "conditional inference"\n')
				}

				if(method=="class"){
					test_confs=confusionMatrix(test_pred,test[,y])
					R2HTML::HTML(R2HTML::as.title("Confusion Matrix"),HR=3,file=stdout()) 
					cft <- test_confs$table
					colnames(cft) <- rownames(cft) <- paste0('Class: ',colnames(cft))
					R2HTML::HTML(cft,file=stdout(),innerBorder = 1,align="left",caption='<div style="text-align:left"> Row : Predicted class <br> Column : True class')
				
					R2HTML::HTML(R2HTML::as.title("Overall Statistics"),HR=3,file=stdout()) 
					ov <- t(as.matrix(test_confs$overall))
					ov <- ov[,c('Accuracy','AccuracyLower','AccuracyUpper','AccuracyNull','AccuracyPValue','Kappa','McnemarPValue'),drop=F]
					colnames(ov) <- c("ACC","Lower bound of 95% CI for Acc","Upper bound of 95% CI for Acc","NIR","P-value [Acc > NIR]","Kappa","P-value for Mcnemar's test")
					R2HTML::HTML(round(ov,4),file=stdout(),innerBorder = 1,align="left",digits=4,row.names=F)

					R2HTML::HTML(R2HTML::as.title("Statistics by Class"),HR=3,file=stdout()) 
					R2HTML::HTML(test_confs$byClass,file=stdout(),innerBorder = 1,align="left",digits=4)
				}else{
					REx_ANA_PLOT()
					if(method=="anova"){
						plot(test_pred,test[,y],xlab="Predicted",ylab="Observed",main="Comparison between predicted value and observed value")
						abline(0,1)
					}else{
						plot(test_pred,test[,y],xlab="Predicted",ylab="Observed",main="Comparison between predicted value and observed value")
					}
					REx_ANA_PLOT_OFF("")
				}
			}
			
						

			## predict with new dataset
			if(new_predict){ 
				tryCatch(
					{if(Prune==FALSE){ 
						new_dataset[complete.cases(new_dataset),"Tree_newpred"]=predict(fit,newdata=new_dataset[complete.cases(new_dataset),],type=type)
					}
					if(method=="class"|method=="anova"){
						if(Prune==TRUE){ 
							new_dataset[complete.cases(new_dataset),"Tree_newpred"]=predict(pfit_tree,newdata=new_dataset[complete.cases(new_dataset),],type=type)
						}
					}else{
						if(Prune==TRUE){ 
							warn.msg4 <- 'Please choose correct tree method:\nPruning is not available for "conditional inference"'
						}
					}
					out_dataset=new_dataset[rowSums(is.na(new_dataset))!=ncol(new_dataset),]
					Ot=data.frame(Tree_newpred=out_dataset[,"Tree_newpred"])},		
					error=function(cond){ 
						R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
						warn.msg9="\a Error : The variables in prediction should be similar with variables in the model fitting. Please check the variables again."	
						R2HTML::HTML(warn.msg9,file=stdout())
					},
					warning=function(cond){ 
						R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
						warn.msg9="\a Error : The variables in prediction should be similar with variables in the model fitting. Please check the variables again."	
						R2HTML::HTML(warn.msg9,file=stdout())
					}
				)
			}
		}
		## Time Track for Analysis 
		R2HTML::HTMLhr(file=stdout()) ;
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Decision Tree Analysis",sep=""),file=stdout()) ;
		R2HTML::HTMLhr(file=stdout()) ;
	})

	if(exists('Ot')) {
		return(list(html=html.output,Output=Ot))
	} else {
		return(html.output)
	}
}

REx_SVM=function(dataset, y, quan_x=NULL, qual_x=NULL, scale= FALSE, train_percent=0.7, type="C-classification", kernel="radial", degree=NULL, gamma=NULL, cost=1, new_predict=FALSE, new_quan_x=NULL, new_qual_x=NULL){


	###################
	#### Arguments ####
	###################
	##Require four inputs: 
	
	# dataset:  the data which will be used in analysis. (데이터셋 이름) 
	# y: name of quantitative response variable. (종속변수:회귀분석시 연속형 변수 입력,분류분석시 범주형 변수 입력, 1개만 입력 가능)
	# quan_x: names of quantitative attribute. (설명변수: 질적변수) 
	# qual_x: names of qualitative attribute. (설명변수: 양적변수)
	# scale: logical, if FAlSE, tells the svm funciton not to scale each feature to have mean zero or standard deviation one. (설명젼수 scale(평균0,표준편차1)여부)  
	# type: svm can be used as a classsification or regression machine. y is factor, run "C-classification", y is not a factor then run "eps-regression". (분석 type: 종속변수 이산형이면 classsification, 종속변수 연속형이면 regression.) 
	# kernel: the kernel used in training and predicting, default is radial. Available options : "liner", "polynomial", "radial", "sigmoid". (SVM 분석 kernel 종류선택)
	# train_percent: percentage of training data to divide, default=0.7 (range of 0~1). (전체데이터중 몇%를 훈련데이터로 사용할지)  
	# degree: parameter needed for kernel of type polynomial (default:3). (다항식 커널 계수 지정: default 3) 
	# gamma: parameter needed for all kernels except linear (default: 1/(dimension of data without y)). (감마 파라미터; kernel 에서 "linear"옵션 제외) 
	# cost: cost of constraints violation (default:1) --it is the "C"-constant of the regulrizaiton term in the Lagrange formulation.(cost 파라미터: default로 1지정)   
	# new_predict: logical, if TRUE, predict new dataset with fitted model, default is FALSE. (새로운 데이터로 예측 여부) 
	# new_quan_x : (새로 입력받는 값 중 quan_x에 대응되는 새로운 변수)
	# new_qual_x : (새로 입력받는 값 중 qual_x에 대응되는 새로운 변수)


	###################
	#### Packages #####
	###################
	load.pkg(c("e1071", "caret", "R2HTML"))

	###############################
	#### Functions in Analysis ####
	###############################
	
	html.output <- capture.output({ 


	## Generage HTML File
		R2HTML::HTML(R2HTML::as.title("Support Vector Machine (SVM)"),HR=1,file=stdout(),append=FALSE) 

	## Warnings
		if(!is.null(quan_x)) {
			is.nom <- sapply(quan_x,function(i) !is.numeric(dataset[,i]))
			if(any(is.nom)){
				warn.msg1 <- paste0("\a Warning : The type of variable '",paste(quan_x[is.nom],collapse=', '),"' is not numeric but selected as the quantitative variable. It was excluded from the analysis.")
				quan_x <- quan_x[!is.nom]
				new_quan_x <- new_quan_x[!is.nom]
				if(length(quan_x)==0) quan_x <- new_quan_x <- NULL
			}
		}
		if(is.null(quan_x)&is.null(qual_x)) warn.msg2 <- '\a Error : At least 1 independent variable should be selected. Analysis has been stopped.'
		if(!is.null(qual_x)) {
			is.num <- sapply(qual_x,function(i) is.numeric(dataset[,i]))
			if(any(is.num)){
				warn.msg3 <- paste0("\a Warning : The type of variable '",paste(qual_x[is.num],collapse=', '),"' is numeric but selected as the qualitative variable. It was coreced into character.")
				for(i in qual_x[is.num]) dataset[,i] <- factor(dataset[,i])
				for(i in new_qual_x[is.num]) dataset[,i] <- factor(dataset[,i])
			}
		}
		if (type=="eps-regression" & !is.numeric(dataset[,y])) warn.msg4 <- "\a Error : Inappropriate type of response variable. (EPS-regression is supported only for numeric variable.) Analysis has been stopped."
		if (cost<=0 | is.character(cost)) warn.msg5<- "\a Error : The cost should be numbers and greater than 0. Analysis has been stopped."
		if (gamma<=0 | is.character(gamma)) warn.msg6<- "\a Error : The gamma should be numbers and greater than 0. Analysis has been stopped."
		if (is.numeric(gamma) & kernel == "linear") warn.msg7 <- "\a Error : The gamma value is not available when kernel is linear. Analysis has been stopped."   
		if (is.numeric(degree) & kernel != "polynomial") warn.msg8 <- "\a Error : The degree is available only kernel is polynomial and should be greater than 0. The default value of degree is 3. Analysis has been stopped." 
		
		if(exists('warn.msg2')|exists('warn.msg4')|exists('warn.msg5')|exists('warn.msg6')|exists('warn.msg7')){
			R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
			if(exists('warn.msg2')){
				if(exists('warn.msg1')) R2HTML::HTML(warn.msg1,file=stdout())
				R2HTML::HTML(warn.msg2,file=stdout())
				}
			if (exists('warn.msg4')) R2HTML::HTML(warn.msg4,file=stdout())
			if (exists('warn.msg5')) R2HTML::HTML(warn.msg5,file=stdout())
			if (exists('warn.msg6')) R2HTML::HTML(warn.msg6,file=stdout())
			if (exists('warn.msg7')) R2HTML::HTML(warn.msg7,file=stdout())
			
		} else {
			if(is.numeric(degree) & kernel == "polynomial") {
				degree=degree
			}else if(is.numeric(degree) & kernel!="polynomial") { 
				R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
				if (exists('warn.msg8')) R2HTML::HTML(warn.msg8,file=stdout())
			}else if(!is.numeric(degree) & kernel=="polynomial"){ 
				degree=3
			}else{ 
				degree=NULL
			}

			if(is.null(gamma) & kernel != "linear") gamma <- 1/(length(quan_x)+length(qual_x))
			tuning <- F
			if(length(gamma)>1 | length(cost)>1) tuning <- T

			D <- dataset
			x <- c(quan_x,qual_x)
			dataset <- D[,c(y,x),drop=F]
			dataset <- dataset[complete.cases(dataset),,drop=F]
			
			if(new_predict){
				new_x <- c(new_quan_x,new_qual_x)
				new_dataset <- D[,c(new_x),drop=F]
				colnames(new_dataset) <- x
			}

			## If type is classification set y as factor	
			if (type=="C-classification"){ 
				dataset[,y]=factor(dataset[,y])
			}
			if (type=="eps-regression") {
				dataset[,y]=as.numeric(dataset[,y])
			}
		
			## Generate model formula
			x.vars=paste(x,collapse="+")
			model=as.formula(paste(y,x.vars,sep="~"))

			## Data Slicing (Training & Testing data)
			y.val=dataset[,y]
			#idx=createDataPartition(y=y.val,p=train_percent,list=F) 
			Wmes1=c("Some classes have a single record")
			tryCatch(
				{idx=createDataPartition(y=y.val,p=train_percent,list=F)},		
				warning=function(cond){ 
					R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
					mes=capture.output(cond)
					warn.msg11="\a WARNING : Please select correct type of response variable." 
					if (grepl(Wmes1,mes)) R2HTML::HTML(warn.msg11,file=stdout())	
				}
			)
			
			if(exists("idx")){
				train=dataset[idx,]
				test=dataset[-idx,]

				## Data Structure
				R2HTML::HTML(R2HTML::as.title("Data Structure"),HR=2,file=stdout())
				DS <- matrix(c('Number of observations' ,nrow(D)),ncol=2)
				DS <- rbind(DS,c('Number of observations for model fitting',nrow(dataset)))
				if(new_predict) DS <- rbind(DS,c('Number of observations for prediction',nrow(new_dataset[complete.cases(new_dataset),])))
				DS <- rbind(DS,c('Number of variables',length(x)+1))
				R2HTML::HTML(DS,file=stdout(), innerBorder = 1,align="left")
			
				## Varibale list
				R2HTML::HTML(R2HTML::as.title("Variable List"),HR=2,file=stdout())
				varlist <- matrix(c('Response variable',y,'Attributes for model fitting',paste(x,collapse=', ')),ncol=2,byrow=T)
				if(new_predict) varlist <- rbind(varlist,c('Attributes for prediction',paste(new_x,collapse=', ')))
				R2HTML::HTML(varlist,file=stdout(), innerBorder = 1,align="left")
				if(exists('warn.msg1')) R2HTML::HTML(warn.msg1,file=stdout())
				if(exists('warn.msg3')) R2HTML::HTML(warn.msg3,file=stdout())

				## Analysis Description
				R2HTML::HTML(R2HTML::as.title("Analysis Description"),HR=2,file=stdout())
				methd <- switch(kernel,linear="Linear", polynomial="Polynomial", radial="Radial", sigmoid="Sigmoid")
				AD <- matrix(c('Fitted model',paste(y,x.vars,sep="~"),'% of data that goes to training',paste0(train_percent*100,'%'),'SVM-Type',type,'SVM-Kernel ',methd,'Tuning hyperparametrs using a grid search',tuning),ncol=2,byrow=T)
				if(kernel=="polynomial") AD <- rbind(AD,c('Hyperparameter - Degree',degree))
				if(kernel%in%c('radial','polynomial','sigmoid')) AD <- rbind(AD,c('Hyperparameter - gamma',paste(gamma,collapse=', ')))
				AD <- rbind(AD,c('Hyperparameter - cost', paste(cost,collapse=', ')))
				R2HTML::HTML(AD,file=stdout(), innerBorder = 1,align="left")

				## Tuning parameters using grid search 
				Wmes=c("reaching max number of iterations")
				if(tuning){ 
					#fit.tune=tune.svm(model, data=train, degree=degree, gamma=gamma, cost=cost, na.action=na.omit) 
					warn.msg9="\a WARNING : Reaching max number of iterations. Please reset parameter options." 
					tryCatch(
						{fit.tune=tune.svm(model, data=train, degree=degree, gamma=gamma, cost=cost, na.action=na.omit)},		
						warning=function(cond){ 
							R2HTML::HTML(warn.msg9,file=stdout())	
						}
					)
					tune_lst=data.frame(t(list(Sample_Method=fit.tune$sampling,'Best_Parameters(gamma,cost)'=fit.tune$best.parameters, Best_Performance=fit.tune$best.performance))) # best performance is smallest error rate 
					colnames(tune_lst)=c("Sample Method","Best Parameters(gamma,cost)","Best Performance")
					R2HTML::HTML(R2HTML::as.title("Parameter Tuning of Function using Grid Search"),HR=2,file=stdout()) 	
					R2HTML::HTML(R2HTML::as.title("Information of Tuning Result"),HR=3,file=stdout()) 
					R2HTML::HTML(tune_lst,file=stdout(),innerBorder = 1,align="left",row.names=F)
				} else {
					## Fit model with trainning dataset 
					#fit=svm(model,data=train, scale=scale, type=type, kernel=kernel, degree=degree gamma=gamma, cost=cost, na.action=na.omit) 
					tryCatch(
						{fit=svm(model,data=train, scale=scale, type=type, kernel=kernel, degree=degree, gamma=gamma, cost=cost, na.action=na.omit)},		
						warning=function(cond){ 
							R2HTML::HTML(warn.msg9,file=stdout())	
						}
					)
					
					R2HTML::HTML(R2HTML::as.title(paste("Model fit with training datasets (",methd,")",sep="")),HR=2,file=stdout()) 
					R2HTML::HTML(R2HTML::as.title("Information of Result"),HR=3,file=stdout()) 
					res_lst=data.frame(t(list(Kernel=methd, SV_Total_Number=fit$tot.nSV, SV_Number_Each_Level=fit$nSV, Levels=fit$levels)))
					colnames(res_lst)=c("Kernel","Total number of SV","SV number of each level","Levels")
					R2HTML::HTML(res_lst,file=stdout(),innerBorder = 1,align="left",row.names=F)
				}


				## visualize (classes by color, sv by crosses) 
				if(length(x)<=2){ 
					if(tuning==FALSE){ 
						R2HTML::HTML(R2HTML::as.title("SVM Plot"),HR=2,file=stdout())
						REx_ANA_PLOT()
						plot(fit,train[,c(y,quan_x)])
						REx_ANA_PLOT_OFF("")
					}else{
						R2HTML::HTML(R2HTML::as.title("SVM Plot"),HR=2,file=stdout())
						REx_ANA_PLOT()
						plot(fit.tune$best.model,train[,c(y,quan_x)])
						REx_ANA_PLOT_OFF("")
					}
				}else{
					if(tuning==FALSE){
						R2HTML::HTML(R2HTML::as.title("SVM Plot after CMD Scaling"),HR=2,file=stdout())
						REx_ANA_PLOT()
						plot(cmdscale(dist(dataset[, quan_x])), xlab="CMDScale_1", ylab="CMDScale_2" , col=as.integer(dataset[,y]),pch=c("o","+")[1:nrow(dataset) %in% fit$index +1],main="Visualization of Data after CMD Scaling (2 Dimension)")
						REx_ANA_PLOT_OFF("")
					}else{
						R2HTML::HTML(R2HTML::as.title("SVM Plot after CMD Scaling"),HR=2,file=stdout())
						REx_ANA_PLOT()
						plot(cmdscale(dist(dataset[, quan_x])), xlab="CMDScale_1", ylab="CMDScale_2" , col=as.integer(dataset[,y]),pch=c("o","+")[1:nrow(dataset) %in% fit.tune$best.model$index+1],main="Visualization of Data after CMD Scaling (2 Dimension)")
						REx_ANA_PLOT_OFF("")
					}
				}
			
				## Predict with test dataset
				if(nrow(train)<nrow(dataset)){ 
					R2HTML::HTML(R2HTML::as.title("Predicted Result with Test Datasets"),HR=2,file=stdout()) 

					if(tuning) test[complete.cases(test),"test_pred"]=predict(fit.tune$best.model, newdata=test[complete.cases(test),])
					if(!tuning) test[complete.cases(test),"test_pred"]=predict(fit,newdata=test[complete.cases(test),])

					if(type=="C-classification"){
						test_confs=confusionMatrix(test[complete.cases(test),"test_pred"],test[complete.cases(test),y])	

						R2HTML::HTML(R2HTML::as.title("Confusion Matrix"),HR=3,file=stdout()) 
						cft <- test_confs$table
						colnames(cft) <- rownames(cft) <- paste0('Class: ',colnames(cft))
						R2HTML::HTML(cft,file=stdout(),innerBorder = 1,align="left",caption='<div style="text-align:left"> Row : Predicted class <br> Column : True class')
						
						R2HTML::HTML(R2HTML::as.title("Overall Statistics"),HR=3,file=stdout()) 
						ov <- t(as.matrix(test_confs$overall))
						ov <- ov[,c('Accuracy','AccuracyLower','AccuracyUpper','AccuracyNull','AccuracyPValue','Kappa','McnemarPValue'),drop=F]
						colnames(ov) <- c("ACC","Lower bound of 95% CI for Acc","Upper bound of 95% CI for Acc","NIR","P-value [Acc > NIR]","Kappa","P-value for Mcnemar's test")
						R2HTML::HTML(round(ov,4),file=stdout(),innerBorder = 1,align="left",digits=4,row.names=F)

						R2HTML::HTML(R2HTML::as.title("Statistics by Class"),HR=3,file=stdout()) 
						R2HTML::HTML(test_confs$byClass,file=stdout(),innerBorder = 1,align="left",digits=4)
					} else {
						REx_ANA_PLOT()
						plot(test[complete.cases(test),"test_pred"],test[complete.cases(test),y],xlab="Predicted",ylab="Observed",main="Comparison between predicted value and observed value")
						abline(0,1)
						REx_ANA_PLOT_OFF("")
					}					
				}

				## predict with new dataset
				if(new_predict){ 	
					tryCatch(
						{if(!tuning) new_dataset[complete.cases(new_dataset),"SVM_newpred"]=predict(fit,newdata=new_dataset[complete.cases(new_dataset),])		
						 if(tuning) new_dataset[complete.cases(new_dataset),"SVM_newpred"]=predict(fit.tune$best.model,newdata=new_dataset[complete.cases(new_dataset),])
						 out_dataset=new_dataset[rowSums(is.na(new_dataset))!=ncol(new_dataset),]
						 Ot=data.frame(SVM_newpred=out_dataset[,"SVM_newpred"])},		
						error=function(cond){ 
							R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
							warn.msg10="\a Error : The variables in prediction should be similar with variables in the model fitting. Please check the variables again."	
							R2HTML::HTML(warn.msg10,file=stdout())
						},
						warning=function(cond){ 
							R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
							warn.msg10="\a Error : The variables in prediction should be similar with variables in the model fitting. Please check the variables again."	
							R2HTML::HTML(warn.msg10,file=stdout())
						}
					)
				}
			}else{
				R2HTML::HTML("Please check the input variables again.",file=stdout())
			}
		}
		## Time Track for Analysis 
		R2HTML::HTMLhr(file=stdout()) 
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Support Vector Machine (SVM)",sep=""),file=stdout()) 
		R2HTML::HTMLhr(file=stdout()) 
	
	})

	if(exists('Ot')) {
		return(list(html=html.output,Output=Ot))
	} else {
		return(html.output)
	}
}

REx_DIST <- function(dataset,conti_var=NULL,binary_var=NULL,conti_method="euclidean",power=3,binary_method="JACCARD S3",
                     by="var",unit="log",zscore=FALSE,filepath=getwd()) {
  
  #### 변수 설명 ################################################################################
  ## dataset : 분석에 이용할 데이터셋
  ## conti_var : 분석셋내 분석할 연속형변수. c('','','') 형식으로 입력. 필수는 아니며, 입력되면 두개 변수 이상 입력받아야 함.
  ## binary_var : 분석셋내 분석할 이분형변수. c('','','') 형식으로 입력. 필수는 아니며, 입력되면 두개 변수 이상 입력받아야 함.
  ## conti_var와 binary_var 둘 중 하나는 무조건 입력받아야 함.
  ## conti_method : conti_var가 입력되면 필수로 입력되어야 함. 아래의 방법 중 선택할 수 있으며 중복선택 가능. c('','')의 형식으로 입력. 
  ##	"euclidean", "manhattan", "minkowski", "chebyshev", "sorensen", "gower", "soergel", "kulczynski_d", "canberra", "lorentzian", "intersection", "non-intersection", 
  ##	"wavehedges", "czekanowski", "motyka", "kulczynski_s", "tanimoto", "ruzicka", "inner_product", "harmonic_mean", 
  ##	"cosine", "hassebrook", "jaccard", "dice", "fidelity", "bhattacharyya", "squared_chord",  "squared_euclidean", 
  ##	"pearson", "neyman", "squared_chi", "prob_symm", "divergence", "clark", "additive_symm", 
  ##	"jeffreys", "topsoe", "jensen-shannon", "jensen_difference", "taneja", "kumar-johnson", "avg"
  ## power : minkowski가 선택되면 활성화. 3이상의 값만 가질 수 있음. 디폴트는 3.
  ## binary_method : binary_var가 입력되면 필수로 입력되어야 함. 아래의 방법 중 선택할 수 있으며 중복선택 가능. c('','')의 형식으로 입력. 
  ##	"JACCARD S3", "SOCKAL & MICHENER S4", "SOCKAL & SNEATH S5", 
  ##	"ROGERS & TANIMOTO S6", "CZEKANOWSKI S7", "GOWER & LEGENDRE S9", 
  ##	"OCHIAI S12", "SOKAL & SNEATH S13", "Phi of PEARSON S14", 
  ##	"GOWER & LEGENDRE S2"
  ## zscore : zscore 변환; true, false
  ## by : 변수,케이스간 거리계산 선택; "var"(변수간), "obs"(관측치간)
  ## power : minkowski distance에서 쓰일 power
  ## unit : the logarithm unit ; "log", "log2", "log10"
  ## filepath : distance by obs의 결과를 .csv 파일 형식으로 저장할 경로. getwd()함수를 이용해서 default 경로 설정. 
  ###############################################################################################
  
  #######################################################
  #### Required packages : R2HTML, philentropy, ade4 ####
  #######################################################
  load.pkg(c("R2HTML", "philentropy", "ade4"))
  
  ### method list 
  ## "dist.binary" function
  db_methods <-  c("JACCARD S3", "SOCKAL & MICHENER S4", "SOCKAL & SNEATH S5", 
                   "ROGERS & TANIMOTO S6", "CZEKANOWSKI S7", "GOWER & LEGENDRE S9", 
                   "OCHIAI S12", "SOKAL & SNEATH S13", "Phi of PEARSON S14", "GOWER & LEGENDRE S2")
  
  html.output <- capture.output({
    ### Analysis Name
    R2HTML::HTML(R2HTML::as.title("Distance Measurement"),HR=1,file=stdout(),append=FALSE) 
    
    ### Filepath Check
    if(file.exists(filepath)) {
      setwd(filepath)
    } else {
      dir.create(filepath)
      setwd(filepath)
    }
    
    ### Warnings
    if(!is.null(conti_var)) {
      is.nom <- sapply(conti_var,function(i) !is.numeric(dataset[,i]))
      if(any(is.nom)){
        isare1 <- ifelse(sum(any(is.nom))>1,'are','is')
        isare2 <- ifelse(sum(any(is.nom))>1,'They were','It was')
        warn.msg1 <- paste0("\a Warning : The type of variable '",paste(conti_var[is.nom],collapse=', '),"' ",isare1," not numeric but selected as the continuous variable. ",isare2," excluded from the analysis.")
        conti_var <- conti_var[!is.nom]
        if(length(conti_var)==0) conti_var <- NULL
      }
    }
    if(!is.null(binary_var)) {
      not.bin <- sapply(binary_var,function(i) length(levels(factor(dataset[,i])))!=2)
      if(any(not.bin)){
        isare1 <- ifelse(sum(any(not.bin))>1,'are','is')
        isare2 <- ifelse(sum(any(not.bin))>1,'They were','It was')
        warn.msg2 <- paste0("\a Warning : The variable '",paste(binary_var[not.bin],collapse=', '),"' ",isare1," not binary. ",isare2," excluded from the analysis.")
        binary_var <- binary_var[!not.bin]
        if(length(binary_var)==0) binary_var <- NULL
      }
    }
    if(length(conti_var)<2 & length(binary_var)<2)	warn.msg3 <- "\a Error : Analysis can be conducted with at least 2 continuous variables (or 2 binary variables). Analysis has been stopped."
    
    if(exists('warn.msg3')) {
      if(exists('warn.msg1')) HTML(warn.msg1,file=stdout())
      if(exists('warn.msg2')) HTML(warn.msg2,file=stdout())
      HTML(warn.msg3,file=stdout())
    } else {
      dataset <- dataset[,c(conti_var,binary_var),drop=F]
      
      ## Data structure
      HTML(as.title("Data Structure"),HR=2,file=stdout())
      DS <- matrix(c("Number of Observations",nrow(dataset),"Number of Variables",length(c(conti_var,binary_var))),ncol=2,byrow=T)
      HTML(DS,file=stdout(), innerBorder = 1,align="left")
      
      # Varibale list
      R2HTML::HTML(R2HTML::as.title("Variable List"),HR=2,file=stdout())
      varlist <- c()
      if(!is.null(conti_var)) varlist <- rbind(varlist,c('Continuous variables',paste(conti_var,collapse=', ')))
      if(!is.null(binary_var)) varlist <- rbind(varlist,c('Binary variables',paste(binary_var,collapse=', ')))
      R2HTML::HTML(varlist,file=stdout(), innerBorder = 1,align="left")
      if(exists('warn.msg1')) R2HTML::HTML(warn.msg1,file=stdout())
      if(exists('warn.msg2')) R2HTML::HTML(warn.msg2,file=stdout())
      
      # Analysis Description
      R2HTML::HTML(R2HTML::as.title("Analysis Description"),HR=2,file=stdout())
      By <- switch(by,var='Variables',obs='Observations')
      AD <- matrix(c('Measuring distances between',By),ncol=2,byrow=T)
      if(!is.null(conti_var)) {
        AD <- rbind(AD,matrix(c('Standardization for continuous variable',zscore,'Method (Continuous variable)',paste(conti_method,collapse=', ')),ncol=2,byrow=T))
        if('minkowski'%in%conti_method) AD <- rbind(AD,c('Power for minkowski measure',power))
      }
      if(!is.null(binary_var)) AD <- rbind(AD,c('Method (Binary variable)',paste(binary_method,collapse=', ')))
      Base <- switch(unit,log='e',log2='2',log10='10')
      AD <- rbind(AD,c('Base of log-transformation',Base))
      R2HTML::HTML(AD,file=stdout(),innerBorder = 1,align="left")
      
      ### 케이스 처리 요약
      HTML(as.title("Case Proccess Summary"),HR=2,file=stdout()) 
      NV <- nrow(dataset[complete.cases(dataset),])
      NM <- nrow(dataset)-nrow(dataset[complete.cases(dataset),])
      NT <- nrow(dataset)
      CPS <- matrix(c("Valid(%)","Missing(%)","Total(%)", paste0(NV," (",round(NV/NT*100,1),")"),paste0(NM," (",round(NM/NT*100,1),")"),paste0(NT," (100)")),nrow=2,byrow=TRUE)
      HTML(CPS,file=stdout(), innerBorder = 1,align="left")
      
      if(!is.null(conti_var)) {
        conti_data <- dataset[complete.cases(dataset[,conti_var]),conti_var]
        if(zscore) for(i in conti_var) conti_data[,i] <- (conti_data[,i]-mean(conti_data[,i]))/sd(conti_data[,i])
      }
      if(!is.null(binary_var)) {
        binary_data <- dataset[,binary_var]
        for(i in binary_var) {
          temp <- factor(binary_data[,i])
          temp.1 <- replace(binary_data[,i],binary_data[,i]==levels(temp)[1],0)
          temp.2 <- replace(temp.1,temp.1==levels(temp)[2],1)
          binary_data[,i] <- as.numeric(temp.2)
        }
      }
      
      ### Analysis
      if (by=="var") {
        HTML(as.title("Distances between Variables"),HR=2,file=stdout()) 
        if(!is.null(conti_var)){
          HTML(as.title("Continuous Variables"),HR=3,file=stdout())
          for(i in conti_method){
            HTML(as.title(paste("Method:",i,"distance")),HR=4,file=stdout())
            if(i=="minkowski") {
              result <- distance(t(conti_data),method=i,unit=unit,p=power)
            } else {
              result <- distance(t(conti_data),method=i,unit=unit)
            }
            if (length(conti_var)>2) {
              colnames(result) <- rownames(result) <- conti_var
            }
            HTML(round(result,4),file=stdout(),innerBorder = 1,align="left",digits=4)
          }
        }
        
        if(!is.null(binary_var)){
          HTML(as.title("Binary Variables"),HR=3,file=stdout())
          for(i in match(binary_method,db_methods)){
            HTML(as.title(paste("Method:",db_methods[i],"distance")),HR=4,file=stdout())
            result <- as.matrix(dist.binary(t(binary_data),method=i,diag = TRUE, upper = TRUE))
            colnames(result) <- rownames(result) <- binary_var
            HTML(round(result,4),file=stdout(),innerBorder = 1,align="left",digits=4)
          }
        }
      }
      
      if (by=="obs") {
        ## Check Filepath. if it doesn't exist then create it.
        if(file.exists(filepath)) {
          setwd(filepath)
        } else {
          dir.create(filepath)
          setwd(filepath)
        }
        HTML(as.title("Distances between Observations"),HR=2,file=stdout()) 
        HTML(paste0("\a The results of distances between observations are saved in '",filepath,"' with row names of observation ID and column names of 'method_obserbation ID.'"),file=stdout())
        if(!is.null(conti_var)){
          conti_O <- data.frame(obs_conti=rownames(conti_data))
          for(i in conti_method){
            if(i=="minkowski") {
              result <- distance(conti_data,method=i,unit=unit,p=power)
            } else {
              result <- distance(conti_data,method=i,unit=unit)
            }
            colnames(result) <- paste0(i,'_',rownames(conti_data))
            write.csv(result,paste0(filepath,"/continous_distance_",i,"_result_by_obs.csv"))
            conti_O <- data.frame(conti_O,result)
          }
        }
        
        if(!is.null(binary_var)){
          binary_O <- data.frame(obs_binary=rownames(binary_data))
          for(i in match(binary_method,db_methods)){
            result <- as.matrix(dist.binary(binary_data,method=i,diag = TRUE, upper = TRUE))
            colnames(result) <- paste0(db_methods[i],'_',rownames(binary_data))
            write.csv(result,paste0(filepath,"/binary_distance_",binary_method[i],"_result_by_obs.csv"))
            binary_O <- data.frame(binary_O,result)
          }
        }
      }
    }
    
    ## Time Track for Analysis 
    R2HTML::HTMLhr(file=stdout()) ;
    R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Distance Measurement",sep=""),file=stdout()) ;
    R2HTML::HTMLhr(file=stdout()) ;
  }) #html.output over
  return(html.output)
}

REx_DA=function(dataset, y, quan_x=NULL, qual_x=NULL, method="LDA", train_percent=0.7, Prior="RF" ,prior = NULL, prior_y=NULL, new_predict=FALSE, new_quan_x = NULL, new_qual_x = NULL){

	###################
	#### Arguments ####
	###################
	##Require four inputs: 

	# dataset: the data which will be used in analysis. (데이터셋 이름) 
	# y: name of target variable (or true classification). (종속변수) 
	# quan_x: names of quantitative attribute. (설명변수: 질적변수) 
	# qual_x: names of qualitative attribute. (설명변수: 양적변수)
	# method : 'LDA','QDA' or 'MDA'. (분석방법 선택: "LDA" 혹은 "QDA" 혹은 "MDA")
	# train_percent: percentage of training data to divide, default=0.7 (range of 0.5~1). (훈련데이터셋을 전체 데이터의 몇% 비율로 지정) 
	# Prior : "RF","Equal","Custom". (사전확률 지정) 
	# prior : Prior가 "Custom"일 때 사용자가 입력한 값 사용. 합은 1이어야 함. 그 외에는 NULL입력.
	# prior_y : Prior가 "Custom"일 때 prior에 대응되는 y값,즉 사용자 지정. 그 외에는 NULL. 
	# new_predict: logical, if TRUE, predict new dataset with fitted model, default is FALSE. 
	# new_quan_x : 새로 입력받는 값 중 quan_x에 대응되는 새로운 변수.
	# new_qual_x : 새로 입력받는 값 중 qual_x에 대응되는 새로운 변수.


	###################
	#### Packages #####
	###################
	load.pkg(c("klaR", "caret", "R2HTML", "e1071", "MASS"))

	###############################
	#### Functions in Analysis ####
	###############################
	
	html.output <- capture.output({

		## Generage HTML File
		R2HTML::HTML(R2HTML::as.title("Discriminant Analysis"),HR=1,file=stdout(),append=FALSE) 

		## Warnings
		if(!is.null(quan_x)) {
			is.nom <- sapply(quan_x,function(i) !is.numeric(dataset[,i]))
			if(any(is.nom)){
				warn.msg1 <- paste0("\a Warning : The type of variable '",paste(quan_x[is.nom],collapse=', '),"' is not numeric but selected as the quantitative variable. It was excluded from the analysis.")
				quan_x <- quan_x[!is.nom]
				new_quan_x <- new_quan_x[!is.nom]
				if(length(quan_x)==0) quan_x <- new_quan_x <- NULL
			}
		}
		if(is.null(quan_x)&is.null(qual_x)) warn.msg2 <- '\a Error : At least 1 independent variable should be selected. Analysis has been stopped.'
		if(!is.null(qual_x)) {
			is.num <- sapply(qual_x,function(i) is.numeric(dataset[,i]))
			if(any(is.num)){
				warn.msg3 <- paste0("\a Warning : The type of variable '",paste(qual_x[is.num],collapse=', '),"' is numeric but selected as the qualitative variable. It was coreced into character.")
				for(i in qual_x[is.num]) dataset[,i] <- factor(dataset[,i])
				for(i in new_qual_x[is.num]) dataset[,i] <- factor(dataset[,i])
			}
		}
		if(!Prior%in%c('RF','Equal')) if(sum(prior)!=1) warn.msg4 <- "\a Error : Sum of prior probabilities should be 1. Analysis has been stopped."

		if(exists('warn.msg2')|exists('warn.msg4')){
			R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
			if(exists('warn.msg2')){
				if(exists('warn.msg1')) R2HTML::HTML(warn.msg1,file=stdout())
				R2HTML::HTML(warn.msg2,file=stdout())
			}
			if(exists('warn.msg4')) R2HTML::HTML(warn.msg2,file=stdout())
		} else {
			D <- dataset
			x <- c(quan_x,qual_x)
			dataset <- D[,c(y,x),drop=F]
			dataset <- dataset[complete.cases(dataset),,drop=F]
			
			if(new_predict){
				new_x <- c(new_quan_x,new_qual_x)
				new_dataset <- D[,c(new_x),drop=F]
				colnames(new_dataset) <- x
			}

			## Generate model formula
			x.vars=paste(x,collapse="+")
			model=as.formula(paste(y,x.vars,sep="~"))

			## Data Slicing (Training & Testing data)
			dataset[,y] <- factor(dataset[,y])
			y.val=dataset[,y]
			train_idx=createDataPartition(y=y.val,p=train_percent,list=F) 
			train=dataset[train_idx,,drop=F]
			test=dataset[-train_idx,,drop=F]

			## Prior probabilities
			if(Prior=="Custom") {
				prior <- prior[match(levels(y.val),prior_y)]
			} else if(Prior=="Equal") {
				prior <- rep(1,length(levels(y.val)))/length(levels(y.val))
			}

			# Data Structure
			R2HTML::HTML(R2HTML::as.title("Data Structure"),HR=2,file=stdout())
			DS <- matrix(c('Number of observations' ,nrow(D)),ncol=2)
			DS <- rbind(DS,c('Number of observations for model fitting',nrow(dataset)))
			if(new_predict) DS <- rbind(DS,c('Number of observations for prediction',nrow(new_dataset[complete.cases(new_dataset),])))
			DS <- rbind(DS,c('Number of variables',length(x)+1))
			R2HTML::HTML(DS,file=stdout(), innerBorder = 1,align="left")

			# Varibale list
			R2HTML::HTML(R2HTML::as.title("Variable List"),HR=2,file=stdout())
			varlist <- matrix(c('Response variable',y,'Attributes for model fitting',paste(x,collapse=', ')),ncol=2,byrow=T)
			if(new_predict) varlist <- rbind(varlist,c('Attributes for prediction',paste(new_x,collapse=', ')))
			R2HTML::HTML(varlist,file=stdout(), innerBorder = 1,align="left")
			if(exists('warn.msg1')) R2HTML::HTML(warn.msg1,file=stdout())
			if(exists('warn.msg3')) R2HTML::HTML(warn.msg3,file=stdout())

			# Analysis Description
			R2HTML::HTML(R2HTML::as.title("Analysis Description"),HR=2,file=stdout())
			methd <- switch(method,LDA='Linear Discriminant Analysis (LDA)',QDA='Quadratic Discriminant Analysis (QDA)',MDA='Mixture Discriminant Analysis (MDA)')
			PRIOR <- switch(Prior,RF='Relative frequancies of groups in training data',Equal='Equal probabilities across the group in training data',Custom='Input by user')
			AD <- matrix(c('Fitted model',paste(y,x.vars,sep=" ~ "),'% of data that goes to training',paste0(train_percent*100,'%'),'Method',methd,'Prior probabilities',PRIOR),ncol=2,byrow=T)
			R2HTML::HTML(AD,file=stdout(), innerBorder = 1,align="left")

			## Fit model 
			if(method=="LDA") {
				if(Prior=="RF") {
					fit=lda(model,data=train,na.action="na.omit")
				} else {
					fit=lda(model,data=train,na.action="na.omit",prior=prior)
				}
			} else if(method=="QDA") {
				if(Prior=="RF") {
					Emes="rank deficiency in group 1"
					Emes1="some group is too small for 'qda'"
					tryCatch(
						{fit=qda(model,data=train,na.action="na.omit")},		
						 error=function(cond){
							R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
							mes=capture.output(cond)
							warn.msg5="\a WARNING : There is a rank deficiency, some variables are collinear. Please check the  explanatory variables." 
							warn.msg6="\a WARNING : The factor explanatory variable is not available in 'qda'. Please select correct variables."  
							if (grepl(Emes,mes)) R2HTML::HTML(warn.msg5,file=stdout())
							if (grepl(Emes1,mes)) R2HTML::HTML(warn.msg6,file=stdout()) 
						}
					)
				} else {
					Emes="rank deficiency in group 1"
					Emes1="some group is too small for 'qda'"
					tryCatch(
						{fit=qda(model,data=train,na.action="na.omit",prior=prior)},		
						 error=function(cond){
							R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
							mes=capture.output(cond)
							warn.msg5="\a WARNING : There is a rank deficiency, some variables are collinear. Please check the  explanatory variables." 
							warn.msg6="\a WARNING : The factor explanatory variable is not available in 'qda'. Please select correct variables."  
							if (grepl(Emes,mes)) R2HTML::HTML(warn.msg5,file=stdout())
							if (grepl(Emes1,mes)) R2HTML::HTML(warn.msg6,file=stdout()) 
						}
					)				
				}
			} else if(method=="MDA") {
				if(Prior=="RF") {
					fit=mda(model,data=train,na.action="na.omit")
				} else {
					fit=mda(model,data=train,na.action="na.omit",weights=prior)
				}
			}

			R2HTML::HTML(R2HTML::as.title(paste("Model Fit with Training Datasets (",method,")",sep="")),HR=2,file=stdout()) 
			# prior prob.
			R2HTML::HTML(R2HTML::as.title("Prior Probabilities of Groups"),HR=3,file=stdout())
			if(method=="MDA"){
				pp <- matrix(fit$prior,nrow=1); colnames(pp) <- levels(train[,y])
			} else {
				pp <- t(as.matrix(fit$prior))
			}
			R2HTML::HTML(round(pp,4), file=stdout(),innerBorder = 1,align="left",row.names=F)

			if(method=="LDA"){
				# coefficient of LD
				R2HTML::HTML(R2HTML::as.title(paste0("Coefficients of ",switch(method,LDA="Linear",QDA="Quadratic",MDA="Mixture")," Discriminants")),HR=3,file=stdout())
				R2HTML::HTML(round(fit$scaling,4), file=stdout(),innerBorder = 1,align="left",row.names=F)

				# Singluar value (https://www.r-bloggers.com/computing-and-visualizing-lda-in-r/)
				R2HTML::HTML(R2HTML::as.title("Proportion of the Trace"),HR=3,file=stdout())
				prop <- matrix(with(fit,svd^2/sum(svd^2)),nrow=1)
				colnames(prop) <- paste0("LD",1:ncol(prop))
				msg <- paste0("About ",round(prop[,1]*100,2),"% of the between-group variance is explained by first LD.")
				R2HTML::HTML(round(prop,4), file=stdout(),innerBorder = 1,align="left",row.names=F)
				R2HTML::HTML(msg, file=stdout())

				# LDA histogram
				R2HTML::HTML(R2HTML::as.title(paste0("Stacked Histogram of the ",method," Values")),HR=3,file=stdout())
				lda.values <- predict(fit)
				for(i in 1:ncol(lda.values$x)){
					R2HTML::HTML(R2HTML::as.title(paste0("LD",i)),HR=4,file=stdout())
					REx_ANA_PLOT()
					ldahist(data=lda.values$x[,i],g=train[,y],col=3)
					REx_ANA_PLOT_OFF("")
				}
			}
			# confusion matrix
			R2HTML::HTML(R2HTML::as.title("Confusion Matrix"),HR=3,file=stdout()) 
			train_pred <- predict(fit)
			if(method!="MDA") train_pred <- train_pred$class
			train_confs=confusionMatrix(train_pred,train[,y])
			cft <- train_confs$table
			colnames(cft) <- rownames(cft) <- paste0('Class: ',colnames(cft))
			R2HTML::HTML(cft,file=stdout(),innerBorder = 1,align="left",caption='<div style="text-align:left"> Row : Predicted class <br> Column : True class')

			## Predict with test dataset
			if(nrow(train)<nrow(dataset)){ 
				R2HTML::HTML(R2HTML::as.title("Predicted Result with Test Datasets"),HR=2,file=stdout()) 

				test_pred=predict(fit, newdata=test)
				if(method!="MDA") test_pred <- test_pred$class
				test_confs=confusionMatrix(test_pred,test[,y])
				
				R2HTML::HTML(R2HTML::as.title("Confusion Matrix"),HR=3,file=stdout()) 
				cft <- test_confs$table
				colnames(cft) <- rownames(cft) <- paste0('Class: ',colnames(cft))
				R2HTML::HTML(cft,file=stdout(),innerBorder = 1,align="left",caption='<div style="text-align:left"> Row : Predicted class <br> Column : True class')

				R2HTML::HTML(R2HTML::as.title("Overall Statistics"),HR=3,file=stdout()) 
				ov <- t(as.matrix(test_confs$overall))
				ov <- ov[,c('Accuracy','AccuracyLower','AccuracyUpper','AccuracyNull','AccuracyPValue','Kappa','McnemarPValue'),drop=F]
				colnames(ov) <- c("ACC","Lower bound of 95% CI for Acc","Upper bound of 95% CI for Acc","NIR","P-value [Acc > NIR]","Kappa","P-value for Mcnemar's test")
				R2HTML::HTML(round(ov,4),file=stdout(),innerBorder = 1,align="left",digits=4,row.names=F)
				
				R2HTML::HTML(R2HTML::as.title("Statistics by Class"),HR=3,file=stdout()) 
				R2HTML::HTML(test_confs$byClass,file=stdout(),innerBorder = 1,align="left",digits=4)
			}

			## predict with new dataset 
			if(new_predict){ 	
				tryCatch(
					{new_dataset[complete.cases(new_dataset),"DA_newpred"]=predict(fit, newdata=new_dataset[complete.cases(new_dataset),])$class
					#if(method!="MDA") new_pred <- new_pred$class
					out_dataset=new_dataset[rowSums(is.na(new_dataset))!=ncol(new_dataset),]
					Ot=data.frame(DA_newpred=out_dataset[,"DA_newpred"])},		
					error=function(cond){ 
						R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
						warn.msg6="\a Error : The variables in prediction should be similar with variables in the model fitting. Please check the variables again."	
						R2HTML::HTML(warn.msg6,file=stdout())
					},
					warning=function(cond){ 
						R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
						warn.msg6="\a Error : The variables in prediction should be similar with variables in the model fitting. Please check the variables again."	
						R2HTML::HTML(warn.msg6,file=stdout())
					}
				)
			}
		
		}
		## Time Track for Analysis 
		R2HTML::HTMLhr(file=stdout()) ;
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Discriminant Analysis",sep=""),file=stdout()) ;
		R2HTML::HTMLhr(file=stdout()) ;
	})

	if(exists('Ot')) {
		return(list(html=html.output,Output=Ot))
	} else {
		return(html.output)
	}
}

REx_Longi <- function(dataset,res_var,time_var, id_var, qual_var=NULL, quan_var=NULL, vars=NULL, noint=FALSE,
		      time=NULL, time_code=NULL, rand_int=FALSE, rand_time=FALSE, rand_corr=FALSE,var_method="ML",
		      resid_plot=FALSE, profile_plot=FALSE, Predict=FALSE, Resid=FALSE){

	#### 변수 설명 ####
	## dataset : 데이터셋이름
	## res_var : 종속변수, 필수로 입력되어야 하며 2개 이상 입력가능
	## time_var : 변수설정탭의 시간변수 (필수)
	## id_var : 변수설정탭의 관측단위 변수 (필수)
	## qual_var : 질적설명변수(independent variables) (2개 이상의 변수가 있는 경우 c(,,)로 구분, 필수사항아님)
	## quan_var : 양적설명변수(independent variables) (2개 이상의 변수가 있는 경우 c(,,)로 구분, 필수사항아님)
	## vars : 최종모형
	## noint : 모형에 절편한 포함 여부. TRUE or FALSE.
	## time, time_code : UI 설명 참고
	## rand_int : 분석옵션탭의 '임의효과 중 관측 단위에 따라 절편이 다름을 가정'이 선택되면 TRUE, 그렇지 않으면 FALSE
	## rand_time : 분석옵션탭의 '관측 단위에 따라 시간에 따른 종속변수의 변화율이  다름을 가정'이 선택되면 TRUE, 그렇지 않으면 FALSE
	## rand_corr : 분석옵션탭의 '임의효과 사이에 상관성이 있음을 가정'이 선택되면 TRUE, 그렇지 않으면 FALSE
	## var_method : 
	## resid_plot
	## profile_plot
	## 
	###################
	#### Required packages : R2HTML, car, lattice, lme4 ####
	###################
	load.pkg(c("R2HTML", "car", "lme4", "lattice"))

	html.output <- capture.output({
		# Title
		R2HTML::HTML(R2HTML::as.title("Longitudinal Data Analysis"),HR=1,file=stdout(),append=FALSE)

		## Warnings
		if(!is.numeric(dataset[,res_var])) warn.msg1 <- '\a Error : Response variable should be numeric. Analysis has been stopped.'

		if(is.null(vars)) {
			Vars <- c()
		} else {
			Vars <-  unique(unlist(strsplit(vars,":")))
		}
		qual_var <- qual_var[qual_var%in%Vars]
		quan_var <- c(quan_var[quan_var%in%Vars],res_var)

		if(!is.null(quan_var)) {
			is.nom <- sapply(quan_var,function(i) !is.numeric(dataset[,i]))
			if(any(is.nom)){
				warn.msg3 <- paste0("\a Warning : The type of variable '",paste(quan_var[is.nom],collapse=', '),"' is not numeric but selected as the quantitative variable. It was excluded from the analysis.")

				# explanatory variables
				vars <- vars[!vars%in%ec]
				if(length(vars)==0) vars <- NULL
			}
		}
		if(!is.null(qual_var)) {
			is.num <- sapply(qual_var,function(i) is.numeric(dataset[,i]))
			if(any(is.num)) warn.msg5 <- paste0("\a Warning : The type of variable '",paste(qual_var[is.num],collapse=', '),"' is continuous but selected as the qualitative variable. It was coreced into character.")
			for(i in qual_var) dataset[,i] <- factor(dataset[,i])
		}
		if(!rand_int & !rand_time) warn.msg6 <- "\a Error : No random effects term was specified by user. Analysis has been stopped."

		if(exists('warn.msg1')|exists('warn.msg6')){
			R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
			if(exists('warn.msg1')) R2HTML::HTML(warn.msg1,file=stdout())
			if(exists('warn.msg6')) R2HTML::HTML(warn.msg6,file=stdout())
		} else {

			## Model formula
			if(rand_int){
				if(rand_time) {
					if(rand_corr){
						re <- paste0('(1 + ',time_var,' | ',id_var,')')
					} else {
						re <- paste0('(1 | ',id_var,') + (0 + ',time_var,' | ',id_var,')')
					}
				} else {
					re <- paste0('(1 | ',id_var,')')
				}
			} else {
				if(rand_time) {
					re <- paste0('(',paste(c(-1,time_var),collapse=' + '),' | ',id_var,')')
				}
			}
			fe <- paste(vars,collapse=' + ')
			if(fe=="") fe <- NULL
			form.0 <- paste0(res_var,' ~ ',paste(c(time_var,fe,re),collapse=' + '))
			if(noint) form.0 <- paste(form.0,-1)
			form.1 <- as.formula(form.0)

			## Data processing
			# id var
			dataset[,id_var] <- factor(dataset[,id_var])
			# time var
			temp <- rep(NA,nrow(dataset))
			for(i in 1:length(time)) temp[dataset[,time_var]==time[i]] <- time_code[i]
			dataset[,time_var] <- temp

			## Model Fitting
			if(var_method=="ML") fit <- suppressWarnings(try(lmer(form.1,data=dataset,REML=0)))
			if(var_method=="REML") fit <- suppressWarnings(try(lmer(form.1,data=dataset)))

			## Data Structure
			R2HTML::HTML(R2HTML::as.title("Data Structure"),HR=2,file=stdout())
			total.var <- ncol(dataset)
			used.var <- ifelse(is.null(vars),0,length(unique(unlist(strsplit(vars,":")))))
			DS <- matrix(c('Number of observations','Number of total variables','Number of used variables',nrow(dataset),total.var,used.var+1),ncol=2)
			R2HTML::HTML(DS,file=stdout(), innerBorder = 1,align="left")

			## Varibale list
			R2HTML::HTML(R2HTML::as.title("Variable List"),HR=2,file=stdout())
			if(is.null(vars)) {
				Vars <- c()
			} else {
				Vars <-  unique(unlist(strsplit(vars,":")))
			}
			qual <- qual_var[qual_var%in%Vars]
			quan <- c(quan_var[quan_var%in%Vars],res_var)
			varlist <- matrix(c('Quantitative variable',paste(quan,collapse=', ')),ncol=2,byrow=T)
			if(length(qual)!=0) varlist <- rbind(varlist,c('Qualitative variable',paste(qual,collapse=', ')))
			R2HTML::HTML(varlist,file=stdout(), innerBorder = 1,align="left")
			if(exists('warn.msg3')) R2HTML::HTML(warn.msg3,file=stdout())
			if(exists('warn.msg5')) R2HTML::HTML(warn.msg5,file=stdout())

			## Analysis Description
			R2HTML::HTML(R2HTML::as.title("Analysis Description"),HR=2,file=stdout())
			AD <- matrix(c('Response variables',res_var,'Explanatory variable',paste(vars,collapse=', '),'Time variable',time_var,'ID variable',id_var,'Intercept included',!noint,'Random effect on intercept',rand_int,'Random effect on time',rand_time),ncol=2,byrow=T)
			if(rand_int&rand_time) AD <- rbind(AD,c('Correlation between random effect',rand_corr))
			varmethod <- switch(var_method,ML="Maximum likelihood",REML="Restricted maximum likelihood")
			AD <- rbind(AD,matrix(c('Distribution of response variable',family(fit)$family,'Link function',family(fit)$link,'Method for estimating variance',varmethod),ncol=2,byrow=T))
			R2HTML::HTML(AD,file=stdout(), innerBorder = 1,align="left")

			## Profile plot
			if(profile_plot){
				R2HTML::HTML(R2HTML::as.title("Profile Plot"),HR=2,file=stdout())
				n.obs <- length(levels(dataset[,id_var]))
				wh <- ifelse(n.obs<=16,500,ifelse(n.obs>16 & n.obs<=25,550,ifelse(n.obs>25 & n.obs<=36,600,ifelse(n.obs>36 & n.obs<=49,650,"\a Profile pictures are only available for observations of less than 50 observations."))))
				if(is.numeric(wh)){
					REx_ANA_PLOT(w=wh,h=wh)
					xyplot(as.formula(paste(res_var,'~',time_var,'|',id_var)),type=c("p","r"),data=dataset, grid=T,col.line="dark orange",pch=16)
					REx_ANA_PLOT_OFF('')
				} else {
					R2HTML::HTML(wh,file=stdout())
				}
			}

			### Model Fitting
			R2HTML::HTML(R2HTML::as.title(paste0("Linear Mixed Model Fit by ",varmethod)),HR=2,file=stdout())
			
			## Check convergence
			CC <- capture.output(summary(fit))
			is.conv <- grep("Model failed to converge",CC)
			if(length(is.conv)==0){
				## Time variable
				R2HTML::HTML(R2HTML::as.title("Transformed Time Variable"),HR=3,file=stdout())
				TV <- data.frame(time,time_code)
				colnames(TV) <- c(time_var,'Coded value')
				R2HTML::HTML(TV,file=stdout(), innerBorder = 1,align="left",row.names=F)
				
				# Fixed effect
				R2HTML::HTML(R2HTML::as.title("Fixed effects"),HR=3,file=stdout())
				CE <- data.frame(summary(fit)$coef)
				CE <- cbind(CE,pt(abs(CE[,3]),df=1,lower.tail=F))
				colnames(CE) <- c('Estimate','SE','T-value','P-value')
				CE[,ncol(CE)] <- format(CE[,ncol(CE)],scientific=T,digits=4)
				R2HTML::HTML(CE,file=stdout(), innerBorder = 1,align="left",digits=4)

				# Random effect
				R2HTML::HTML(R2HTML::as.title("Random effects"),HR=3,file=stdout())
				RE <- as.data.frame(summary(fit)$varcor)
				if(rand_corr){
					RE <- RE[,-1]
					RE[nrow(RE),1] <- "Residual"
					colnames(RE) <- c('Variable 1','Variable 2','Variance (Covariance)','SD (Correlation)')
					RE[is.na(RE)] <- ""
					R2HTML::HTML(RE,file=stdout(), innerBorder = 1,align="left",digits=4,row.names=F)

				} else {
					RE <- RE[,-c(1,3)]
					RE[nrow(RE),1] <- "Residual"
					colnames(RE) <- c('Variable','Variance','SD')
					RE[is.na(RE)] <- ""
					R2HTML::HTML(RE,file=stdout(), innerBorder = 1,align="left",digits=4,row.names=F)
				}

				# Measure of GOF
				if(var_method=="ML"){
					R2HTML::HTML(R2HTML::as.title("Measures of Goodness of Fit"),HR=3,file=stdout())
					GOF <- as.data.frame(t(as.matrix(summary(fit)$AIC)))
					R2HTML::HTML(GOF,file=stdout(), innerBorder = 1,align="left",digits=4,row.names=F)
				}

				# Residual Plot
				if(resid_plot){
					R2HTML::HTML(R2HTML::as.title("Residual Plot"),HR=3,file=stdout())
					REx_ANA_PLOT()
					plot(fit,pch=16,ylab="Residual",xlab="Fitted Value")
					REx_ANA_PLOT_OFF('')
				}

				# 다음의 값들은 엑셀 시트에 저장
				if(Predict) {
					temp.0 <- data.frame(Fitted_Longi=fitted(fit))
					temp <- data.frame(Fitted_Longi=rep(NA,nrow(dataset)))
					row.dataset <- rownames(dataset)
					row.output <- rownames(temp.0)
					for(i in row.dataset) if(i%in%row.output) temp[row.dataset==i,1] <- temp.0[row.output==i,1]
					O <- temp
				}
				if(Resid) {
					temp.0 <- data.frame(Resid_Longi=resid(fit))
					temp <- data.frame(Resid_Longi=rep(NA,nrow(dataset)))
					row.dataset <- rownames(dataset)
					row.output <- rownames(temp.0)
					for(i in row.dataset) if(i%in%row.output) temp[row.dataset==i,1] <- temp.0[row.output==i,1]
					if(exists('O')) {
						O <- cbind(O,temp)
					} else { 
						O <- temp
					}
				}
			} else {
				conv.msg <- paste('\a',CC[is.conv])
				R2HTML::HTML(conv.msg,file=stdout())
			}
		}
		R2HTML::HTMLhr(file=stdout())
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Multivariate Linear Regression",sep=""),stdout())
		R2HTML::HTMLhr(file=stdout())
	})

	if(exists('O')){
		return(list(html=html.output,Output=O))
	} else {
		return(html.output)
	}
}

REx_MultinomGLM	<- function(dataset, res_var, quan_var=NULL, qual_var=NULL, vars=NULL, res_type=c("nominal", "ordinal"),  baseline, 
			link=c("logit", "probit", "cloglog", "cauchit", "foldsqrt", "logc"), prop_odds=TRUE, ord_levels, 
			test_type=c("Chisq", "F"), ss_type=c(2, "2", "II", 3, "3", "III"), CI_level=0.95,
			exp_estim=FALSE, testing=TRUE, covar_mat=FALSE, corr_mat=FALSE,
			fitted_value=FALSE, resid_orig=FALSE, resid_pearson=FALSE, linear_pred=FALSE, hat_value=FALSE, num_digit=3){
	#### Arguments ####
	## res_var : GUI의 '반응변수', 필수로 1개 필요
	## quan_var : 양적변수(2개 이상의 변수가 있는 경우 +로 구분, 필수아님. 아무변수도 선택되지 않으면 1로 코딩, 질적변수와 중복되어 사용될 수 없음)
	## qual_var : 질적변수(2개 이상의 변수가 있는 경우 +로 구분, 필수아님. 아무변수도 선택되지 않으면 1로 코딩, 양적변수와 중복되어 사용될 수 없음)
	## inter_var : 교호작용 변수
	## link: link function / GUI의 '연결함수'와 대응
	## res_type: the type of the response variable / GUI의 '명목형 (nominal)' 및 '순서형 (ordinal)'과 대응
	## baseline: baseline category if the response is nominal / GUI의 ''과 대응
	## prop_odds: GUI의 'proportional odds 모형 사용'과 대응
	## ord_levels: order of response variable categories: GUI의 '변수 순서 지정'과 대응, 반드시 낮은 것부터 와야 함
	## testing: testing whether each parameter is equal to 0 / GUI의 '변수적합성 검정'에 대응
	## ss_type: type of chi-squared statistics / GUI의 '통계량 유형'과 대응
	## test_type: type of hypothesis testing / GUI의 '검정 방법'과 대응
	## CI_level: confidence level of the confidence interval / GUI의 '신뢰수준'에 대응
	## exp_estim: parameter estimates & exponential mode / GUI의 '모수 추정값' 및 '지수 모수 추정값'에 대응
	## covar_mat: corr_mat: covariance(correlation) matrix of parameter estimates / GUI의 '추정값 공분산 행렬' 및 '추정값 상관행렬'에 대응
	## fiited_value: print/save fitted values if TRUE / GUI의 '적합값'에 대응
	## resid_orig: print/save response residuals if TRUE / GUI의 '기본 확률 잔차'에 대응
	## resid_pearson: print/save Pearson residuals if TRUE / GUI의 'Pearson 잔차'에 대응
	## linear_pred: print/save linear predictors for each observations / GUI의 '선형 예측값'에 대응
	## hat_value: diagonals of hat matrix / GUI의 'Hat 행렬 대각값'에 대응
	## num_digit: number of output digits
	###################

	###################
	#### Required packages : R2HTML, VGAM, car, AICcmodavg ####
	###################
	load.pkg(c("R2HTML", "VGAM", "car", "AICcmodavg"))

	html.output <- capture.output({	
		R2HTML::HTML(R2HTML::as.title("Generalized Linear Regression for Multinomial Data"), HR=1, file=stdout(), append=FALSE) ;
		getOutput <- FALSE
		
		inter_var	<- vars[sapply(vars, function(s) grepl(":", s))]

		## Warnings
		if(any(table(dataset[,res_var]) <= 1)) warn.msg1 <- "\a Warning : Response variable seems not to be categorized."
		if(!is.null(quan_var)) {
			is.nom <- sapply(quan_var, function(i) !is.numeric(dataset[,i]))
			if(any(is.nom)){
				warn.msg3 <- paste0("\a Warning : The type of variable '",paste(quan[is.nom],collapse=', '),"' is not numeric but selected as the quantitative variable. It was excluded from the analysis.")
				# interaction
				ec_inter <- unique(do.call(c, lapply(quan_var[is.nom], function(i) grep(i, inter_var))))
				if(length(ec_inter) != 0) inter_var <- inter_var[-ec_inter]

				# main effect
				quan_var <- quan[!is.nom] ;
			}
		}
		if((!prop_odds) & (length(qual_var) != 0)) warn.msg4 <- "\a Error : Non-proportional odds model is not supported for qualitative explanatory variable. Analysis has been stopped."
		if(!is.null(qual_var)) {
			is.num <- sapply(qual_var, function(i) is.numeric(dataset[,i]))
			if(any(is.num)) warn.msg5 <- paste0("\a Warning : The type of variable '", paste(qual_var[is.num],collapse=', '), "' is numeric but selected as the qualitative variable. It was coreced into character.")
			for(i in qual_var) dataset[,i] <- factor(dataset[,i])
		}

		if(exists("warn.msg4")){
			R2HTML::HTML(R2HTML::as.title("Warnings"), HR=2, file=stdout()) ;
			R2HTML::HTML(warn.msg4, HR=2, file=stdout())
		} else {
			ss_type	<- as.character(ss_type) ;

			link		<- match.arg(link) ;
			test_type	<- match.arg(test_type) ;
			ss_type		<- match.arg(ss_type) ;

			form.0 <- ifelse(is.null(c(quan_var,qual_var)), paste0(res_var,' ~ 1'),  
				paste0(res_var, ' ~ ', paste(c(quan_var, qual_var, inter_var), collapse=' + '))) ;
			form.1 <- as.formula(form.0)
			if (missing(dataset)) dataset <- environment(form.1) ;

			dist	<- ifelse(res_type == "nominal", multinomial, cumulative) ;
			
			if (!missing(ord_levels)) {
				dataset[[res_var]]	<- factor(dataset[[res_var]], levels=ord_levels, ordered=TRUE) ;
			} else if (res_type == "ordinal") {
				dataset[[res_var]]	<- factor(dataset[[res_var]], ordered=TRUE) ;
			} else {
				dataset[[res_var]]	<- factor(dataset[[res_var]], ordered=FALSE) ;
			}

			alllevel	<- levels(dataset[[res_var]]) ;
			
			if (res_type == "nominal") {
				if(missing(baseline)) baseline	<- alllevel[length(alllevel)] ;
				fit_glm	<- suppressWarnings(vglm(form.1, data=dataset, family=dist(refLevel = baseline)))
				
				printlevel	<- alllevel[alllevel != baseline] ;
			} else {
				link_actual	<- link ;
				command_str	<- paste0("fit_glm<- vglm(form.1, data=dataset, family=dist(link=", link_actual, ", parallel=prop_odds))") ;
				suppressWarnings(eval(parse(text=command_str)))
				
				printlevel	<- alllevel[-length(alllevel)] ;
			}
			
			
			variables	<- as.character(form.1) ;
			
			var_list	<- unlist(strsplit(variables[3], c(" \\* | \\+ |:"))) ;
			var_list	<- unique(var_list)
			if(variables[3] == "1") {
				summ_var	<- 1 ;
				var_lists	<- variables[2] ;
			} else {
				var_list	<- var_list[var_list != "1"] ;
				var_lists	<- c(variables[2], var_list) ;
				summ_var	<- length(var_lists) ;
			}
			data_str	<- matrix(c("Number of observations", "Number of variables", nobs(fit_glm), summ_var), ncol=2) ;

			R2HTML::HTML(R2HTML::as.title("Data Structure"), HR=2, file=stdout()) ;
			R2HTML::HTML(data_str, file=stdout(), row.names=FALSE, col.names=TRUE, innerBorder = 1, align="left") ;

			## Variable list
			level <- paste(alllevel, collapse=ifelse(res_type=='nominal', ', ', ' < '))
			if(res_type=='nominal') level <- gsub(baseline,paste0(baseline,' (baseline)'), level, fixed=TRUE)
			VL <- matrix(c('Response variable',res_var,'Type of response variable',res_type,'Levels of response variable',level,
					'Explanatory variable (Qualitative)',ifelse(is.null(qual_var),'None',paste(qual_var, collapse=", ")),
					'Explanatory variable (Quantitative)',ifelse(is.null(quan_var),'None',paste(quan_var, collapse=", "))),ncol=2,byrow=T)
			R2HTML::HTML(R2HTML::as.title("Variable List"), HR=2, file=stdout()) ;
			R2HTML::HTML(VL, file=stdout(), row.names=FALSE, col.names=TRUE, innerBorder = 1, align="left") ;
			if(exists('warn.msg1')) R2HTML::HTML(warn.msg1, file=stdout(), row.names=FALSE)
			if(exists('warn.msg3')) R2HTML::HTML(warn.msg3, file=stdout(), row.names=FALSE)
			if(exists('warn.msg5')) R2HTML::HTML(warn.msg5, file=stdout(), row.names=FALSE)

			if (res_type == "nominal") {
				fit_print	<- matrix(c(paste0("log(P(", variables[2], "=y)/P(", variables[2], "=", baseline, ")) = ", variables[3]), "Logit"))
			} else  {
				link_print	<- switch(link,
								logit = "Logit",
								probit = "Probit",
								cloglog = "Complementary log-log", 
								cauchit = "Cauchit", 
								foldsqrt = "Folded squareroot", 
								logc = "Complementary log") ;
				fit_print	<- matrix(c(paste0(link, "(P(", variables[2], "<=j)) = ", variables[3]), link_print))
			}
			fit_print	<- cbind(c("Model", "Link function"), fit_print) ;
			colnames(fit_print)	<- c("Category", "Value") ;
				
			R2HTML::HTML(R2HTML::as.title("Model Information"), HR=2, file=stdout()) ;
			R2HTML::HTML(fit_print, file=stdout(), row.names=FALSE, col.names=TRUE, innerBorder = 1, align="left") ;

			# Summary stat	
			summ_print <- summary(dataset[, var_lists, drop=F]) ;
			summ_print <- data.frame(unclass(summ_print), stringsAsFactors=FALSE) ;
			colnames(summ_print)[1] <- variables[2] ; colnames(summ_print)[-1] <- var_list ;
			colnames(summ_print) <- gsub(" ", "", colnames(summ_print), fixed=TRUE) ;
			R2HTML::HTML(R2HTML::as.title("Summary Statistics"), HR=2, file=stdout()) ;
			R2HTML::HTML(summ_print, file=stdout(),row.names=FALSE,na="", innerBorder = 1, align="left") ;

			### Analysis results
			R2HTML::HTML(R2HTML::as.title("Results of Generalized Linear Regression for Multinomial Data"), HR=2, file=stdout()) ;
			
			## Goodness of fit
			fitness_print	<- data.frame(numeric(0), numeric(0)) ;
			fitness_print	<- rbind(fitness_print, c(deviance(fit_glm), fit_glm@df.residual)) ;

			chisq_pearson	<- sum(residuals(fit_glm, type="pearson")^2) ;
			fitness_print	<- rbind(fitness_print, c(chisq_pearson, fit_glm@df.residual)) ;

			fitness_print	<- rbind(fitness_print, c(as.numeric(logLik(fit_glm)), fit_glm@df.total - fit_glm@df.residual)) ;
				
			fitness_print	<- rbind(fitness_print, c(AIC(fit_glm), NA)) ;
			fitness_print	<- rbind(fitness_print, c(AICc(fit_glm), NA)) ;
			fitness_print	<- rbind(fitness_print, c(BIC(fit_glm), NA)) ;
			names(fitness_print)	<- c("Value", "DF") ;
			row.names(fitness_print)	<- c("Deviance", "Pearson's chi-square", "log-likelihood", "AIC", "AICc", "BIC") ;
			fitness_print[1]	<- format(round(fitness_print[1], num_digit), nsmall=num_digit) ;

			R2HTML::HTML(R2HTML::as.title("Model Fitness Measurements"), HR=3, file=stdout()) ;
			R2HTML::HTML(fitness_print, file=stdout(), na="", innerBorder = 1, align="left") ;

			coef_print	<- data.frame(coef(summary(fit_glm))) ;
			coef_print	<- cbind(coef_print, confint.default(fit_glm, level=CI_level)) ;
			if(exp_estim){
				coef_print <- data.frame(exp(coef_print[,1]),coef_print)
				coef_print	<- cbind(coef_print, exp(confint.default(fit_glm, level=CI_level)))
				colnames(coef_print)	<- c("exp(Estimate)","Estimate","Std.Error","Z-value","P-value",paste0(100*CI_level, "% CI of Est. ", c("(lower)", "(upper)")),paste0(100*CI_level, "% CI exp(Est.) ", c("(lower)", "(upper)")))
			} else {
				colnames(coef_print)	<- c("Estimate","Std.Error","Z-value","P-value",paste0(100*CI_level, "% CI of Est. ", c("(lower)", "(upper)")))
			}
			for (i in seq(nrow(coef_print))) {
				rown_tmp	<- unlist(strsplit(row.names(coef_print)[i], ":", fixed=TRUE)) ;
				if (rown_tmp[1] == "(Intercept)") rown_tmp[1] <- "Intercept" ;
				if (length(rown_tmp) < 2) next ;
				rown_new	<- paste(rown_tmp, collapse=": ") ;
				row.names(coef_print)[i]	<- rown_new ;
			}
				
			R2HTML::HTML(R2HTML::as.title("Regression Coefficient Estimates and Confidence Intervals"), HR=3, file=stdout()) ;
			R2HTML::HTML(format(round(coef_print, num_digit), nsmall=num_digit), file=stdout(), na="", innerBorder = 1, align="left") ;	
			
			if (testing) {
				R2HTML::HTML(R2HTML::as.title(paste0("Testing Regression Coefficients using Type ",ss_type," SS")), HR=3, file=stdout()) ;

				if (variables[3]=="1"){
					R2HTML::HTML("Testing is not available when the model contains only an intercept.", file=stdout(), na="")
				} else {
						table_tmp	<- data.frame(Anova(fit_glm, test=test_type, type=ss_type)) ;
						colnames(table_tmp)	<- c("DF", ifelse(test_type == "Chisq", "Chi-squared", "F"), "P-value") ;
							
						if (ss_type %in% c(3, "3", "III")) {
							test_print	<- table_tmp[-1,] ;
						} else {
							test_print	<- table_tmp ;
						}
						
						test_print[-1]	<- format(round(test_print[-1], num_digit), nsmall=num_digit) ;
				}
				if (exists("test_print")) R2HTML::HTML(as.data.frame(test_print), file=stdout(), na="", innerBorder = 1, align="left") ;
				
			}

			if (covar_mat) {
				covar_print	<- vcov(fit_glm) ;
				
				for (i in seq(nrow(covar_print))) {
					rown_tmp	<- unlist(strsplit(row.names(covar_print)[i], ":", fixed=TRUE)) ;
					if (rown_tmp[1] == "(Intercept)") rown_tmp[1] <- "Intercept" ;
					if (length(rown_tmp) < 2) next ;
					rown_new	<- paste(rown_tmp, collapse=": ") ;
					row.names(covar_print)[i]	<- colnames(covar_print)[i]	<- rown_new ;
				}
				
				R2HTML::HTML(R2HTML::as.title("Covariance Matrix for Estimates of Regression Coefficients"), HR=3, file=stdout()) ;
				R2HTML::HTML(format(round(covar_print, num_digit), nsmall=num_digit), file=stdout(), na="", innerBorder = 1, align="left") ;
			}
			
			if (corr_mat) {
				corr_print	<- summary(fit_glm, correlation=TRUE)@correlation ;
				
				for (i in seq(nrow(corr_print))) {
					rown_tmp	<- unlist(strsplit(row.names(corr_print)[i], ":", fixed=TRUE)) ;
					if (rown_tmp[1] == "(Intercept)") rown_tmp[1] <- "Intercept" ;
					if (length(rown_tmp) < 2) next ;
					rown_new	<- paste(rown_tmp, collapse=": ") ;
					row.names(corr_print)[i]	<- colnames(corr_print)[i]	<- rown_new ;
				}
				
				R2HTML::HTML(R2HTML::as.title("Correlation Matrix for Estimates of Regression Coefficients"), HR=3, file=stdout()) ;
				R2HTML::HTML(format(round(corr_print, num_digit), nsmall=num_digit), file=stdout(), na="", innerBorder = 1, align="left") ;
			}

			if (fitted_value) {
				MultinomGLM_output	<- data.frame(fitted(fit_glm)) ;
				colnames(MultinomGLM_output)	<- paste0("MultinomGLM_Fitted_", colnames(fit_glm@y)) ;

				getOutput <- TRUE
			}
			
			if (resid_orig) {
				if (exists("MultinomGLM_output")) {
					MultinomGLM_output	<- cbind(MultinomGLM_output, residuals(fit_glm, "response")) ;
					colnames(MultinomGLM_output)[ncol(MultinomGLM_output)-(length(colnames(fit_glm@y))-1):0]	<- paste0("MultinomGLM_ResidOriginal_", alllevel) ;
				} else {
					MultinomGLM_output	<- data.frame(residuals(fit_glm, "response")) ;
					colnames(MultinomGLM_output)	<- paste0("MultinomGLM_ResidOriginal", colnames(fit_glm@y)) ;
				}
				getOutput <- TRUE
			}
			
			if (resid_pearson) {
				if (exists("MultinomGLM_output")) {
					MultinomGLM_output	<- cbind(MultinomGLM_output, residuals(fit_glm, "pearson")) ;
				} else {
					MultinomGLM_output	<- residuals(fit_glm, "pearson") ;
				}
				
				names_bef	<- colnames(MultinomGLM_output)[ncol(MultinomGLM_output)-(length(colnames(fit_glm@y))-2):0] ;
				names_aft	<- paste0("MultinomGLM_ResidPearson_", names_bef) ;
					
				for (i in seq(length(colnames(fit_glm@y)), 1, -1)) {
					if (res_type == "nominal") {
						names_aft	<- gsub(paste0(",", as.character(i)), colnames(fit_glm@y)[i], names_aft, fixed=TRUE) ;
					} else {
						names_aft	<- gsub(paste0("Y<=", as.character(i)), paste0(res_var, "<=", colnames(fit_glm@y)[i]), names_aft, fixed=TRUE) ;
					}
				}	
				colnames(MultinomGLM_output)[ncol(MultinomGLM_output)-(length(colnames(fit_glm@y))-2):0]	<- names_aft ;
				
				getOutput <- TRUE
			}
			
			if (linear_pred) {
				if (exists("MultinomGLM_output")) {
					MultinomGLM_output	<- cbind(MultinomGLM_output, fit_glm@predictors) ;
				} else {
					MultinomGLM_output	<- fit_glm@predictors ;
				}
				
				names_bef	<- colnames(MultinomGLM_output)[ncol(MultinomGLM_output)-(length(colnames(fit_glm@y))-2):0] ;
				names_aft	<- paste0("MultinomGLM_LinearPred_", names_bef) ;
					
				for (i in seq(length(colnames(fit_glm@y)), 1, -1)) {
					if (res_type == "nominal") {
						names_aft	<- gsub(paste0(",", as.character(i)), colnames(fit_glm@y)[i], names_aft, fixed=TRUE) ;
					} else {
						names_aft	<- gsub(paste0("Y<=", as.character(i)), paste0(res_var, "<=", colnames(fit_glm@y)[i]), names_aft, fixed=TRUE) ;
					}
				}	
				colnames(MultinomGLM_output)[ncol(MultinomGLM_output)-(length(colnames(fit_glm@y))-2):0]	<- names_aft ;
				
				getOutput <- TRUE
			}
			
			if (hat_value) {
				if (exists("MultinomGLM_output")) {
					MultinomGLM_output	<- cbind(MultinomGLM_output, hatvalues(fit_glm)) ;
				} else {
					MultinomGLM_output	<- hatvalues(fit_glm) ;
					attr(MultinomGLM_output, "predictors.names")	<- NULL ;
					attr(MultinomGLM_output, "ncol.X.vlm")			<- NULL ;
				}
				
				names_bef	<- colnames(MultinomGLM_output)[ncol(MultinomGLM_output)-(length(colnames(fit_glm@y))-2):0] ;
				names_aft	<- paste0("MultinomGLM_HatValue_", names_bef) ;
					
				for (i in seq(length(colnames(fit_glm@y)), 1, -1)) {
					if (res_type == "nominal") {
						names_aft	<- gsub(paste0(",", as.character(i)), colnames(fit_glm@y)[i], names_aft, fixed=TRUE) ;
					} else {
						names_aft	<- gsub(paste0("Y<=", as.character(i)), paste0(res_var, "<=", colnames(fit_glm@y)[i]), names_aft, fixed=TRUE) ;
					}
				}		
				colnames(MultinomGLM_output)[ncol(MultinomGLM_output)-(length(colnames(fit_glm@y))-2):0]	<- names_aft ;
				
				getOutput <- TRUE
			}
		}
		R2HTML::HTMLhr(file=stdout()) ;
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Generalized Linear Regression for Multinomial Data",sep=""),file=stdout()) ;
		R2HTML::HTMLhr(file=stdout()) ;
	})
	
	if(getOutput) {
		row.dataset <- rownames(dataset)
		row.output <- rownames(MultinomGLM_output)
		O <- matrix(NA,nrow=nrow(dataset),ncol=ncol(MultinomGLM_output))
		for(i in row.dataset) if(i%in%row.output) O[row.dataset==i,] <- as.numeric(MultinomGLM_output[row.output==i,])
		O <- as.data.frame(O)
		colnames(O) <- colnames(MultinomGLM_output)
		return(list(html=html.output, Output=O))
	} else {
		return(html.output)
	}
}

REx_LGM <- function(dataset,dep_var,dep_ref=NULL,indep_numeric_var=NULL,indep_cat_var=NULL,vars=NULL,indep_ref=NULL,indep_baseline=NULL,noint=FALSE,link='logit',
			CI=FALSE,confint.level=0.95,odds=FALSE,Validation=FALSE,k=10,ANOVA=TRUE,ss="I",GOF=TRUE,ng=10,R2=TRUE,classtab=TRUE,cutpoint=.5,VIP=FALSE,corInform=FALSE,
			Predict_prob=FALSE,Predict_g=FALSE,Predict_g_criteria=.5,Resid=FALSE,stdResid=FALSE,
			select=FALSE,keep_var=NULL,direct="forward",Plot=TRUE){		
		
	#### 변수 설명 ####
	## dep_var : 종속변수(dependent variables),필수변수(무조건 한개변수만 선택할수있음,분석적에 자동 factor 처리 필요,
	## dep_ref : 종속변수 기저범주
	## indep_numeric_var : 양적설명변수(independent variables) (2개 이상의 변수가 있는 경우 +로 구분, 필수사항아님,단,상수항옵션이(noint=TRUE)인경우 돌아갈수없음)
	## indep_cat_var : 질적설명변수(independent variables) (자동factor처리,2개 이상의 변수가 있는 경우 +로 구분, 필수사항아님,단,상수항옵션이(noint=TRUE)인경우 돌아갈수없음)
	## vars : 선택된 변수들 (2개 이상의 변수가 있는 경우 c(,,)로 구분, 필수사항아님)
	## indep_ref : 선택된 범주형 설명변수 중 기저범주 바꾸고 싶은 범주형 설명변수 (필수아님,2개이상의 변수가 있는 경우 c(,,)로 구분)
	## indep_baseline : indef_ref에서 선택된 범주형 설명변수별로 바꾸고 싶은 기저값 (indep_ref=NULL이 아닌경우 무조건 필수 아니면 indep_baseline=NULL, 2개이상의 값이 있는 경우 c(,,)로 구분)
	## indep_ref=NULL :기본설정된 순서로 기저범주 설정됨
	## indep_baseline=NULL: 기본설정된 순서로 기저범주 설정됨
	## noint : 모델설정 탭의 '상수항 포함하지 않음' 옵션
	## link : 'logit', 'probit', 'cauchit', 'log', 'cloglog'
	## CI : 출력옵션 탭의 신뢰구간
	## confint.level : 통계 설정 탭의 '신뢰수준' 옵션 (신뢰수준범위 변경가능,0.95 default)
	## odds : exp(b) 즉 오즈비를 출력할것인지 말것인지 결정하는 옵션 (TRUE/FALSE, default=FALSE)
	## Validation : k-fold cross validation 선택 옵션 (TRUE/FALSE)
	## k : k-fold cross validation에서 분할의 개수k (default, k=10)
	## ANOVA : 출력옵션 탭의 'ANOVA Table' 옵션 (TRUE/FALSE)
	## ss : 출력옵션 탭의 "SS type". "I","II","III" 중 하나를 선택할 수 있음.
	## GOF : 적합도 검정(goodness of fit test), Hosmer-Lemeshow Test (TRUE/FALSE)
	## ng : 적합도 검정(GOF)할때 number of groups for C and H statistic.(defalut, ng=10)
	## R2 : pseudoR2옵션(TRUE/FALSE,defalt=TRUE)
	## classtab : outcome의 실제값과 예측값의 분류표 출력 옵션 (TRUE/FALSE,default=TRUE)
 	## cutpoint : 예측값의 0/1의 분류 기준 (default=05, 기준점 변경가능)
	## VIP : 변수의 중요도 (variable importance) (TRUE/FALSE)
	## corInform : 추정값들의 상관계수 출력 옵션(TRUE/FALSE)
	## Predict_prob : 출력옵션 탭의 '예측값' 옵션에서 확률 (TRUE/FALSE)
	## Predict_g : 출력옵션 탭의 '예측값' 옵션에서 소속집단 (TRUE/FALSE)
	## Predict_g_criteria : 출력옵션 탭의 '예측값' 옵션에서 소속집단의 기준 (default, p=.5)
	## Resid : 출력옵션 탭의 '잔차' 옵션
	## stdResid: 출력옵션 탭의 '표준화잔차' 옵션
	## select : 통계설정 탭의 '변수선택' 옵션(2개 이상의 변수가 있는 경우 +로 구분)
	## keep_var : 질적변수 or 양적변수 상관없이 모든 설명변수가 대상임, 만약에 설명변수에 변수가 설정이 안되있으면 활성화되어도 설명변수창에 아무것도 뜨지 않게해야됨
	## 통계설정 탭의 '고정변수선택' 옵션 (설정없으면(select=FALSE,default) "" or 설정(select=TRUE),2개 이상의 변수가 있는 경우 +로 구분)
	## direct : 통계설정 탭의 '변수선택의 세부항목' 옵션("forward","backward","both",가능)
	## Plot : 출력옵션 탭의 '도표' 옵션(TRUE/FALSE,TRUE:default)##
	###################
	#### Required packages : R2HTML , rms, MASS, caret,pscl, ResourceSelection, MKmisc, lmtest,popbio ####
	###################
	load.pkg(c("R2HTML", "rms", "MASS", "caret", "pscl", "ResourceSelection", "MKmisc", "lmtest", "VGAM", "popbio"))

	html.output <- capture.output({
		# Title
		R2HTML::HTML(R2HTML::as.title("Generalized Linear Regression for Binomial Data"),HR=1,file=stdout(),append=FALSE)

		## Warnings
		# Response variable type
		dataset[,dep_var] <- as.factor(dataset[,dep_var])
		if(length(levels(factor(dataset[,dep_var])))!=2) warn.msg1 <- '\a Error : Response variable should be binary. Analysis has been stopped.'		
		# explanatory variable type
		if(!is.null(indep_numeric_var)) {
			is.nom <- sapply(indep_numeric_var,function(i) !is.numeric(dataset[,i]))
			if(any(is.nom)){
				warn.msg2 <- paste0("\a Warning : The type of variable '",paste(indep_numeric_var[is.nom],collapse=', '),"' is not numeric but selected as the quantitative variable. It was excluded from the analysis.")
				ec <- indep_numeric_var[is.nom]
				vars <- vars[-grep(paste0('[',paste(ec,collapse=','),']'),vars)]
				if(length(vars)==0) vars <- NULL

			}
		}
		if(!is.null(indep_cat_var)) {
			is.num <- sapply(indep_cat_var,function(i) is.numeric(dataset[,i]))
			if(any(is.num)) warn.msg3 <- paste0("\a Warning : The type of variable '",paste(indep_cat_var[is.num],collapse=', '),"' is numeric but selected as the qualitative variable. It was coreced into character.")
			for(i in indep_cat_var[is.num]) dataset[,i] <- as.factor(dataset[,i])
		}

		# no intercept & no explanatory variable : error
		if(noint&is.null(vars)) warn.msg4 <- '\a Error : With no intercept, at least 1 independent variable should be selected. Analysis has been stopped.'


		# warnings
		if(exists('warn.msg1')|exists('warn.msg4')){
			R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
			if(exists('warn.msg1')) R2HTML::HTML(warn.msg1,file=stdout())
			if(exists('warn.msg4')){
				if(exists('warn.msg2')) R2HTML::HTML(warn.msg2,file=stdout())
				R2HTML::HTML(warn.msg4,file=stdout())
			}
		} else {
			## Data Processing
			var_info <- c(dep_var,indep_numeric_var, indep_cat_var)
			is.na <- sapply(var_info,function(i) is.na(dataset[,i]))
			raw.dat <- dataset
			oo <- which(rowSums(is.na)>0)
			if(length(oo)>0) dataset <- dataset[-oo,]
			## To define baseline for dep/indep cat
			if(!is.null(dep_ref))dataset[,dep_var] <- relevel(factor(dataset[,dep_var]),ref=dep_ref)
			if(!is.null(indep_ref)) for(i in seq_along(indep_ref)) dataset[,indep_ref[i]] <- relevel(factor(dataset[,indep_ref[i]]),ref=indep_baseline[i])
			## model formula
			form.0 <- ifelse(is.null(vars),paste0(dep_var,' ~ 1'),paste0(dep_var,' ~ ',paste(vars,collapse=' + ')))
			if(noint) form.0 <- paste(form.0,-1)
			form.1 <- as.formula(form.0)
			
			## Fitting LGM
			command_str	<- paste0("res_LGM_1<- glm(formula(form.1), data=dataset, family=binomial(",link,"))")
			eval(parse(text=command_str)) ;
			res_LGM <- res_LGM_1
			# res_LGM_1 <- glm(formula(form),data=dataset,family=binomial(link))
			
			## Variable selection
			if(select){
				## 고정변수 선택여부에따른 분석방법(선택안한경우)
				f = file()
				sink(file=f)
				if(is.null(keep_var)){
					res_LGM2 <- stepAIC(res_LGM_1,direction=direct)
					form2 <- formula(res_LGM2)
					command_str	<- paste0("res_LGM <- glm(form2, data=dataset, family=binomial(",link,"))")
					eval(parse(text=command_str)) ;
				## 고정변수 선택여부에따른 분석방법(선택한경우)	
				}else{
					form1 <- formula(paste('~',paste(keep_var,collapse='+'),sep=''))
					res_LGM2 <- stepAIC(res_LGM_1,direction=direct,scope=list(lower=form1))
					form2 <- formula(res_LGM2)
					command_str	<- paste0("res_LGM <- glm(form2, data=dataset, family=binomial(",link,"))")
					eval(parse(text=command_str)) ;
				}
				sink()
				close(f)
			}

			# 다음의 값들은 엑셀 시트에 저장
			if(Predict_prob) {
				temp.0 <- data.frame(Fitted_LGM=fitted(res_LGM))
				temp <- data.frame(Fitted_LGM=rep(NA,nrow(raw.dat)))
				row.dataset <- rownames(raw.dat)
				row.output <- rownames(temp.0)
				for(i in seq_along(row.output)) if(i%in%seq_along(row.output)) temp[row.dataset==row.output[i],1] <- temp.0[row.output[i],1]				
				O <- temp
			}
			if(Predict_g) {
				temp.0 <- data.frame(Fitted_LGM_g=ifelse(fitted(res_LGM)<Predict_g_criteria,levels(dataset[,dep_var])[1],levels(dataset[,dep_var])[2]))
				temp.0[,1] <- as.character(temp.0[,1])
				temp <- data.frame(Fitted_LGM_g=rep(NA,nrow(raw.dat)))
				row.dataset <- rownames(raw.dat)
				row.output <- rownames(temp.0)
				for(i in seq_along(row.output)) if(i%in%seq_along(row.output)) temp[row.dataset==row.output[i],1] <- temp.0[row.output[i],1]
				if(exists('O')) {
					O <- cbind(O,temp)
				} else {
					O <- temp
				}
			}
			if(Resid) {
				temp.0 <- data.frame(Resid_LGM=resid(res_LGM))
				temp <- data.frame(Resid_LGM=rep(NA,nrow(raw.dat)))
				row.dataset <- rownames(raw.dat)
				row.output <- rownames(temp.0)
				for(i in seq_along(row.output)) if(i%in%seq_along(row.output)) temp[row.dataset==row.output[i],1] <- temp.0[row.output[i],1]
				if(exists('O')) {
					O <- cbind(O,temp)
				} else { 
					O <- temp
				}
			}
 			if(stdResid) {
				temp.0 <- data.frame(stdResid_LGM=studres(res_LGM))
				temp <- data.frame(stdResid_LGM=rep(NA,nrow(raw.dat)))
				row.dataset <- rownames(raw.dat)
				row.output <- rownames(temp.0)
				for(i in seq_along(row.output)) if(i%in%seq_along(row.output)) temp[row.dataset==row.output[i],1] <- temp.0[row.output[i],1]

				if(exists('O')) {
					O <- cbind(O,temp)
				} else { 
					O <- temp
				}
			}

			### Default output
			# Data Structure
			R2HTML::HTML(R2HTML::as.title("Data Structure"),HR=2,file=stdout())
			total.var <- ncol(dataset)
			used.var <- ifelse(is.null(vars),0,length(unique(unlist(strsplit(vars,":")))))
			DS <- matrix(c('Number of observations','Number of total variables','Number of used variables',nrow(dataset),total.var,used.var+1),ncol=2)
			R2HTML::HTML(DS,file=stdout(), innerBorder = 1,align="left")

			# Varibale list
	
			R2HTML::HTML(R2HTML::as.title("Variable List"),HR=2,file=stdout())
			if(is.null(vars)) {
				Vars <- c()
			} else {
				Vars <-  unique(unlist(strsplit(vars,":")))
			}

			qual <- c(indep_cat_var[indep_cat_var%in%Vars],dep_var)
			quan <- indep_numeric_var[indep_numeric_var%in%Vars]
			varlist <- matrix(c('Quantitative variable',paste(quan,collapse=', ')),ncol=2,byrow=T)
			if(length(qual)!=0) varlist <- rbind(varlist,c('Qualitative variable',paste(qual,collapse=', ')))
			R2HTML::HTML(varlist,file=stdout(), innerBorder = 1,align="left")
			if(exists('warn.msg2')) R2HTML::HTML(warn.msg2,file=stdout())
			if(exists('warn.msg3')) R2HTML::HTML(warn.msg3,file=stdout())

			# Analysis Description
			#R2HTML::HTML(R2HTML::as.title("Analysis Description"),HR=2,file=stdout())
			R2HTML::HTML(R2HTML::as.title("Analysis Description"),HR=2,file=stdout())
			AD <- matrix(c('Link function',link,'Response variable',dep_var),ncol=2,byrow=T)
			if(!is.null(vars)) AD <- rbind(AD,c('Explanatory variable',paste(vars,collapse=', ')))
			AD <- rbind(AD,matrix(c('Intercept included',!noint,'Odds Rtio (OR)',odds),ncol=2,byrow=T))

			if(CI) AD <- rbind(AD,c('Significance level of CI',confint.level))
			if(ANOVA)AD <- rbind(AD,matrix(c('Print ANOVA table',paste(ANOVA,"(Type ",ss," SS)",sep="")),ncol=2,byrow=T))
			AD <- rbind(AD,matrix(c('Variable selection',select),ncol=2,byrow=T))			
			#variable selection
			if(select) {
				direction <- switch(direct,forward="Forward selection",backward="Backward elimination",both="Stepwise regression")
				AD <- rbind(AD,matrix(c('Variable selection method',direction),ncol=2,byrow=T))
			}
			if(Validation) AD <- rbind(AD,c('Validation',paste(Validation,"(k=",k,")",sep="")))
			if(GOF) AD <- rbind(AD,c('Goodness of Fit Test',paste(GOF,"(n=",ng,")",sep='')))
			if(classtab) AD <- rbind(AD,c('Classfication Table',paste(classtab,"(cut point=",cutpoint,")",sep='')))
			if(corInform) AD <- rbind(AD,c('Covariance Matrix',corInform))
			if(R2) AD <- rbind(AD,c('pseudo-R2',R2))
			if(VIP) AD <- rbind(AD,c('Variable Importance',VIP))
			R2HTML::HTML(AD,file=stdout(), innerBorder = 1,align="left")

			#### Results
			R2HTML::HTML(R2HTML::as.title("Results of Generalized Linear Regression"),HR=2,file=stdout())
			# Coefficient Estimate
			# odds=TRUE;confint.level=.95
			R2HTML::HTML(R2HTML::as.title("Coefficient Estimates"),HR=3,file=stdout())
			VS1 <- capture.output(summary(res_LGM))
			blank1 <- which(VS1=="")
			IF1 <- paste(VS1[blank1[length(blank1)]-1],collapse="")
			IF1 <- gsub("[ ]{2,}"," ",IF1)
			IF2 <- paste(VS1[blank1[length(blank1)]-3],collapse="")
			IF2 <- gsub("[ ]{2,}"," ",IF2)
			IF3 <- paste(VS1[blank1[length(blank1)]-4],collapse="")
			IF3 <- gsub("[ ]{2,}"," ",IF3)
			IF4 <- paste(VS1[blank1[length(blank1)]-5],collapse="")
			IF4 <- gsub("[ ]{2,}"," ",IF4)
			IF5 <- paste(VS1[blank1[length(blank1)]-7],collapse="")
			IF5 <- gsub("[ ]{2,}"," ",IF5)

			IF <- rbind(IF5,IF4,IF3,IF2,IF1)
			rownames(IF) <- c("","","","","")
			
			CE <- as.data.frame(summary(res_LGM)$coef)
			colnames(CE) <- c('Estimate','SE','Z-value','P-value')
			if(CI) {CE <- cbind(CE,matrix(confint(res_LGM, level=confint.level),ncol=2)); colnames(CE)[(ncol(CE)-1):ncol(CE)] <- paste0(c('Lower','Upper'),' bound of ',confint.level*100,'% CI')}
			if(odds) {ORs <- cbind(exp(summary(res_LGM)$coef[,1]),exp(matrix(confint(res_LGM, level=confint.level),ncol=2))); colnames(ORs) <- c('exp(Estimate)',paste0(c('Lower','Upper'),' bound of ',confint.level*100,'% CI'));CE <- cbind(CE,ORs)}
			R2HTML::HTML(CE,file=stdout(), innerBorder = 1,align="left",digits=4)
			R2HTML::HTML(IF,file=stdout(), align="left",digits=4)
		
			#### Anova table, Goodness of fit, R-squared and classfication table
			# ANOVA=TRUE;GOF=TRUE;ng=10;R2=TRUE;classtab=TRUE;cutpoint=.5
			if(ANOVA==TRUE){
				R2HTML::HTML(R2HTML::as.title(paste0("ANOVA Table with Type ",ss," SS")),HR=3,file=stdout())
				if(length(vars)!=0){
					if(ss=='I'){
						AT <- try(anova(res_LGM,test="Chisq"),s=T)
						rowNames <- rownames(AT)
						if(!noint)rowNames[1] <- ""
						colNames <- colnames(AT)
						if(class(AT)[1]!='try-error'){
							AT <- data.frame(cbind(AT[,1],round(AT[,2],1),AT[,3],round(AT[,4],1),AT[,5]))
							rownames(AT) <- rowNames
							colnames(AT) <- colNames


							if(any(AT[,5]<0.001)){
								oo3 <- which(AT[,5]<0.001)
								AT[oo3,5] <- "<0.001"
							}
							AT[,ncol(AT)] <- format(AT[,ncol(AT)],scientific=T)
							AT[1,c(1:2,5)] <- "" 
							R2HTML::HTML(AT,file=stdout(), innerBorder = 1,align="left",digits=4)
						} else {
							warn.msg6 <- "\a Error : Fail to fit the Model."
							R2HTML::HTML(warn.msg6,file=stdout())
						}
					}

					if(ss %in% c('II','III')){
						if(ss=='II') res<-try(car::Anova(res_LGM,type='II'))
						if(ss=='III'){ options(contrasts=c("contr.sum", "contr.poly"))
									res<-try(car::Anova(res_LGM,type='III'))}
						if(class(res)[1]!='try-error'){
							AT <- data.frame(res)
							AT[,ncol(AT)] <- format(AT[,ncol(AT)],scientific=T,digits=4)
							colnames(AT) <- c('LR_Chisq','DF','Pr(>Chisq)')
							R2HTML::HTML(AT,file=stdout(), innerBorder = 1,align="left",digits=4)
						} else {
							warn.msg7 <- "\a Error : Fail to fit the Model."
							R2HTML::HTML(warn.msg7,file=stdout())
						}
					}
				} else {
					warn.msg8 <- "\a Warning : ANOVA table is not supported for intercept only model."
					R2HTML::HTML(warn.msg8,file=stdout())
				}
			}

			if(GOF==TRUE){
				R2HTML::HTML(R2HTML::as.title("Goodness of Fit Test"),HR=3,file=stdout())
				if(length(vars)!=0){

					## loglikelhodd ratio
					#R2HTML::HTML(R2HTML::as.title("Analysis of Deviance Table (Likelhood Ratio Test)"),HR=3,file=stdout())
					GOF1 <- anova(glm(formula(paste(dep_var,'~1',sep='')),data=dataset,family="binomial"), res_LGM, test ="Chisq")
					#VS3 <- capture.output(GOF1)
					#IF10 <- paste(VS3[1],collapse="")
					#IF11 <- paste(VS3[2],collapse="")
					#IF12 <- paste(VS3[3],collapse="")
					#IF13 <- paste(VS3[4:(length(VS3)-5)],collapse="")
					#IF13 <- gsub("[ ]{2,}"," ",IF13)
					#IF14 <- rbind(IF10,IF11,IF12,IF13)
					#rownames(IF14) <- rep("",nrow(IF14))
					
					A3 <- as.matrix(GOF1)
					A3[1,] <- as.numeric(A3[1,])
					A3[2,] <- as.numeric(A3[2,])
					A3[,c(2,4)] <- round(A3[,c(2,4)],2)
					A3[,ncol(A3)] <- format(A3[,ncol(A3)],scientific=T,digits=4)
					A3[1,3:5] <- ""
					rownames(A3) <- c("Null Model","Proposed Model")
					R2HTML::HTML(R2HTML::as.title("Analysis of Deviance Table (Likelhood Ratio Test)"),HR=4,file=stdout())
					R2HTML::HTML(A3,file=stdout(),innerBorder = 1,align="left",digits=4)

					##Hosmer-Lemeshow Test
					GOF2 <- hoslem.test(res_LGM$y,fitted(res_LGM), g=ng)
					#HLgof.test(fit = fitted(res_LGM), obs = as.numeric(dataset[,dep_var])-1,ngr=ng)

					A4 <- data.frame(round(GOF2$statistic,2),GOF2$parameter,round(GOF2$p.value,4))
					colnames(A4) <- c("X-squared","df","p-value")
					rownames(A4) <- ""
					A4[1,] <- as.character(A4[1,])
					A5 <- cbind(GOF2$observed,GOF2$expected)
					colnames(A5) <-paste(c("observed_","observed_","expected_","expected_"),c(rep(levels(dataset[,dep_var]),2)),sep="") 
										
					R2HTML::HTML(R2HTML::as.title("Hosmer-Lemeshow test"),HR=4,file=stdout())
					R2HTML::HTML(A4,file=stdout(), align="left", innerBorder = 1,row.names=F,digits=4)
					R2HTML::HTML(R2HTML::as.title("Contingency Table for Hosmer-Lemeshow Test"),HR=4,file=stdout())
					R2HTML::HTML(A5,file=stdout(), align="left",innerBorder = 1,digits=4)
				} else {
					warn.msg12 <- "\a Warning : Goodness of fit test is not supported for intercept only model."
					R2HTML::HTML(warn.msg12,file=stdout())
				}
			}

			if(R2==TRUE){
				R2HTML::HTML(R2HTML::as.title("Pseudo R-squared Measures"),HR=3,file=stdout())		
				# look for 'McFadden'
				A6 <- t(data.frame(round(pR2(res_LGM),3)))
				rownames(A6) <- ""
				IF15 <- "llh     : The log-likelihood from the fitted model"
				IF16 <- "llhNull : The log-likelihood from the intercept-only restricted model"
				IF17 <- "G2      : Minus two times the difference in the log-likelihoods"
				IF18 <- "McFadden: McFadden's pseudo r-squared"
				IF19 <- "r2ML    : Maximum likelihood pseudo r-squared"
				IF20 <- "r2CU    : Cragg and Uhler's pseudo r-squared"
				IF21 <- rbind(IF15,IF16,IF17,IF18,IF19,IF20)
				rownames(IF21) <- rep("",nrow(IF21))

				R2HTML::HTML(data.frame(A6),file=stdout(), innerBorder = 1,align="left",row.names=F,digits=4)
				# R2HTML::HTML(IF21,file=stdout(),align='left',digits=4) # Manual
			}

			# classification table
			# classtab=TRUE;cutpoint=.5
			if(classtab==TRUE){
				R2HTML::HTML(R2HTML::as.title(paste0("Classification Table - cut point : ",cutpoint)),HR=3,file=stdout())		
				Observed <- res_LGM$y
				level <- levels(dataset[,dep_var])
				Predicted <- as.factor(ifelse(fitted(res_LGM)<cutpoint,level[1],level[2]))
				Predicted <- factor(Predicted,levels=level)
				A7 <- table(Observed,Predicted)
				rownames(A7) <- level
				Total <- colSums(A7)
				A7 <- rbind(A7,Total)
				Total <- rowSums(A7)
				A7 <- cbind(A7,Total)						
				colnames(A7) <- paste("Predicted_",colnames(A7),sep="")
				rownames(A7) <- paste("Observed_",rownames(A7),sep="")
				R2HTML::HTML(A7,file=stdout(), innerBorder = 1,align="left",digits=4)
			}

			#### Variable importance(VIF), covariance matrix,and validation
			# VIP=TRUE;corInform=TRUE;Validation=TRUE;K=10
			if(VIP==TRUE){
				if(length(vars)!=0){
					R2HTML::HTML(R2HTML::as.title("Variable Importance"),HR=3,file=stdout())		
					vip <- round(varImp(res_LGM),2)
					R2HTML::HTML(vip,file=stdout(), innerBorder = 1,align="left",digits=4)
				} else {
					warn.msg9 <- "\a Warning : Variable importance table is not supported for intercept only model."
					R2HTML::HTML(warn.msg9,file=stdout())
				}
			}

			if(corInform==TRUE){
				R2HTML::HTML(R2HTML::as.title("The Estimated Correlations of the Estimated Coefficients"),HR=3,file=stdout())		
				Corr1 <- round(summary(res_LGM)$cov.unscaled,3)
				#Corr2 <- round(summary(res_LGM)$cov.scaled,3)
				#R2HTML::HTML(R2HTML::as.title("Unscaled (Dispersion = 1)"),HR=4,file=stdout())		
				R2HTML::HTML(Corr1,file=stdout(), innerBorder = 1,align="left",digits=4)
				#R2HTML::HTML(R2HTML::as.title("Scaled by Dispersion"),HR=4,file=stdout())		
				#R2HTML::HTML(Corr2,file=stdout(), innerBorder = 1,align="left")
			}

			if(Validation==TRUE){
				R2HTML::HTML(R2HTML::as.title(paste0(k,"-Fold Cross Validation")),HR=3,file=stdout())
				if(length(vars)!=0){
					ctrl <- trainControl(method = "repeatedcv", number = k, savePredictions = TRUE)
					newdat <- dataset[,c(dep_var,indep_numeric_var,indep_cat_var),drop=F]
					newdat <- newdat[complete.cases(newdat),,drop=F]
					command_str	<- paste0("mod_fit <- train(res_LGM$formula,  data=newdat, method='glm',family=binomial(",link,"),trControl = ctrl, tuneLength = 5)")
					eval(parse(text=command_str)) ;
					A8 <- mod_fit$resample
					A8 <- A8[,-3]
					A8 <- rbind(A8,colMeans(A8))
					rownames(A8) <- c(paste('Fold',1:k),'Mean')
					R2HTML::HTML(A8,file=stdout(),innerBorder = 1,align="left",digits=4)
				} else {
					warn.msg10 <- "\a Warning : K-fold cross validation is not supported for intercept only model."
					R2HTML::HTML(warn.msg10,file=stdout())
				}
			}

			# Variable Selection
			if(select){
				R2HTML::HTML(R2HTML::as.title("Variable Selction"),HR=3,file=stdout())
				if(length(vars)!=0){
					VS <- capture.output(res_LGM2$anova)
					blank <- which(VS=="")
					IM <- paste(VS[(blank[1]+2):(blank[2]-1)],collapse="")
					IM <- gsub("[ ]{2,}"," ",IM)
					FM <- paste(VS[(blank[2]+2):(blank[3]-1)],collapse="")
					FM <- gsub("[ ]{2,}"," ",FM)

					if(noint){
						IM <- paste0(gsub(" - 1","",IM)," (Intercept is not included)")	# Initial Model
						FM <- paste0(gsub(" - 1","",FM)," (Intercept is not included)")	# Final Model
					}
					if(is.null(keep_var)){
						resVS <- matrix(c('Method',direction,'Initial Model',IM,'Final Model',FM),byrow=T,ncol=2)
					} else {
						resVS <- matrix(c('Method',direction,'Fixed Variable',paste(keep_var,collapse=', '),'Initial Model',IM,'Final Model',FM),byrow=T,ncol=2)
					}
					R2HTML::HTML(resVS,file=stdout(), innerBorder = 1,align="left",digits=4)
				} else {
					warn.msg11 <- "\a Warning : Variable selection is not supported for intercept only model."
					R2HTML::HTML(warn.msg11,file=stdout())
				}
			}
			### plot
			if(Plot) {
				R2HTML::HTML(R2HTML::as.title(paste("Generalized Linear Regression plot",", Link function=",link,sep="")),HR=2,file=stdout(),append=TRUE)
				REx_ANA_PLOT(800,800)
				par(mfrow=c(2,2))
				plot(res_LGM)
				REx_ANA_PLOT_OFF("")
			}
		}
		# Analysis end time
		R2HTML::HTMLhr(file=stdout())
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Generalized Linear Regression.",sep=""),file=stdout())
		R2HTML::HTMLhr(file=stdout())
	})

	if(exists('O')){
		return(list(html=html.output,Output=O))
	} else {
		return(html.output)
	}
}

REx_WLM <- function(dataset,dep_var,indep_numeric_var=NULL,indep_cat_var=NULL,vars=NULL,weight_var=NULL,weight_option=FALSE,noint=FALSE,CI=TRUE,confint.level=0.95,ANOVA=TRUE,Predict=FALSE,Resid=FALSE,stdResid=FALSE,weightVar=FALSE,Plot=TRUE){
	#### 변수 설명 ####
	## dep_var : 종속변수(dependent variables),필수변수(무조건 한개변수만 선택할수있음)
	## indep_numeric_var : 질적설명변수(independent variables) (2개 이상의 변수가 있는 경우 c(,,)로 구분, 필수사항아님,단,상수항옵션이(noint=TRUE)인경우 돌아갈수없음)
	## indep_cat_var : 양적설명변수(independent variables) (2개 이상의 변수가 있는 경우 c(,,)로 구분, 필수사항아님,단,상수항옵션이(noint=TRUE)인경우 돌아갈수없음)
	## vars : 선택된 변수들 (2개 이상의 변수가 있는 경우 c(,,)로 구분, 필수사항아님)
	## weight_var : 등분산성을 만족하지 않을때, 사용하는 방법으로 가중변수를 만들기위해 다음의 식을 적용하기 위해서( 1/(weight_var)^2 )계산되기전 weight값의 변수를 "weight_var"변수로 설정한다.
	## wts_option : 모형추정방법 선택 (TRUE/FALSE, FALSE:default)
	## noint : 모델설정 탭의 '상수항 포함하지 않음' 옵션(TRUE/FALSE,FALSE:default)
	## CI : 출력옵션 탭의 신뢰구간(TRUE/FALSE,TRUE:default)
	## confint.level : 출력옵션 탭의 '신뢰수준' 옵션 (신뢰수준범위 변경가능,0.95 default)
	## ANOVA : 출력옵션 탭의 'ANOVA Table' 옵션(TRUE/FALSE,TRUE:default)
	## Predict : 출력옵션 탭의 '예측값' 옵션(TRUE/FALSE,FALSE:default)
	## Resid : 출력옵션 탭의 '잔차' 옵션(TRUE/FALSE,FALSE:default)
	## stdResid: 출력옵션 탭의 '표준화잔차' 옵션(TRUE/FALSE,FALSE:default)
	## weightVar: 출력옵션 탭의 '(최종)가중변수' 옵션(TRUE/FALSE,FALSE:default)
	## Plot : 출력옵션 탭의 '도표' 옵션(TRUE/FALSE,TRUE:default)
	###################
	#### Required packages : R2HTML , rms, MASS,markdown ####
	###################
	load.pkg(c("R2HTML", "rms", "MASS"))
		
	html.output <- capture.output({

		# Title
		R2HTML::HTML(R2HTML::as.title("Weighted Linear Regression"),HR=1,file=stdout(),append=FALSE)

		## Warnings
		# Response variable type
		if(!is.numeric(dataset[,dep_var])) warn.msg1 <- '\a Error : Response variable should be numeric. Analysis has been stopped.'
		## weight_variable
		if(is.null(weight_var)) warn.msg5 <- '\a Error : Weight variable should be selected. Analysis has been stopped.'

		# explanatory variable type
		if(!is.null(indep_numeric_var)) {
			is.nom <- sapply(indep_numeric_var,function(i) !is.numeric(dataset[,i]))
			if(any(is.nom)){
				warn.msg2 <- paste0("\a Warning : The type of variable '",paste(indep_numeric_var[is.nom],collapse=', '),"' is not numeric but selected as the quantitative variable. It was excluded from the analysis.")
				ec <- indep_numeric_var[is.nom]
				vars <- vars[-grep(paste0('[',paste(ec,collapse=','),']'),vars)]
				if(length(vars)==0) vars <- NULL

			}
		}
		if(!is.null(indep_cat_var)) {
			is.num <- sapply(indep_cat_var,function(i) is.numeric(dataset[,i]))
			if(any(is.num)) warn.msg3 <- paste0("\a Warning : The type of variable '",paste(indep_cat_var[is.num],collapse=', '),"' is numeric but selected as the qualitative variable. It was coreced into character.")
			for(i in indep_cat_var[is.num]) dataset[,i] <- as.factor(dataset[,i])
		}

		# no intercept & no explanatory variable : error
		if(noint&is.null(vars)) warn.msg4 <- '\a Error : With no intercept, at least 1 independent variable should be selected. Analysis has been stopped.'

		# warnings
		if(exists('warn.msg1')|exists('warn.msg4')|exists('warn.msg5')){
			R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
			if(exists('warn.msg1')) R2HTML::HTML(warn.msg1,file=stdout())
			if(exists('warn.msg4')){
				if(exists('warn.msg2')) R2HTML::HTML(warn.msg2,file=stdout())
				R2HTML::HTML(warn.msg4,file=stdout())
			}
			if(exists('warn.msg5')) R2HTML::HTML(warn.msg5,file=stdout())
			
		} else {
			## model formula
			var_info <- c(dep_var,indep_numeric_var, indep_cat_var)
			is.na <- sapply(var_info,function(i) is.na(dataset[,i]))
			oo <- which(rowSums(is.na)>0)
			raw.dat <- dataset
			if(length(oo)>0) dataset <- dataset[-oo,]

			form.0 <- ifelse(is.null(vars),paste0(dep_var,' ~ 1'),paste0(dep_var,' ~ ',paste(vars,collapse=' + ')))
			if(noint) form.0 <- paste(form.0,-1)
			form.1 <- as.formula(form.0)

			## Fitting WLS
			wts <- dataset[,weight_var]
			wts1 <- wts2 <- wts3 <- wts4 <- rep(NA,nrow(dataset))
			# 가중변수가 연속형변수일때
			if(is.numeric(wts)){	
				if(!weight_option){
					#weight_option=FALSE
					wts1 <- 1/wts^2
					res_LM_1 <- lm(formula(form.1),data=dataset,weight=wts1)
					WTS <- wts1
				}else{
					#weight_option=TRUE
					wts1 <- 1/wts^2
					model.1 <- lm(formula(form.1),data=dataset,weight=wts1)
					wts2 <- 1/fitted(lm(abs(residuals(model.1)) ~ fitted(model.1)))^2
					res_LM_1 <- lm(formula(form.1),data=dataset,weight=wts2)
					WTS <- wts2
				}
			}else{
			# 가중변수가 범주형변수일때
					if(!weight_option){
					#weight_option=FALSE
					model.1 <- lm(formula(form.1),data=dataset)
					VV <- tapply(residuals(model.1), wts, var)
					ll <- levels(wts)
					for(i in seq_along(ll)){
						wts3[wts==ll[i]] <- 1/VV[i]
					}
					res_LM_1 <- lm(formula(form.1),data=dataset,weight=wts3)
					WTS <- wts3
				}else{
					#weight_option=TRUE
					model.1 <- lm(formula(form.1),data=dataset)
					VV <- tapply(residuals(model.1), wts, var)
					ll <- levels(wts)
					for(i in seq_along(ll)){
						wts3[wts==ll[i]] <- 1/VV[i]
					}
					model.2 <- lm(formula(form.1),data=dataset,weight=wts3)
					wts4 <- 1/fitted(lm(abs(residuals(model.2)) ~ fitted(model.2)))^2
					res_LM_1 <- lm(formula(form.1),data=dataset,weight=wts4)
					WTS <- wts4
				}
			}


			### final model:res_LM_1
			### weight_var = WTS

			### Options
			# 다음의 값들은 엑셀 시트에 저장
if(Predict) {
				temp.0 <- data.frame(Fitted_LM=fitted(res_LM_1))
				temp <- data.frame(Fitted_LM=rep(NA,nrow(raw.dat)))
				row.dataset <- rownames(raw.dat)
				row.output <- rownames(temp.0)
				for(i in seq_along(row.output)) if(i%in%seq_along(row.output)) temp[row.dataset==row.output[i],1] <- temp.0[row.output[i],1]
				O <- temp
			}

			if(Resid) {
				temp.0 <- data.frame(Resid_LM=resid(res_LM_1))
				temp <- data.frame(Resid_LM=rep(NA,nrow(raw.dat)))
				row.dataset <- rownames(raw.dat)
				row.output <- rownames(temp.0)
				for(i in seq_along(row.output)) if(i%in%seq_along(row.output)) temp[row.dataset==row.output[i],1] <- temp.0[row.output[i],1]
				if(exists('O')){
					O <- cbind(O,temp)
				} else {
					O <- temp
				}
			}

			if(stdResid) {
				temp.0 <- data.frame(stdResid_LM=studres(res_LM_1))
				is.singular <- class(try(temp.0[,1],silent=T))=="try-error"
				if(is.singular){
					warn.msg6 <- "\a Error : Cannnot calculate standardized residuals because of the (near) linear dependences in explatory variables. Please check multicollinearity problem."
					R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
					R2HTML::HTML(warn.msg6,file=stdout())
				} else {

					temp <- data.frame(stdResid_LM=rep(NA,nrow(raw.dat)))
					row.dataset <- rownames(raw.dat)
					row.output <- rownames(temp.0)
					for(i in seq_along(row.output)) if(i%in%seq_along(row.output)) temp[row.dataset==row.output[i],1] <- temp.0[row.output[i],1]
				
					if(exists('O')){
						O <- cbind(O,temp)
					} else {
						O <- temp
					}
				}
			}

			if(weightVar) {
				temp.0 <- data.frame(weightVar=WTS)
				temp <- data.frame(weightVar=rep(NA,nrow(raw.dat)))
				row.dataset <- rownames(raw.dat)
				row.output <- rownames(temp.0)
				for(i in seq_along(row.output)) if(i%in%seq_along(row.output)) temp[row.dataset==row.output[i],1] <- temp.0[row.output[i],1]
				if(exists('O')){
					O <- cbind(O,temp)
				} else {
					O <- temp
				}
			}

			# Data Structure
			R2HTML::HTML(R2HTML::as.title("Data Structure"),HR=2,file=stdout())
			total.var <- ncol(dataset)
			used.var <- ifelse(is.null(vars),0,length(unique(unlist(strsplit(vars,":")))))
			weight.var <- weight_var
			DS <- matrix(c('Number of observations','Number of total variables','Number of used variables','Number of weight variables',nrow(dataset),total.var,used.var+1,1),ncol=2)
			R2HTML::HTML(DS,file=stdout(), innerBorder = 1,align="left")

			# Varibale list
			R2HTML::HTML(R2HTML::as.title("Variable List"),HR=2,file=stdout())
			if(is.null(vars)) {
				Vars <- c()
			} else {
				Vars <-  unique(unlist(strsplit(vars,":")))
			}

			qual <- indep_cat_var[indep_cat_var%in%Vars]
			quan <- c(indep_numeric_var[indep_numeric_var%in%Vars],dep_var)
			varlist <- matrix(c('Quantitative variable',paste(quan,collapse=', ')),ncol=2,byrow=T)
			if(length(qual)!=0) varlist <- rbind(varlist,c('Qualitative variable',paste(qual,collapse=', ')))
			varlist <- rbind(varlist,c('Weight variable',weight_var))
			R2HTML::HTML(varlist,file=stdout(), innerBorder = 1,align="left")
			if(exists('warn.msg2')) R2HTML::HTML(warn.msg2,file=stdout())
			if(exists('warn.msg3')) R2HTML::HTML(warn.msg3,file=stdout())


			# Analysis Description
			#R2HTML::HTML(R2HTML::as.title("Analysis Description"),HR=2,file=stdout())
			R2HTML::HTML(R2HTML::as.title("Analysis Description"),HR=2,file=stdout())
			AD <- matrix(c('Response variable',dep_var),ncol=2,byrow=T)
			if(!is.null(vars)) AD <- rbind(AD,c('Explanatory variable',paste(vars,collapse=', ')))
			AD <- rbind(AD,matrix(c('Intercept included',!noint),ncol=2,byrow=T))
			if(CI) AD <- rbind(AD,c('Significance level of CI',confint.level))
			AD <- rbind(AD,matrix(c('Print ANOVA table',ANOVA,'Adjusted weight variable',weight_option),ncol=2,byrow=T))
			R2HTML::HTML(AD,file=stdout(), innerBorder = 1,align="left")



			# Coefficient Estimate
			R2HTML::HTML(R2HTML::as.title("Weighted Linear Regression Results"),HR=2,file=stdout())
			R2HTML::HTML(R2HTML::as.title("Coefficients"),HR=3,file=stdout())
			CE <- summary(res_LM_1)$coef
			if(CI){
				CE <- cbind(CE,confint(res_LM_1, level=confint.level))
				colnames(CE)[(ncol(CE)-1):ncol(CE)] <- paste0(c('Lower','Upper'),' bound of ',confint.level*100,'% CI')
			}
			R2HTML::HTML(round(CE,4),file=stdout(), innerBorder = 1,align="left")

			# Anova table and R-squared
			if(ANOVA){
				R2HTML::HTML(R2HTML::as.title("ANOVA Table"),HR=3,file=stdout())
				Anova <- as.matrix(anova(lm(formula(paste(dep_var,'~1',sep='')),x=TRUE,data=dataset),res_LM_1))
				A1 <- round(rbind(Anova[2,4:3], Anova[2:1,c(2,1)]),5)
				MS <- c(round(A1[1:2,1]/A1[1:2,2],5)," ")
				A2 <- data.frame(round(A1,5),MS=MS,F=c(round(Anova[2,5],5)," "," "),Pvalue=c(signif(Anova[2,6],5)," "," "),R2=c(round(summary(res_LM_1)$r.squared,5),"","") ,adj.R2=c(round(summary(res_LM_1)$adj.r.squared,5),"",""))
				colnames(A2)[1] <- 'SS'
				rownames(A2) <- c('SSR','SSE','SST')
				R2HTML::HTML(as.matrix(A2),file=stdout(), innerBorder = 1,align="left")
			}
			### plot
			if(Plot) {
				R2HTML::HTML(R2HTML::as.title("Weighted Regression plot"),HR=2,file=stdout(),append=TRUE)
				REx_ANA_PLOT(800,800)
				par(mfrow=c(2,2))
				plot(res_LM_1)
				REx_ANA_PLOT_OFF("")
			}

		}
		# Analysis end time
		R2HTML::HTMLhr(file=stdout())
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Weighted Linear Regression.",sep=""),file=stdout())
		R2HTML::HTMLhr(file=stdout())
	})

	if(any(Predict,Resid,stdResid,weightVar)){
		return(list(html=html.output,Output=O))
	} else {
		return(html.output)
	}
}

REx_LM <- function(dataset,dep_var,indep_numeric_var=NULL,indep_cat_var=NULL,vars=NULL,noint=FALSE,CI=TRUE,confint.level=0.95,VIF=FALSE,select=FALSE,keep_var=NULL,direct="forward",ANOVA=TRUE,Predict=FALSE,Resid=FALSE,stdResid=FALSE,Plot=TRUE){
	#### 변수 설명 ####
	## dep_var : 종속변수(dependent variables),필수변수(무조건 한개변수만 선택할수있음)
	## indep_numeric_var : 질적설명변수(independent variables) (2개 이상의 변수가 있는 경우 c(,,)로 구분, 필수사항아님,단,상수항옵션이(noint=TRUE)인경우 돌아갈수없음)
	## indep_cat_var : 양적설명변수(independent variables) (2개 이상의 변수가 있는 경우 c(,,)로 구분, 필수사항아님,단,상수항옵션이(noint=TRUE)인경우 돌아갈수없음)
	## vars : 선택된 변수들 (2개 이상의 변수가 있는 경우 c(,,)로 구분, 필수사항아님)
	## noint : 모델설정 탭의 '상수항 포함하지 않음' 옵션 (FALSE:default)
	## CI : 출력옵션 탭의 신뢰구간 (TRUE/FALSE, TRUE:default)
	## confint.level : 출력옵션 탭의 '신뢰수준' 옵션 (신뢰수준범위 변경가능,0.95 default)
	## VIF : 출력옵션 탭의 '다중공선성' 옵션 (TRUE/FALSE,FALSE:default)
	## ANOVA : 출력옵션 탭의 'ANOVA Table' 옵션 (TRUE/FALSE, TRUE:default)
	## Predict : 출력옵션 탭의 '예측값' 옵션 (TRUE/FALSE,FALSE:default)
	## Resid : 출력옵션 탭의 '잔차' 옵션 TRUE/FALSE,FALSE:default)
	## stdResid: 출력옵션 탭의 '표준화잔차' 옵션 TRUE/FALSE,FALSE:default)
	## select : 변수선택 탭의 '변수선택' 옵션(TRUE/FALSE,FALSE:default))
	## keep_var : 질적변수 or 양적변수 상관없이 모든 설명변수가 대상임, 만약에 설명변수에 변수가 설정이 안되있으면 활성화되어도 설명변수창에 아무것도 뜨지 않게해야됨
	## 변수선택 탭의 '고정변수선택' 옵션 (설정없으면(select=FALSE,default) NULL or 설정(select=TRUE),2개 이상의 변수가 있는 경우 c(,,)로 구분)
	## direct : 통계설정 탭의 '변수선택의 세부항목' 옵션("forward","backward","both",가능)
	## Plot : 회귀분석 그림출력 옵션 (TRUE/FALSE, TRUE:default)
	###################
	#### Required packages : R2HTML , rms, MASS, markdown ####
	###################
	load.pkg(c("R2HTML", "rms", "MASS"))
	
	html.output <- capture.output({

		# Title
		#R2HTML::HTML(R2HTML::as.title("Linear Regression"),HR=1,file="./test.html",append=FALSE)
		R2HTML::HTML(R2HTML::as.title("Linear Regression"),HR=1,file=stdout(),append=FALSE)

		## Warnings
		# Response variable type
		if(!is.numeric(dataset[,dep_var])) warn.msg1 <- '\a Error : Response variable should be numeric. Analysis has been stopped.'
		# explanatory variable type
		if(!is.null(indep_numeric_var)) {
			is.nom <- sapply(indep_numeric_var,function(i) !is.numeric(dataset[,i]))
			if(any(is.nom)){
				warn.msg2 <- paste0("\a Warning : The type of variable '",paste(indep_numeric_var[is.nom],collapse=', '),"' is not numeric but selected as the quantitative variable. It was excluded from the analysis.")
				ec <- indep_numeric_var[is.nom]
				vars <- vars[-grep(paste0('[',paste(ec,collapse=','),']'),vars)]
				if(length(vars)==0) vars <- NULL

			}
		}
		if(!is.null(indep_cat_var)) {
			is.num <- sapply(indep_cat_var,function(i) is.numeric(dataset[,i]))
			if(any(is.num)) warn.msg3 <- paste0("\a Warning : The type of variable '",paste(indep_cat_var[is.num],collapse=', '),"' is numeric but selected as the qualitative variable. It was coreced into character.")
			for(i in indep_cat_var[is.num]) dataset[,i] <- as.factor(dataset[,i])
		}

		# no intercept & no explanatory variable : error
		if(noint&is.null(vars)) warn.msg4 <- '\a Error : With no intercept, at least 1 independent variable should be selected. Analysis has been stopped.'

		# warnings
		if(exists('warn.msg1')|exists('warn.msg4')){
			R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
			if(exists('warn.msg1')) R2HTML::HTML(warn.msg1,file=stdout())
			if(exists('warn.msg4')){
				if(exists('warn.msg2')) R2HTML::HTML(warn.msg2,file=stdout())
				R2HTML::HTML(warn.msg4,file=stdout())
			}
		} else {
			## model formula
			var_info <- c(dep_var,indep_numeric_var, indep_cat_var)
			is.na <- sapply(var_info,function(i) is.na(dataset[,i]))
			oo <- which(rowSums(is.na)>0)
			raw.dat <- dataset
			if(length(oo)>0) dataset <- dataset[-oo,]

			form.0 <- ifelse(is.null(vars),paste0(dep_var,' ~ 1'),paste0(dep_var,' ~ ',paste(vars,collapse=' + ')))
			if(noint) form.0 <- paste(form.0,-1)
			form.1 <- as.formula(form.0)

			## Fitting WLS
			
			res_LM_1 <- lm(formula(form.1),data=dataset)

			## Variable selection
			if(select){
				## 고정변수 선택여부에따른 분석방법(선택안한경우)
				f = file()
				sink(file=f)
				if(is.null(keep_var)){
					res_LM2 <- stepAIC(res_LM_1,direction=direct)
					form2 <- formula(res_LM2)
					res_LM <- lm(form2,data=dataset)
				## 고정변수 선택여부에따른 분석방법(선택한경우)	
				}else{
					form1 <- formula(paste('~',paste(keep_var,collapse='+'),sep=''))
					res_LM2 <- stepAIC(res_LM_1,direction=direct,scope=list(lower=form1))
					form2 <- formula(res_LM2)
					res_LM <- lm(form2,data=dataset)
				}
				sink()
				close(f)
			}
			
			### final model
			if(select) res_LM_1 <- res_LM

			### Options
			# 다음의 값들은 엑셀 시트에 저장
			if(Predict) {
				temp.0 <- data.frame(Fitted_LM=fitted(res_LM_1))
				temp <- data.frame(Fitted_LM=rep(NA,nrow(raw.dat)))
				row.dataset <- rownames(raw.dat)
				row.output <- rownames(temp.0)
				for(i in seq_along(row.output)) if(i%in%seq_along(row.output)) temp[row.dataset==row.output[i],1] <- temp.0[row.output[i],1]
				O <- temp
			}

			if(Resid) {
				temp.0 <- data.frame(Resid_LM=resid(res_LM_1))
				temp <- data.frame(Resid_LM=rep(NA,nrow(raw.dat)))
				row.dataset <- rownames(raw.dat)
				row.output <- rownames(temp.0)
				for(i in seq_along(row.output)) if(i%in%seq_along(row.output)) temp[row.dataset==row.output[i],1] <- temp.0[row.output[i],1]
				if(exists('O')){
					O <- cbind(O,temp)
				} else {
					O <- temp
				}
			}

			if(stdResid) {
				temp.0 <- data.frame(stdResid_LM=studres(res_LM_1))
				is.singular <- class(try(temp.0[,1],silent=T))=="try-error"
				if(is.singular){
					warn.msg6 <- "\a Error : Cannnot calculate standardized residuals because of the (near) linear dependences in explatory variables. Please check multicollinearity problem."
					R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
					R2HTML::HTML(warn.msg6,file=stdout())
				} else {

					temp <- data.frame(stdResid_LM=rep(NA,nrow(raw.dat)))
					row.dataset <- rownames(raw.dat)
					row.output <- rownames(temp.0)
					for(i in seq_along(row.output)) if(i%in%seq_along(row.output)) temp[row.dataset==row.output[i],1] <- temp.0[row.output[i],1]
				
					if(exists('O')){
						O <- cbind(O,temp)
					} else {
						O <- temp
					}
				}
			}

			# Data Structure
			R2HTML::HTML(R2HTML::as.title("Data Structure"),HR=2,file=stdout())
			total.var <- ncol(dataset)
			used.var <- ifelse(is.null(vars),0,length(unique(unlist(strsplit(vars,":")))))
			DS <- matrix(c('Number of observations','Number of total variables','Number of used variables',nrow(dataset),total.var,used.var+1),ncol=2)
			R2HTML::HTML(DS,file=stdout(), innerBorder = 1,align="left")

			# Varibale list
			R2HTML::HTML(R2HTML::as.title("Variable List"),HR=2,file=stdout())
			if(is.null(vars)) {
				Vars <- c()
			} else {
				Vars <-  unique(unlist(strsplit(vars,":")))
			}

			qual <- indep_cat_var[indep_cat_var%in%Vars]
			quan <- c(indep_numeric_var[indep_numeric_var%in%Vars],dep_var)
			varlist <- matrix(c('Quantitative variable',paste(quan,collapse=', ')),ncol=2,byrow=T)
			if(length(qual)!=0) varlist <- rbind(varlist,c('Qualitative variable',paste(qual,collapse=', ')))
			R2HTML::HTML(varlist,file=stdout(), innerBorder = 1,align="left")
			if(exists('warn.msg2')) R2HTML::HTML(warn.msg2,file=stdout())
			if(exists('warn.msg3')) R2HTML::HTML(warn.msg3,file=stdout())


			# Analysis Description
			R2HTML::HTML(R2HTML::as.title("Analysis Description"),HR=2,file=stdout())
			#R2HTML::HTML(R2HTML::as.title("Analysis Description"),HR=2,file="./test.html")
			AD <- matrix(c('Response variable',dep_var),ncol=2,byrow=T)
			if(!is.null(vars)) AD <- rbind(AD,c('Explanatory variable',paste(vars,collapse=', ')))
			AD <- rbind(AD,matrix(c('Intercept included',!noint),ncol=2,byrow=T))

			if(CI) AD <- rbind(AD,c('Significance level of CI',confint.level))
			AD <- rbind(AD,matrix(c('Print VIF',VIF,'Print ANOVA table',ANOVA,'Variable selection',select),ncol=2,byrow=T))
			#variable selection
			if(select) {
				direction <- switch(direct,forward="Forward selection",backward="Backward elimination",both="Stepwise regression")
				AD <- rbind(AD,matrix(c('Variable selection method',direction),ncol=2,byrow=T))
			}
			if(!is.null(keep_var)) AD <- rbind(AD,c('Keep variable in variable selection',paste(keep_var,collapse=', ')))
			R2HTML::HTML(AD,file=stdout(), innerBorder = 1,align="left")



			# Coefficient Estimate
			R2HTML::HTML(R2HTML::as.title("Linear Regression Results"),HR=2,file=stdout())
			R2HTML::HTML(R2HTML::as.title("Coefficients"),HR=3,file=stdout())
			CE <- summary(res_LM_1)$coef
			if(CI){
				CE <- cbind(CE,confint(res_LM_1, level=confint.level))
				colnames(CE)[(ncol(CE)-1):ncol(CE)] <- paste0(c('Lower','Upper'),' bound of ',confint.level*100,'% CI')
			}
			if(VIF)	{
				VIF.1 <- vif(res_LM_1)
				if(!noint) VIF.1 <- c(NA,VIF.1)
				CE <- cbind(CE,VIF=VIF.1)
			}
			R2HTML::HTML(round(CE,4),file=stdout(), innerBorder = 1,align="left")

			# Anova table and R-squared
			if(ANOVA){
				R2HTML::HTML(R2HTML::as.title("ANOVA Table"),HR=3,file=stdout())
				Anova <- as.matrix(anova(lm(formula(paste(dep_var,'~1',sep='')),x=TRUE,data=dataset),res_LM_1))
				A1 <- round(rbind(Anova[2,4:3], Anova[2:1,c(2,1)]),5)
				MS <- c(round(A1[1:2,1]/A1[1:2,2],5)," ")
				A2 <- data.frame(round(A1,5),MS=MS,F=c(round(Anova[2,5],5)," "," "),Pvalue=c(signif(Anova[2,6],5)," "," "),R2=c(round(summary(res_LM_1)$r.squared,5),"","") ,adj.R2=c(round(summary(res_LM_1)$adj.r.squared,5),"",""))
				colnames(A2)[1] <- 'SS'
				rownames(A2) <- c('SSR','SSE','SST')
				R2HTML::HTML(as.matrix(A2),file=stdout(), innerBorder = 1,align="left",digits=4)
			}
			### plot
			if(Plot) {
				R2HTML::HTML(R2HTML::as.title("Regression plot"),HR=2,file=stdout(),append=TRUE)
				REx_ANA_PLOT(800,800)
				par(mfrow=c(2,2))
				plot(res_LM_1)
				REx_ANA_PLOT_OFF("")
			}

		}
		# Analysis end time
		R2HTML::HTMLhr(file=stdout())
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Linear Regression.",sep=""),file=stdout())
		R2HTML::HTMLhr(file=stdout())
	})

	if(any(Predict,Resid,stdResid)){
		return(list(html=html.output,Output=O))
	} else {
		return(html.output)
	}
}

REx_loglinear <- function(data, freq.form = T, freq.var=NULL, variables,
                          table = T, cond.tab = T, margin.tab = F,
                          fitness = T, estim = TRUE, exp_estim = FALSE, CI = TRUE, CI_mode = c("likelihood", "wald", "both"), CI_level = 0.95, 
                          testing = T,chisq_mode = c("likelihood", "wald", "both"),  ss_type = c("I","II","III")) {
  ###################
  #### Variables
  ###################
  
  ### [data tab]
  ## data                : data frame object, input data
  ## freq.form           : a logical indicating whether input data has a frequecy variable
  ##                       (입력 데이터에 빈도 변수가 있는 경우, 즉, UI에서 빈도변수가 1개 선택된 경우)
  ## freq.var            : frequency variable
  ##                       (UI의 모델설정 탭에서 선택된 빈도 변수)
  ## variables           : defined variables
  ##                       (UI의 모델설정 탭에서 선택된 설명 변수)
  ## table               : contigency table
  ##                       (UI의 출력옵션 탭에서 '빈도표'와 대응)
  ## fitness             : model fitness measures such as deviance
  ##                       (UI의 출력옵션 탭에서 '적합성 측도'와 대응)
  ## estim, exp_estim    : parameter estimates & exponential mode
  ##                       (UI의 출력옵션 탭에서 '모수 추정값' 및 '지수 모수 추정값'에 대응)
  ## chisq_mode, CI_mode : Chi-squared statistics & confidence intervals mode
  ##                       (UI의 출력옵션 탭에서 '카이제곱 통계량', '신뢰구간'에 대응)
  ## CI                  : confidence interval
  ##                       (UI의 출력옵션 탭에서 '신뢰구간'이 선택되면 TRUE)
  ## CI_level            : confidence level of the confidence interval
  ##                       (UI의 출력옵션 탭에서 '신뢰구간'의 입력값에 대응)
  ## testing             : testing whether each parameter is equal to 0 
  ##                       (UI의 출력옵션 탭에서 '가설 검정'에 대응)
  ## ss_type             : type of chi-squared statistics
  ##                       (UI의 출력옵션 탭에서 '통계량 유형'과 대응, type1,2,3에 따라 "I","II","III"로 입력)
  
  
  ###################
  #### Required packages
  ###################
  
  ### Required packages : R2HTML, car, MASS, AICcmodavg
  ## - R2HTML     : to output in the form of HTML
  ## - car        : to make anova table using 'Anova' function
  ## - MASS       : to perform log linear analysis using 'loglm' function 
  ## - AICcmodavg : to estimate dispersion for Poisson GLM
  load.pkg(c("R2HTML", "car", "MASS"))
  
  #####################################################################################
  ################### Function(myTryCatch) : Catch and save warning and error message
  #####################################################################################
  
  myTryCatch <- function(expr) {
    warn <- err <- NULL
    value <- withCallingHandlers(
      tryCatch(expr, error = function(e) {
        err <<- e
        NULL
      }), warning = function(w) {
        warn <<- w
        invokeRestart("muffleWarning")
      })
    
    return(list(warning = warn, error = err))
  }
  
  
  #####################################################################################
  ################### Function (combVarLevel, tableName) : Format 2*2 contingency table
  #####################################################################################
  
  #### Combine Variable name and levels
  combVarLevel <- function(table, i) { # if row, i = 1, if column i = 2
    dimname <- dimnames(table)
    
    var.name <- names(dimname)[i]
    level.name <- dimname[[i]]
    var.lev.name <- paste0(var.name, "(", level.name, ")")   
    
    return(var.lev.name)
  }
  
  #### Assign names to a table
  tableName <- function(table.dim2) {
    x.table.dim2 <- as.data.frame.matrix(table.dim2)
    rownames(x.table.dim2) <- combVarLevel(table.dim2, 1)
    colnames(x.table.dim2) <- combVarLevel(table.dim2, 2)
    
    return(x.table.dim2)
  }
  
  #### Calculate row-wise sum and column-wise sum
  new.rowSums <- function(mat) apply(mat,1,function(x) sum(as.numeric(x)))
  new.colSums <- function(mat) apply(mat,2,function(x) sum(as.numeric(x)))
  
  #### Add total sum to a table
  sumTl <- function(x.table.dim2) {
    tab <- tableName(x.table.dim2)
    tab <- cbind(tab, Total = new.rowSums(tab))
    tab <- rbind(tab, Total = new.colSums(tab))
    
    return(tab)
  }
  
  
  
  #####################################################################################
  ################### Function (cbind.fill) : column bind multiple data frames when some of them are empty
  #####################################################################################
  
  #### cbind with empty data frame
  cbind.fill <- function(...){
    nm <- list(...) 
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
      rbind(x, matrix(, n-nrow(x), ncol(x))))) 
  }
  
  
  
  #####################################################################################
  ################### Analysis and Output
  #####################################################################################
  
  html.output <- capture.output({
    
    ss_type	<- as.character(ss_type)
    ss_type	<- match.arg(ss_type)
    chisq_mode	<- match.arg(chisq_mode)
    CI_mode	<- match.arg(CI_mode)
    num_digit <- 3
    
    
    #### HTML file title
    R2HTML::HTML(R2HTML::as.title("Log Linear Analysis"), HR = 1, file = stdout(), append = F)
    
    
    #### Create table
    if (!is.null(freq.var)) {
      model <- paste(freq.var, "~", paste(variables, collapse = " + "))
    } else model <- paste("Frequency ~", paste(variables, collapse = "+"))
    model.x <- strsplit(model, "~")[[1]][2]
    covariates <- unique(unlist(strsplit(variables,":")))

    if (any(!(covariates %in% colnames(data))) | is.null(covariates)) {
      R2HTML::HTML("Unable to find covariates in the data.", file=stdout())
      stop("Unable to find covariates in the data.")
    }
    
    if (!is.null(freq.var)) {
      x.table <- xtabs(formula = paste(freq.var, "~", paste(covariates, collapse = " + ")), data = data)
    } else if (is.null(freq.var)) {
      x.table <- xtabs(paste("~", paste(covariates, collapse = " + ")), data = data)
    }
    
    n <- sum(x.table)
    p <- length(dimnames(x.table))
    
    
    #### Log linear analysis
    loglm.1 <- loglm(formula(model), data = x.table)
    loglm.2 <- glm(formula = paste("Freq ~", strsplit(model, "~")[[1]][2]), 
                   data = x.table, family = poisson)
    
    
    #### Basic Information
    
    ## Data structure
    data_str	<- matrix(c("Number of observations", "Number of variables", n, p), ncol = 2)
    if(!freq.form) rbind(data_str, c("Variable list", paste(covariates, collapse=", ")))
    R2HTML::HTML(R2HTML::as.title("Data Structure"), HR = 2, file = stdout())
    R2HTML::HTML(data_str, file=stdout(), row.names = FALSE, col.names = TRUE, innerBorder = 1, align = "left") ;
    
    ## Variable list
    if(freq.form){
      R2HTML::HTML(R2HTML::as.title("Variable List"), HR = 2, file = stdout())
      VL <- matrix(c('Explanatory variables',paste(covariates, collapse = ", "),'Frequency variable', freq.var), ncol=2, byrow=T)
      R2HTML::HTML(VL, file = stdout(), row.names = FALSE, col.names = TRUE, innerBorder = 1, align = "left")
    }
    
    ## Check variable
    check.var <- sapply(data[,covariates,drop=F],is.numeric)
    numvars <- covariates[check.var]
    war.msg <- ifelse(length(numvars)==0,NA, paste0("Warning message : The type of variable '",paste(numvars,collapse=', '),"' was numeric, but was coerced into character. <br>(Log linear analysis is supported only for categorical variables.)"))
    if (!is.na(war.msg)) R2HTML::HTML(war.msg, file = stdout(), align = "left")
    
    ## Analysis description
    R2HTML::HTML(R2HTML::as.title("Analysis Description"), HR = 2, file = stdout())
    TCI <- switch(chisq_mode,wald = "wald", likelihood = "likelihood", both = "wald, likelihood")
    AD <- matrix(c('Model', model), ncol = 2)
    if(CI) {
      SCI <- switch(CI_mode, wald="wald", likelihood = "likelihood", both = "wald, likelihood")
      AD <- rbind(AD, matrix(c('Statistics for CI', SCI, 'Significance level for CI', CI_level), ncol = 2, byrow = T))
    }
    if(testing){
      TCI <- switch(chisq_mode, wald="wald", likelihood="likelihood", both="wald, likelihood")    
      AD <- rbind(AD, matrix(c('Chi-square statistics for testing', TCI, 'Type of statistics', paste0('type ', ss_type)), ncol = 2, byrow = T))
    }
    R2HTML::HTML(AD, file = stdout(), row.names = FALSE, col.names = TRUE, innerBorder = 1, align = "left")
    
    
    
    #### Contingency Table
    
    ## when the dimension is 2
    if (table) {
      R2HTML::HTML(R2HTML::as.title("Contingency Table"), HR = 2, file = stdout())
      
      if (length(dim(x.table)) == 1) {
        cont.result <- x.table
        R2HTML::HTML(cont.result, innerBorder = 1, file = stdout(), row.names = T, align = "left")
      } else if (length(dim(x.table)) == 2) {
        cont.result <- sumTl(x.table)
        R2HTML::HTML(cont.result, innerBorder = 1, file = stdout(), row.names = T, align = "left")
        
        ## when the dimension is 3
      } else if (length(dim(x.table)) == 3) {
        ## when conditional contigency table is selected
        if (cond.tab == T) {
          R2HTML::HTML(R2HTML::as.title("Conditional Contingency Table"), HR = 3, file = stdout())
          
          # For 1st variable
          R2HTML::HTML(R2HTML::as.title(paste0("Control variable - ",names(dimnames(x.table))[1])), HR = 4, file = stdout())
          for(i in 1:dim(x.table)[1]) {
            x.table.dim2 <- sumTl(x.table[i,,])
            R2HTML::HTML(R2HTML::as.title(paste(names(dimnames(x.table))[1], "=", dimnames(x.table)[[1]][i])), 
                         HR = 5, file = stdout())
            R2HTML::HTML(x.table.dim2, innerBorder = 1, file = stdout(), row.names = T, align = "left")
          }
          
          # For 2nd variable
          R2HTML::HTML(R2HTML::as.title(paste0("Control variable - ",names(dimnames(x.table))[2])), HR = 4, file = stdout())
          for(i in 1:dim(x.table)[2]) {
            x.table.dim2 <- sumTl(x.table[,i,])
            R2HTML::HTML(R2HTML::as.title(paste(names(dimnames(x.table))[2], "=", dimnames(x.table)[[2]][i])), 
                         HR = 5, file = stdout())
            R2HTML::HTML(x.table.dim2, innerBorder = 1, file = stdout(), row.names = T, align = "left")
          }
          
          # For 3rd variable
          R2HTML::HTML(R2HTML::as.title(paste0("Control variable - ",names(dimnames(x.table))[3])), HR = 4, file = stdout())
          for(i in 1:dim(x.table)[3]) {
            x.table.dim2 <- sumTl(x.table[,,i])
            R2HTML::HTML(R2HTML::as.title(paste(names(dimnames(x.table))[3], "=", dimnames(x.table)[[3]][i])), 
                         HR = 5, file = stdout())
            R2HTML::HTML(x.table.dim2, innerBorder = 1, file = stdout(), row.names = T, align = "left")
          }
        }
        
        ## when marginal contigency table is selected
        if (margin.tab == T) {
          R2HTML::HTML(R2HTML::as.title("Marginal Contingency Table"), HR = 3, file = stdout())
          getMargin <- function(i,j){
            R2HTML::HTML(R2HTML::as.title(paste(names(dimnames(x.table))[i], "-", names(dimnames(x.table))[j])), HR = 4, file = stdout())
            mar <- margin.table(x.table, c(i, j))
            R2HTML::HTML(sumTl(mar), innerBorder = 1, file = stdout(), row.names = T, align = "left")
          }
          getMargin(1,2); getMargin(1,3); getMargin(2,3)
        }
        
        ## when the dimension over 3
      } else if (length(dim(x.table)) > 3) {
        R2HTML::HTML("Warning messages : Contingency table is not supported for data over 3 dimesions.", file = stdout(), align = "left")
      }
    }
    
    
    #### Goodness-of-fit Tests
    if (fitness) {
      fitness_print	<- data.frame(numeric(0), numeric(0))
      fitness_print	<- rbind(fitness_print, c(deviance(loglm.2), loglm.2$df.residual))
      
      chisq_pearson	<- sum(residuals(loglm.2, type="pearson")^2)
      fitness_print	<- rbind(fitness_print, c(chisq_pearson, summary(loglm.2)$df.residual))
      fitness_print	<- rbind(fitness_print, c(as.numeric(logLik(loglm.2)), attr(logLik(loglm.2), "df")))
      fitness_print	<- rbind(fitness_print, c(AIC(loglm.2), NA))
      fitness_print	<- rbind(fitness_print, c(BIC(loglm.2), NA))
      
      names(fitness_print)	<- c("Value", "DF")
      row.names(fitness_print)	<- c("Deviance", "Pearson's chi-square", "log-likelihood", "AIC", "BIC")
      fitness_print[1]	<- format(round(fitness_print[1], num_digit), nsmall=num_digit)
      
      R2HTML::HTML(R2HTML::as.title("Model Fitness Measurements"), HR=2, file=stdout())
      R2HTML::HTML(fitness_print, file=stdout(), na="", innerBorder = 1, align="left")
    }
    
    
    #### Estimates
    if (estim) {
      R2HTML::HTML(R2HTML::as.title("Coefficient Estimates and Confidence Intervals"), HR=2, file=stdout()) ;
      
      estimate	<- data.frame(coef(summary(loglm.2)))
      colnames(estimate)	<- c("Estimate", "Std.Error", "Z.value","P.value")
      
      if (CI_mode %in% c("wald", "both")) {
        estimate	<- cbind(estimate, confint.default(loglm.2, level = CI_level))
        names(estimate)[ncol(estimate)-1:0]	<- paste0(100*CI_level, "% CI ", c("lower", "upper"), " (Wald)")
      }
      
      if (CI_mode %in% c("likelihood", "both")) {
        fit.msg <- myTryCatch(suppressMessages(matrix(exp(confint(loglm.2, level=CI_level)),ncol=2)))
        # fit.warn <- sapply(fit.msg$warning, function(x) x$message)
        fit.error <- fit.msg$error$message
        
        if (!is.null(fit.error)) {
          R2HTML::HTML(paste("\a Error :", fit.error), file=stdout())
          coef_exp_print <- cbind(estimate, data.frame(NA, NA))
          names(estimate)[ncol(estimate)-1:0]	<- paste0(100*CI_level, "% CI ", c("lower", "upper"), " (Profile likelihood)")
        } else if (is.null(fit.error)) {
          lCI <- suppressMessages(matrix(exp(confint(loglm.2, level=CI_level)),ncol=2))
          rownames(lCI) <- rownames(estimate)
          estimate	<- cbind(estimate,lCI)
          names(estimate)[ncol(estimate)-1:0]	<- paste0(100*CI_level, "% CI ", c("lower", "upper"), " (Profile likelihood)")
        }
     }
    R2HTML::HTML(format(round(estimate, num_digit), nsmall=num_digit), file=stdout(), na="", innerBorder = 1, align="left") ;
      
      if (estim & exp_estim) {
        R2HTML::HTML(R2HTML::as.title("Exponential Coefficient Estimates and Confidence Intervals"), HR=2, file=stdout()) 
        
        estimate	<- data.frame(coef(summary(loglm.2)))
        colnames(estimate)	<- c("Estimate","Std.Error","Z.value","Pvalue")
        exp_Est <- exp(estimate[,1])
        coef_exp_print <- data.frame(exp.Estimate = exp_Est, estimate)
        
        if (CI_mode %in% c("wald", "both")) {
          coef_exp_print	<- cbind(coef_exp_print, exp(confint.default(loglm.2, level = CI_level)))
          names(coef_exp_print)[ncol(coef_exp_print)-1:0]	<- paste0(100*CI_level, "% CI ", c("lower", "upper"), " (Wald)") 
        }
        if (CI_mode %in% c("likelihood", "both")) {
          fit.msg <- myTryCatch(suppressMessages(matrix(exp(confint(loglm.2, level=CI_level)),ncol=2)))
          # fit.warn <- sapply(fit.msg$warning, function(x) x$message)
          fit.error <- fit.msg$error$message
          
          if (!is.null(fit.error)) {
            R2HTML::HTML(paste("\a Error :", fit.error), file=stdout())
            coef_exp_print <- cbind(coef_exp_print, data.frame(NA, NA))
            names(coef_exp_print)[ncol(coef_exp_print)-1:0]	<- paste0(100*CI_level, "% CI ", c("lower", "upper"), " (Profile likelihood)")
          } else if (is.null(fit.error)) {
            lCI <- suppressMessages(matrix(exp(confint(loglm.2, level=CI_level)),ncol=2))
            rownames(lCI) <- rownames(coef_exp_print)
            coef_exp_print	<- cbind(coef_exp_print,lCI) ;
            names(coef_exp_print)[ncol(coef_exp_print)-1:0]	<- paste0(100*CI_level, "% CI ", c("lower", "upper"), " (Profile likelihood)") 
          }
        }
        R2HTML::HTML(format(round(coef_exp_print, num_digit), nsmall=num_digit), file=stdout(), na="", innerBorder = 1, align="left") 
      }
    }
    
    
    #### Testing
    if (testing) {
      nowald_flag	<- FALSE
      
      R2HTML::HTML(R2HTML::as.title("Testing Regression Coefficients"), HR=2, file=stdout()) ;
      
      if (ss_type %in% c("II","III")) {
        
        if (chisq_mode %in% c("wald", "both")) {
          
          err.msg <- myTryCatch(Anova(loglm.2, test="Wald", type=ss_type))$error$message
          if (length(err.msg) > 0) {
            catch.error <- grepl("residual df = 0", err.msg, ignore.case = T)
          } else catch.error <- FALSE
          
          if (catch.error) {
            test_print <- data.frame()
            R2HTML::HTML("Error message : Wald test is not supported for the saturated model.", file = stdout(), align = "left")
          } else {
            anov <- Anova(loglm.2, test="Wald", type=ss_type)
            table_tmp	<- data.frame(anov)
            colnames(table_tmp)	<- c("DF", "Chi-squared (Wald)", "P-value (Wald)")
            
            if (ss_type %in% c("III")) {
              test_print	<- table_tmp[-1,]
            } else {
              test_print	<- table_tmp
            }
          }
        }
        
        if (chisq_mode %in% c("likelihood", "both")) {
          table_tmp	<- data.frame(Anova(loglm.2, test="LR", type=ss_type))
          colnames(table_tmp)	<- c("Chi-squared (LR)", "DF", "P-value (LR)")
          table_tmp	<- table_tmp[c(2,1,3)]
          
          if (chisq_mode == "likelihood") { 
            test_print	<- table_tmp
          } else {
            test_print	<- cbind.fill(test_print, table_tmp[-1])
          }
        }
        
        if (length(test_print) > 0) {
          test_print[2] <- sapply(test_print[[2]], function(x) sprintf(paste0("%.", num_digit, "f"), round(as.numeric(x), num_digit)))
          test_print[3] <-  sapply(test_print[[3]], function(x) {
            if (x < 0.001) {
              form.x <- sprintf(paste0("%.", num_digit, "e"), x) 
            } else form.x <- sprintf(paste0("%.", num_digit, "f"), round(as.numeric(x), num_digit)) 
            return(form.x)
          })
        }
        
      } else if (ss_type %in% c("I")) {
        if (chisq_mode %in% c("wald", "both")) {
          test_print <- data.frame()
          nowald_flag	<- TRUE
        }
        
        if (chisq_mode %in% c("likelihood", "both")) {
          table_tmp	<- data.frame(anova(loglm.2, test="LRT"))[-1,-c(3,4)]
          colnames(table_tmp)	<- c("DF", "Chi-squared (LR)", "P-value (LR)")
          
          test_print	<- table_tmp
          
          test_print[2] <- sapply(test_print[[2]], function(x) sprintf(paste0("%.", num_digit, "f"), round(as.numeric(x), num_digit)))
          test_print[3] <-  sapply(test_print[[3]], function(x) {
            if (x < 0.001) {
              form.x <- sprintf(paste0("%.", num_digit, "e"), x) 
            } else form.x <- sprintf(paste0("%.", num_digit, "f"), round(as.numeric(x), num_digit)) 
            return(form.x)
          })
        }
      }
      
      if (exists("test_print") & length(test_print) > 0) R2HTML::HTML(test_print, file=stdout(), na="", innerBorder = 1, align="left")
      if (nowald_flag) R2HTML::HTML("Error message : Wald test for type I SS is not available.", file=stdout(), na="")
    }

    #### End of output
    R2HTML::HTMLhr(file=stdout())
    R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Log linear Analysis",sep=""),file=stdout())
    R2HTML::HTMLhr(file=stdout())
  })
  
  return(html.output)
}

REx_ContinTabAnal <- function(data, freq.table = T,
                              row.sum = F, col.sum = F, row.var.name = T, col.var.name = T,
                              x = NULL, y = NULL, adj.var = NULL, freq.var = NULL,
                              chsq.test = T, likelihood.test = F, fisher.exact.test = F,
                              lbl.test = F, lbl.x.score = NULL, lbl.y.score = NULL, 
                              whitney.wilcox.test = F,
                              nominal.order.test = F, cmh.y.score = NULL,
                              cmh.test = F,
                              contingency.coefficient = F, phi = F, cramerV = F, lambda = F, uncertainty.coefficient = F,
                              pearson.cor = F, spearman.cor = F,
                              gamma = F, somerD = F, kendall.taub = F, kendall.tauc = F,
                              eta = F, eta.num = NULL,
                              kappa = F, mcnemar = F,
                              freq.obs = T, freq.exp = F, freq.dec = 0,
                              percent.row = F, percent.col = F, percent.tl = F, percent.dec= 2,
                              res.unstnd = F, res.stnd = F, res.adj = F, res.dec = 2) {
  ###################
  #### Variables
  ###################
  
  ### [data tab]
  ## data         : input data
  ## freq.table   : a logical indicating whether input is frequecy table
  ##              (UI 데이터 탭에서 데이터를 선택한 경우 F, 빈도표를 선택한 경우 T)
  ## row.sum      : a logical indicating whether frequecy table contains a row-wise sum varable
  ##              (UI 데이터 탭에서 빈도표를 선택한 경우, 빈도표에 행별 합계가 포함된 경우 T)
  ## col.sum      : a logical indicating whether frequecy table contains a column-wise sum varable
  ##              (UI 데이터 탭에서 빈도표를 선택한 경우, 빈도표에 열별 합계가 포함된 경우 T)
  ## row.var.name : a logical indicating whether frequecy table contains a name of row variables
  ##              (UI 데이터 탭에서 빈도표를 선택한 경우, 빈도표에 행 이름이 포함된 경우 T)
  ## col.var.name : a logical indicating whether frequecy table contains a name of column variables
  ##              (UI 데이터 탭에서 빈도표를 선택한 경우, 빈도표에 열 이름이 포함된 경우 T)
  ## x            : two-dimensional contingency table(matrix) or factor object
  ##              (UI 데이터 탭에서 데이터를 선택하였을 때, 첫번째 변수)
  ## y            : factor object, ignored if x is a contingency table
  ##              (UI 데이터 탭에서 데이터를 선택하였을 때, 두번째 변수)
  ## adj.var      : adjust variables
  ##              (UI 데이터 탭에서 데이터를 선택하였을 때, 제어 변수)
  ## freq.var     : numeric object, frequency for each combination of x and y
  ##              (UI 데이터 탭에서 데이터를 선택하였을 때, 빈도 변수)
  
  ### [statistics tab - Independent test]
  ## chsq.test           : a logical indicating whether to conduct chisquare test
  ##                       (UI 통계량 탭에서 독립성 검정 중 피어슨의 카이제곱 검정을 택하였을 경우 T)
  ## likelihood.test     : a logical indicating whether to conduct likelihood ratio test
  ##                       (UI 통계량 탭에서 독립성 검정 중 우도비 검정을 택하였을 경우 T)
  ## fisher.exact.test   : a logical indicating whether to conduct Fisher's exact test
  ##                       (UI 통계량 탭에서 독립성 검정 중 Fisher의 정확 검정을 택하였을 경우 T)
  ## nominal.order.test  : a logical indicating whether to conduct independent test for nominal*order data 
  ##                       (UI 통계량 탭에서 독립성 검정 중 명목형*순서형 독립성 검정을 택하였을 경우 T)
  ## cmh.y.score         : a numeric vector indicating score for Y variable 
  ##                       (UI 통계량 탭에서 독립성 검정 중 명목형*순서형 독립성 검정을 택하였을 때, y.score의 입력 값, 숫자형 벡터)
  ## lbl.test            : a logical indicating whether to conduct independent test for order*order data 
  ##                       (UI 통계량 탭에서 독립성 검정 중 순서형*순서형 독립성 검정을 택하였을 경우 T)
  ## lbl.x.score         : a numeric vector indicating score for X variable doing lbl.test
  ##                       (UI 통계량 탭에서 독립성 검정 중 순서형*순서형 독립성 검정을 택하였을 때, x.score의 입력 값, 숫자형 벡터)
  ## lbl.y.score         : a numeric vector indicating score for Y variable doing lbl.test
  ##                       (UI 통계량 탭에서 독립성 검정 중 순서형*순서형 독립성 검정을 택하였을 때, y.score의 입력 값, 숫자형 벡터)
  ## whitney.wilcox.test : a logical indicating whether to conduct Whitney Wilcox test
  ##                       (UI 통계량 탭에서 독립성 검정 중 Whitney wilcox 검정을 택하였을 경우 T)
  
  ### [statistics tab - Nominal data]
  ## contingency.coefficient : a logical indicating whether to calculate contingency coefficient
  ##                           (UI 통계량 탭에서 명목형 데이터 중 분할계수를 택하였을 경우 T)
  ## phi                     : a logical indicating whether to calculate Phi
  ##                           (UI 통계량 탭에서 명목형 데이터 중 파이를 택하였을 경우 T)
  ## caramerV                : a logical indicating whether to calculate Cramer's V
  ##                           (UI 통계량 탭에서 명목형 데이터 중 Cramer의 V를 택하였을 경우 T)
  ## lambda                  : a logical indicating whether to calculate lambda
  ##                           (UI 통계량 탭에서 명목형 데이터 중 람다를 택하였을 경우 T)
  ## uncertainty.coefficient : a logical indicating whether to calculate an uncertainty coefficient
  ##                           (UI 통계량 탭에서 명목형 데이터 중 불확실성 계수를 택하였을 경우 T)
  
  ### [statistics tab - Ordinal data]
  ## gamma        : a logical indicating whether to calculate a Gamma
  ##                (UI 통계량 탭에서 순서형 데이터 중 감마를 택하였을 경우 T)
  ## somerD       : a logical indicating whether to calculate a Somer's D
  ##                (UI 통계량 탭에서 순서형 데이터 중 Somer의 D를 택하였을 경우 T)
  ## kendall.taub : a logical indicating whether to calculate Kendall's tau-b
  ##                (UI 통계량 탭에서 순서형 데이터 중 Kendall의 tau-b를 택하였을 경우 T)
  ## kendall.tauc : a logical indicating whether to calculate Kendall's tau-c
  ##                (UI 통계량 탭에서 순서형 데이터 중 Kendall의 tau-c를 택하였을 경우 T)
  
  ### [statistics tab - Interval data]
  ## eta     : a logical indicating whether to calculate an Eta
  ##           (UI 통계량 탭에서 명목:구간데이터 중 에타를 택하였을 경우 T)
  ## eta.num : numbers indicating each group for a nominal variable 
  
  ### [statistics tab ]
  ## kappa    : a logical indicating whether to calclate Kappa value
  ##            (UI 통계량 탭에서 Kappa를 택하였을 경우 T)
  ## McNemar  : a logical indicating whether to calculate McNemar
  ##            (UI 통계량 탭에서 McNemar를 택하였을 경우 T)
  ## cmh.test : a logical indicating whether to conduct Cochran-Mantel-Haenzel Test
  ##            (UI 통계량 탭에서 Cochran 및 Mantel-Haenszel 통계량을 택하였을 경우 T)
  
  ### [output tab - Frequency]
  ## freq.obs : a logical indicating whether to output observed frequency
  ##            (UI 출력옵션 탭에서 빈도 중 관측 빈도를 택하였을 경우 T)
  ## freq.exp : a logical indicating whether to output expected frequency
  ##            (UI 출력옵션 탭에서 빈도 중 기대 빈도를 택하였을 경우 T)
  ## freq.dec : integer indicaitng the number of decimal places when freq.obs or freq.exp is T
  ##            (UI 출력옵션 탭에서 빈도 중 소숫점 자릿수 값, 숫자형 벡터)
  
  ### [output tab - Percentile]
  ## percent.row : a logical indicating whether to output percentile by row
  ##               (UI 출력옵션 탭에서 백분율 중 행을 택하였을 경우 T)
  ## percent.col : a logical indicating whether to output percentile by column
  ##               (UI 출력옵션 탭에서 백분율 중 열을 택하였을 경우 T)
  ## percent.tl  : a logical indicating whether to output percentile by total
  ##               (UI 출력옵션 탭에서 백분율 중 총계를 택하였을 경우 T)
  ## percent.dec : integer indicaitng the number of decimal places when percent.row, percent.col or percent.tl is T
  ##               (UI 출력옵션 탭에서 백분율 중 소숫점 자릿수 값, 숫자형 벡터)
  
  ### [output tab - res.unstnd]
  ## res.unstnd : a logical indicating whether to calculate unstandardized residuals
  ##              (UI 출력옵션 탭에서 잔차 중 비표준화 잔차를 택하였을 경우 T)
  ## res.stnd   : a logical indicating whether to calculate standardized residuals
  ##              (UI 출력옵션 탭에서 잔차 중 표준화 잔차를 택하였을 경우 T)
  ## res.adj    : a logical indicating whether to calculate adjusted standardized residuals
  ##              (UI 출력옵션 탭에서 잔차 중 수정된 표준화 잔차를 택하였을 경우 T)
  ## res.dec    : integer indicaitng the number of decimal places when res.unstnd, res.stnd or res.adj is T
  ##              (UI 출력옵션 탭에서 잔차 중 소숫점 자릿수 값, 숫자형 벡터)
  
  ###################
  #### Required packages
  #### - R2HTML    : to output in the form of HTML
  ####               and to calculate uncertainty coefficient using 'UncertCoef' function
  #### - oii      : to measure association using 'assocstats' function
  #### - coin     : to test linear by linear independence using 'lbl_test', 'wilcox_test' function
  #### - vcd      : for Multi-way Contingency Tables
  #### - vcdExtra : to conduct CMH test using 'CMHtest' function
  ###################
  load.pkg(c("R2HTML", "DescTools", "vcd", "oii", "coin", "vcdExtra", "lsr"))
  
  
  #####################################################################################
  ################### Function : Catch and save warning and error message
  #####################################################################################
  
  myTryCatch <- function(expr) {
    warn <- err <- NULL
    value <- withCallingHandlers(
      tryCatch(expr, error = function(e) {
        err <<- e
        NULL
      }), warning = function(w) {
        warn <<- w
        invokeRestart("muffleWarning")
      })
    
    return(list(warning = warn, error = err))
  }
 
  #####################################################################################
  ################### Function : Formatting small vale
  #####################################################################################
  
  numForm <- function(x) {
    x <- as.numeric(x)
    
    form.x <- ifelse(x == 0, 0, 
                     ifelse(x %% 1 == 0, sprintf("%.0f", x),
                            ifelse(abs(x) < 0.0001, sprintf("%.2e", x), 
                                   ifelse(abs(x) >= 10000, sprintf("%.2e", x), sprintf("%.4f", x)))))
    return(form.x)
  }
  
  
  #####################################################################################
  ################### Function : Drop unused levels from table
  #####################################################################################
  
  dropUnusedLevel <- function(table.dim2) {
    unused1 <- apply(table.dim2, 1, function(x) all(x == 0))
    filter.tab1 <- table.dim2[!unused1, ]
    
    unused2 <- apply(filter.tab1, 2, function(x) all(x == 0))
    filter.tab2 <- filter.tab1[, !unused2]
    
    return(as.table(filter.tab2))
  }
  
  #####################################################################################
  ################### Function : Create contingency tables
  #####################################################################################
  
  #### Combine Variable name and levels
  combVarLevel <- function(table, i) { # if row, i = 1, if column i = 2
    dimname <- dimnames(table)
    
    var.name <- names(dimname)[i]
    level.name <- dimname[[i]]
    var.lev.name <- paste0(var.name, "(", level.name, ")")   
    
    return(var.lev.name)
  }
  
  #### Assign names to a table
  tableName <- function(table.dim2) {
    x.table.dim2 <- as.data.frame.matrix(table.dim2)
    rownames(x.table.dim2) <- combVarLevel(table.dim2, 1)
    colnames(x.table.dim2) <- combVarLevel(table.dim2, 2)
    
    return(x.table.dim2)
  }
  
  #### Calculate row-wise sum and column-wise sum
  new.rowSums <- function(mat) apply(mat,1,function(x) sum(as.numeric(x)))
  new.colSums <- function(mat) apply(mat,2,function(x) sum(as.numeric(x)))
  
  #### Add total sum to a table
  sumTl <- function(x.table.dim2) {
    tab <- tableName(x.table.dim2)
    tab <- cbind(tab, Total = new.rowSums(tab))
    tab <- rbind(tab, Total = new.colSums(tab))
    
    return(tab)
  }
  
  #### Add row sum to a table
  sumRowTl <- function(x.table.dim2) {
    tab <- tableName(x.table.dim2)
    tab <- cbind(tab, Total = new.rowSums(tab))
    
    return(tab)
  }
  
  #### Add column sum to a table
  sumColTl <- function(x.table.dim2) {
    tab <- tableName(x.table.dim2)
    tab <- rbind(tab, Total = new.colSums(tab))
    
    return(tab)
  }
  
  
  #### Create a 2 dimensional table
  table_2dim <- function(x.dim2, subtitle = NULL, current.hr = 2, ...) {
    #### output contingency table in the desired format
    
    if (is.null(subtitle)) {
      subtitle <- ""
    } else subtitle <- paste(":", subtitle)
    
    hr <- current.hr + 1
    
    # calculate frequency and proportion
    cpr <- prop.table(x.dim2, 1)
    cpc <- prop.table(x.dim2, 2)
    cpt <- prop.table(x.dim2)
    rs <- rowSums(x.dim2)
    cs <- colSums(x.dim2)
    ts <- sum(x.dim2)
    
    cst <- suppressWarnings(chisq.test(x.dim2, correct = F))
    
    # the observed freqeucy
    if (freq.obs) {
      freq.obs.print <- format(round(sumTl(cst$observed), freq.dec), nsmall = freq.dec)
      R2HTML::HTML(R2HTML::as.title(paste("Observed Frequency", subtitle)), HR = hr, file=stdout())
      R2HTML::HTML(freq.obs.print, innerBorder = 1, file=stdout(), row.names = T, align = "left")
    }
    
    # the expected freqeucy
    if (freq.exp) {
      freq.exp.print <- format(round(sumTl(cst$expected), freq.dec), nsmall = freq.dec)
      R2HTML::HTML(R2HTML::as.title(paste("Expected Frequency", subtitle)), HR = hr, file=stdout(), stringsAsFactors = F, fix.empty.names = T)
      R2HTML::HTML(freq.exp.print, innerBorder = 1, file=stdout(), row.names = T, align = "left")
    }
    
    # the percent of observations by row
    if (percent.row) {
      percent.row.print <- format(round(sumRowTl(cpr*100), percent.dec), nsmall = percent.dec)
      R2HTML::HTML(R2HTML::as.title(paste("Row Percentage", subtitle)), HR = hr, file=stdout())
      R2HTML::HTML(percent.row.print, innerBorder = 1, file=stdout(), row.names = T, align = "left")
    }
    
    # the percent of observations by column
    if (percent.col) {
      percent.col.print <- format(round(sumColTl(cpc*100), percent.dec), nsmall = percent.dec)
      R2HTML::HTML(R2HTML::as.title(paste("Column Percentage", subtitle)), HR = hr, file=stdout())
      R2HTML::HTML(percent.col.print, innerBorder = 1, file=stdout(), row.names = T, align = "left")
    }
    
    # the percent of observations by total
    if (percent.tl) {
      percent.tl.print <- format(round(sumTl(cpt*100), percent.dec), nsmall = percent.dec)
      R2HTML::HTML(R2HTML::as.title(paste("Total Percentage", subtitle)), HR = hr, file=stdout())
      R2HTML::HTML(percent.tl.print, innerBorder = 1, file=stdout(), row.names = T, align = "left")
    }
    
    # unstandardized residuals
    if (res.unstnd) {
      res.unstnd.print <- format(round(tableName(cst$observed - cst$expected), res.dec), nsmall = res.dec)
      R2HTML::HTML(R2HTML::as.title(paste("Unstandardized Residuals", subtitle)), HR = hr, file=stdout())
      R2HTML::HTML(res.unstnd.print, innerBorder = 1, file=stdout(), row.names = T, align = "left")                        
    }
    
    # standardized residuals
    if (res.stnd) {
      res.stnd.print <- format(round(tableName(cst$residuals), res.dec), nsmall = res.dec)
      R2HTML::HTML(R2HTML::as.title(paste("Standardized Residuals", subtitle)), HR = hr, file=stdout())
      R2HTML::HTML(res.stnd.print, innerBorder = 1, file=stdout(), row.names = T, align = "left")
    }
    
    # adjusted residuals
    if (res.adj) {
      res.adj.print <- format(round(tableName((cst$observed - cst$expected)/sqrt(cst$expected * ((1-rs/ts) %*% t(1-cs/ts)))), res.dec), nsmall = res.dec)
      R2HTML::HTML(R2HTML::as.title(paste("Adjusted Residuals", subtitle)), HR = hr, file=stdout())
      R2HTML::HTML(res.adj.print, innerBorder = 1, file=stdout(), row.names = T, align = "left")
    }
  }
  
  
  #####################################################################################
  ################### Function : Likelihood ratio Test 
  ################### version by Peter Hurd at the University of Alberta and modified it to make the Williams correction (often recommended) the default
  #####################################################################################
  
  g.test <- function(x, y = NULL, correct="williams",
                     p = rep(1/length(x), length(x)), simulate.p.value = FALSE, B = 2000)
    #can also use correct="none" or correct="yates"
  {
    DNAME <- deparse(substitute(x))
    if (is.data.frame(x)) x <- as.matrix(x)
    if (is.matrix(x)) {
      if (min(dim(x)) == 1) 
        x <- as.vector(x)
    }
    if (!is.matrix(x) && !is.null(y)) {
      if (length(x) != length(y)) 
        stop("x and y must have the same length")
      DNAME <- paste(DNAME, "and", deparse(substitute(y)))
      OK <- complete.cases(x, y)
      x <- as.factor(x[OK])
      y <- as.factor(y[OK])
      if ((nlevels(x) < 2) || (nlevels(y) < 2)) 
        stop("x and y must have at least 2 levels")
      x <- table(x, y)
    }
    if (any(x < 0) || any(is.na(x))) 
      stop("all entries of x must be nonnegative and finite")
    if ((n <- sum(x)) == 0) 
      stop("at least one entry of x must be positive")
    #If x is matrix, do test of independence
    if (is.matrix(x)) {
      #Test of Independence
      nrows<-nrow(x)
      ncols<-ncol(x)
      if (correct=="yates"){ # Do Yates' correction?
        if(dim(x)[1]!=2 || dim(x)[2]!=2) # check for 2x2 matrix
          stop("Yates' correction requires a 2 x 2 matrix")
        if((x[1,1]*x[2,2])-(x[1,2]*x[2,1]) > 0)
        {
          x[1,1] <- x[1,1] - 0.5
          x[2,2] <- x[2,2] - 0.5
          x[1,2] <- x[1,2] + 0.5
          x[2,1] <- x[2,1] + 0.5
        }
        else
        {
          x[1,1] <- x[1,1] + 0.5
          x[2,2] <- x[2,2] + 0.5
          x[1,2] <- x[1,2] - 0.5
          x[2,1] <- x[2,1] - 0.5
        }
      }
      
      sr <- apply(x,1,sum)
      sc <- apply(x,2,sum)
      E <- outer(sr,sc, "*")/n
      # are we doing a monte-carlo?
      # no monte carlo GOF?
      if (simulate.p.value){
        METHOD <- paste("Log likelihood ratio (G-test) test of independence\n\t with simulated p-value based on", B, "replicates")
        tmp <- .C("gtestsim", as.integer(nrows), as.integer(ncols),
                  as.integer(sr), as.integer(sc), as.integer(n), as.integer(B),
                  as.double(E), integer(nrows * ncols), double(n+1),
                  integer(ncols), results=double(B), PACKAGE= "ctest")
        g <- 0
        for (i in 1:nrows){
          for (j in 1:ncols){
            if (x[i,j] != 0) g <- g + x[i,j] * log(x[i,j]/E[i,j])
          }
        }
        STATISTIC <- G <- 2 * g
        PARAMETER <- NA
        PVAL <- sum(tmp$results >= STATISTIC)/B
      }
      else {
        # no monte-carlo
        # calculate G
        g <- 0
        for (i in 1:nrows){
          for (j in 1:ncols){
            if (x[i,j] != 0) g <- g + x[i,j] * log(x[i,j]/E[i,j])
          }
        }
        q <- 1
        if (correct=="williams"){ # Do Williams' correction
          row.tot <- col.tot <- 0    
          for (i in 1:nrows){ row.tot <- row.tot + 1/(sum(x[i,])) }
          for (j in 1:ncols){ col.tot <- col.tot + 1/(sum(x[,j])) }
          q <- 1+ ((n*row.tot-1)*(n*col.tot-1))/(6*n*(ncols-1)*(nrows-1))
        }
        STATISTIC <- G <- 2 * g / q
        PARAMETER <- (nrow(x)-1)*(ncol(x)-1)
        PVAL <- 1-pchisq(STATISTIC,df=PARAMETER)
        if(correct=="none")
          METHOD <- "Log likelihood ratio (G-test) test of independence without correction"
        if(correct=="williams")
          METHOD <- "Log likelihood ratio (G-test) test of independence with Williams' correction"
        if(correct=="yates")
          METHOD <- "Log likelihood ratio (G-test) test of independence with Yates' correction"
      }
    }
    else {
      # x is not a matrix, so we do Goodness of Fit
      METHOD <- "Log likelihood ratio (G-test) goodness of fit test"
      if (length(x) == 1) 
        stop("x must at least have 2 elements")
      if (length(x) != length(p)) 
        stop("x and p must have the same number of elements")
      E <- n * p
      
      if (correct=="yates"){ # Do Yates' correction
        if(length(x)!=2)
          stop("Yates' correction requires 2 data values")
        if ( (x[1]-E[1]) > 0.25) {
          x[1] <- x[1]-0.5
          x[2] <- x[2]+0.5
        }
        else if ( (E[1]-x[1]) > 0.25){
          x[1] <- x[1]+0.5
          x[2] <- x[2]-0.5
        }
      }
      names(E) <- names(x)
      g <- 0
      for (i in 1:length(x)){
        if (x[i] != 0) g <- g + x[i] * log(x[i]/E[i])
      }
      q <- 1
      if (correct=="williams"){ # Do Williams' correction
        q <- 1+(length(x)+1)/(6*n)
      }
      STATISTIC <- G <- 2*g/q
      PARAMETER <- length(x) - 1
      PVAL <- pchisq(STATISTIC, PARAMETER, lower = FALSE)
    }
    names(STATISTIC) <- "Log likelihood ratio statistic (G)"
    names(PARAMETER) <- "X-squared df"
    names(PVAL) <- "p.value"
    structure(list(statistic = STATISTIC,
                   parameter = PARAMETER,
                   p.value = PVAL,
                   method = METHOD,
                   data.name = DNAME, 
                   observed = x, 
                   expected = E),
              class = "htest")
  }
  
  
  #####################################################################################
  ################### Function : Uncertainty Coefficient
  #####################################################################################
  
  calc_UC <- function(x) {
    x <- matrix(as.numeric(x), dim(x))
    
    SumR <- rowSums(x)
    SumC <- colSums(x)
    n <- sum(x)
    
    HY <- -sum((SumC / n) * log(SumC / n))
    HX <- -sum((SumR / n) * log(SumR / n))
    HXY <- -sum((x / n) * log(x / n))
    
    UC.RC <- (HX + HY - HXY) / HX
    UC.CR <- (HY + HX - HXY) / HY
    UC.S <- 2 * (HX + HY - HXY) / (HX + HY)
    
    UClist <- list(UC.RC, UC.CR, UC.S)
    names(UClist) <- c("UC.RC", "UC.CR", "UC.S")
    
    return(UClist)
  }
  
  
  #####################################################################################
  ################### Function : Lambda
  #####################################################################################
  
  calc_lambda <- function(x) {
    x <- matrix(as.numeric(x), dim(x))
    
    SumRmax <- sum(apply(x, 1, max))
    SumCmax <- sum(apply(x, 2, max))
    MaxCSum <- max(colSums(x))
    MaxRSum <- max(rowSums(x))
    n <- sum(x)
    
    L.CR <- (SumRmax - MaxCSum) / (n - MaxCSum)
    L.RC <- (SumCmax - max(rowSums(x))) / (n - MaxRSum)
    L.S <- (SumRmax + SumCmax - MaxCSum - MaxRSum) /
      ((2 * n) - MaxCSum - MaxRSum)
    
    Llist <- list(L.CR, L.RC, L.S)
    names(Llist) <- c("L.CR", "L.RC", "L.S")
    
    return(Llist)
  }
  
    
  
  
  #####################################################################################
  ################### Function : Independent Test 
  #####################################################################################
  
  ind_2dim <- function(x.dim2, ...) {
    #### create result format
    result.ind <- data.frame(statistics = character(0), df = numeric(0), value = numeric(0), p.value = numeric(0))
    
    
    #### Analysis
    
    ## nominal - nomial
    # Chi-square test
    if (chsq.test == T) {
      exr.chisq <- suppressWarnings(chisq.test(x.dim2, correct = F))
      if (!is.null(myTryCatch(chisq.test(x.dim2, correct = F))$warning$message)) {
        war.msg1 <- paste("Warning messages :", myTryCatch(chisq.test(x.dim2, correct = F))$warning$message)
      } else war.msg1 <- NULL
          
      exr.chisq.result <- data.frame(statistics = "Pearson's Chi-square",
                                     df = exr.chisq$parameter, 
                                     value = exr.chisq$statistic, 
                                     p.value = exr.chisq$p.value, 
                                     row.names = NULL, stringsAsFactors = F)
      result.ind <- rbind(result.ind, exr.chisq.result)
      
      exr.chisq.con <- suppressWarnings(chisq.test(x.dim2, correct = T))
      if (!is.null( myTryCatch(chisq.test(x.dim2, correct = T))$warning$message)) {
        war.msg2 <- paste("Warning messages :", myTryCatch(chisq.test(x.dim2, correct = T))$warning$message)
      } else war.msg2 <- NULL
      
      exr.chisq.con.result <- data.frame(statistics = "Pearson's Chi-square with Yates' continuity correction",
                                         df = exr.chisq.con$parameter, 
                                         value = exr.chisq.con$statistic, 
                                         p.value = exr.chisq.con$p.value, 
                                         row.names = NULL, stringsAsFactors = F)
      result.ind <- rbind(result.ind, exr.chisq.con.result)
    }
    
    # Likelihood ratio test
    if (likelihood.test == T) {
      exr.gtest <- GTest(x.dim2)
      if (!is.null(myTryCatch(GTest(x.dim2))$warning$message)) {
        war.msg1 <- paste("Warning messages :", myTryCatch(GTest(x.dim2))$warning$message)
        war.msg2 <- war.msg1
      } else war.msg2 <- war.msg1 <- NULL

      exr.gtest.result <- data.frame(statistics = "Likelihood Ratio Chi-square",
                                     df = exr.gtest$parameter,
                                     value = exr.gtest$statistic,
                                     p.value = exr.gtest$p.value,
                                     row.names = NULL, stringsAsFactors = F)
      result.ind <- rbind(result.ind, exr.gtest.result)
    }
    
    # Fisher's exact test
    if (fisher.exact.test == T) {
      exr.fisher.exact <- suppressWarnings(fisher.test(x.dim2, hybrid = T))
      if (!is.null(myTryCatch(fisher.test(x.dim2, hybrid = T))$warning$message)) {
        war.msg1 <- paste("Warning messages :", myTryCatch(fisher.test(x.dim2, hybrid = T))$warning$message)
        war.msg2 <- war.msg1
      } else war.msg2 <- war.msg1 <- NULL
      exr.fisher.exact.result <- data.frame(statistics = "Fisher's exact - odds ratio",
                                            df = NA,
                                            value = if(!is.null(exr.fisher.exact$estimate)) {
                                              exr.fisher.exact$estimate
                                            } else NA,
                                            p.value = exr.fisher.exact$p.value,
                                            row.names = NULL, stringsAsFactors = F)
      result.ind <- rbind(result.ind, exr.fisher.exact.result)
    }
    
    ## nomial - ordinal
    if (nominal.order.test) {
      if(is.null(cmh.y.score)) {
        cmh.y.score1 <- 1:ncol(x.dim2)
      } else {
        y.score.pos <- match(dimnames(x.dim2)[[2]], names(cmh.y.score))
        cmh.y.score1 <- unlist(cmh.y.score)[y.score.pos]
      }
      exr.cmh.r <- CMHtest(x.dim2, rscores = NULL, cscores = cmh.y.score1)
      if (!is.null(myTryCatch(CMHtest(x.dim2, rscores = NULL, cscores = cmh.y.score1))$warning$message)) {
        war.msg1 <- paste("Warning messages :", myTryCatch(CMHtest(x.dim2, rscores = NULL, cscores = cmh.y.score1))$warning$message)
        war.msg2 <- war.msg1
      } else war.msg2 <- war.msg1 <- NULL
      exr.cmh.r.result <- data.frame(statistics = "Cochran-Mantel-Haenszel Statistics",
                                     df = exr.cmh.r$table["rmeans", "Df"],
                                     value = exr.cmh.r$table["rmeans", "Chisq"],
                                     p.value = exr.cmh.r$table["rmeans", "Prob"],
                                     row.names = NULL, stringsAsFactors = F)
      result.ind <- rbind(result.ind, exr.cmh.r.result)
    }
    
    ## ordinal - ordinal
    # Ordinal Chi-square test
    if (lbl.test) {
      dimnames(x.dim2) <- list("X" = dimnames(x.dim2)[[1]], "Y" = dimnames(x.dim2)[[2]])
      
      if(is.null(lbl.x.score)) {
        lbl.x.score1 <- 1:nrow(x.dim2)
      } else {
        lbl.x.score.pos <- match(dimnames(x.dim2)[[2]], names(lbl.x.score))
        lbl.x.score1 <- unlist(lbl.x.score)[lbl.x.score.pos]
      }
      
      if(is.null(lbl.y.score)) {
        lbl.y.score1 <- 1:ncol(x.dim2)
      } else {
        lbl.y.score.pos <- match(dimnames(x.dim2)[[2]], names(lbl.y.score))
        lbl.y.score1 <- unlist(lbl.y.score)[lbl.y.score.pos]
      }
      
      exr.lbl <- lbl_test(x.dim2, scores = list("X" = lbl.x.score1, "Y" = lbl.y.score1))
      if (!is.null(myTryCatch(lbl_test(x.dim2, scores = list("X" = lbl.x.score1, "Y" = lbl.y.score1)))$warning$message)) {
        war.msg1 <- paste("Warning messages :", myTryCatch(lbl_test(x.dim2, scores = list("X" = lbl.x.score1, "Y" = lbl.y.score1)))$warning$message)
        war.msg2 <- war.msg1
      } else war.msg2 <- war.msg1 <- NULL
      exr.lbl.result <- data.frame(statistics = "Linear-by-Linear Association Test",
                                   df = 1,
                                   value = exr.lbl@statistic@teststatistic,
                                   p.value = (1-pnorm(abs(exr.lbl@statistic@teststatistic)))*2,
                                   row.names = NULL, stringsAsFactors = F)
      result.ind <- rbind(result.ind, exr.lbl.result)
    }
    
    # Wilcoxon-Mann-Whitney test
    if (whitney.wilcox.test) {
      dimnames(x.dim2) <- list("X" = dimnames(x.dim2)[[1]], "Y" = dimnames(x.dim2)[[2]])
      x.expand <- expand.table(x.dim2)
      exr.wilcox <- wilcox_test(as.numeric(ordered(Y)) ~ X, data = x.expand)
      if (!is.null(myTryCatch(wilcox_test(as.numeric(ordered(Y)) ~ X, data = x.expand))$warning$message)) {
        war.msg1 <- paste("Warning messages :", myTryCatch(wilcox_test(as.numeric(ordered(Y)) ~ X, data = x.expand))$warning$message)
        war.msg2 <- war.msg1
      } else war.msg2 <- war.msg1 <- NULL
      exr.wilcox.result <- data.frame(statistics = "Asymptotic Wilcoxon-Mann-Whitney Test",
                                      df = 1,
                                      value = exr.wilcox@statistic@teststatistic,
                                      p.value = (1-pnorm(abs(exr.wilcox@statistic@teststatistic)))*2,
                                      row.names = NULL, stringsAsFactors = F)
      result.ind <- rbind(result.ind, exr.wilcox.result)
      
    }
    
    
    #### output independence test result in the desired format
    ## adjust the digits
    result.ind$value <- sprintf("%.4f", result.ind$value)
    result.ind$p.value <- sprintf("%.4f", result.ind$p.value)
    
    ## warning message
    if (any(!is.null("war.msg1"), !is.null("war.msg2"))) {
      war <- unique(c(war.msg1, war.msg2))
     }
    
    return(list(result = result.ind, message = war))
  }
  

  #####################################################################################
  ################### Function : Association Test  
  #####################################################################################
  
  asso_2dim <- function(x.dim2, ...) {
    #### create result format
    result.assoc <- data.frame(statistics = character(0), value = numeric(0))
    
    
    #### Analysis
    
    ## nominal
    
    if (any(contingency.coefficient, phi, cramerV,
            gamma, somerD, kendall.taub, kendall.tauc)) {
      exr.assoc <- association.measures(x.dim2)
    }
    
    # contingency coefficient
    if (contingency.coefficient) {
      exr.cont.coef.result <- data.frame(statistics = "Contingency Coefficient",
                                         value = exr.assoc$contingency,
                                         row.names = NULL, stringsAsFactors = F)
      result.assoc <- rbind(result.assoc, exr.cont.coef.result)
    }
    
    # Phi coefficient
    if (phi) {
      exr.phi.result <- data.frame(statistics = "Phi-coefficient",
                                   value = exr.assoc$phi,
                                   row.names = NULL, stringsAsFactors = F)
      result.assoc <- rbind(result.assoc, exr.phi.result)
    }
    
    # Cramer's V
    if (cramerV) {
      exr.cramerV.result <- data.frame(statistics = "Cramer's V",
                                       value = exr.assoc$cramer,
                                       row.names = NULL, stringsAsFactors = F)
      result.assoc <- rbind(result.assoc, exr.cramerV.result)
    }
    
    # Lambda
    if (lambda) {
      exr.lambda <- calc_lambda(x.dim2)
      exr.lambda.sym <- exr.lambda$L.S
      exr.lambda.row <- exr.lambda$L.RC
      exr.lambda.col <- exr.lambda$L.CR
      
      exr.lambda.result <- data.frame(statistics = c("Lambda (Symmetric)", "Lambda (Row|Column)", "Lambda (Column|Row)"),
                                      value = c(exr.lambda.sym, exr.lambda.row, exr.lambda.col),
                                      row.names = NULL, stringsAsFactors = F)
      result.assoc <- rbind(result.assoc, exr.lambda.result)
    }
    
    # Uncertainty coefficient
    if (uncertainty.coefficient) {
      exr.uncert.coef <- calc_UC(x.dim2)
      exr.uncert.coef.sym <- exr.uncert.coef$UC.S
      exr.uncert.coef.row <- exr.uncert.coef$UC.RC
      exr.uncert.coef.col <- exr.uncert.coef$UC.CR
      
      exr.uncert.coef.result <- data.frame(statistics = c("Uncertainty Coefficient (Symmetric)", "Uncertainty Coefficient (Row|Column)", "Uncertainty Coefficient (Column|Row)"),
                                           value = c(exr.uncert.coef.sym, exr.uncert.coef.row, exr.uncert.coef.col),
                                           row.names = NULL, stringsAsFactors = F)
      result.assoc <- rbind(result.assoc, exr.uncert.coef.result)
    }
    
    
    ## Ordinal
    # Pearson Correlation
    if (pearson.cor) {
      x.expand <- apply(expand.table(x.dim2), 2, function(x) as.numeric(as.factor(x)))
      exr.pearson <- cor.test(x.expand[,1], x.expand[,2], method = "pearson", alternative = "two.sided")
      exr.pearson.result <- data.frame(statistics = "Pearson Correlation",
                                       value = exr.pearson$statistic,
                                       row.names = NULL, stringsAsFactors = F)
      result.assoc <- rbind(result.assoc, exr.pearson.result)
    }
    
    # Spearman Correlation
    if (spearman.cor) {
      x.expand <- apply(expand.table(x.dim2), 2, function(x) as.numeric(as.factor(x)))
      exr.spearman <- cor.test(x.expand[,1], x.expand[,2], method = "spearman", alternative = "two.sided")
      exr.spearman.result <- data.frame(statistics = "Spearman Correlation",
                                        value = exr.spearman$statistic,
                                        row.names = NULL, stringsAsFactors = F)
      result.assoc <- rbind(result.assoc, exr.spearman.result)
    }
    
    # Goodman-Kruskal Gamma
    if (gamma) {
      exr.gamma.result <- data.frame(statistics = "Goodman-Kruskal Gamma",
                                     value = exr.assoc$gamma,
                                     row.names = NULL, stringsAsFactors = F)
      result.assoc <- rbind(result.assoc, exr.gamma.result)
      
    }
    
    # Somer's D
    if (somerD) {
      exr.somerD.result <- data.frame(statistics = "Somer's D",
                                      value = exr.assoc$somersd,
                                      row.names = NULL, stringsAsFactors = F)
      result.assoc <- rbind(result.assoc, exr.somerD.result)
    }
    
    # Kendall's tau-b
    if (kendall.taub) {
      exr.taub.result <- data.frame(statistics = "Kendall's tau-b",
                                    value = exr.assoc$taub,
                                    row.names = NULL, stringsAsFactors = F)
      result.assoc <- rbind(result.assoc, exr.taub.result)
    }
    
    # Kendall's tau-c
    if (kendall.tauc) {
      exr.tauc.result <- data.frame(statistics = "Kendall's tau-c",
                                    value = exr.assoc$tauc,
                                    row.names = NULL, stringsAsFactors = F)
      result.assoc <- rbind(result.assoc, exr.tauc.result)
    }
    
    # Eta
    if (eta) {
      x.eta <- x.dim2
      x.eta.expand <- expand.table(x.eta)
      
      if (is.null(eta.num)) {
        x.eta.expand[,2] <- as.numeric(as.factor(x.eta.expand[,2]))
      } else {
        x.eta.expand[, 2] <- as.character(x.eta.expand[,2])
        
        level <- names(eta.num)
        num <- unlist(eta.num)
        for(i in 1:length(level)) {
          pos <- which(x.eta.expand[,2] == level[i])
          x.eta.expand[pos, 2] <- num[i] 
        }
        x.eta.expand[, 2] <- as.numeric(x.eta.expand[, 2])
      }
      
      x.aov <- aov(x.eta.expand[,2] ~ x.eta.expand[, 1])
      
      exr.eta <- etaSquared(x.aov)
      exr.eta.result <- data.frame(statistics = "Eta-squared",
                                   value = exr.eta[, "eta.sq"],
                                   row.names = NULL, stringsAsFactors = F)
      result.assoc <- rbind(result.assoc, exr.eta.result)
    }
    
    
    #### output independence test result in the desired format
    
    ## adjust digits
    result.assoc$value <- numForm(result.assoc$value)
    
    return(result.assoc)
  }
  
  
  
  
  #####################################################################################
  ################### Function : Agreement for paired data 
  #####################################################################################
  
  #### Kappa
  kappa_2dim <- function(x.dim2, ...) {
    exr.kappa <- Kappa(x.dim2)
    exr.kappa.result <- data.frame(value = numForm(c(exr.kappa$Unweighted["value"], exr.kappa$Weighted["value"])),
                                   ASE = numForm(c(exr.kappa$Unweighted["ASE"], exr.kappa$Weighted["ASE"])),
                                   Z = numForm(c(exr.kappa$Unweighted["value"]/exr.kappa$Unweighted["ASE"], 
                                         exr.kappa$Weighted["value"]/exr.kappa$Weighted["ASE"])),
                                   p.value = numForm(c((1-pnorm(abs(exr.kappa$Unweighted["value"]/exr.kappa$Unweighted["ASE"])))*2, 
                                               (1-pnorm(abs(exr.kappa$Weighted["value"]/exr.kappa$Weighted["ASE"])))*2)))
    rownames(exr.kappa.result) <- c("Unweighted", "Weighted")
    
    return(exr.kappa.result)
  }
  
  #### McNemar
  mcnemar_2dim <- function(x.dim2, ...) {
    exr.mcnemar <- mcnemar.test(x.dim2)
    exr.mcnemar.result <- data.frame(statistic = exr.mcnemar$method,
                                     df = numForm(exr.mcnemar$parameter),
                                     value = numForm(exr.mcnemar$statistic),
                                     p.value = numForm(exr.mcnemar$p.value),
                                     row.names = NULL, stringsAsFactors = F)
    
    return(exr.mcnemar.result)
  }
  
  
  
  #####################################################################################
  ################### Analysis
  #####################################################################################
  
  html.output <- capture.output({
    #### HTML file title
    R2HTML::HTML(R2HTML::as.title("Contingency Table Analysis"), HR = 1, file=stdout(), append = F)
    
    #### create contingency table
    if (freq.table) {
      x.table <- as.table(as.matrix(data[1:2, 1:2]))
    } else if (is.null(adj.var) & is.null(freq.var)) {
      x.table <- xtabs(~ data[, x] + data[,y ], drop.unused.levels = F)
      names(dimnames(x.table)) <- c(x, y)
    } else if (is.null(adj.var) & !is.null(freq.var)) {
      x.table <- xtabs(data[, freq.var] ~ data[, x] + data[, y], drop.unused.levels = F)
      names(dimnames(x.table)) <- c(x, y)
    } else if (!is.null(adj.var) & is.null(freq.var)) {
      x.table <- xtabs(~ data[, x] + data[, y] + data[, adj.var], drop.unused.levels = F)
      names(dimnames(x.table)) <- c(x, y, adj.var)
    } else if (!is.null(adj.var) & !is.null(freq.var)) {
      x.table <- xtabs(data[, freq.var] ~ data[, x] + data[, y] + data[, adj.var], drop.unused.levels = F)
      names(dimnames(x.table)) <- c(x, y, adj.var)
    }
    
    n <- sum(x.table)
    p <- ifelse(freq.table == T, 2, 2 + !is.null(adj.var))

    
    #### Data Structure
    R2HTML::HTML(R2HTML::as.title("Data Structure"), HR = 2, file=stdout())
    DS <- matrix(c("Number of observations",n,"Number of variables",p),byrow=T,ncol=2)
    R2HTML::HTML(DS, file=stdout(), row.names = T, align = "left", innerBorder = 1)
    
    
    #### Variable List
    R2HTML::HTML(R2HTML::as.title("Variable List"), HR = 2, file=stdout())
    varlist <- matrix(c("Row Variable (levels)", paste(names(dimnames(x.table))[1], "(", paste(dimnames(x.table)[[1]], collapse = ", "), ")"),
                        "Column Variable (levels)", paste(names(dimnames(x.table))[2], "(", paste(dimnames(x.table)[[2]], collapse = ", "), ")")),
                      byrow = T, ncol = 2)
    if (!is.null(adj.var)) varlist <- rbind(varlist, c("Group Variable (levels)", 
                                                       paste(names(dimnames(x.table))[3], "(", paste(dimnames(x.table)[[3]], collapse = ", "), ")")))
    if (!is.null(freq.var)) varlist <- rbind(varlist,c("Frequency Variable", freq.var))
    R2HTML::HTML(varlist, file=stdout(), row.names = T, align = "left", innerBorder = 1)
    
    ## Check variable
    check.var <- sapply(names(dimnames(x.table)), function(i) {
      if (is.numeric(data[, i])) {
        war.msg <- paste0("Warning message : The type of variable '", i, "' was numeric, but was coerced into character.")
      } else war.msg <- NULL
      
      return(war.msg)
    })
    
    if (!is.null(unlist(check.var))) {
      war <- paste0(unlist(check.var), "<br>")
      R2HTML::HTML(war, file=stdout(), align = "left")
    }
    
    
    #### Analysis Description
    R2HTML::HTML(R2HTML::as.title("Analysis Description"), HR = 2, file=stdout())
    ind.test <- data.frame(factor = c(cmh.test, chsq.test, likelihood.test, fisher.exact.test,
                                      nominal.order.test, lbl.test, whitney.wilcox.test),
                           name = c("Cochran-Mantel-Heanszel Statistic","Pearson Chi-square Test", "Likelihood Ratio Test", "Fisher's Exact Test",
                                    "Nominal*Ordinal Independent Test", "Ordinal*Ordinal Test", "Whitney Wilcox Test"),
                           stringsAsFactors = F)
    assoc.test <- data.frame(factor = c(contingency.coefficient, phi, cramerV, lambda, uncertainty.coefficient, 
                                        gamma, somerD, kendall.taub, kendall.tauc, eta),
                             name = c("Contingency Coefficient", "Phi", "Cramer's V", "Lambda", "Uncertainty Coefficient",
                                      "Gamma", "Somer's D", "Kendall's tau-b", "Kendall's tau-c", "Eta"),
                             stringsAsFactors = F)
    agr.test <- data.frame(factor = c(kappa, mcnemar),
                           name = c("Kappa", "McNemar"),
                           stringsAsFactors = F)
    
    anal.desc <- list(`Independent Test` = ind.test[ind.test$factor == "TRUE", "name"],
                      `Association Test` = assoc.test[assoc.test$factor == "TRUE", "name"],
                      `Agreement Test` = agr.test[agr.test$factor == "TRUE", "name"])
    attributes(anal.desc) <- list(names = names(anal.desc),
                                  row.names = 1:max(length(anal.desc[[1]]), length(anal.desc[[2]]), length(anal.desc[[3]])), 
                                  class = 'data.frame')
    suppressWarnings(R2HTML::HTML(anal.desc, file=stdout(), row.names = F, align = "left", innerBorder = 1))
    
    
    
    #### when the dimension is 2
    if (length(dim(x.table)) == 2) {
      
      # contingency table
      if(any(freq.obs, freq.exp, percent.row, percent.col, percent.tl, percent.dec, res.unstnd, res.stnd, res.adj)) {
        R2HTML::HTMLhr(file=stdout())
        R2HTML::HTML(R2HTML::as.title("Contingency Table"), HR = 2, file=stdout())
        table_2dim(x.dim2 = x.table, current.hr = 2)
      }
      
      level.check <- (any(dim(dropUnusedLevel(x.table)) <= 1) | (length(dim(dropUnusedLevel(x.table))) <= 1))
      if (level.check) {
        err.msg <- "Error message : After dropping unsed levels, only one level has left. Analysis does not work.<br>"
        R2HTML::HTML(err.msg, file=stdout(), align = "left")
     
      } else {
        # independent test
        if(any(chsq.test, likelihood.test, fisher.exact.test, nominal.order.test, lbl.test, whitney.wilcox.test)) {
          ind.result <- ind_2dim(x.dim2 = x.table)
          R2HTML::HTMLhr(file=stdout())
          R2HTML::HTML(R2HTML::as.title("Independence Test"), HR = 2, file=stdout())
          if (length(ind.result$message) > 0) {
            R2HTML::HTML(paste(ind.result$message, collapse = "<br>"), file=stdout(), align = "left")
          }
          R2HTML::HTML(ind.result$result, innerBorder = 1, file=stdout(), row.names = F, align = "left")
          
        }
        
        # association test 
        if (any(contingency.coefficient, phi, cramerV, lambda, uncertainty.coefficient,
                gamma, somerD, kendall.taub, kendall.tauc, eta)) {
          assoc.result <- asso_2dim(x.dim2 = x.table)
          R2HTML::HTMLhr(file=stdout())
          R2HTML::HTML(R2HTML::as.title("Association Test"), HR = 2, file=stdout())
          R2HTML::HTML(assoc.result, innerBorder = 1, file=stdout(), row.names = F, align = "left")
        }
        
        # kappa
        if (kappa) {
          R2HTML::HTMLhr(file=stdout())
          R2HTML::HTML(R2HTML::as.title("Cohen's Kappa Test"), HR = 2, file=stdout())
          if(dim(x.table)[1] == dim(x.table)[2]) {
            kappa.result <- kappa_2dim(x.dim2 = x.table)
            R2HTML::HTML(kappa.result, innerBorder = 1, file=stdout(), row.names = T, align = "left")
          } else {
            err.msg <- "Error message : Kappa can be applied for the confusion matrix.<br>"
            R2HTML::HTML(err.msg, file=stdout(), align = "left")
          }
        }
        
        # mcnemar
        if (mcnemar) {
          R2HTML::HTML(R2HTML::as.title("McNemar's Chi-squared Test"), HR = 2, file=stdout())
          
          if(dim(x.table)[1] == dim(x.table)[2]) {
            mcnemar.result <- mcnemar_2dim(x.dim2 = x.table)
            R2HTML::HTML(mcnemar.result, innerBorder = 1, file=stdout(), row.names = F, align = "left")
          } else {
            err.msg <- "Error message : McNemar's Chi-squared test can be applied for a square matrix.<br>"
            R2HTML::HTML(err.msg, file=stdout(), align = "left")
          }
        }
      }
      
      #### when the dimension is over 2
    } else if (length(dim(x.table)) >= 3) {
      
      # contingency table
      if(any(freq.obs, freq.exp, percent.row, percent.col, percent.tl, res.unstnd, res.stnd, res.adj)) {
        R2HTML::HTMLhr(file=stdout())
        R2HTML::HTML(R2HTML::as.title("Contingency Table"), HR = 2, file=stdout())
        
        # conditional contingency table
        if (sum(freq.obs, freq.exp, percent.row, percent.col, percent.tl, percent.dec, res.unstnd, res.stnd, res.adj) > 0) {
          for(i in 1:dim(x.table)[3]) {
            x.table.dim2 <- x.table[, , i]
            subtt <- paste("(", names(dimnames(x.table))[3], "=", dimnames(x.table)[[3]][i], ")")
            
            R2HTML::HTML(R2HTML::as.title(paste("Conditional Contingency Table", subtt)), HR = 3, file=stdout())
            table_2dim(x.dim2 = x.table.dim2, current.hr = 3)
          }
          
          # marginal contingency table
          R2HTML::HTML(R2HTML::as.title("Marginal Contingency Table"), HR = 3, file=stdout())
          mar.table <- margin.table(x.table, c(1, 2))
          table_2dim(x.dim2 = mar.table, current.hr = 3)
        } 
      }

      level.check <- lapply(1:dim(x.table)[3], function(i) {
        (any(dim(dropUnusedLevel(x.table[,,i])) <= 1) | (length(dim(dropUnusedLevel(x.table[,,i]))) <= 1))
      })
    
     if (do.call(sum, level.check) >= 1) {
        err.msg1 <- "Error message : After dropping unsed levels, only one level has left. Analysis does not work.<br>"
        R2HTML::HTML(err.msg1, file=stdout(), align = "left")
        
      } else {
        
        x.table2 <- lapply(1:dim(x.table)[3], function(i) {
          tmp <- dropUnusedLevel(x.table[,,i])
          return(tmp)
        })
        
        level.check2 <- lapply(1:dim(x.table)[3], function(i) {
          identical(dim(x.table2[[i]]), dim(x.table[,,i]))
        })
        
        ### CMH Test
        if (cmh.test) {
          R2HTML::HTMLhr(file=stdout())
          R2HTML::HTML(R2HTML::as.title("Cochran-Mantel-Haenszel Test"), HR = 2, file=stdout())
          
          if(any(!unlist(level.check2))) {
            err.msg2 <- "Error message : The conditional contingency table contains zero margins. CMH test does not work.<br>"
            R2HTML::HTML(err.msg2, file=stdout(), align = "left")
            
          } else {
            
            exr.cmh <- CMHtest(x.table, strata = adj.var, overall = T)
            
            for(i in 1:length(exr.cmh)) {
              exr.cmh.result <- exr.cmh[i]
              R2HTML::HTML(R2HTML::as.title(paste(paste(exr.cmh.result[[1]]$names, collapse = " by "), "in stratum", names(exr.cmh.result))), 
                           HR = 3, file=stdout())
              rname <- sapply(rownames(exr.cmh.result[[1]]$table), function(x) {
                switch(x, cor = "Nonzero correlation", rmeans = "Row mean scores differ",
                       cmeans = "Col mean scores differ", general = "General association")
              })
              
              exr.cmh.print <- data.frame('Alternative Hypothesis' = rname,
                                          Chisq = numForm(exr.cmh.result[[1]]$table[, "Chisq"]),
                                          DF = numForm(exr.cmh.result[[1]]$table[, "Df"]),
                                          Pvalue = numForm(exr.cmh.result[[1]]$table[, "Prob"]))
              R2HTML::HTML(exr.cmh.print, innerBorder = 1, file=stdout(), row.names = F, align = "left")
            }
          }
        }
        
        # independent test
        if(any(chsq.test, likelihood.test, fisher.exact.test, nominal.order.test, lbl.test, whitney.wilcox.test)) {
          R2HTML::HTMLhr(file=stdout())
          R2HTML::HTML(R2HTML::as.title("Independence Test"), HR = 2, file=stdout())
          
          for(i in 1:dim(x.table)[3]) {
            x.table.dim2 <- dropUnusedLevel(x.table[, , i]) 
            subtt <- paste(names(dimnames(x.table))[1], "by", names(dimnames(x.table))[2], "in stratum", names(dimnames(x.table))[3], ":", dimnames(x.table)[[3]][i])
            
            R2HTML::HTML(R2HTML::as.title(subtt), HR = 3, file=stdout())
            ind.result <- ind_2dim(x.dim2 = x.table.dim2)
            if (length(ind.result$message) > 0) {
              R2HTML::HTML(paste(ind.result$message, collapse = "<br>"), file=stdout(), align = "left")
            }
            R2HTML::HTML(ind.result$result, innerBorder = 1, file=stdout(), row.names = F, align = "left")
          }
        }
        
        # association test 
        if(any(contingency.coefficient, phi, cramerV, lambda, uncertainty.coefficient,
               gamma, somerD, kendall.taub, kendall.tauc, eta)) {
          R2HTML::HTMLhr(file=stdout())
          R2HTML::HTML(R2HTML::as.title("Association Test"), HR = 2, file=stdout())
          
          for(i in 1:dim(x.table)[3]) {
            x.table.dim2 <- dropUnusedLevel(x.table[, , i])
            subtt <- paste(names(dimnames(x.table))[1], "by", names(dimnames(x.table))[2], "in", names(dimnames(x.table))[3], ":", dimnames(x.table)[[3]][i])
            
            R2HTML::HTML(R2HTML::as.title(subtt), HR = 3, file=stdout())
            asso.result.tmp <- asso_2dim(x.dim2 = x.table.dim2)
            asso.result <- asso.result.tmp[complete.cases(asso.result.tmp), ]
            R2HTML::HTML(asso.result, innerBorder = 1, file=stdout(), row.names = F, align = "left")
          }
        }
        
        # Kappa test 
        if(kappa) {
          R2HTML::HTMLhr(file=stdout())
          R2HTML::HTML(R2HTML::as.title("Cohen's Kappa Test"), HR = 2, file=stdout())
          
          for(i in 1:dim(x.table)[3]) {
            x.table.dim2 <- dropUnusedLevel(x.table[, , i])
            subtt <- paste(names(dimnames(x.table))[1], "by", names(dimnames(x.table))[2], "in", names(dimnames(x.table))[3], ":", dimnames(x.table)[[3]][i])
            
            R2HTML::HTML(R2HTML::as.title(subtt), HR = 3, file=stdout())
            if(dim(x.table.dim2)[1] == dim(x.table.dim2)[2]) {
              kappa.result <- kappa_2dim(x.dim2 = x.table.dim2)
              R2HTML::HTML(kappa.result, innerBorder = 1, file=stdout(), row.names = T, align = "left")
            } else {
              err.msg <-"Error message : Kappa can be applied for the confusion matrix.<br>"
              R2HTML::HTML(err.msg, file=stdout(), align = "left")
            }
          }
        }
        
        # Mcnemar
        if(mcnemar) {
          R2HTML::HTMLhr(file=stdout())
          R2HTML::HTML(R2HTML::as.title("McNemar's Chi-squared Test"), HR = 2, file=stdout())
          
          for(i in 1:dim(x.table)[3]) {
            x.table.dim2 <- dropUnusedLevel(x.table[, , i])
            subtt <- paste(names(dimnames(x.table))[1], "by", names(dimnames(x.table))[2], "in", names(dimnames(x.table))[3], ":", dimnames(x.table)[[3]][i])
            
            R2HTML::HTML(R2HTML::as.title(subtt), HR = 3, file=stdout())
            
            if(dim(x.table.dim2)[1] == dim(x.table.dim2)[2]) {
              mcnemar.result <- mcnemar_2dim(x.dim2 = x.table.dim2)
              R2HTML::HTML(mcnemar.result, innerBorder = 1, file=stdout(), row.names = F, align = "left")
            } else {
              err.msg <- "Error message : McNemar's Chi-squared test can be applied for a square matrix.<br>"
              R2HTML::HTML(err.msg, file=stdout(), align = "left")
            }
          }
        }
      }
    }
    
    R2HTML::HTMLhr(file=stdout())
    R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Contingency Table Analysis",sep=""), file=stdout())
    R2HTML::HTMLhr(file=stdout())
  })
  return(html.output)
}

REx_PLS <- function(dataset,dep_numeric_var=NULL,dep_cat_var=NULL,dep_cat_baseline=NULL,
indep_numeric_var=NULL,indep_cat_var=NULL,indep_cat_baseline=NULL,dimnum=2,crosval=TRUE,
Xscore=FALSE,Yscore=FALSE,Pred=FALSE,Resid=FALSE){		
			
	#### 변수 설명 ####
	## 종속변수는 missing값이 있으면 안됨(dep_numeric_var,dep_car_var)
	## dep_numeric_var : 질적종속변수(dependent variables)(2개 이상의 변수가 있을수 있음)
	## dep_cat_var : 양적종속변수(dependent variables)(2개 이상의 변수가 있을수 있음)
	## dep_cat_baseline : 기저범주
	## indep_numeric_var : 질적설명변수(independent variables) (2개 이상의 변수가 있는 경우 +로 구분, 필수사항아님,단,상수항옵션이(noint=TRUE)인경우 돌아갈수없음)
	## indep_cat_var : 양적설명변수(independent variables) (자동factor처리,2개 이상의 변수가 있는 경우 +로 구분, 필수사항아님,단,상수항옵션이(noint=TRUE)인경우 돌아갈수없음)
	## indep_cat_baseline : 기저범주
	## dimnum : 출력옵션에서 dimension 수 선택 최소 2
	## Xscore : 출력옵션(X-scroed (T-componets), default=FALSE)
	## Yscore : 출력옵션(Y-scroed (U-componets), default=FALSE)
	## Pred : 출력옵션(Y-predicted, default=FALSE)
	## Resid : 출력옵션(Residulas, default=FALSE)
	## 도표관련 옵션은 함수안에 코드로 구성되어있지만 UI구현 안함(Plot=FALSE이면 활성화 안됨)
	## Plot : 출력옵션 (Plot, default=FALSE)
	## plot_option : 출력옵션 (Plot=TRUE인경우, "observations" 또는 "variables"(default)중 선택할수 있음)
	## comp          : 출력옵션 (Plot=TRUE인경우,  X vs Y 축에 표현할 dimension의 숫자 선택가능 (숫자는 dimnum이상이 될 수 없음) )
	## comp (comp_x) : x축의 dimension의 숫자 (default=1<=dimnum)
	## comp (comp_y) : y축의 dimension의 숫자 (default=2<=dimnum)
	## where      	    : 출력옵션 (Plot=TRUE인경우,  X vs Y 축에 표현할 components의 조합)
	## where (wh_options) :  Possible options are: c("t","u") for using xy components
     	##                                     c("t","t"), for using x components (default)
      	##                                     c("u","u") for using y components. 
	## show.names (sn_options) : 출력옵션 (Plot=TRUE인경우, observation name으로 점 출력할 것인가? TRUE(default)/FALSE)
	###################
	#### Required packages : R2HTML , rms, MASS, plsdepot  ####
	###################
	#load.pkg(c("R2HTML", "rms", "MASS", "plsdepot"))
	library("R2HTML")
	library("rms")
	library("MASS")
	library("plsdepot")

	html.output <- capture.output({
		# Title
		R2HTML::HTML(R2HTML::as.title("Partial Least Squares (PLS) Data Analysis Methods"),HR=1,file=stdout(),append=FALSE)

		## Warnings


		## dependent variables
		## numeric
		if(!is.null(dep_numeric_var)) {
			is.nom <- sapply(dep_numeric_var,function(i) !is.numeric(dataset[,i]))
			if(any(is.nom)){
				warn.msg2 <- paste0("\a Warning : The type of variable '",paste(dep_numeric_var[is.nom],collapse=', '),"' is not continuous but selected as the quantitative variable. It was excluded from the analysis.")
				ec <- dep_numeric_var[is.nom]
				dep_numeric_var <- dep_numeric_var[!dep_numeric_var%in%ec]
			}
		}
		if(length(dep_numeric_var)==0) indep_numeric_var <- NULL
		##category
		if(!is.null(dep_cat_var)) {
			is.num <- sapply(dep_cat_var,function(i) is.numeric(dataset[,i]))
			if(any(is.num)) warn.msg3 <- paste0("\a Warning : The type of variable '",paste(dep_cat_var[is.num],collapse=', '),"' is continuous but selected as the qualitative variable. It was coreced into character.")
			for(i in dep_cat_var) dataset[,i] <- factor(dataset[,i])
		}
		if(is.null(dep_numeric_var)&is.null(dep_cat_var)) warn.msg1 <- '\a Error : At least 1 independent variable should be selected. Analysis has been stopped.'
		## independent variables
		## numeric
		if(!is.null(indep_numeric_var)) {
			is.nom <- sapply(indep_numeric_var,function(i) !is.numeric(dataset[,i]))
			if(any(is.nom)){
				warn.msg5 <- paste0("\a Warning : The type of variable '",paste(indep_numeric_var[is.nom],collapse=', '),"' is not continuous but selected as the quantitative variable. It was excluded from the analysis.")
				ec <- indep_numeric_var[is.nom]
				indep_numeric_var <- indep_numeric_var[!indep_numeric_var%in%ec]
			}
		}
		if(length(indep_numeric_var)==0) indep_numeric_var <- NULL
		## category
		if(!is.null(indep_cat_var)) {
			is.num <- sapply(indep_cat_var,function(i) is.numeric(dataset[,i]))
			if(any(is.num)) warn.msg6 <- paste0("\a Warning : The type of variable '",paste(indep_cat_var[is.num],collapse=', '),"' is continuous but selected as the qualitative variable. It was coreced into character.")
			for(i in indep_cat_var) dataset[,i] <- factor(dataset[,i])
		}

		if(is.null(indep_numeric_var)&is.null(indep_cat_var)) warn.msg4 <- '\a Error : At least 1 independent variable should be selected. Analysis has been stopped.'

		if(!crosval&is.null(dimnum)) warn.msg7 <- '\a Error : TRUE must be selected for cross validation when number of dimension is NULL. Analysis has been stopped.'
		
		# warnings
		if(exists('warn.msg1')|exists('warn.msg4')|exists('warn.msg7')|exists('warn.msg8')){
			R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
			if(exists('warn.msg1')) R2HTML::HTML(warn.msg1,file=stdout())
			if(exists('warn.msg2')) R2HTML::HTML(warn.msg2,file=stdout())
			if(exists('warn.msg4')) R2HTML::HTML(warn.msg4,file=stdout())
			if(exists('warn.msg5')) R2HTML::HTML(warn.msg5,file=stdout())
			if(exists('warn.msg7')) R2HTML::HTML(warn.msg7,file=stdout())
			if(exists('warn.msg8')) R2HTML::HTML(warn.msg8,file=stdout())
		} else {
			## warnings
			if(exists('warn.msg2')|exists('warn.msg3')|exists('warn.msg5')|exists('warn.msg6'))R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
			if(exists('warn.msg2')) R2HTML::HTML(warn.msg2,file=stdout())
			if(exists('warn.msg3')) R2HTML::HTML(warn.msg3,file=stdout())
			if(exists('warn.msg5')) R2HTML::HTML(warn.msg5,file=stdout())
			if(exists('warn.msg6')) R2HTML::HTML(warn.msg6,file=stdout())

			## Data Processing ## 
			var_info <- c(dep_numeric_var,dep_cat_var,indep_numeric_var, indep_cat_var)
			is.na <- sapply(var_info,function(i) is.na(dataset[,i]))
			raw.dat <- dataset
			oo <- which(rowSums(is.na)>0)
			if(length(oo)>0) dataset <- dataset[-oo,]


			if(!is.null(dep_cat_var)) {
				dep_o1 <- gsub("::.+","",dep_cat_baseline)
				dep_o2 <- gsub(".+::","",dep_cat_baseline)
				names(dep_o2) <- dep_o1
				sapply(dep_o1,function(i) relevel(factor(dataset[,i]),ref=dep_o2[i]))
			}
			if(!is.null(indep_cat_var)){
				indep_o1 <- gsub("::.+","",indep_cat_baseline)
				indep_o2 <- gsub(".+::","",indep_cat_baseline)
				names(indep_o2) <- indep_o1
				sapply(indep_o1,function(i) relevel(factor(dataset[,i]),ref=indep_o2[i]))
			}
			### 분석을 위해서 종속변수/독립변수 행렬 및 벡터는 무조건 numeric 형태여야함
			### 범주형변수는 numeric 0/1벡터로 바꿔서 범주가 n개이면 새로생성되는 벡터(변수)는 n-1개
			f <- function(i,dataset,cat,cat_baseline){
				level <- levels(dataset[,cat[i]])
				nn <- length(level)
				resout <- as.data.frame(matrix(0,nrow(dataset),nn))
				colnames(resout) <- paste(cat[i],"_",level,sep="")
				oo <- which(level==cat_baseline[i])		
				for(j in 1:nn){
					oo1 <- which(dataset[,cat[i]]==level[j])
					resout[oo1,j] <- 1
				}
				if(nn>2){
					resout <- resout[,-oo]
				}else{	
					resout <- as.data.frame(resout[,-oo])
					colnames(resout) <- paste(cat[i],"_",level[-oo],sep="")
				}
				resout
			}

			## 종속변수에 범주형 없을때
			if(is.null(dep_cat_var)){
				## 독립변수에 범주형 없을때
				if(is.null(indep_cat_var)){
					responses  <- as.data.frame(dataset[,dep_numeric_var])
					colnames(responses) <- dep_numeric_var
					predictors <- as.data.frame(dataset[,indep_numeric_var])				
					colnames(predictors) <- indep_numeric_var

				## 독립변수에 범주형 있을때
				}else{
					responses <- as.data.frame(dataset[,dep_numeric_var])
					colnames(responses) <- dep_numeric_var

					nn1 <- length(indep_cat_var)					
					if(nn1==1){
						indep_cat_new <- f(1,dataset,indep_o1,indep_o2)
					}else{
						indep_cat_new <- f(1,dataset,indep_o1,indep_o2)
						for(k in 2:nn1){
							indep_cat <- f(k,dataset,indep_o1,indep_o2)
							indep_cat_new <- cbind(indep_cat_new,indep_cat)
						}
					}
				
					## 독립변수에 연속형변수 여부에따라서 predictors 구성
					if(!is.null(indep_numeric_var)){
						predictors <- cbind(dataset[,indep_numeric_var],indep_cat_new)
						colnames(predictors) <- c(indep_numeric_var,colnames(indep_cat_new))
					}else{
						predictors <- indep_cat_new
					}
				}	
			## 종속변수에 범주형 있을때
			}else{
				nn2 <- length(dep_cat_var)
				if(nn2==1){
						dep_cat_new <- f(1,dataset,dep_o1,dep_o2)
				}else{
						dep_cat_new <- f(1,dataset,dep_o1,dep_o2)
						for(k in 2:nn2){
							dep_cat <- f(k,dataset,dep_o1,dep_o2)
							dep_cat_new <- cbind(dep_cat_new,dep_cat)
						}
				}
				## 종속변수에 연속형변수 여부에따라서 responses 구성
				if(!is.null(indep_numeric_var)){
					responses <- cbind(dataset[,dep_numeric_var],dep_cat_new)
					colnames(responses) <- c(dep_numeric_var,colnames(dep_cat_new))
				}else{
					responses <- dep_cat_new
				}

				## 독립변수에 범주형 없을때
				if(is.null(indep_cat_var)){
					predictors <- as.data.frame(dataset[,indep_numeric_var])				
					colnames(predictors) <- indep_numeric_var					
				## 독립변수에 범주형 있을때
				}else{
					nn3 <- length(indep_cat_var)					
					if(nn3==1){
						indep_cat_new <- f(1,dataset,indep_o1,indep_o2)
					}else{
						indep_cat_new <- f(1,dataset,indep_o1,indep_o2)
						for(k in 2:nn3){
							indep_cat <- f(k,dataset,indep_o1,indep_o2)
							indep_cat_new <- cbind(indep_cat_new,indep_cat)
						}
					}

				## 독립변수에 연속형변수 여부에따라서 predictors 구성
					if(!is.null(indep_numeric_var)){
						predictors <- cbind(dataset[,indep_numeric_var],indep_cat_new)
						colnames(predictors) <- c(indep_numeric_var,colnames(indep_cat_new))
					}else{
						predictors <- indep_cat_new
					}
				}
			}


			## responses : 종속변수 dataframe
			## predictors : 독립변수 dataframe

			## model
			form.0 <- paste0(paste(colnames(responses),collapse=' , '))
			#form <- ifelse(noint, paste0(form.0,' - 1'),form.0)
			form.1 <- paste0(paste(colnames(predictors),collapse=' , '))

			## Fitting PLS
			if(ncol(responses)==1){
				command_str	<- paste0("res_PLS_1<- plsreg1(predictors,responses,comps=dimnum,crosval=crosval)")
			}else{
				command_str	<- paste0("res_PLS_1<- plsreg2(predictors,responses,comps=dimnum,crosval=crosval)")
			}
			eval(parse(text=command_str)) ;
			#res_PLS_1 <- plsreg1(predictors,responses,comps=dimnum,crosval=crosval)
			#res_PLS_1 <- plsreg2(predictors,responses,comps=dimnum,crosval=crosval)

			### Default output
			# Data Structure
			R2HTML::HTML(R2HTML::as.title("Data Structure"),HR=2,file=stdout())
			DS <- matrix(c('Number of observations','Number of Response variables','Number of Predictor variables','Number of Response dimensios','Number of Predictor dimensions',nrow(dataset),length(c(dep_numeric_var,dep_cat_var)),length(c(indep_numeric_var,indep_cat_var)),ncol(responses),ncol(predictors)),ncol=2)
			R2HTML::HTML(DS,file=stdout(), innerBorder = 1,align="left")
			R2HTML::HTML(R2HTML::as.title("Model Information"),HR=2,file=stdout())
			DS1 <- matrix(c('Reponses','Predictors',form.0,form.1),ncol=2)
			R2HTML::HTML(DS1,file=stdout(), innerBorder = 1,align="left")
			#if(exists('warn.msg3')) R2HTML::HTML(warn.msg3,file=stdout())
			#### Results
			R2HTML::HTML(R2HTML::as.title("Results of Partial Least Squares Data Analysis Methods"),HR=2,file=stdout())

			#################### response의 열이 하나인경우 출력결과(plsreg1) ##############
			if(ncol(responses)==1){
				####### Results of Models
				Xscores <- as.data.frame(res_PLS_1$x.scores)
				Xload   <- as.data.frame(res_PLS_1$x.loads)
				Yscores <- as.data.frame(res_PLS_1$y.scores)
				rownames(Yscores) <- rownames(Xscores)
				Yload   <- as.data.frame(res_PLS_1$y.loads)

				corxyt <- as.data.frame(res_PLS_1$cor.xyt)
			
				raw_weight <- as.data.frame(res_PLS_1$raw.wgs)
				mod_weight <- as.data.frame(res_PLS_1$mod.wgs)

				std <- as.data.frame(res_PLS_1$std.coefs)
				regcoeff <- as.data.frame(res_PLS_1$reg.coefs)	

				pred <- as.data.frame(res_PLS_1$y.pred)
				rownames(pred) <- rownames(Xscores)
				resid <- as.data.frame(res_PLS_1$resid)
				rownames(resid) <- rownames(Xscores)
				R2 <- as.data.frame(res_PLS_1$R2)
				R2xy <- as.data.frame(res_PLS_1$R2Xy)
				T2 <- as.data.frame(res_PLS_1$T2)
				Q2 <- as.data.frame(res_PLS_1$Q2)

				R2HTML::HTML(R2HTML::as.title("X-Loadings"),HR=3,file=stdout())
				R2HTML::HTML(Xload,file=stdout(), innerBorder = 1,align="left")

				R2HTML::HTML(R2HTML::as.title("Y-Loadings"),HR=3,file=stdout())
				R2HTML::HTML(Yload,file=stdout(), innerBorder = 1,align="left")	

				R2HTML::HTML(R2HTML::as.title("Score correlation"),HR=3,file=stdout())
				R2HTML::HTML(corxyt,file=stdout(), innerBorder = 1,align="left")

				R2HTML::HTML(R2HTML::as.title("Raw Weights"),HR=3,file=stdout())
				R2HTML::HTML(raw_weight,file=stdout(), innerBorder = 1,align="left")

				R2HTML::HTML(R2HTML::as.title("Modified Weights"),HR=3,file=stdout())
				R2HTML::HTML(mod_weight,file=stdout(), innerBorder = 1,align="left")	
	
				R2HTML::HTML(R2HTML::as.title("Regular Coefficients"),HR=3,file=stdout())
				R2HTML::HTML(regcoeff,file=stdout(), innerBorder = 1,align="left")

				R2HTML::HTML(R2HTML::as.title("Standard Coefficients"),HR=3,file=stdout())
				R2HTML::HTML(std,file=stdout(), innerBorder = 1,align="left")

				R2HTML::HTML(R2HTML::as.title("R-squared"),HR=3,file=stdout())
				R2HTML::HTML(R2,file=stdout(), innerBorder = 1,align="left")

				R2HTML::HTML(R2HTML::as.title("Explained Variance of X-Y by T"),HR=3,file=stdout())
				R2HTML::HTML(R2xy,file=stdout(), innerBorder = 1,align="left")

				R2HTML::HTML(R2HTML::as.title("T2 hotelling"),HR=3,file=stdout())
				R2HTML::HTML(T2,file=stdout(), innerBorder = 1,align="left")

				if(crosval==TRUE){
					R2HTML::HTML(R2HTML::as.title("Q2 cross validation"),HR=3,file=stdout())
					R2HTML::HTML(Q2,file=stdout(), innerBorder = 1,align="left")
				}


			#################### response의 열이 하나인경우 출력결과(plsreg2) ##############
			}else{			
				####### Results of Models
				Xscores <- as.data.frame(res_PLS_1$x.scores)
				Xload   <- as.data.frame(res_PLS_1$x.loads)
				Yscores <- as.data.frame(res_PLS_1$y.scores)
				rownames(Yscores) <- rownames(Xscores)
				Yload   <- as.data.frame(res_PLS_1$y.loads)

				corxt <- as.data.frame(res_PLS_1$cor.xt)
				coryt <- as.data.frame(res_PLS_1$cor.yt)
				corxu <- as.data.frame(res_PLS_1$cor.xu)
				coryu <- as.data.frame(res_PLS_1$cor.yu)
				cortu <- as.data.frame(res_PLS_1$cor.tu)			

				raw_weight <- as.data.frame(res_PLS_1$raw.wgs)
				mod_weight <- as.data.frame(res_PLS_1$mod.wgs)

				std <- as.data.frame(res_PLS_1$std.coefs)
				regcoeff <- as.data.frame(res_PLS_1$reg.coefs)	

				pred <- as.data.frame(res_PLS_1$y.pred)
				rownames(pred) <- rownames(Xscores)
				resid <- as.data.frame(res_PLS_1$resid)
				rownames(resid) <- rownames(Xscores)			
				exp_var <- as.data.frame(res_PLS_1$expvar)
				vip <- as.data.frame(res_PLS_1$VIP)
				q2 <- as.data.frame(res_PLS_1$Q2)
				q2sum <- as.data.frame(res_PLS_1$Q2cum)

				R2HTML::HTML(R2HTML::as.title("X-Loadings"),HR=3,file=stdout())
				R2HTML::HTML(Xload,file=stdout(), innerBorder = 1,align="left")

				R2HTML::HTML(R2HTML::as.title("Y-Loadings"),HR=3,file=stdout())
				R2HTML::HTML(Yload,file=stdout(), innerBorder = 1,align="left")	

				R2HTML::HTML(R2HTML::as.title("(X,T) correlation"),HR=3,file=stdout())
				R2HTML::HTML(corxt,file=stdout(), innerBorder = 1,align="left")
	
				R2HTML::HTML(R2HTML::as.title("(Y,T) correlation"),HR=3,file=stdout())
				R2HTML::HTML(coryt,file=stdout(), innerBorder = 1,align="left")	
	
				R2HTML::HTML(R2HTML::as.title("(X,U) correlation"),HR=3,file=stdout())
				R2HTML::HTML(corxu,file=stdout(), innerBorder = 1,align="left")

				R2HTML::HTML(R2HTML::as.title("(Y,U) correlation"),HR=3,file=stdout())
				R2HTML::HTML(coryu,file=stdout(), innerBorder = 1,align="left")

				R2HTML::HTML(R2HTML::as.title("(T,U) correlation"),HR=3,file=stdout())
				R2HTML::HTML(cortu,file=stdout(), innerBorder = 1,align="left")

				R2HTML::HTML(R2HTML::as.title("Raw Weights"),HR=3,file=stdout())
				R2HTML::HTML(raw_weight,file=stdout(), innerBorder = 1,align="left")

				R2HTML::HTML(R2HTML::as.title("Modified Weights"),HR=3,file=stdout())
				R2HTML::HTML(mod_weight,file=stdout(), innerBorder = 1,align="left")	
	
				R2HTML::HTML(R2HTML::as.title("Regular Coefficients"),HR=3,file=stdout())
				R2HTML::HTML(regcoeff,file=stdout(), innerBorder = 1,align="left")

				R2HTML::HTML(R2HTML::as.title("Standard Coefficients"),HR=3,file=stdout())
				R2HTML::HTML(std,file=stdout(), innerBorder = 1,align="left")

				R2HTML::HTML(R2HTML::as.title("Explained Variance"),HR=3,file=stdout())
				colnames(exp_var) <- c("X variance","Cumulative X variance","Y variance","Cummulative Y variance (R-square)")
				R2HTML::HTML(exp_var,file=stdout(), innerBorder = 1,align="left")

				R2HTML::HTML(R2HTML::as.title("Variable Importance for Projection(VIP)"),HR=3,file=stdout())
				R2HTML::HTML(vip,file=stdout(), innerBorder = 1,align="left")
				if(crosval==TRUE){
					R2HTML::HTML(R2HTML::as.title("Q2 Index"),HR=3,file=stdout())
					colnames(q2) <- c(colnames(responses),"Total")
					R2HTML::HTML(q2,file=stdout(), innerBorder = 1,align="left")
					R2HTML::HTML(R2HTML::as.title("Cummulated Q2"),HR=3,file=stdout())
					colnames(q2sum) <- c(colnames(responses),"Total")
					R2HTML::HTML(q2sum,file=stdout(), innerBorder = 1,align="left")
				}

			}

			# 다음의 값들은 엑셀 시트에 저장
			if(Xscore) {
				temp.0 <- Xscores
				temp <- as.data.frame(matrix(NA,nrow(raw.dat),ncol(Xscores)))
				colnames(temp) <- paste("Xscore_",colnames(temp.0),se="")
				row.dataset <- rownames(raw.dat)
				row.output <- rownames(temp.0)
				for(i in seq_along(row.output)) if(i%in%seq_along(row.output)) temp[row.dataset==row.output[i],] <- temp.0[row.output[i],]
				if(exists('O')) {
					O <- cbind(O,temp)
				} else {
					O <- temp
				}
			}

			if(Yscore) {
				temp.0 <- Yscores
				temp <- as.data.frame(matrix(NA,nrow(raw.dat),ncol(Yscores)))
				colnames(temp) <- paste("Yscore_",colnames(temp.0),sep="")
				row.dataset <- rownames(raw.dat)
				row.output <- rownames(temp.0)
				for(i in seq_along(row.output)) if(i%in%seq_along(row.output)) temp[row.dataset==row.output[i],] <- temp.0[row.output[i],]
				if(exists('O')) {
					O <- cbind(O,temp)
				} else {
					O <- temp
				}
			}

			if(Pred) {
				temp.0 <- pred
				temp <- as.data.frame(matrix(NA,nrow(raw.dat),ncol(pred)))
				colnames(temp) <- paste("Pred_",colnames(responses),sep="")
				row.dataset <- rownames(raw.dat)
				row.output <- rownames(temp.0)
				for(i in seq_along(row.output)) if(i%in%seq_along(row.output)) temp[row.dataset==row.output[i],] <- temp.0[row.output[i],]
				if(exists('O')) {
					O <- cbind(O,temp)
				} else {
					O <- temp
	
				}
			}
			if(Resid) {
				temp.0 <- resid
				temp <- as.data.frame(matrix(NA,nrow(raw.dat),ncol(resid)))
				colnames(temp) <- paste("Resid_",colnames(responses),sep="")
				row.dataset <- rownames(raw.dat)
				row.output <- rownames(temp.0)
				for(i in seq_along(row.output)) if(i%in%seq_along(row.output)) temp[row.dataset==row.output[i],] <- temp.0[row.output[i],]
				if(exists('O')) {
					O <- cbind(O,temp)
				} else {
					O <- temp
	
				}
			}
		}
		# Analysis end time
		R2HTML::HTMLhr(file=stdout())
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Partial Least Squares Analysis.",sep=""),file=stdout())
		R2HTML::HTMLhr(file=stdout())
	})

	if(exists('O')){
		return(list(html=html.output,Output=O))
	} else {
		return(html.output)
	}
}

REx_TSE <- function(dataset,ts_var,tsdisplay=TRUE,season=FALSE,freq=2,trans=c("notrans","Log","Sqrt","bc"),
                    lambda=NULL,diff=FALSE,ndiff=1,sdiff=FALSE,nsdiff=1,bc.lambda=F,lagp=F,l=9,
                    seasonp=FALSE,tsdiag=FALSE,uni=FALSE,extract=FALSE,varname="Trans_TSdata")
{
  ######
  ## dataset : dataset 이름
  ## ts_var : 변수선택 탭의 시계열 변수: 결측값이 없는 numeric 변수 하나필수 선택.
  ## season : 변수선택 탭의 계절형 자료 여부
  ## freq : 변수선택 탭의 계절길이
  ## trans : 변수선택 탭의 변환방법 "notrans","Log","Sqrt","bc" 4가지중 선택
  ## lambda : 변수선택 탭의 boxcox변환을 선택했을 때 lambda값
  ## diff : 변수선택 탭의 차분 선택되면 ndiffs 가 반드시 입력이 되야함
  ## ndiff : 1이상 정수
  ## sdiff : 변수선택 탭의 계절 차분 선택되면 sndiffs가 반드시 입력이 되야함
  ## sndiff : 1이상 정수
  ## tsdisplay : 시계열 acf pacf 도표
  ## seasonp : 출력옵션 탭의 계절별 도표
  ## lagp : 지연 시차 도표
  ## l : 최대 지연시차 값
  ## tsdiag : 백색잡음 검정
  ## uni : 정상성을 위한 차분 횟수 검정
  ## extract : 변환한 데이터 출력
  ## varname : 변환한 데이터 이름
  
  
  
  ######
  ## required packages : forecast, ggplot2, R2HTML, ggfortify
  load.pkg(c("forecast", "R2HTML", "ggplot2", "ggfortify"))
  
  
  html.output <- capture.output({
    ### HTML Output
    ## Title
    HTML(as.title('Time Series Data Exploration'),HR=1,file=stdout(),append=FALSE)
    ## Warnings
    if (!(is.numeric(dataset[,ts_var]))) warn.msg1 <- "\a Error : Time series variable should be numeric. Analysis has been stopped."
    if (sum(is.na(dataset[,ts_var]))!=0 ) warn.msg2 <- "\a Error : Missing values encountered in Time series variable. Analysis has been stopped."
    data.positive <- (min(dataset[,ts_var]) > 0)
    if (trans=="bc" & is.null(lambda)==F) {if(lambda<=0 & !data.positive) warn.msg3<- "\a Error : Inappropriate transformation for data with negative or zero values. Analysis has been stopped."}
    if (trans=="Log" & !data.positive) warn.msg3<- "\a Error : Inappropriate transformation for data with negative or zero values. Analysis has been stopped."
    data.positive2 <- (min(dataset[,ts_var]) >= 0)
    if (trans=="Sqrt" & !data.positive2) warn.msg3<- "\a Error : Inappropriate transformation for data with negative or zero values. Analysis has been stopped."
    if (diff){if(ndiff > dim(dataset)[1]-2) warn.msg4<-  "\a Error :  Inappropriate number of differencing is used. Analysis has been stopped."}
    if (sdiff){if(nsdiff > (dim(dataset)[1]-ndiff)/freq-2) warn.msg4<-  "\a Error :  Inappropriate number of differencing is used. Analysis has been stopped."}
    if(exists("warn.msg1")|exists("warn.msg2")|exists("warn.msg3")|exists("warn.msg4")) {
      R2HTML::HTML(R2HTML::as.title('Warnings'),HR=2,file=stdout())
      if(exists("warn.msg1")) R2HTML::HTML(warn.msg1,file=stdout())
      if(exists("warn.msg2")) R2HTML::HTML(warn.msg2,file=stdout())
      if(exists("warn.msg3")) R2HTML::HTML(warn.msg3,file=stdout())
      if(exists("warn.msg4")) R2HTML::HTML(warn.msg4,file=stdout())
    }
    else {
      
      R2HTML::HTML(as.title('Data Structure'),HR=2,file=stdout())
      DS <- matrix(c('Number of observations',dim(dataset)[1]),ncol=2,byrow=T)
      if(season){ DS<-rbind(DS,matrix(c('Seasonal period',freq),ncol=2,byrow=T))}
      if(trans!="notrans" & is.null(lambda)) {DS<-rbind(DS,matrix(c('Transformation',trans),ncol=2,byrow=T))}
      if(trans!="notrans" & !is.null(lambda)) {DS<-rbind(DS,matrix(c('Box-Cox transformation',ifelse(lambda==0,paste0('f(y)=log(y)'),paste0('f(y)=(y^(',lambda,')-1)/(',lambda,')'))),ncol=2,byrow=T))}
      
      HTML(DS,file=stdout(), innerBorder = 1,align="left")
      
      ## Analysis Description
      HTML(as.title('Analysis Description'),HR=2,file=stdout())
      AD <- matrix(c("Method","Data exploration","Variable name",ts_var),ncol=2,byrow=T)
      HTML(AD,file=stdout(), innerBorder = 1,align="left")
      
      
      t<-dim(dataset)[1]
      ## seasonal data
      if(season) {dataset[,ts_var] <- ts(dataset[,ts_var],freq=freq)}
      else {dataset[,ts_var]<-ts(dataset[,ts_var]) }
      
      ## trans
      if(trans!="notrans")
      {
        if(trans=="Log") dataset[,ts_var]<-log(dataset[,ts_var])
        if(trans=="Sqrt") dataset[,ts_var]<-sqrt(dataset[,ts_var])
        if(trans=="bc")  dataset[,ts_var]<-BoxCox(dataset[,ts_var],lambda)
      }
      
      ##diffrences
      if(diff | sdiff){
        y<-dataset[,ts_var]
        
        if(diff)
        {
          y<-diff(y,1,ndiff)
        }
        
        if(sdiff)
        {
          y<-diff(y,freq,nsdiff)
        }
        
        dataset<-as.data.frame(y)
        colnames(dataset)<-ts_var
      }
      
      #time series display
      if(tsdisplay)
      {
        R2HTML::HTML(R2HTML::as.title(paste0('Time series Plot: ',ts_var)),HR=3,file=stdout(),append=TRUE)
        REx_ANA_PLOT()
        print(autoplot(dataset[,ts_var],main=paste0("Series: ",ts_var),ylab=ts_var))
        REx_ANA_PLOT_OFF("")
        REx_ANA_PLOT()
        print(ggAcf(dataset[,ts_var],main=NULL))
        REx_ANA_PLOT_OFF("")
        REx_ANA_PLOT()
        print(ggPacf(dataset[,ts_var],main=NULL))
        REx_ANA_PLOT_OFF("")
        
      }
      
      ##seasonal plot
      if(seasonp)
      {
        R2HTML::HTML(R2HTML::as.title('Seasonal Plot'),HR=3,file=stdout(),append=TRUE)
        REx_ANA_PLOT()
        print(ggseasonplot(dataset[,ts_var],main=paste0("Series: ",ts_var),ylab=NULL))
        REx_ANA_PLOT_OFF("")
        R2HTML::HTML(R2HTML::as.title('Monthly Plot'),HR=3,file=stdout(),append=TRUE)
        REx_ANA_PLOT()
        print(ggmonthplot(dataset[,ts_var],main=paste0("Series: ",ts_var),ylab=NULL))
        REx_ANA_PLOT_OFF("")
      }
      
      
      if(lagp)
      {
        R2HTML::HTML(R2HTML::as.title('Lag Plot'),HR=3,file=stdout(),append=TRUE)
        REx_ANA_PLOT()
        print(forecast:::gglagplot(dataset[,ts_var],main=paste0("Series: ",ts_var),l))
        REx_ANA_PLOT_OFF("")
        
      }
      
      #white noise test
      if(tsdiag)
      {
        REx_ANA_PLOT()
        ##
        conf.int = TRUE; conf.int.colour = "#0000FF"; 
        conf.int.linetype = "dashed"; conf.int.fill = NULL; conf.int.alpha = 0.3; 
        ad.colour = "#888888"; ad.linetype = "dashed"; ad.size = 0.2;
        gof.lag<-10
        nlag <- gof.lag
        pval <- numeric(nlag)
        for (i in 1L:nlag) pval[i] <- max(0.0001,stats::Box.test(dataset[,ts_var], i, type = "Ljung-Box")$p.value)
        lb.df <- data.frame(Lag = 1L:nlag, `p value` = pval, lower = -0.05, 
                            upper = 0.05)
        colnames(lb.df) <- c("Lag", "p value", "lower", "upper")
        p.lb <- ggplot2::ggplot(data = lb.df, mapping = ggplot2::aes_string(x = "Lag")) + 
          ggplot2::geom_point(mapping = ggplot2::aes_string(y = "`p value`")) + 
          ggplot2::scale_y_continuous(limits = c(0, 1)) + ggplot2::ggtitle("p values for Ljung-Box statistics") +
                                                         theme(plot.title = element_text(size = 20))
        p.lb <- ggfortify:::plot_confint(p = p.lb, data = lb.df, conf.int = conf.int, 
                                         conf.int.colour = conf.int.colour, conf.int.linetype = conf.int.linetype, 
                                         conf.int.fill = conf.int.fill, conf.int.alpha = conf.int.alpha)
        
        
        print(p.lb)
        ##
        REx_ANA_PLOT_OFF("Blue dashed line means alpha=0.05")
        
      }
      
      #boxcox lambda
      if(bc.lambda)
      {
        HTML(R2HTML::as.title('Estimate Box-Cox lambda'),HR=3,file=stdout(),append=TRUE)
        VL <- matrix(c("Lambda",round(BoxCox.lambda(dataset[,ts_var]),4)),nrow=2,ncol=1)
        HTML(VL,file=stdout(), innerBorder = 1,align="left")
      }
      
      #uniroot test
      if(uni)
      {
        HTML(R2HTML::as.title('Stationarity Test'),HR=3,file=stdout(),append=TRUE)
        VL <- matrix(c("N. of differences estimated",ndiffs(dataset[,ts_var])),nrow=2,ncol=1)
        if(season) VL<- rbind(VL,matrix(c("N. of seasonal differences estimated",nsdiffs(dataset[,ts_var])),nrow=2,ncol=1))
        HTML(VL,file=stdout(), innerBorder = 1,align="left")
      }
      
      ##add transformation variable
      if(extract)
      {
        
        O<-as.data.frame(c(dataset[,ts_var],rep(NA,t-dim(dataset)[1])))
        colnames(O)<-varname
      }
      
    }  
    
    
    ### TIME 
    R2HTML::HTMLhr(file=stdout())
    R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Time Series Data Exploration",sep=""),file=stdout())
    R2HTML::HTMLhr(file=stdout())  
  })
  
  if(exists('O')){
    return(list(html=html.output,Output=O))
  } else {
    return(html.output)
  }
}

REx_ARIMA <- function(dataset,ts_var,ind_var=NULL,season=F,freq=NULL,trans=F,lambda=NULL,i.m=F,i.d=F,auto=F,p=0,d=0,q=0,P=0,
                      D=0,Q=0,timeplot=T,diag=T,pred=F,pred.plot=F,h=NULL,conf=0.95,fitted=F,resid=F)
{
  
  ######
  ## dataset : dataset 이름
  ## ts_var : 변수선택 탭의 시계열 변수: 결측값이 없는 numeric 변수 하나필수 선택.
  ## ind_var : 변수선택 탭의 독립 변수 : 결측값이 없는 numeric 변수 (필수아님)
  ## season : 변수선택 탭의 계절형 자료 : 체크되면 T
  ## freq : 변수선택 탭의 계절길이
  ## trans : 분석옵션 탭의 Box-cox 변환 : 체크되면 T
  ## lambda : 분석옵션 탭의 lambda값 
  ## p,d,q : 분석옵션 탭의 비계절 차수 차례대로 (AR, 차분, MA)
  ## P,D,Q : 분석옵션 탭의 계절 차수 차례대로(AR, 차분, MA)
  ## auto : 분석옵션 탭의 비계절 및 계절 차수 자동선택 : 체크되면 T
  ## i.m : 분석옵션 탭의 절편 포함 : 체크되면 T
  ## i.d : 분석옵션 탭의 1차 추세 포함 : 체크되면 T
  ## timeplot : 출력옵션 탭의 시계열 도표 : 체크되면 T
  ## diag : 출력옵션 탭의 잔차진단 도표 : 체크되면 T
  ## pred : 출력옵션 탭의 예측값 : 체크되면 T
  ## pred.plot : 출력옵션 탭의 예측 도표 체크되면 T
  ## h : 출력옵션 탭의 예측시차 (1이상의 정수만 입력 가능)
  ## conf : 출력옵션 탭의 신뢰구간 ( 0부터 1사이의 값)
  ## fitted : 출력옵션 탭의 적합값 : 체크되면 T
  ## resid :  출력옵션 탭의 잔차 : 체크되면 T
  
  
  ######
  ## required packages : forecast, R2HTML
  
  load.pkg(c("forecast", "R2HTML"))
  html.output <- capture.output({
    ### HTML Output
    ## Title
    HTML(as.title('ARIMA Model'),HR=1,file=stdout(),append=FALSE)
    if(!trans) lambda<-NULL
    if(season) {dataset[,ts_var]<-ts(dataset[,ts_var],freq=freq)}
    
    
    
    if(auto==T & is.null(ind_var)==T)  model<-tryCatch(auto.arima(dataset[,ts_var],lambda=lambda),error=function(e) e)
    if(auto==F & is.null(ind_var)==T)  model<-tryCatch(Arima(dataset[,ts_var],order=c(p,d,q),seasonal=c(P,D,Q),include.mean=i.m,include.drift=i.d,lambda=lambda),error=function(e) e)
    if(auto==T & is.null(ind_var)==F)  model<-tryCatch(auto.arima(dataset[,ts_var],lambda=lambda,xreg=dataset[,ind_var]),error=function(e) e)
    if(auto==F & is.null(ind_var)==F)  model<-tryCatch(Arima(dataset[,ts_var],order=c(p,d,q),seasonal=c(P,D,Q),include.mean=i.m,include.drift=i.d,lambda=lambda,xreg=dataset[,ind_var]),error=function(e) e)
    
    
    
    if (!(is.numeric(dataset[,ts_var]))) warn.msg1 <- "\a Error : Time series variable should be numeric. Analysis has been stopped."
    if (sum(is.na(dataset[,ts_var]))!=0 ) warn.msg2 <- "\a Error : Missing values encountered in Time series variable. Analysis has been stopped."
    data.positive <- (min(dataset[,ts_var]) > 0)
    if(trans)
    {
      if(data.positive==F & lambda<=0) warn.msg3<- "\a Error : Inappropriate transformation for data with negative or zero values. Analysis has been stopped."
    }
    
    if(is.null(model$message)==F | exists("warn.msg1")|exists("warn.msg2")|exists("warn.msg3")) {
      
      R2HTML::HTML(R2HTML::as.title('Warnings'),HR=2,file=stdout())
      if(is.null(model$message)==F) R2HTML::HTML(paste("\a",model$message),file=stdout())
      if(exists("warn.msg1")) R2HTML::HTML(warn.msg1,file=stdout())
      if(exists("warn.msg2")) R2HTML::HTML(warn.msg2,file=stdout())
      if(exists("warn.msg3")) R2HTML::HTML(warn.msg3,file=stdout())
      
    } else 
    {
      
      
      ## Data Structure
      R2HTML::HTML(as.title('Data Structure'),HR=2,file=stdout())
      DS <- matrix(c('Number of observations',dim(dataset)[1]),ncol=2,byrow=T)
      if(season) DS<-rbind(DS,c('Seasonal period',freq))
      HTML(DS,file=stdout(), innerBorder = 1,align="left")
      
      ## Analysis Description
      HTML(as.title('Analysis Description'),HR=2,file=stdout())
      p<-model$arma[1];d<-model$arma[6];q<-model$arma[2];P<-model$arma[3];D<-model$arma[7];Q<-model$arma[4]
      if(season==T & is.null(ind_var)==T) {method<-paste0("ARIMA(",p,",",d,",",q,")(",P,",",D,",",Q,")[",freq,"]")}
      if(season==F & is.null(ind_var)==T) {method<-paste0("ARIMA(",p,",",d,",",q,")")}
      if(season==T & is.null(ind_var)==F) {method<-paste0("Regression with ARIMA(",p,",",d,",",q,")(",P,",",D,",",Q,")[",freq,"] errors")}
      if(season==F & is.null(ind_var)==F) {method<-paste0("Regression with ARIMA(",p,",",d,",",q,") errors")}
      
      AD <- matrix(c("Method",method,"Variable name",ts_var),ncol=2,byrow=T)
      if(trans) AD <- rbind(AD,matrix(c('Box-Cox transformation',ifelse(lambda==0,paste0('f(y)=log(y)'),paste0('f(y)=(y^(',lambda,')-1)/(',lambda,')'))),ncol=2,byrow=T))
      R2HTML::HTML(AD,HR=2,file=stdout(),append=TRUE,innerBorder = 1,align="left")
      
      ### Default output	
      
      HTML(R2HTML::as.title('Estimates of ARIMA Parameters'),HR=3,file=stdout(),append=TRUE) 	
      
      VL <- matrix(0,nrow=1,ncol=3)
      if(sum(p,q,P,Q)!=0){
        for(i in 1:length(model$coef))
        {
          VL <- rbind(VL,c(attributes(model$coef)$names[i],round(model$coef[i], 4),round(sqrt(model$var.coef[i,i]),4)))
        }
        VL <- VL[-1,,drop=F]
      }
      colnames(VL) <- c('Parameter','Estimate','s.e')
      HTML(data.frame(VL),file=stdout(), innerBorder = 1,align="left",row.names=F)
      
      R2HTML::HTML(R2HTML::as.title('Goodness of Fit Measures'),HR=3,file=stdout(),append=TRUE)
      fit1<-matrix(round(c(model$aic,model$aicc,model$bic),4),ncol=3)
      colnames(fit1)<-c("AIC","AICc","BIC"); rownames(fit1) <- 'Criterion'
      HTML(fit1,file=stdout(), innerBorder = 1,align="left")
      
      R2HTML::HTML(R2HTML::as.title('Training Set Accuracy Measures'),HR=3,file=stdout(),append=TRUE)
      fit2<-round(accuracy(model),4)
      rownames(fit2)<-'Accuracy'
      HTML(fit2,file=stdout(), innerBorder = 1,align="left")
      
      ### Option_plot : Time-series Plot
      if(timeplot) 
      {
        R2HTML::HTML(R2HTML::as.title('Time Series Plot'),HR=3,file=stdout(),append=TRUE)
        REx_ANA_PLOT()
        plot(dataset[,ts_var],ylab=ts_var,xlab="Time",type="l")
        REx_ANA_PLOT_OFF("")
      }
      
      ### Option_plot : Diagnostic Plot
      if(diag)
      {
        R2HTML::HTML(R2HTML::as.title('Diagnostic Plots'),HR=3,file=stdout(),append=TRUE)
        REx_ANA_PLOT()
        tsdiag(model)
        REx_ANA_PLOT_OFF("")
      }
      
      ### Option : Prediction
      if(pred==T & is.null(ind_var)==T){
        step<-as.character(1:h)
        fp<-forecast(model,h=h,level=conf)
        f.mean<-fp$mean;f.lo<-fp$lower;f.up<-fp$upper
        f.table<-matrix(c(round(c(f.mean,f.lo,f.up),4)),ncol=3)
        f.table<-data.frame(step,f.table)
        colnames(f.table)<-c("Step","Point Forecast",paste("Lo",conf*100),paste("Hi",conf*100))
        R2HTML::HTML(R2HTML::as.title('Forecast Intervals'),HR=3,file=stdout(),append=TRUE)
        HTML(f.table,file=stdout(), innerBorder = 1,align="left")
      }
      
      if(pred==T & is.null(ind_var)==F){
        newdata<-NULL
        for(i in 1:length(ind_var)){
          xmodel<-auto.arima(dataset[,ind_var[i]])
          predx<-forecast(xmodel,h=h)$mean
          newdata<-cbind(newdata,predx)
        }
        newdata<-as.data.frame(newdata)
        colnames(newdata)<-ind_var
        step<-as.character(1:h)
        fp<-forecast(model,h=h,level=conf,xreg=newdata)
        f.mean<-fp$mean;f.lo<-fp$lower;f.up<-fp$upper
        f.table<-matrix(c(round(c(f.mean,f.lo,f.up),4)),ncol=3)
        f.table<-data.frame(step,f.table)
        colnames(f.table)<-c("Step","Point Forecast",paste("Lo",conf*100),paste("Hi",conf*100))
        R2HTML::HTML(R2HTML::as.title('Forecast Intervals'),HR=3,file=stdout(),append=TRUE)
        HTML(f.table,file=stdout(), innerBorder = 1,align="left",row.names=F)
      }
      
      ### Option_plot : Forecast Plot
      if(pred.plot)
      {
        R2HTML::HTML(R2HTML::as.title('Forecast Plot'),HR=3,file=stdout(),append=TRUE)
        REx_ANA_PLOT()
        plot(fp,ylab=ts_var,xlab="Time")
        REx_ANA_PLOT_OFF("")
      }
      
      
      if(fitted | resid)
      {
        O<-data.frame(x=1:length(dataset[,ts_var]))
        if(fitted) O<-data.frame(O,Fitted_ets=as.matrix(model$fitted))
        if(resid) O<-data.frame(O,Resid_ets=as.matrix(model$residuals))
        O<-O[,-1,drop=F]
      }	
    }
    
    
    ### TIME 
    R2HTML::HTMLhr(file=stdout())
    R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : ARIMA",sep=""),file=stdout())
    R2HTML::HTMLhr(file=stdout())
    
  })
  
  if(exists('O')){
    return(list(html=html.output,Output=O))
  } else {
    return(html.output)
  }
}

REx_GARCH<-function(dataset,ts_var,ind_var=NULL,method="sGARCH",arma.p,arma.q,garch.p,garch.q,
                    arfima=F,archm=F,archpow=1,inc.mean=T,dist="norm",robust=F,nyblom=F,
                    sign=F,goodness=F,tsplot=F,acfplot=F,standplot=F,NIC=F,pred=F,h=NULL,predplot=F,
                    boot=F,resid=F,Sresid=F,condSD=F)
{
  ######
  ## dataset : dataset 이름
  ## ts_var : 변수선택 탭의 시계열 변수: 결측값이 없는 numeric 변수 하나필수 선택.
  ## ind_var : 변수선택 탭의 독립 변수 : 결측값이 없는 numeric 변수 (필수아님)
  ## arma.p : 분석옵션 탭의 평균모형 ARMA Order p (0이상의 정수)
  ## arma.q : 분석옵션 탭의 평균모형 ARMA Order q (0이상의 정수)
  ## arfima : fractional differencing이 체크, 아니면 F
  ## inc.mean : 절편포함이 체크,T 아니면 F
  ## archm : ARCH-in-mean 포함 체크,아니면 F
  ## archpow : ARCH-in-mean 포함이 체크, SD가 선택되면 1 VAR이 선택되면 2
  ## 
  ## method : 분석옵션 탭의 분산모형 모형의 종류.
  ##            standard GARCH : "sGARCH"
  ##            integrated GARCH : "iGARCH"
  ##            exponential GARCH : "eGARCH"
  ##            GJR GARCH : "gjrGARCH"
  ##            Threshold GARCH : "fGARCH"
  ##            asymmetric power ARCH : "apARCH"
  ##            component standard GARCH : "csGARCH"
  ## garch.p : 분석옵션 탭의 분산모형 GARCH p  (0이상의 정수)
  ## garch.q : 분석옵션 탭의 분산모형 GARCH q (0이상의 정수)
  ## dist : 분석옵션 탭의 오차의 분포
  ##          normal : "norm"
  ##          student-t : "std"
  ##          generalized error : "ged"
  ##          skewed normal : "snorm"
  ##          skewed student : "sstd"
  ##          skewed generalized error : "sged"
  ##          normal inverse gaussian : "nig"
  ##   tsplot : 출력옵션 탭의 진단그림 Original time series plots
  ##   acfplot : 출력옵션 탭의 진단그림 Original time series ACFs
  ##   standplot : 출력옵션 탭의 진단그림 Standardized residual plot
  ##   NIC : 출력옵션 탭의 진단그림 News Impact curve
  ##   pred : 출력옵션 탭의 예측 체크 되면 T 아니면 F
  ##   h : 출력옵션 탭의 예측 예측시차 (이상의 정수)
  ##   predplot 출력옵션 탭의 예측 예측도표 체크되면 T 아니면 F
  ##   boot : 출력옵션 탭의 예측 GARCH Bootstrap 체크되면 T 아니면 F
  ##   resid : 출력옵션 탭의 저장 Residuals 체크되면 T 아니면 F
  ##   Sresid : 출력옵션 탭의 저장 Standardized Residuals 체크되면 T 아니면 F
  ##   condSD : 출력옵션 탭의 저장 Conditional SD(volatility) 체크되면 T 아니면 F
  ##   robust : 출력옵션 탭의 검정 Robust Standard Errors 체크되면 T 아니면 F
  ##   nyblom : 출력옵션 탭의 검정 Nyblom stability test 체크되면 T 아니면 F
  ##   sign : 출력옵션 탭의 검정 Sign Bias Test 체크되면 T 아니면 F
  ##   goodness : 출력옵션 탭의 검정 Adjusted Pearson Goodness of Fit test 체크되면 T 아니면 F
  
  ######
  ## required packages : rugarch, R2HTML
  load.pkg(c("rugarch", "R2HTML"))
  options("scipen"=999,"digits"=4)
  
  if(is.null(ind_var)==F) {X<-as.matrix(dataset[,ind_var]);colnames(X)<-ind_var;} 
  if(is.null(ind_var)==T) {X<-NULL}
  if(method=="fGARCH") {submodel<-"TGARCH"} else{submodel<-NULL}
  
  html.output <- capture.output({
    
    ### HTML Output
    ## Title
    HTML(as.title('GARCH Model'),HR=1,file=stdout(),append=FALSE)
    
    model <- tryCatch(ugarchspec(variance.model = list(model = method, garchOrder = c(garch.p, garch.q), 
                                                       submodel = submodel, external.regressors = NULL, variance.targeting = FALSE), 
                                 mean.model = list(armaOrder = c(arma.p, arma.q), include.mean = inc.mean, archm = archm, 
                                                   archpow = archpow, arfima = arfima, external.regressors = X, archex = FALSE), 
                                 distribution.model = dist),error=function(e) e)
    
    fit<-tryCatch(ugarchfit(data=dataset[,ts_var], model),error=function(e) e)
    
    if (!(is.numeric(dataset[,ts_var]))) warn.msg1 <- "\a Error : Time series variable should be numeric. Analysis has been stopped."
    if (sum(is.na(dataset[,ts_var]))!=0 ) warn.msg2 <- "\a Error : Missing values encountered in Time series variable. Analysis has been stopped."
    if (sum(is.na(dataset[,ind_var]))!=0 ) warn.msg8 <- "\a Error : Missing values encountered in Independent variable. Analysis has been stopped."
    if (!all(sapply(dataset[,ind_var], is.numeric))) warn.msg3 <- "\a Error : Independent variable variable should be numeric. Analysis has been stopped."
    if (class(model)[1]!="uGARCHspec") warn.msg4<-paste("\a ",model$message)
    if (class(fit)[1]=="uGARCHfit"){if(fit@fit$convergence!=0) warn.msg5<-paste("\a Error : Convergence Problem <br> \a Solver Message : ", fit@fit$message)}
    if (class(fit)[1]!="uGARCHfit") warn.msg6<-paste("\a Error : ",fit$message)
    if (garch.p+garch.q==0) warn.msg7 <- "\a Error : GARCH p+q = 0 Analysis has been stopped."
    
    if(exists("warn.msg1")|exists("warn.msg2")|exists("warn.msg3")|exists("warn.msg4")|
       exists("warn.msg5")|exists("warn.msg6")|exists("warn.msg7")|exists("warn.msg8"))
    {
      
      R2HTML::HTML(R2HTML::as.title('Warnings'),HR=2,file=stdout())
      if(exists("warn.msg1")) R2HTML::HTML(warn.msg1,file=stdout())
      if(exists("warn.msg2")) R2HTML::HTML(warn.msg2,file=stdout())
      if(exists("warn.msg3")) R2HTML::HTML(warn.msg3,file=stdout())
      if(exists("warn.msg4")) R2HTML::HTML(warn.msg4,file=stdout())
      if(exists("warn.msg5")) R2HTML::HTML(warn.msg5,file=stdout())
      if(exists("warn.msg6")) R2HTML::HTML(warn.msg6,file=stdout())
      if(exists("warn.msg7")) R2HTML::HTML(warn.msg7,file=stdout())
      if(exists("warn.msg8")) R2HTML::HTML(warn.msg8,file=stdout())
    } else {
      
      vmodel = fit@model$modeldesc$vmodel
      model = fit@model
      modelinc = fit@model$modelinc
      stdresid = fit@fit$residuals/fit@fit$sigma
      
      ##plot1
      
      .plot.garchfit.1 = function(x, ...)
      {
        vmodel  = x@model$modeldesc$vmodel
        T = x@model$modeldata$T
        insample = 1:T
        xdates  = x@model$modeldata$index[insample]
        xseries = x@model$modeldata$data[insample]
        xsigma  = x@fit$sigma
        ci = 2
        plot(xdates, xseries, type = "l",  col = "steelblue", ylab = "Returns", xlab="Time",
             main = "Series with 2 Conditional SD Superimposed",
             cex.main = 1.6,cex.axis = 1, cex.lab = 0.8,cex=2)
        lines(xdates, +ci* xsigma , col = "tomato1")
        
        lines(xdates, -ci* xsigma, col = "tomato1")
        mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray",cex=1)
        if(vmodel == "fGARCH"){
          mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray",cex=1)
        }
        abline(h = 0, col = "grey", lty = 3)
        grid()
      }
      
      ## plot2
      .plot.garchfit.2 = function(x, ...)
      {
        vmodel  = x@model$modeldesc$vmodel
        T = x@model$modeldata$T
        insample = 1:T
        xseries = x@model$modeldata$data[insample]
        xdates  = x@model$modeldata$index[insample]
        xsigma 	= x@fit$sigma
        distribution = x@model$modeldesc$distribution
        xcmu = fitted(x)
        idx = x@model$pidx
        pars  = x@fit$ipars[,1]
        skew  = pars[idx["skew",1]]
        shape = pars[idx["shape",1]]
        if(distribution == "ghst") ghlambda = -shape/2 else ghlambda = pars[idx["ghlambda",1]]
        z1 	= 0.01
        z2 	= 0.99
        
        q01 	= fitted(x) + sigma(x)* qdist(distribution, z1, 0, 1, lambda = ghlambda, skew, shape)
        q99 	= fitted(x) + sigma(x)* qdist(distribution, z2, 0, 1, lambda = ghlambda, skew, shape)
        plot(xdates, xseries, type = "l", col = "steelblue", ylab = "Returns", xlab="Time", 
             main = "Series with with 1% VaR Limits", 
             cex.main = 1.6,cex.axis = 1, cex.lab = 0.8)
        lines(xdates, q01, col = "tomato1")
        lines(xdates, q99, col = "green")
        mtext(paste("GARCH model :", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 1)
        if(vmodel == "fGARCH"){
          mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 1)
        }
        abline(h = 0, col = "grey", lty = 3)
        grid()
      }
      
      ##plot3
      .plot.garchfit.3 = function(x, ...)
      {
        vmodel  = x@model$modeldesc$vmodel
        T = x@model$modeldata$T
        insample = 1:T
        xseries = abs(x@model$modeldata$data[insample])
        xdates  = x@model$modeldata$index[insample]
        xsigma 	= x@fit$sigma
        plot(xdates, xseries, type = "l", col = "lightgrey", ylab = "Volatility", 
             xlab="Time", main = "Conditional SD (vs |returns|)", 
             cex.main = 1.6,cex.axis = 1, cex.lab = 0.8)
        lines(xdates, xsigma, col = "steelblue")
        mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 1)
        if(vmodel == "fGARCH"){
          mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 1)
        }
        abline(h = 0, col = "grey", lty = 3)
        grid()
      }
      
      ##plot4
      .plot.garchfit.4 = function(x, ...)
      {
        vmodel  = x@model$modeldesc$vmodel
        T = x@model$modeldata$T
        insample = 1:T
        xseries = x@model$modeldata$data[insample]
        lag.max = as.integer(10*log10(T))
        acfx 	= acf(xseries, lag.max = lag.max, plot = FALSE)
        clim0	= qnorm((1 + 0.95)/2)/sqrt(acfx$n.used)
        ylim 	= range(c(-clim0, clim0, as.numeric(acfx$acf)[-1]))
        clx 	= vector(mode = "character", length = lag.max)
        clx[which(as.numeric(acfx$acf)[-1]>=0)] = "steelblue"
        clx[which(as.numeric(acfx$acf)[-1]<0)] = "orange"
        barplot(height = as.numeric(acfx$acf)[-1], names.arg = as.numeric(acfx$lag)[-1], ylim = 1.2*ylim, col = clx,
                ylab = "ACF", xlab="lag", main = "ACF of Observations", 
                cex.main = 1.6,cex.axis = 1, cex.lab = 0.8)
        abline(h = c(clim0, -clim0), col = "tomato1", lty = 2)
        abline(h = 0, col = "black", lty = 1)
        box()
        mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray")
        if(vmodel == "fGARCH"){
          mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray")
        }
        grid()
      }
      
      ##plot5
      .plot.garchfit.5 = function(x, ...)
      {
        vmodel  = x@model$modeldesc$vmodel
        T = x@model$modeldata$T
        insample = 1:T
        xseries = x@model$modeldata$data[insample]
        lag.max = as.integer(10*log10(T))
        acfx 	= acf(xseries^2, lag.max = lag.max, plot = FALSE)
        clim0	= qnorm((1 + 0.95)/2)/sqrt(acfx$n.used)
        ylim 	= range(c(-clim0, clim0, as.numeric(acfx$acf)[-1]))
        clx 	= vector(mode="character",length=lag.max)
        clx[which(as.numeric(acfx$acf)[-1]>=0)] = "steelblue"
        clx[which(as.numeric(acfx$acf)[-1]<0)] = "orange"
        barplot(height = as.numeric(acfx$acf)[-1], names.arg = as.numeric(acfx$lag)[-1], ylim = 1.2*ylim, col = clx,
                ylab = "ACF", xlab="lag", main = "ACF of Squared Observations", 
                cex.main = 1.6,cex.axis = 1, cex.lab = 0.8)
        abline(h = c(clim0, -clim0), col = "tomato1", lty = 2)
        abline(h = 0, col = "black", lty = 1)
        box()
        mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 1)
        if(vmodel == "fGARCH"){
          mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 1)
        }
        grid()
      }
      
      # ACF of absolute observations
      .plot.garchfit.6 = function(x, ...)
      {
        vmodel  = x@model$modeldesc$vmodel
        T = x@model$modeldata$T
        insample = 1:T
        xseries = x@model$modeldata$data[insample]
        lag.max = as.integer(10*log10(T))
        acfx 	= acf(abs(xseries), lag.max = lag.max, plot = FALSE)
        clim0	= qnorm((1 + 0.95)/2)/sqrt(acfx$n.used)
        ylim 	= range(c(-clim0, clim0, as.numeric(acfx$acf)[-1]))
        clx 	= vector(mode="character",length=lag.max)
        clx[which(as.numeric(acfx$acf)[-1]>=0)] = "steelblue"
        clx[which(as.numeric(acfx$acf)[-1]<0)] = "orange"
        barplot(height=as.numeric(acfx$acf)[-1], names.arg = as.numeric(acfx$lag)[-1], ylim = 1.2*ylim, col = clx,
                ylab = "ACF", xlab = "lag", main = "ACF of Absolute Observations",
                cex.main = 1.6,cex.axis = 1, cex.lab = 0.8)
        abline(h = c(clim0, -clim0), col = "tomato1", lty = 2)
        abline(h = 0, col = "black", lty = 1)
        box()
        mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 1)
        if(vmodel == "fGARCH"){
          mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 1)
        }
        grid()
      }
      
      # Cross Correlations
      .plot.garchfit.7 = function(x, ...)
      {
        vmodel  = x@model$modeldesc$vmodel
        T = x@model$modeldata$T
        insample = 1:T
        xseries = x@model$modeldata$data[insample]	
        lag.max = as.integer(10*log10(T))
        ccfx	= ccf(xseries^2, xseries, lag.max = lag.max, plot = FALSE)
        clim0	= qnorm((1 + 0.95)/2)/sqrt(ccfx$n.used)
        ylim 	= range(c(-clim0, clim0, as.numeric(ccfx$acf)))
        clx 	= vector(mode="character",length=(2*lag.max)+1)
        clx[which(as.numeric(ccfx$acf)>=0)]="steelblue"
        clx[which(as.numeric(ccfx$acf)<0)]="orange"
        barplot(height=as.numeric(ccfx$acf), names.arg=as.numeric(ccfx$lag),ylim=1.2*ylim,col=clx,
                ylab = "ACF", xlab="lag", main = "Cross-Correlations of \nSquared vs Actual Observations", 
                cex.main = 1.6,cex.axis = 1, cex.lab = 0.8)
        abline(h = c(clim0, -clim0), col = "tomato1", lty = 2)
        abline(h = 0, col = "black")
        box()
        mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 1)
        if(vmodel == "fGARCH"){
          mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 1)
        }
        grid()
      }
      
      # Standardized Residuals Empirical Denisty
      .plot.garchfit.8 = function(x, ...)
      {
        vmodel  = x@model$modeldesc$vmodel
        zseries = as.numeric(residuals(x, standardize=TRUE))
        distribution = x@model$modeldesc$distribution
        idx = x@model$pidx
        pars  = x@fit$ipars[,1]
        skew  = pars[idx["skew",1]]
        shape = pars[idx["shape",1]]
        if(distribution == "ghst") ghlambda = -shape/2 else ghlambda = pars[idx["ghlambda",1]]	
        xmean 	= mean(zseries)
        xmedian = median(zseries)
        xsd 	= sd(zseries)
        xlim 	= c(min(zseries), max(zseries))
        result 	= hist(x = zseries, col = "grey", border = "white",
                       breaks = "Scott", main = "Empirical Density of Standardized Residuals", xlim = xlim, ylim = c(0,0.6),
                       probability = TRUE, ylab="Probability",
                       cex.main = 1.6,cex.axis = 1, cex.lab = 0.8, ...)
        box()
        
        s = result$breaks
        y	= ddist(distribution, s, lambda = ghlambda, skew = skew, shape = shape)
        lines(s, dnorm(s, 0, 1), lwd = 2, col = "blue")
        lines(s, y, lwd = 2, col = "orange")
        abline(v = xmean, lwd = 2, col = "red")
        abline(v = xmedian, lwd = 2, col = "darkgreen")
        Text = paste("Median: ", round(xmedian, 2), "| Mean: ",signif(xmean, 3))
        mtext(Text, side = 3, adj = 0, col = "darkgrey", cex = 1.4)
        mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 1)
        if(vmodel == "fGARCH"){
          mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj = 1.5, col = "gray", cex = 1)
        }
        abline(h = 0, col = "grey")
        lg.txt = c("normal Density", paste(distribution," (0,1) Fitted Density",sep=""))
        legend("topleft", legend = lg.txt, col = c("blue","orange"), pch = 1, cex = 1.6, bty = "n")
        grid()
      }
      
      .qqDist<-function (y, dist = "norm", ylim = NULL, main = paste(dist, "- QQ Plot"), 
                         xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", 
                         doplot = TRUE, datax = FALSE, cex.main = 1.6, ...) 
      {
        y = as.vector(y)
        if (has.na <- any(ina <- is.na(y))) {
          yN = y
          y = y[!ina]
        }
        if (0 == (n <- length(y))) 
          stop("y is empty or has only NAs")
        x = qdist(distribution = dist, p = ppoints(n), ...)[order(order(y))]
        if (has.na) {
          y = x
          x = yN
          x[!ina] = y
          y = yN
        }
        if (doplot) {
          if (is.null(ylim)) 
            ylim = range(y)
          if (datax) {
            plot(y, x, main = main, xlab = ylab, ylab = xlab, 
                 xlim = ylim, col = "steelblue", cex = 0.7,
                 cex.main = 1.6,cex.axis = 1, cex.lab = 0.8)
          }
          else {
            plot(x, y, main = main, xlab = xlab, ylab = ylab, 
                 ylim = ylim, col = "steelblue", cex = 0.7,
                 cex.main = 1.6,cex.axis = 1, cex.lab = 0.8)
          }
          rugarch:::.qqLine(y = y, dist = dist, datax = datax, ...)
          grid()
        }
        invisible(if (datax) list(x = y, y = x) else list(x = x, 
                                                          y = y))
      }
      
      # QQ-Plot of Standardized Residuals
      .plot.garchfit.9 = function(x, ...)
      {
        vmodel  = x@model$modeldesc$vmodel
        zseries = as.numeric(residuals(x, standardize=TRUE))
        distribution = x@model$modeldesc$distribution
        idx = x@model$pidx
        pars  = x@fit$ipars[,1]
        skew  = pars[idx["skew",1]]
        shape = pars[idx["shape",1]]
        if(distribution == "ghst") ghlambda = -shape/2 else ghlambda = pars[idx["ghlambda",1]]
        .qqDist(y = zseries, dist = distribution, lambda = ghlambda, skew = skew, shape = shape)
        mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 1)
        if(vmodel == "fGARCH"){
          mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj = 1.5, col = "gray", cex = 1)
        }
      }
      
      # ACF of standardized residuals
      .plot.garchfit.10 = function(x, ...)
      {
        vmodel  = x@model$modeldesc$vmodel
        zseries = as.numeric(residuals(x, standardize=TRUE))
        zseries[is.na(zseries)] = 0
        n 		= length(zseries)
        lag.max = as.integer(10*log10(n))
        acfx 	= acf(zseries, lag.max = lag.max, plot = FALSE)
        clim0	= qnorm((1 + 0.95)/2)/sqrt(acfx$n.used)
        ylim 	= range(c(-clim0, clim0, as.numeric(acfx$acf)[-1]))
        clx 	= vector(mode = "character", length = lag.max)
        clx[which(as.numeric(acfx$acf)[-1]>=0)] = "steelblue"
        clx[which(as.numeric(acfx$acf)[-1]<0)] = "orange"
        barplot(height = as.numeric(acfx$acf)[-1], names.arg = as.numeric(acfx$lag)[-1], ylim=1.2*ylim, col = clx,
                ylab = "ACF", xlab="lag", main = "ACF of Standardized Residuals",
                cex.main = 1.6,cex.axis = 1, cex.lab = 0.8)
        abline(h = c(clim0, -clim0), col = "tomato1", lty = 2)
        abline(h = 0, col = "black", lty = 1)
        box()
        mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 1)
        if(vmodel == "fGARCH"){
          mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj = 1.5, col = "gray", cex = 1)
        }
        grid()
      }
      
      # ACF of squared standardized residuals
      .plot.garchfit.11 = function(x, ...)
      {
        vmodel  = x@model$modeldesc$vmodel
        zseries = as.numeric(residuals(x, standardize=TRUE))
        zseries[is.na(zseries)] = 0
        n 		= length(zseries)
        lag.max = as.integer(10*log10(n))
        acfx 	= acf(zseries^2, lag.max = lag.max, plot = FALSE)
        clim0	= qnorm((1 + 0.95)/2)/sqrt(acfx$n.used)
        ylim 	= range(c(-clim0, clim0, as.numeric(acfx$acf)[-1]))
        clx 	= vector(mode="character",length=lag.max)
        clx[which(as.numeric(acfx$acf)[-1]>=0)]="steelblue"
        clx[which(as.numeric(acfx$acf)[-1]<0)]="orange"
        barplot(height = as.numeric(acfx$acf)[-1], names.arg = as.numeric(acfx$lag)[-1], ylim=1.2*ylim, col = clx,
                ylab = "ACF", xlab="lag", main = "ACF of Squared Standardized Residuals", 
                cex.main = 1.6,cex.axis = 1, cex.lab = 0.8)
        abline(h = c(clim0, -clim0), col = "tomato1", lty = 2)
        abline(h = 0, col = "black", lty = 1)
        box()
        mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 1)
        if(vmodel == "fGARCH"){
          mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj = 1.5, col = "gray", cex = 1)
        }
        grid()
      }
      
      # News Impact Curve
      .plot.garchfit.12 = function(x, ...)
      {
        vmodel  = x@model$modeldesc$vmodel
        if(vmodel == "iGARCH"){
          warning(paste("\nplot-->: iGARCH newsimpact not available"))
        } else{
          ni = newsimpact(z = NULL, x)
          ni.y = ni$zy
          ni.x = ni$zx
          xf = ni$xexpr
          yf  = ni$yexpr
          plot( ni.x, ni.y, ylab = yf, xlab = xf, type = "l", lwd = 2, col = "steelblue", main = "News Impact Curve",
                cex.main = 1.6,cex.axis = 1, cex.lab = 0.8)
          
          mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 1)
          if(vmodel == "fGARCH"){
            mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj = 1.5, col = "gray", cex = 1)
          }
          grid()
        }
      }
      # Time Series forecast Plot
      .plot.garchforecast.1 = function(x, n.roll = 0, ...)
      {
        vmodel = x@model$modeldesc$vmodel
        # 1. Time Series:
        nr = x@forecast$n.roll
        if(n.roll > nr) stop("plot-->error: n.roll choice is invalid", call. = FALSE)
        n = x@forecast$n.ahead
        N = x@forecast$N - x@forecast$n.start
        forseries = x@forecast$seriesFor[,n.roll+1]
        forsigma = x@forecast$sigmaFor[,n.roll+1]
        xdates = x@model$modeldata$index[(N+n.roll-min(N,100)):(N+n.roll)]
        fdates = seq(tail(xdates,1), by = x@model$modeldata$period, length.out = n+1)[-1]
        series = x@model$modeldata$data[(N+n.roll-min(N,100)):(N+n.roll)]
        
        xforseries = c(series, forseries)
        series = c(series, rep(NA, n))
        ylim=c(0.95*min(xforseries,na.rm=TRUE), 1.05*max(xforseries,na.rm=TRUE))
        plot(c(xdates, fdates), as.numeric(xforseries), type="l", col="steelblue", 
             main = paste("Forecast Series\n w/th unconditional 1-Sigma bands", sep = "") ,
             ylab="Series",xlab="Time/Horizon", ylim = ylim, cex.main = 1.4, cex.axis = 1, cex.lab = 0.8)
        abline(h = 0, col = "grey", lty = 3)
        Zup = forseries+1*forsigma
        Zdn = forseries-1*forsigma
        for(i in 2:n) rect(fdates[i-1], Zdn[i-1], fdates[i], Zup[i], col = colors()[142], border=NA)
        lines(c(xdates, fdates), series, col="steelblue")
        lines(fdates, forseries, col="tomato1")
        mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
        if(vmodel=="fGARCH"){
          mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
          mtext(paste("Horizon: ",n,sep=""), side = 3, adj = 0, col = "gray", cex = 0.5)
        } else{
          mtext(paste("Horizon: ",n,sep=""), side = 3, adj = 0, col = "gray", cex = 0.5)
        }
        lg.txt = c("Actual","Forecast")
        legend("topleft", legend = lg.txt, col = c("steelblue", "tomato1"), y.intersp = 1.5, pch = 21, cex = 0.7, bty="n")
        box()
        grid()
      }
      
      .plot.garchforecast.3 = function(x, n.roll = 0, ...)
      {
        vmodel = x@model$modeldesc$vmodel
        nr = x@forecast$n.roll
        if(n.roll > nr) stop("plot-->error: n.roll choice is invalid", call. = FALSE)
        n = x@forecast$n.ahead
        N = x@forecast$N - x@forecast$n.start
        forsigma = x@forecast$sigmaFor[,n.roll+1]
        xdates = x@model$modeldata$index[(N+n.roll-min(N,25)):(N+n.roll)]	
        fdates = seq(tail(xdates,1), by = x@model$modeldata$period, length.out = n+1)[-1]	
        sigma = x@model$modeldata$sigma[(N+n.roll-min(N,25)):(N+n.roll)]
        forsigma = c(sigma, forsigma)
        sigma = c(sigma, rep(NA,n))
        ylim=c(0.95*min(forsigma,na.rm=TRUE),1.05*max(forsigma,na.rm=TRUE))
        plot(c(xdates,fdates), forsigma, type = "l", col = "black", 
             main = paste("Forecast Unconditional Sigma\n (n.roll = ", n.roll,")", sep = "") , 
             ylab = "Sigma", xlab = "Time/Horizon", ylim = ylim, cex.main = 1.4, 
             cex.axis = 1, cex.lab = 1.1)
        abline(h = 0, col = "grey", lty = 3)
        lines(c(xdates,fdates), forsigma, , col = "tomato1")
        lines(c(xdates,fdates), sigma, col = "steelblue")
        mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
        if(vmodel == "fGARCH"){
          mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
          mtext(paste("Horizon: ",n,sep=""), side = 3, adj = 0, col = "gray", cex = 0.5)
        } else{
          mtext(paste("Horizon: ",n,sep=""), side = 3, adj = 0, col = "gray", cex = 0.5)
        }
        lg.txt<-c("Actual","Forecast")
        legend("topleft", legend=lg.txt,col=c("steelblue","tomato1"), 
               y.intersp=1.5,pch=21,cex=0.7, bty="n")
        box()
        grid()
      }
      
      .plot.garchboot.2 = function(x, ...)
      {
        vmodel = x@model$modeldesc$vmodel
        vsubmodel = x@model$modeldesc$vsubmodel
        n.ahead = x@model$n.ahead
        seriesfor = fitted(x@forc)
        xser = as.data.frame(x, which = "series")
        N = dim(xser)[1]
        serp = as.data.frame(x, which = "series", type = "q", qtile = c(0.05, 0.25, 0.75, 0.95))
        miny = min(serp[1,])
        maxy = max(serp[4,])
        meanser = apply(xser, 2, FUN = function(x) mean(x))
        plot(seriesfor, type = "l", col = "red", ylim = c(miny, maxy), main = "Series Forecast
             with Bootstrap Error Bands\n (q: 5%, 25%, 75%, 95%)", cex.main = 1.2,cex.axis = 1, cex.lab = 1.1,
             ylab = "returns", xlab = "T+i")
        lines(as.numeric(meanser), col = "black")
        #lines(as.numeric(rx), col = "green")
        points(as.numeric(serp[1,]), col = "steelblue1", pch = 19, cex = 0.5)
        points(as.numeric(serp[2,]), col = "steelblue2", pch = 19, cex = 0.5)
        points(as.numeric(serp[3,]), col = "steelblue3", pch = 19, cex = 0.5)
        points(as.numeric(serp[4,]), col = "steelblue4", pch = 19, cex = 0.5)
        n.overlays = round(0.1*n.ahead)
        if(n.overlays == 0) n.overlays = 1
        dens1 = apply(xser, 2, FUN = function(x) density(x, kernel = "gaussian", 
                                                         n = max(512, round(0.01*N))))
        rugarch:::.densityoverlay(as.numeric(meanser), dens1, n.overlays = n.overlays)
        mtext(paste("GARCH model :",vmodel), side = 4, adj = 0, padj=0, col = "gray", 
              cex = 0.5)
        if(vmodel == "fGARCH"){
          mtext(paste("fGARCH submodel: ", vsubmodel,sep=""), side = 4, adj = 0, 
                padj = 1.5, col = "gray", cex = 0.5)
        }
        legend.txt = c("forecast", "bootstrapped")
        legend("bottomleft", legend = legend.txt, fill = c("red", "black"), col = c("red", "black"),
               bg = "n", bty = "n", cex = 0.6)
      }
      
      .plot.garchboot.3 = function(x, ...)
      {
        vmodel = x@model$modeldesc$vmodel
        vsubmodel = x@model$modeldesc$vsubmodel
        n.ahead = x@model$n.ahead
        fs = rep(NA, n.ahead)
        colr = NULL
        namr = NULL
        if(!is.null(x@model$modeldata$filtered.s)){
          n = length(x@model$modeldata$filtered.s)
          if(n.ahead>n) fs[1:n] = x@model$modeldata$filtered.s else fs = x@model$modeldata$filtered.s[1:n.ahead]
          colr = "green"
          namr = "filtered"
        }
        sigmafor = sigma(x@forc)
        sigdist = vector(mode = "list", length = n.ahead)
        sigp = as.data.frame(x, which = "sigma", type = "q", qtile = c(0.05, 0.25, 0.75, 0.95))
        miny = min(sigp[1,])
        maxy = max(sigp[4,])
        meansig = apply(x@fsigma, 2, FUN = function(x) mean(x))
        plot(sigmafor, type = "l", col = "red", ylim = c(miny, maxy), main = "Sigma Forecast
             with Bootstrap Error Bands\n (q: 5%, 25%, 75%, 95%)", cex.main = 1.2,cex.axis = 1, cex.lab = 1.1,
             ylab = "sigma", xlab = "T+i")
        lines(as.numeric(meansig), col = "black")
        lines(as.numeric(fs), col = colr)
        points(as.numeric(sigp[1,]), col = "steelblue1", pch = 19, cex = 0.5)
        points(as.numeric(sigp[2,]), col = "steelblue2", pch = 19, cex = 0.5)
        points(as.numeric(sigp[3,]), col = "steelblue3", pch = 19, cex = 0.5)
        points(as.numeric(sigp[4,]), col = "steelblue4", pch = 19, cex = 0.5)
        mtext(paste("GARCH model :",vmodel), side = 4, adj = 0, padj=0, col = "gray", 
              cex = 0.5)
        if(vmodel == "fGARCH"){
          mtext(paste("fGARCH submodel: ", vsubmodel,sep=""), side = 4, adj = 0, 
                padj = 1.5, col = "gray", cex = 0.5)
        }
        legend.txt = c("forecast", "bootstrapped", namr)
        legend("topleft", legend = legend.txt, fill = c("red", "black", colr), col = c("red", "black", colr),
               bg = "n", bty = "n", cex = 0.6)
      }
      
      ## Data Structure
      HTML(as.title('Data Structure'),HR=2,file=stdout())
      DS <- matrix(c('Number of observations',dim(dataset)[1]),ncol=2,byrow=T)
      HTML(DS,file=stdout(), innerBorder = 1,align="left")
      
      ## Analysis Description
      HTML(as.title('Analysis Description'),HR=2,file=stdout())
      AD <- matrix(c("GARCH Model",paste0(vmodel,"(", modelinc[8], ",", modelinc[9], ")", sep="")),ncol=2,byrow=T)
      if(method=="fGARCH")  AD <- rbind(AD,c("fGARCH Sub-Model", model$modeldesc$vsubmodel))
      AD <- rbind(AD,c("Mean Model",paste0("ARFIMA(", modelinc[2],",",ifelse(modelinc[4]>0, "d", 0),",",modelinc[3],")", sep = "")))
      AD <- rbind(AD,c("Distribution", model$modeldesc$distribution))            
      AD <- rbind(AD,c("Variable name",ts_var))
      HTML(AD,file=stdout(), innerBorder = 1,align="left")
      
      
      ## Estimates of Model Parameters
      HTML(R2HTML::as.title('Estimates of Model Parameters'),HR=3,file=stdout(),append=TRUE) 
      VL<-round(fit@fit$matcoef,4)
      HTML(VL,file=stdout(), innerBorder = 1,align="left",row.names=F)
      
      
      ## Robust Standard Errors
      if(robust)
      {
        HTML(R2HTML::as.title('Robust Standard Errors'),HR=3,file=stdout(),append=TRUE) 
        VL<-round(fit@fit$robust.matcoef,4)
        if(!is.null(fit@fit$hessian.message)){
          VL<- fit@fit$hessian.message
        }
        HTML(VL,file=stdout(), innerBorder = 1,align="left",row.names=F)
      }
      
      ##LogLikelihood
      LR<-paste("LogLikelihood :", round(fit@fit$LLH,4))
      HTML(LR,file=stdout(), innerBorder = 1,align="left",row.names=F)  
      
      ##Information Criteria
      HTML(R2HTML::as.title('Information Criteria'),HR=3,file=stdout(),append=TRUE)
      itestm = infocriteria(fit)
      colnames(itestm)<-NULL
      HTML(round(itestm,4),file=stdout(), innerBorder = 1,align="left",row.names=F)
      
      
      ##Weighted Ljung-Box Test on Standardized Residuals
      HTML(R2HTML::as.title('Weighted Ljung-Box Test on Standardized Residuals'),HR=3,file=stdout(),append=TRUE) 
      tmp1 = round(rugarch:::.weightedBoxTest(stdresid, p = 1, df = sum(modelinc[2:3])),4)
      HTML(tmp1,file=stdout(), innerBorder = 1,align="left",row.names=F)  
      HTML(paste0("d.o.f=", sum(modelinc[2:3]),"<br>","H0 : No serial correlation"),file=stdout())
      
      
      ##Weighted Ljung-Box Test on Standardized Squared Residuals
      HTML(R2HTML::as.title('Weighted Ljung-Box Test on Standardized Squared Residuals'),HR=3,file=stdout(),append=TRUE) 
      tmp2 = round(rugarch:::.weightedBoxTest(stdresid, p = 2, df = sum(modelinc[8:9])),4)
      HTML(tmp2,file=stdout(), innerBorder = 1,align="left",row.names=F)  
      HTML(paste0("d.o.f=", sum(modelinc[8:9])),file=stdout())
      
      
      ##Weighted ARCH LM Tests
      HTML(R2HTML::as.title('Weighted ARCH LM Tests'),HR=3,file=stdout(),append=TRUE)
      gdf = sum(modelinc[8:9])
      L2 = rugarch:::.weightedarchlmtest(residuals(fit), sigma(fit), lags  = gdf+1, fitdf=gdf)
      L5 = rugarch:::.weightedarchlmtest(residuals(fit), sigma(fit), lags  = gdf+3, fitdf=gdf)
      L10 = rugarch:::.weightedarchlmtest(residuals(fit), sigma(fit), lags = gdf+5, fitdf=gdf)
      alm = matrix(0, ncol = 4, nrow = 3)
      alm[1,1:4] = as.numeric(c(L2$statistic, L2$parameter, L2$p.value))
      alm[2,1:4] = as.numeric(c(L5$statistic, L5$parameter, L5$p.value))
      alm[3,1:4] = as.numeric(c(L10$statistic, L10$parameter, L10$p.value))
      colnames(alm) = c("Statistic", "Shape", "Scale", "P-Value")
      rownames(alm) = c(paste("ARCH Lag[",gdf+1,"]",sep=""), paste("ARCH Lag[",gdf+3,"]",sep=""), paste("ARCH Lag[",gdf+5,"]",sep=""))
      alm<-round(alm,4)
      HTML(alm,file=stdout(), innerBorder = 1,align="left",row.names=F)
      
      
      ##Nyblom stability test
      if(nyblom)
      {
        
        HTML(R2HTML::as.title('Nyblom stability test'),HR=3,file=stdout(),append=TRUE)
        nyb = rugarch:::.nyblomTest(fit)
        if(is.character(nyb$JointCritical)){
          
          HTML("Joint Statistic:  no.parameters>20 (not available)",file=stdout())
          HTML("Individual Statistics:",file=stdout())
          HTML(round(nyb$IndividualStat, 4),file=stdout(), innerBorder = 1,align="left",row.names=F)
          HTML("Asymptotic Critical Values (10% 5% 1%)",file=stdout())
          HTML(matrix(c("Individual Statistic:\t", round(nyb$IndividualCritical, 2)),ncol=2),file=stdout())
          
        } else{
          
          HTML(paste("Joint Statistic: ", round(nyb$JointStat,4),"<br>","Individual Statistics:"),file=stdout())
          HTML(round(nyb$IndividualStat,4),file=stdout(), innerBorder = 1,align="left",row.names=F)
          
          NY<-matrix(c("Asymptotic Critical Values ","(10%","5%", "1%)","Joint Statistic:     ", round(nyb$JointCritical, 3),"Individual Statistic:\t", 
                       round(nyb$IndividualCritical, 2)),ncol=4,byrow=T)
          HTML(NY,file=stdout(), innerBorder = 1,align="left",row.names=F)
          
          
        }
      }
      
      
      ##Sign Bias Test
      if(sign)
      {
        HTML(R2HTML::as.title('Sign Bias Test'),HR=3,file=stdout(),append=TRUE)
        sgtest = signbias(fit)
        sgtest[,1:2]<-round(sgtest[,1:2],4)
        sgtest[,3]<-ifelse(sgtest[,3]==""," ",sgtest[,3])
        HTML(sgtest,file=stdout(), innerBorder = 1,align="left",row.names=T)
        
      }
      
      ##Adjusted Pearson Goodness-of-Fit Test
      if(goodness)
      {
        HTML(R2HTML::as.title('Adjusted Pearson Goodness-of-Fit Test'),HR=3,file=stdout(),append=TRUE)
        gofm = gof(fit,c(20, 30, 40, 50))
        HTML(round(gofm,4),file=stdout(), innerBorder = 1,align="left",row.names=F)
      }
      
      
      ##Original Time series plots
      if(tsplot)
      {
        R2HTML::HTML(R2HTML::as.title('Original Time series plots'),HR=3,file=stdout(),append=TRUE)
        REx_ANA_PLOT(w=1500)
        par(mfrow=c(1,3))
        .plot.garchfit.1(fit);.plot.garchfit.2(fit);.plot.garchfit.3(fit);
        REx_ANA_PLOT_OFF("")
      }
      
      
      ##Original Time series ACF plots
      if(acfplot)
      {
        R2HTML::HTML(R2HTML::as.title('Original Time series ACF plots'),HR=3,file=stdout(),append=TRUE)
        REx_ANA_PLOT(w=1000,h=1000)
        par(mfrow=c(2,2))
        .plot.garchfit.4(fit);.plot.garchfit.5(fit);.plot.garchfit.6(fit);.plot.garchfit.7(fit);
        REx_ANA_PLOT_OFF("")
      }
      
      
      if(sum(is.na(residuals(fit)))==0)
      {
        
        
        ##News-Impact Curve
        if(NIC)
        {
          R2HTML::HTML(R2HTML::as.title('News-Impact Curve'),HR=3,file=stdout(),append=TRUE)
          REx_ANA_PLOT()
          .plot.garchfit.12(fit);
          REx_ANA_PLOT_OFF("")
        }
        
        ##Standardized residual plots
        if(standplot)
        {
          R2HTML::HTML(R2HTML::as.title('Standardized Residual Plots'),HR=3,file=stdout(),append=TRUE)
          REx_ANA_PLOT(w=1000,h=1000)
          par(mfrow=c(2,2))
          .plot.garchfit.8(fit);
          if(dist!="nig") .plot.garchfit.9(fit);
          .plot.garchfit.10(fit);.plot.garchfit.11(fit);
          REx_ANA_PLOT_OFF("")
        }
        
        
        ##Time-series Forecast
        if(pred)
        {
          fore1<-ugarchforecast(fit,n.ahead=h)
          R2HTML::HTML(R2HTML::as.title('Time Series Forecast'),HR=3,file=stdout(),append=TRUE) 
          fore<-cbind(fore1@forecast$seriesFor,fore1@forecast$sigmaFor)
          colnames(fore)<-c("Series Forecast","Sigma Forecast")
          HTML(fore,file=stdout(), innerBorder = 1,align="left",row.names=F)
          
          if(predplot)
          {
            R2HTML::HTML(R2HTML::as.title('Forecast Plots'),HR=3,file=stdout(),append=TRUE) 
            REx_ANA_PLOT(w=1000,h=500)
            par(mfrow=c(1,2))
            .plot.garchforecast.1(fore1);.plot.garchforecast.3(fore1)
            REx_ANA_PLOT_OFF("")
          }
          
          if(boot)
          {
            R2HTML::HTML(R2HTML::as.title('GARCH Bootstrap'),HR=3,file=stdout(),append=TRUE) 
            fore2<-ugarchboot(fit,method="Partial",n.ahead = h,n.bootpred = 2000)
            seriesfor<-apply(fore2@fseries,2,mean)
            sigmafor<-apply(fore2@fsigma,2,mean)
            fore3<-cbind(seriesfor,sigmafor)
            rownames(fore3)<-paste0("T+",1:h)
            colnames(fore3)<-c("Series Forecast","Sigma Forecast")
            HTML(fore3,file=stdout(), innerBorder = 1,align="left")
            REx_ANA_PLOT(w=1000,h=500)
            par(mfrow=c(1,2))
            .plot.garchboot.2(fore2);.plot.garchboot.3(fore2)
            REx_ANA_PLOT_OFF("")
          }
          
        }
      } else { 
        R2HTML::HTML(R2HTML::as.title('Warnings'),HR=2,file=stdout())
        warn.msg7 <- "\a Error : Model residual is Not a Number. Analysis has been stopped."
        R2HTML::HTML(warn.msg7,file=stdout())
      }
      
      
      if(resid | Sresid | condSD)
      {
        O<-data.frame(x=1:length(dataset[,ts_var]))
        if(resid) O<-data.frame(O,Resid=residuals(fit))
        if(Sresid) O<-data.frame(O,SResid=residuals(fit,standard=T))
        if(condSD) O<-data.frame(O,condSD=sigma(fit))
        rownames(O)<-NULL
        O<-O[,-1,drop=F]
      }
      
    } 
    ### TIME 
    R2HTML::HTMLhr(file=stdout())
    R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : GARCH Model",sep=""),file=stdout())
    R2HTML::HTMLhr(file=stdout())
  })
  
  if(exists('O')){return(list(html=html.output,Output=O))}
  else{return(html.output)}
}