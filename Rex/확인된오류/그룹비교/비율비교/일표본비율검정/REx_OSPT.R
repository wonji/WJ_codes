# 일표본비율검정

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
			tmp.dataset <- dataset[,res_var,drop=F]; tmp.dataset <- tmp.dataset[!is.na(tmp.dataset)]
			x<-length(which(tmp.dataset==event))
			n<-length(which(!is.na(tmp.dataset)))
			DS <- matrix(c('Number of observations',nrow(dataset),'Number of total variables',ncol(dataset),'Number of used variables',1),byrow=T,ncol=2)	# Data structure
			VS <- matrix(c('Dependent variable',res_var,'Event',paste(event,"-Success, ",unique(tmp.dataset)[which(!unique(tmp.dataset)%in% c(event,NA))],"-Failure",sep='')),byrow=T,ncol=2)	# Variable list
			data_desc<-c("Data type : Original data",
					paste("Number of non-missing trials : ",n,sep=''),
					paste("Number of Successes : ",x,sep=''))
		}

		if(class(try(nrow(dataset),s=T))=='try-error'){
			data_desc<-c("Data type : Summarized data",
					paste("Number of non-missing trials : ",n,sep=''),
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
		
		if(exists('DS')){
			R2HTML::HTML(R2HTML::as.title("Data structure"),HR=2,file=stdout())
			R2HTML::HTML(DS,file=stdout(), innerBorder = 1,align="left")
		}

		if(exists('VS')){
			R2HTML::HTML(R2HTML::as.title("Variable list"),HR=2,file=stdout())
			R2HTML::HTML(VS,file=stdout(), innerBorder = 1,align="left")
		}

		R2HTML::HTML(R2HTML::as.title("Analysis description"),HR=2,file=stdout())
		R2HTML::HTML(f_res_OSPT[1:4,],file=stdout(), innerBorder = 1,align="left")

		R2HTML::HTML(R2HTML::as.title("Results of One sample proportion test"),HR=2,file=stdout())
		R2HTML::HTML(f_res_OSPT[5:nrow(f_res_OSPT),],file=stdout(), innerBorder = 1,align="left")

		R2HTML::HTMLhr(file=stdout())
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : One sample proportion test",sep=""),file=stdout())
		R2HTML::HTMLhr(file=stdout())
	})
	options(global_op)	#수정
	return(html.output)
}
