# 이표본비율검정

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
					tmp.dataset <- dataset[,c(res_var1,res_var2),drop=F]
					n1<-length(which(!is.na(dataset[,res_var1])))
					n2<-length(which(!is.na(dataset[,res_var2])))
					x1<-length(which(dataset[,res_var1]==event1))
					x2<-length(which(dataset[,res_var2]==event2))
	
					data_desc<-c("Data type : Summarized data",
					paste("Number of non-missing trials (Sample 1): ",n1,sep=''),
					paste("Number of Successes (Sample 1): ",x1,sep=''),
					paste("Number of non-missing trials (Sample 2): ",n2,sep=''),
					paste("Number of Successes (Sample 2): ",x2,sep=''))


					DS <- matrix(c('Number of observations',nrow(dataset),'Number of total variables',ncol(dataset),'Number of used variables',2),byrow=T,ncol=2)	# Data structure
					VL1 <- matrix(c('Dependent variable (Sample 1)',res_var1,'Event (Sample 1)',paste(event1,"-Success, ",unique(dataset[,res_var1])[which(!unique(dataset[,res_var1])%in% c(NA,event1))],"-Failure",sep='')),byrow=T,ncol=2)	# Varaible list
					VL2 <- matrix(c('Dependent variable (Sample 2)',res_var2,'Event (Sample 2)',paste(event2,"-Success, ",unique(dataset[,res_var2])[which(!unique(dataset[,res_var2])%in% c(NA,event2))],"-Failure",sep='')),byrow=T,ncol=2)	# Varaible list
					VL <- rbind(VL1,VL2)

				    }
				if(class(try(length(dataid),s=T))!='try-error'){
					tmp.dataset <- dataset[,c(dataid,res_var),drop=F]
					n1<-length(which(dataset[,dataid]==names(table(dataset[,dataid]))[1] & !is.na(dataset[,res_var])))
					n2<-length(which(dataset[,dataid]==names(table(dataset[,dataid]))[2] & !is.na(dataset[,res_var])))
					x1<-length(which(dataset[,dataid]==names(table(dataset[,dataid]))[1] & dataset[,res_var]==event))
					x2<-length(which(dataset[,dataid]==names(table(dataset[,dataid]))[2] & dataset[,res_var]==event))

					data_desc<-c("Data type : Summarized data",
					paste("Number of non-missing trials (Sample 1): ",n1,sep=''),
					paste("Number of Successes (Sample 1): ",x1,sep=''),
					paste("Number of non-missing trials (Sample 2): ",n2,sep=''),
					paste("Number of Successes (Sample 2): ",x2,sep=''))

					DS <- matrix(c('Number of observations',nrow(dataset),'Number of total variables',ncol(dataset),'Number of used variables',2),byrow=T,ncol=2)	# Data structure
					VL <- matrix(c('Sample variable',dataid,'Dependent variable',res_var,'Event',paste(event,"-Success, ",unique(dataset[,res_var])[which(!unique(dataset[,res_var])%in% c(NA,event))],"-Failure",sep='')),byrow=T,ncol=2)	# Varaible list

				    }

			}

			if(class(try(nrow(dataset),s=T))=='try-error'){
				data_desc<-c("Data type : Summarized data",
					paste("Number of non-missing trials (Sample 1): ",n1,sep=''),
					paste("Number of Successes (Sample 1): ",x1,sep=''),
					paste("Number of non-missing trials (Sample 2): ",n2,sep=''),
					paste("Number of Successes (Sample 2): ",x2,sep=''))
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

		if(exists('DS')){
			R2HTML::HTML(R2HTML::as.title("Data structure"),HR=2,file=stdout())
			R2HTML::HTML(DS,file=stdout(), innerBorder = 1,align="left")
		}

		if(exists('VL')){
			R2HTML::HTML(R2HTML::as.title("Variable list"),HR=2,file=stdout())
			R2HTML::HTML(VL,file=stdout(), innerBorder = 1,align="left")
		}

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
