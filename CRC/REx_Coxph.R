# Cox비례위험모형
REx_Coxph <- function(dataset, time1, event, cov_set1=NULL, cov_set2=NULL, tdp_time=NULL, tdp_cov=NULL, tdp_var=NULL, vars=NULL, confi=0.95, Select=FALSE, direct= "both", keep_var = NULL, survival_plot=FALSE, a1_survival_plot=FALSE, CHazard_plot=FALSE, log_survival_plot=FALSE, expestim=FALSE, CI=TRUE, CI_level=0.95,out="./test.html") {
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
	load.pkg(c("R2HTML","survival", "KMsurv", "MASS", "ggplot2", "ggfortify", "ggpubr"))

	html.output <- capture.output({
		### Title
		R2HTML::HTML(R2HTML::as.title('Cox Proportional-Hazards Regression'),HR=1,file=out,append=FALSE) 

		### Warnings
		# string time variable : analysis stopped
		Time <- time1
		if(length(tdp_time)!=0) Time <- c(Time,tdp_time)
		
		test.res <- function(i) {
			gg <- !is.numeric(dataset[,i])
			if(!gg) gg <- any(dataset[,i]<0,na.rm=T)
			return(gg)
		}

		ggg <- sapply(Time,test.res)
		if(any(ggg)) {
			warn.msg1 <- paste0('<li> Error : Analysis stopped due to the type of time variable \'',paste(Time[ggg],collapse=', '),'\'. (Time variable should be non-negative numeric.)')
		} 
		
		# invalid status values : analysis stopped
		Status <- event
		if(length(tdp_cov)!=0) Status <- c(Status,tdp_cov)
		is.invalid <- sapply(1:length(Status), function(i) !all(na.omit(unique(dataset[,Status[i]]))%in%c(0,1)))
		if(any(is.invalid)) warn.msg2 <- paste0('<li> Error : Analysis stopped due to the invalid status in the event variable ',paste(Status[is.invalid],collapse=', '),'. (Expected value : 0 - censoring, 1 - event)')

		# explanatory variable type
		if(is.null(vars)) {
			Vars <- c()
		} else {
			Vars <-  unique(unlist(strsplit(vars,":")))
		}
		cov_set2 <- cov_set2[cov_set2%in%Vars]
		if(length(cov_set2)==0) cov_set2 <- NULL
		cov_set1 <- cov_set1[cov_set1%in%Vars]
		if(length(cov_set1)==0) cov_set1 <- NULL

		# character cov_set1 : analysis stopped
		if(length(cov_set1)!=0){
			is.nom <- !apply(dataset[,cov_set1,drop=F],2,is.numeric)
			if(any(is.nom)) {
				warn.msg3 <- paste0("<div style='text-align:left'> <li> Warning : The type of variable '",paste(cov_set1[is.nom],collapse=', '),"' is not numeric but selected as the quantitative variable. It was excluded from the analysis.")
				ec <- cov_set1[is.nom]
				for(jj in ec) vars <- vars[-grep(jj,vars)]
				if(length(vars)==0) {
					vars <- NULL
					warn.msg5 <- "<li> Error: At least one explanatory variable should be included. Analysis has been stopped."
				}
			}
		}

		# character cov_set2
		if(length(cov_set2)!=0){
			is.num <- apply(dataset[,cov_set2,drop=F],2,is.numeric)
			for(i in 1:length(cov_set2)) dataset[,cov_set2[i]] <- as.factor(dataset[,cov_set2[i]])
			if(any(is.num)) warn.msg4 <- paste0("<div style='text-align:left'> <li> Warning : The type of variable '",paste(cov_set2[is.num],collapse=', '),"' is numeric but selected as the qualitative variable. It was coerced into character.")
		}

		# explanatory variable type
		if(is.null(vars)) {
			Vars <- c()
		} else {
			Vars <-  unique(unlist(strsplit(vars,":")))
		}
		cov_set2 <- cov_set2[cov_set2%in%Vars]
		if(length(cov_set2)==0) cov_set2 <- NULL
		cov_set1 <- cov_set1[cov_set1%in%Vars]
		if(length(cov_set1)==0) cov_set1 <- NULL

		if(length(vars)==1 & length(Vars)>1 & all(Vars%in%cov_set2)) warn.qualinte <- '<li> Error: Cox proportional hazard model is not supported when there is only one interaction term consisting of qualitative variables.'

		# Warning in variable selection
		if(Select) {
			if(is.null(vars)) {
				Select <- FALSE
				warn.VS <- '<li> Warning : Variable selection is not supported for intercept only model.'
			}
			keep_var <- keep_var[keep_var%in%vars]
			if(length(keep_var)==0) keep_var <- NULL
		}

		### Warning & analysis stopped
		if(exists("warn.msg1")|exists("warn.msg2")|exists("warn.msg5")|exists("warn.qualinte")) {
			R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=out)
			if(exists("warn.msg1")) R2HTML::HTML(warn.msg1,file=out)
			if(exists("warn.msg2")) R2HTML::HTML(warn.msg2,file=out)
			if(exists("warn.msg5")) {
				if(exists("warn.msg3")) R2HTML::HTML(warn.msg3,file=out)
				R2HTML::HTML(warn.msg5,file=out)
			}
			if(exists("warn.qualinte")) R2HTML::HTML(warn.qualinte,file=out)
		} else { 
			if(length(tdp_var)!=0){
				new.dataset <- na.omit(dataset[,c(time1,event,cov_set1,cov_set2,tdp_cov,tdp_time)])
				# make time-dependent covariates
				cut.point <- lapply(1:length(tdp_var),function(i) unique(with(new.dataset, new.dataset[new.dataset[,tdp_cov[i]]==1,tdp_time[i]])))
				cut.points <- unique(do.call(c,cut.point))
				# cut.points <- unique(with(new.dataset, new.dataset[new.dataset[,tdp_cov]==1,tdp_time]))
				survobj <- Surv(new.dataset[,time1], new.dataset[,event], type="right")
				
				splited <- survSplit(survobj~., zero=0, id="ID", cut=cut.points, data=new.dataset)

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
				dat <- na.omit(dataset[,c(time1,event,Vars)])
			}

			if(length(unique(dat[,1]))==1){
				warn.time <- '<li> Error : Time variable must consist of at least two levels.'
				R2HTML::HTML(warn.time,file=out)
			} else {
				
				fm <- formula(paste0("Surv(",time1,", ",event,") ~ ",paste(vars,collapse='+')))
				null.fm <- formula(paste0("Surv(",time1,", ",event,") ~ 1"))
				fit <- survival::coxph(fm, data=dat)
				null.fit <- survival::coxph(null.fm, data=dat)

				rn <- cov_set1
				if(length(cov_set2)!=0) for(i in cov_set2) rn <- c(rn,paste(i,levels(factor(dataset[,i]))[-1],sep=' - '))
				rn <- c(rn,tdp_var)

				con <- Digits(unname(summary(fit)$concordance[1]))
				con_se <- Digits(unname(summary(fit)$concordance[2]))
				rsq <- Digits(unname(summary(fit)$rsq[1]))
				rsq_max <- Digits(unname(summary(fit)$rsq[2]))
				logtest <- Digits(unname(summary(fit)$logtest[1])); logtest_df <- Digits(unname(summary(fit)$logtest[2])); logtest_pval <- Digits(unname(summary(fit)$logtest[3]))
				waldtest <- Digits(unname(summary(fit)$waldtest[1])); waldtest_df <- Digits(unname(summary(fit)$waldtest[2])); waldtest_pval <- Digits(unname(summary(fit)$waldtest[3]))
				sctest <- Digits(unname(summary(fit)$sctest[1])); sctest_df <- Digits(unname(summary(fit)$sctest[2])); sctest_pval <- Digits(unname(summary(fit)$sctest[3]))
			
				### Description 
				R2HTML::HTML(R2HTML::as.title("Data structure"),HR=2,file=out,append=TRUE) 
				DS <- matrix(c("Number of observations", "Number of total variables", "Number of used variables",nrow(dataset), ncol(dataset),length(c(time1,event,cov_set1,cov_set2,tdp_time,tdp_cov))), ncol=2)
				R2HTML::HTML(DS, file=out,align="left")

				### Variable list
				R2HTML::HTML(R2HTML::as.title("Variable List"),HR=2,file=out,append=TRUE) 
				VL <- matrix(c('Quantitative variable',paste(c(time1,cov_set1,tdp_time),collapse=', '),'Qualitative variable',paste(c(event,cov_set2,tdp_cov),collapse=', ')),byrow=T,ncol=2)
				R2HTML::HTML(VL, file=out,align="left")
				if(exists('warn.msg3')) R2HTML::HTML(warn.msg3,file=out)
				if(exists('warn.msg4')) R2HTML::HTML(warn.msg4,file=out)

				### Analysis description
				R2HTML::HTML(R2HTML::as.title("Analysis Description"),HR=2,file=out,append=TRUE) 
				AD <- matrix(c("Time variable", time1, "Status variable",event,'Model',paste0('Hazard ratio ~ exp(',paste(vars,collapse=" + "),')')),ncol=2,byrow=T)
				if(!is.null(cov_set1)|!is.null(cov_set2)) AD <- rbind(AD,c('Time-independent explanatory variable',paste(c(cov_set1,cov_set2), collapse=", ")))
				if(length(tdp_var)!=0) AD <- rbind(AD,c('Time-dependent explanatory variable', paste(tdp_var, collapse=", ")))
				AD <- rbind(AD,matrix(c('Variable selection',Select),ncol=2,byrow=T))
				#variable selection
				if(Select) {
					direction <- switch(direct,forward="Forward selection",backward="Backward elimination",both="Stepwise regression")
					AD <- rbind(AD,matrix(c('Method for variable selection',direction),ncol=2,byrow=T))
					if(!is.null(keep_var)) AD <- rbind(AD,c('Fixed variable for variable selection',paste(keep_var,collapse=', ')))
				}
				R2HTML::HTML(AD,file=out,align="left") 

				### 분석가정 만족하는지 test
				R2HTML::HTML(R2HTML::as.title("Assessing Test Assumptions"), HR = 2, file = out)

				## proportional hazards assumption
				R2HTML::HTML(R2HTML::as.title("Proportional Hazards Assumptions"), HR = 3, file = out)
				R2HTML::HTML(R2HTML::as.title("Chi-squared Test"), HR = 4, file = out)
				prop_test <- cox.zph(fit)
				test_res <- Digits(prop_test$table)
				colnames(test_res) <- c('Rho','Chisq','P-value')
				R2HTML::HTML(test_res,file=out,align="left")
				
				R2HTML::HTML(R2HTML::as.title("Graphical Checks"), HR = 4, file = out)
				R2HTML::HTML("<li> Note : If the proportional hazards assupmtion is true, beta(t) will be a horizontal line.",file=out,align="left")
		  
				czp <- try(REx_coxzphplot(prop_test, font.main=10, point.col="black"))
				if(class(czp)=='try-error'){
					czp <- try(REx_coxzphplot(prop_test, font.main=10, point.col="black",df=3))
					if(class(czp)=='try-error'){
						czp <- try(REx_coxzphplot(prop_test, font.main=10, point.col="black",df=2))
					}
				}
				
				for(ii in 1:length(czp)){
					REx_ANA_PLOT()
					print(czp[[ii]])
					REx_ANA_PLOT_OFF('')
				}

				## Checking for nonlinearity
				if(!is.null(cov_set1)){
					R2HTML::HTML(R2HTML::as.title("Checking Nonlinearity for Quantitative Variables"), HR = 3, file = out)
					mres.null <- residuals(null.fit)
					for(ii in cov_set1){
						datchecknl <- data.frame(x=dat[,ii], y=mres.null)
						REx_ANA_PLOT()
						print(REx_scatterplot("x", "y", datchecknl, Smooth=T, xlab=ii, ylab='Martingale Residuals of Null Model', Sm.span=0.5))
						REx_ANA_PLOT_OFF('')
					}
				}

				### Default output  
				R2HTML::HTML(R2HTML::as.title("Results of Cox Proportional-Hazards Regression"),HR=2,file=out,append=TRUE) 	

				## Regression coefficient
				R2HTML::HTML(R2HTML::as.title("Coefficients"),HR=3,file=out,append=TRUE) 	
				CE <- summary(fit)$coef
				colnames(CE) <- c('Estimate','exp(Estimate)','SE(Estimate)','Z-value','P-value')
				if(!expestim) CE <- CE[,-2,drop=F]
				if(CI){
					tmp <- merge(CE,confint(fit, level=CI_level),by="row.names",all=TRUE)
					rownames(tmp) <- tmp[,1]
					CE <- tmp[,-1]
					colnames(CE)[(ncol(CE)-1):ncol(CE)] <- paste0(c('Lower','Upper'),' bound of<br>',CI_level*100,'% CI of<br> Estimate')
					if(expestim){
						CE <- cbind(CE,exp(CE[,(ncol(CE)-1):ncol(CE)]))
						colnames(CE)[(ncol(CE)-1):ncol(CE)] <- paste0(c('Lower','Upper'),' bound of<br>',CI_level*100,'% CI of<br> exp(Estimate)')
					}
				}
				
				R2HTML::HTML(Digits(CE), file=out, align="left",digits=15,
					caption=paste('<div style="text-align:left"> Concordance=',con,' (se=',con_se,') <br> Rsquare=',rsq,' (max possible=',rsq_max,') <br> Likelihood ratio test=',logtest, 'on', logtest_df,' degree of freedom,','p-value=',logtest_pval,
					'<br> Wald test=',waldtest, 'on', waldtest_df,' degree of freedom,','p-value=',waldtest_pval,
					'<br> Score (logrank) test=',sctest, 'on', sctest_df,' degree of freedom,','p-value=',sctest_pval))
			
				#### Figures - jhan :ggplot style로 바꾸기
				if(any(survival_plot,a1_survival_plot,CHazard_plot,log_survival_plot)){
					 tfac <- attr(terms(fit$formula), 'factors')
					 if(any(tfac>1)){
						R2HTML::HTML(R2HTML::as.title("Survival Graphs"),HR=3,file=out,append=TRUE)
						plot.msg <- "<li> Warning : Not able to create a survival graphs for models that contain an interaction of a qualitative variable with two more levels."
						R2HTML::HTML(plot.msg,file=out)
					} else {
				
						### survival_plot : Survival function
						if(survival_plot) {
							R2HTML::HTML(R2HTML::as.title("Survival function"),HR=3,file=out,append=TRUE)
							REx_ANA_PLOT()
							print(autoplot(survfit(fit, data=dataset, conf.int=confi), xlab="Time", ylab="Survival probability", ylim=c(0,1)) + theme_bw())
							REx_ANA_PLOT_OFF(ifelse(confi==0,"",paste("<li>",paste0(confi*100,"%"), "confidence interval bounded")))
						}

						### a1_survival_plot : a1-Survival function
						if(a1_survival_plot) {
							R2HTML::HTML(R2HTML::as.title("1-Survival function"),HR=3,file=out,append=TRUE)
							REx_ANA_PLOT()
							print(autoplot(survfit(fit, data=dataset, conf.int=confi), fun="event", xlab="Time", ylab="1-Survival probability", ylim=c(0,1)) + theme_bw())
							REx_ANA_PLOT_OFF(ifelse(confi==0,"",paste("<li>",paste0(confi*100,"%"), "confidence interval bounded")))
						}

						### log_survival_plot : Log Survival function
						if(log_survival_plot) {
							R2HTML::HTML(R2HTML::as.title("Log Survival function"),HR=3,file=out,append=TRUE)
							REx_ANA_PLOT()
							H.hat <- -log(survfit(fit, data=dataset)$surv); H.hat <- c(H.hat, H.hat[length(H.hat)])
							H.hat[!is.finite(H.hat)] <- NA
							h.sort.of <- survfit(fit, data=dataset)$n.event / survfit(fit, data=dataset)$n.risk
							H.tilde <- vector()
							for(i in 1:length(h.sort.of)) H.tilde[i] <- sum(h.sort.of[1:i])
							H.tilde <- c(H.tilde, H.tilde[length(H.tilde)])
							print(autoplot(survfit(fit, data=dataset, conf.int=confi), fun="log", xlab="Time", ylab="log(Survival probability)", ylim=range(-c(na.omit(H.hat), na.omit(H.tilde)))) + theme_bw())
							REx_ANA_PLOT_OFF(ifelse(confi==0,"",paste("<li>",paste0(confi*100,"%"), "confidence interval bounded")))
						}
						
						### CHazard_plot : Cumulative Hazard function
						if(CHazard_plot) {
							R2HTML::HTML(R2HTML::as.title("Cumulative Hazard function"),HR=3,file=out,append=TRUE)
							REx_ANA_PLOT()
							H.hat <- -log(survfit(fit, data=dataset)$surv); H.hat <- c(H.hat, H.hat[length(H.hat)])
							H.hat[!is.finite(H.hat)] <- NA
							h.sort.of <- survfit(fit, data=dataset)$n.event / survfit(fit, data=dataset)$n.risk
							H.tilde <- vector()
							for(i in 1:length(h.sort.of)) H.tilde[i] <- sum(h.sort.of[1:i])
							H.tilde <- c(H.tilde, H.tilde[length(H.tilde)])
							print(autoplot(survfit(fit, data=dataset, conf.int=confi), fun="cumhaz", xlab="Time", ylab="Cumulative hazard", ylim=range(c(na.omit(H.hat), na.omit(H.tilde)))) + theme_bw())
							REx_ANA_PLOT_OFF(ifelse(confi==0,"",paste("<li>",paste0(confi*100,"%"), "confidence interval bounded")))
						}
							
					}
				}
				
				## Variable selection
				if(Select) stepAIC.wj(object=fit,dataset=dat,type='coxph',noint=T,direct=direct,keep_var=keep_var,hr=3,vars=vars,time1=time1,event=event)
			}
		}

		# Used R packages
		R2HTML::HTML(R2HTML::as.title("Used R Packages"),HR=2,file=out)
		pkg.list <- list(list("Main results","Surv, coxph","survival"))
		R2HTML::HTML(used.pkg(pkg.list), file=out)

		### TIME 
		R2HTML::HTMLhr(file=out)
		R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Cox Proportional-Hazards Regression",sep=""),out)
		R2HTML::HTMLhr(file=out)		
	})
	return(html.output)
}

df20186118174 <- data.frame("id"=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,342,343,344,345,346,347,348,349,350,351,352,353,354,355,356,357,358,359,360,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375,376,377,378,379,380,381,382,383,384,385,386,387,388,389,390,391,392,393,394,395,396,397,398,399,400,401,402,403,404,405,406,407,408,409,410,411,412,413,414,415,416,417,418,419,420,421,422,423,424,425,426,427,428,429,430,431,432,433,434,435,436,437,438,439,440,441,442,443,444,445,446,447,448,449,450,451,452,453,454,455,456,457,458,459,460,461,462,463,464,465,466,467,468,469,470,471,472,473,474,475,476,477,478,479,480,481,482,483,484,485,486,487,488,489,490,491,492,493,494,495,496,497,498,499,500),"bweight"=c(2974,3270,2620,3751,3200,3673,3628,3773,3960,3405,4020,2724,3001,3039,3662,3035,3351,3804,3573,3283,2894,1203,3306,3753,2844,3585,3798,3164,3739,1780,4022,3942,2887,2391,3911,3509,3566,3652,3279,3007,3053,3503,3120,3743,3592,3184,3234,2581,3305,3678,3542,2148,3774,3079,3465,2887,4501,3375,3886,2849,2002,2213,2797,3303,3296,2921,2929,3385,1946,3354,3189,3392,2696,2597,3409,3190,3571,4512,2215,3315,3341,3451,3338,2507,3316,4141,4071,2893,3064,3603,3554,4027,2418,3092,2671,3430,3244,3259,3942,3576,3784,2796,3226,3138,3715,3773,2623,2830,3332,2911,2894,3593,2950,3605,3198,3183,2487,2092,3932,3995,3592,4092,1663,1546,3149,3062,3467,2610,3261,2926,3727,3545,3423,2959,3879,2935,3552,2764,4069,3267,3178,628,3288,3233,3985,3676,2579,3461,2990,2995,2922,3449,3292,3621,3636,3471,3360,2338,3567,3007,2740,2428,3027,3261,1999,3117,2751,3824,2762,4057,3501,3503,3001,1791,3582,3575,3117,3260,2704,2545,1019,2609,2127,3126,3606,2257,3051,3146,3696,1402,2539,3550,2252,2842,2679,3526,3184,3180,3561,2698,4131,1880,3807,3005,3448,4340,3518,3296,2505,3188,3109,3134,2992,708,3571,3463,1741,4122,3141,2328,2404,3079,2837,3216,2951,981,3882,3385,2545,3096,3138,3341,2768,3367,3768,3486,2913,4436,3943,3009,3041,3099,3647,2698,3078,3976,2729,3254,3163,3919,2718,3269,693,2879,2590,3311,2831,2991,2699,3695,3250,1570,3554,4205,3694,3505,3023,3723,2482,2622,3558,2090,3446,2938,864,1764,3096,3003,3197,3664,1618,3059,2360,2961,3323,2982,2974,4423,3105,2446,3590,3450,3291,2866,2802,2953,3379,3386,3770,2188,3813,3576,3290,3502,3591,3734,2989,3546,3783,3754,2297,1500,1595,3376,3349,3727,2719,3365,4553,3007,1801,3451,3249,3092,2646,3739,3767,924,3428,3004,3625,2804,2274,3075,2699,2605,1325,2736,3192,3184,2855,3392,3040,3249,3156,2606,2859,2558,2905,3331,3730,3683,3398,3435,4319,4179,1431,3154,2914,3646,3980,3399,3621,2659,2890,2425,2694,3804,2784,1541,3282,3695,2769,3286,2744,3193,2666,3102,3333,2924,1824,3045,2943,3425,2863,3024,4287,3207,3090,3257,3370,3200,4226,3160,2732,2977,2968,3562,3316,3482,2969,4035,3704,3082,3404,3041,3519,3432,3733,2622,2497,2797,3337,2407,3383,3134,3559,3516,3497,3350,2969,3851,2697,2801,2362,2978,3082,3545,1874,3140,3457,3093,2353,3184,3222,2944,3179,3463,1324,3064,3581,2399,3830,3376,2098,3601,3386,3122,3551,2944,4224,2719,4133,2885,3042,4197,3419,3294,4304,3306,3102,4041,4182,2949,3757,3486,3770,3087,4300,3122,3337,3565,3611,4516,3718,1938,2878,2639,2417,2963,2843,2595,3455,3061,3095,2552,3409,3152,3172,3237,3948,3542,3322,3349,2968,2852,3187,3054,3178,2918),"lowbw"=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,1,0,0,1,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),"gestwks"=c(38.52000046,NA,38.15000153,39.79999924,38.88999939,40.97000122,42.13999939,40.20999908,42.02999878,39.33000183,40.91999817,38.97000122,41.02000046,40.04999924,39.22999954,38.95999908,39.34999847,38.99000168,40.74000168,39.45000076,39.27000046,32.79999924,38.72000122,40.18999863,40.13999939,39.54999924,39.38000107,36.52999878,41.65999985,40.45000076,39.70999908,39.18000031,38.50999832,36.40000153,40.45999908,39.47000122,39.11999893,40.93999863,39.34999847,38.79999924,40.66999817,37.86999893,38.50999832,39.04000092,40.36000061,37.43999863,39.58000183,35.70000076,39.18999863,40.63000107,39.43000031,38.43999863,38.63999939,38.43999863,39.58000183,38.90000153,40.02999878,40,40.91999817,37.61999893,36.47999954,37.68000031,36.34999847,39.95000076,38.20999908,39.61999893,38.27000046,41.68000031,37.88999939,39.56999969,37.16999817,37.74000168,34.59000015,34.66999817,40.36999893,39.04000092,38.27000046,38.68000031,38.36999893,NA,38.58000183,38.13999939,40.47999954,38.47999954,41.24000168,39.00999832,38.93000031,37.33000183,38.91999817,40.06000137,NA,40.18000031,36.36999893,40.58000183,34.86999893,39.86000061,42.06000137,38.18000031,39.97000122,39.84999847,39.00999832,38.09000015,39.49000168,39.45999908,41.04000092,37.70000076,NA,37.08000183,39.09999847,38.25999832,39.11000061,41.04999924,40.15000153,37.47000122,36.54000092,40.84000015,36.27000046,37.95999908,38.74000168,40.81999969,40.09999847,39.83000183,35.13999939,37.84999847,41.11000061,40.54999924,39.86000061,39.63999939,39.45000076,38.63000107,39.75,37.40000153,40.11000061,37.25999832,40.72999954,36.91999817,37.88000107,37.65000153,40.5,40.90000153,37.04000092,26.95000076,40.29000092,39.59999847,40.36000061,37.74000168,39.88000107,37.58000183,41.11000061,39.52000046,39.75999832,41.02000046,39.34999847,41.88000107,40.99000168,38.59999847,38.68999863,38.58000183,40.88999939,37.16999817,37.09999847,37.22000122,40.22999954,37.63000107,36.84999847,40.54999924,34.70000076,40.02999878,37.38999939,38.33000183,38.58000183,40.06999969,37.86000061,37.74000168,39.70999908,38.59999847,40.40999985,37.88999939,37.97000122,35.72999954,28.04000092,39.93000031,39.75999832,40.5,36.40999985,38.40999985,39.77000046,41.22000122,41.84000015,32.47000122,39.02000046,40.04000092,38.06000137,37.93999863,41,40.04000092,39.84000015,40.93000031,38.36000061,36.45000076,42.29999924,33.66999817,40.47999954,39.66999817,38.91999817,39.81000137,39.79999924,38.25,35.61999893,36.15999985,39.13999939,38.88000107,37.59000015,27.32999992,37.77999878,40.34999847,32.63999939,38.93999863,41.95999908,34.13999939,36.93999863,39.84999847,39.41999817,39.93999863,38.99000168,27.98999977,40.33000183,38.34000015,36.81999969,38.18000031,39.97999954,39.43999863,39.70000076,38.75,40.41999817,43.15999985,39.04999924,40.79999924,39.52000046,34.63999939,39.11999893,37.56999969,39.5,40.79999924,38.27000046,41.29000092,38.54999924,39.20000076,40.45000076,40.74000168,37.74000168,39.09000015,30.70999908,37.45999908,37.33000183,40.79999924,38.43999863,37.04999924,38.25999832,42.86000061,40.25,34.02000046,38.02000046,39.54000092,40.99000168,38.65999985,40.75,39.91999817,37.54999924,37.13000107,37.49000168,37.22000122,40.56000137,37.45999908,24.69000053,30.64999962,39.90999985,39.84999847,38.99000168,39.84999847,32.81000137,NA,37.41999817,40.74000168,38.50999832,38.99000168,39.20999908,41.79000092,39.65000153,37.66999817,37.83000183,40.81999969,39.68999863,40.00999832,37.33000183,40.34000015,39.61000061,40.47999954,40.97000122,37.43000031,38.47999954,38.56000137,41.16999817,40.63999939,37.90000153,39.65000153,39.04000092,38.90000153,41.70999908,41.18000031,36.93999863,35.27000046,30.52000046,40.63999939,40.70000076,39.81000137,37.84000015,41.09999847,40.06999969,40.09000015,35.68000031,40.34999847,39.49000168,38.22999954,38.77000046,40.79000092,38.68999863,30.85000038,39.75999832,39.50999832,39.61999893,38.04000092,35.40999985,38.15999985,41.38000107,39.22000122,34.66999817,36.24000168,38.72000122,40.72999954,38.06999969,38.86999893,39.63000107,39.95999908,38.54999924,35.5,40.20999908,39.04000092,38.22999954,38.29000092,41.04000092,37.79000092,40.95000076,40.47999954,39.13000107,41.04999924,32.40999985,38.97000122,38.90999985,40.02999878,42.04999924,NA,41.02999878,33.22999954,37.52999878,36.22999954,37.79000092,38.88999939,38.22999954,32.52999878,39.16999817,38.77999878,39.34999847,40.59000015,38.75,39.40000153,36.43999863,38.95000076,37.27000046,38.68999863,NA,37.65000153,37.95999908,39.49000168,38.15000153,40.56999969,39.56999969,39.75999832,39.95999908,38.36000061,40.11999893,39.75,38.04999924,38.25999832,36.99000168,NA,40.38000107,37.95999908,NA,NA,37.81999969,37.90999985,38.36000061,39.63999939,38.83000183,39.77000046,40.45000076,40.36000061,40.13000107,37.11000061,34.06000137,39.70000076,39.90000153,35.95000076,39.38000107,39.49000168,36.27999878,39.74000168,40.15000153,39.59999847,40.27000046,40.65999985,37.54000092,38.66999817,37.36000061,38.68000031,42.20000076,40.34999847,34.18999863,38.86000061,39.90000153,36.36999893,39.56999969,41.99000168,40.18999863,37.86000061,39.74000168,39.43000031,31.29000092,40.93000031,39.90000153,38.45000076,39.93999863,40.18000031,31.37000084,38.97999954,40.61999893,40.77999878,40.81999969,37.93999863,40.20000076,35.97000122,40.61000061,40.02999878,39.43000031,39.18000031,39.06000137,38.75,41.25999832,41.31000137,39.20000076,41.02000046,41.13000107,36.09000015,38.90000153,39.86000061,39.56000137,39.02000046,40.09000015,38.93000031,38.47999954,37.09999847,40.43999863,40.29999924,40.66999817,31.70999908,39.56000137,39.79999924,34.34000015,39.70000076,39.88000107,39.52000046,38.84999847,39.54000092,39.22000122,38.09000015,40.08000183,38.09000015,38.25999832,37.59000015,39.65999985,38.90999985,39.43000031,38.02999878,41.00999832,38.45000076,38.02999878,38.5,39.91999817,37.97000122),"preterm"=c(0,NA,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,NA,0,0,0,0,0,0,0,0,0,0,NA,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,NA,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,1,0,0,0,1,0,0,1,0,0,1,1,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,NA,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,NA,0,1,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,NA,0,0,0,0,0,0,0,0,0,0,0,0,0,1,NA,0,0,NA,NA,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),"matage"=c(34,30,35,31,33,33,29,37,36,39,37,35,38,34,28,38,34,39,40,33,27,39,33,38,33,32,32,38,35,35,33,36,35,29,41,34,33,32,30,32,35,34,38,35,32,29,34,34,36,34,26,29,28,34,39,28,27,36,40,33,37,39,37,27,35,40,36,29,33,41,39,39,33,32,34,35,38,36,35,29,25,39,40,32,31,33,34,38,37,39,38,33,35,37,31,38,37,41,34,33,32,33,38,40,29,32,40,41,39,38,33,33,33,38,40,38,38,31,27,34,33,34,37,34,34,43,34,38,37,31,41,38,40,36,38,34,32,35,39,36,31,32,37,34,37,30,37,31,28,35,38,33,31,36,36,35,38,37,39,40,36,37,35,25,39,40,34,36,34,42,28,34,39,37,32,40,34,39,34,30,31,27,36,33,36,26,35,34,31,38,32,37,39,38,32,32,42,38,37,36,30,27,39,35,39,40,37,34,40,40,38,35,37,27,32,28,41,31,35,35,40,30,34,38,38,29,34,38,35,29,37,32,32,32,35,33,35,27,38,37,40,37,32,26,34,34,41,36,29,39,31,37,34,37,36,37,38,32,33,38,27,28,40,29,36,39,34,34,34,40,37,31,36,34,40,30,38,40,38,33,31,32,29,29,34,28,38,32,39,36,32,32,32,34,31,33,38,32,36,35,27,35,30,37,32,28,30,27,32,40,33,34,33,25,34,29,37,34,34,29,37,37,32,39,40,31,37,35,35,31,31,30,33,37,35,29,36,24,37,34,35,30,31,30,32,27,37,36,33,37,28,35,34,34,37,27,35,30,32,36,31,29,37,29,32,34,38,29,35,31,30,35,30,27,35,30,33,32,30,32,33,36,35,36,37,37,37,35,40,34,30,31,37,26,34,32,38,30,30,33,38,36,39,37,33,25,33,35,38,43,28,30,33,38,30,29,33,30,26,39,32,39,43,30,30,34,36,40,30,32,30,34,35,39,26,29,33,32,25,29,35,32,35,37,31,35,33,23,29,26,27,34,38,36,29,33,33,29,33,34,32,34,38,38,35,27,26,30,31,33,37,34,38,29,38,36,39,33,37,39,33,36,35,30,31,33,36,38,31,32,31,33,35,33,34,28,38,26,31,31),"hyp"=c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,1,0,0,1,0,0,0,0,0,1,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,1,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0),"sex"=c(2,1,2,1,1,2,2,1,2,1,1,2,2,1,2,1,2,1,2,2,2,1,2,2,2,1,2,1,2,1,2,2,1,1,2,2,1,1,2,1,2,1,1,1,1,1,1,1,1,1,1,2,1,2,1,1,1,1,1,1,2,1,1,2,1,1,2,2,1,2,1,1,1,1,1,2,1,1,2,1,2,1,1,2,1,1,1,2,1,1,1,2,1,2,1,1,2,1,2,2,2,2,2,2,2,1,1,2,2,2,1,2,2,1,2,2,2,2,1,1,2,1,1,2,1,1,2,1,2,1,1,1,1,2,1,2,1,1,1,2,1,2,2,1,1,1,2,2,2,2,2,2,1,2,1,1,2,1,1,2,2,1,2,1,1,1,1,1,2,1,1,2,2,1,1,1,2,2,2,2,1,2,1,2,1,1,2,2,1,2,2,1,2,1,2,2,2,2,1,1,1,2,1,1,2,1,1,1,1,2,2,2,2,2,1,1,2,1,2,2,2,1,2,1,1,1,1,2,1,2,2,2,1,2,2,2,1,1,2,1,1,2,1,2,2,2,1,1,2,2,2,2,1,1,1,2,2,1,2,2,2,2,1,1,1,1,1,1,2,1,2,1,2,2,2,2,2,1,1,1,2,2,2,1,2,1,2,1,2,1,1,1,2,2,2,2,1,2,1,1,2,1,2,1,1,2,1,2,1,2,2,1,1,2,1,1,2,1,1,2,2,1,1,1,1,2,2,2,2,1,1,1,1,2,1,1,2,2,1,2,2,1,1,2,1,1,2,1,2,1,1,1,1,2,1,1,2,2,2,1,1,1,1,1,1,2,2,2,2,2,1,1,2,1,1,2,1,2,1,1,1,2,2,1,1,2,1,1,1,2,1,1,1,1,1,1,1,2,2,1,2,2,1,2,2,1,1,2,2,1,1,2,1,1,2,1,2,1,2,2,1,1,2,2,1,2,2,2,2,1,1,2,2,2,2,2,2,1,2,1,1,2,2,1,1,1,2,2,2,2,1,2,2,2,2,1,1,2,1,2,1,1,1,1,2,1,1,1,2,1,2,1,1,1,1,2,2,1,2,1,2,2,1,1,2,2,1,1,1,1,2,2,2,1,1,2,1,2,2,1),"time_longi"=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25),"obs_longi"=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),"test"=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),check.names=FALSE);
 