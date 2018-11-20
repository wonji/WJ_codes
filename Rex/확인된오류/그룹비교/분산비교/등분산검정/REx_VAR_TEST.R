# 등분산 검정

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
    R2HTML::HTML(R2HTML::as.title("Test for Homogeneity of Variances"),HR=1,file=stdout(),append=FALSE)
    
    # Warning
    if(class(dataset[,res_var])!="numeric") {
      warn.msg1 <- '\a Error : Response variable should be numeric. Analysis has been stopped.'
    }
    
    #group_var
    groupvar_list<-unlist(strsplit(group_var,"\\+"))
    var_num<-which(colnames(dataset) %in% groupvar_list)
    dataset$group_var <- factor(do.call("paste", c(dataset[var_num], sep=" ")))
    is.num <- sapply(groupvar_list,function(i) is.numeric(dataset[,i]))
    if(any(is.num)) warn.msg2 <- paste0("\a Warning : The type of variable '",paste(groupvar_list[is.num],collapse=', '),"' is numeric but selected as the group variable. It was coreced into character.")
    
    # Warning
    if(any(table(dataset[,unlist(strsplit(group_var,'\\+'))])==1)) {
      warn.msg3 <- paste0("\a Error : There must be at least 2 observations in each group. Analysis has been stopped.")
    }
    
    if(exists('warn.msg1')|exists('warn.msg3')){
      R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())
      if(exists('warn.msg1')) R2HTML::HTML(warn.msg1,file=stdout())
      if(exists('warn.msg3')) R2HTML::HTML(warn.msg3,file=stdout())
    } else {
      ## Data structure
      R2HTML::HTML(R2HTML::as.title("Data Structure"),HR=2,file=stdout())
      DS <- matrix(c('Number of observations','Number of total variables','Number of used variables',nrow(dataset),ncol(dataset),length(strsplit(group_var,'\\+')[[1]])+1),ncol=2)
      R2HTML::HTML(DS,file=stdout(), innerBorder = 1,align="left")
      
      ## Analysis description
      R2HTML::HTML(R2HTML::as.title("Analysis Description"),HR=2,file=stdout()) 
      AD <- matrix(c('Method',methods,'Response variable',res_var,'Group variable',gsub('\\+',', ',group_var)),byrow=T,ncol=2)
      
      if(methods=='Bartlett'){
        AD1 <- matrix(c('Null Hypothesis (H0)','Variances in each of the groups are the same','Alternative Hypothesis (H1)','not H0'),byrow=T,ncol=2)
        AD <- rbind(AD,AD1)
        R2HTML::HTML(AD,file=stdout(), innerBorder = 1,align="left")
      } else if(methods=='Levene') {
        AD1 <- matrix(c('Centering Method',centering,'Null Hypothesis (H0)','Variances in each of the groups are the same','Alternative Hypothesis (H1)','not H0'),byrow=T,ncol=2)
        AD <- rbind(AD,AD1)
        R2HTML::HTML(AD,file=stdout(), innerBorder = 1,align="left")
      } else {
        H1 <- ifelse(alternative_hypothesis=="two.sided","Variance ratios between groups are not equal to 1",ifelse(alternative_hypothesis=="less","Variance ratios between groups are less than 1","Variance ratios between groups are greater than 1"))
        AD1 <- matrix(c('Null Hypothesis (H0)','Variances in each of the groups are the same','Alternative Hypothesis (H1)',H1,'Note','F-test is conducted in pairs of groups'),byrow=T,ncol=2)
        AD <- rbind(AD,AD1)
        R2HTML::HTML(AD,file=stdout(), innerBorder = 1,align="left")
      }
      if(exists('warn.msg2')) R2HTML::HTML(warn.msg2,file=stdout())
      
      dataset <- dataset[complete.cases(dataset[,c(res_var,unlist(strsplit(group_var,'\\+')))]),]

      #분석결과
      if(methods=='Bartlett'){
        # Bartlett methods 
        R2HTML::HTML(R2HTML::as.title("Bartlett Test of homogeneity of variances"),HR=2,file=stdout()) 
        
        if(length(strsplit(group_var,split='\\+')[[1]])==1){
          result<-bartlett.test(formula(paste0(res_var, '~ ' ,group_var)), data=dataset)
          results <- matrix(c("Bartlett's K-squared",round(result$statistic,digits),"degree of freedom",result$parameter,"P-value",round(result$p.value,digits)),byrow=T,ncol=2)
        } else{
          interaction.term<-gsub(pattern="\\+", replacement=", ", x=group_var)
          result<- bartlett.test(formula(paste0(res_var, '~interaction(' , interaction.term,')')), data=dataset)
          results <- matrix(c("Bartlett's K-squared",round(result$statistic,digits),"degree of freedom",result$parameter,"P-value",round(result$p.value,digits)),byrow=T,ncol=2)
        }
        R2HTML::HTML(results,file=stdout(), innerBorder = 1,align="left")
      } else if(methods=="Levene"){
        # Levene methods 
        R2HTML::HTML(R2HTML::as.title("Levene's Test for homogeneity of variance"),HR=2,file=stdout()) 
        
        if(length(strsplit(group_var,split='\\+')[[1]])==1){
          result<-leveneTest(formula(paste0(res_var ,'~ as.factor(' ,group_var,')')), data=dataset,center=centering)
          results <- matrix(c("F value",round(result$'F value'[1],digits),"degree of freedom(numerator)",result$Df[1],"degree of freedom(denominator)",result$Df[2],'p-value',round(result$'Pr(>F)'[1],digits)),byrow=T,ncol=2)
        } else {
          result<-leveneTest(formula(paste0(res_var, '~ as.factor(' , gsub(pattern="\\+", replacement="*", x=group_var),')')), data=dataset,center=centering)
          results <- matrix(c("F value",round(result$'F value'[1],digits),"degree of freedom(numerator)",result$Df[1],"degree of freedom(denominator)",result$Df[2],'p-value',round(result$'Pr(>F)'[1],digits)),byrow=T,ncol=2)
        }
        R2HTML::HTML(results,file=stdout(), innerBorder = 1,align="left")
      } else if(methods=="F-test") {
        # F-test
        R2HTML::HTML(R2HTML::as.title("F Test for homogeneity of variance"),HR=2,file=stdout()) 
        
        grp.var0 <- strsplit(group_var,"\\+")[[1]]
	tmp.dataset <- dataset[,grp.var0,drop=F]
	for(ii in grp.var0) tmp.dataset[,ii] <- paste0(ii,'(',tmp.dataset[,ii],')')
        grp.var <- apply(tmp.dataset,1,function(x) paste(x,collapse='&'))
        dataset$grp.var <- grp.var
        
        a <- combn(unique(grp.var),2)
        for(i in 1:ncol(a)){
          dataset1 <- dataset[dataset$grp.var%in%a[,i],]
          result<-var.test(formula(paste0(res_var,"~ grp.var")),alternative=alternative_hypothesis,conf.level=confi.level,data=dataset1) 
          H1 <- ifelse(alternative_hypothesis=="two.sided","Variance ratio is not equal to 1",ifelse(alternative_hypothesis=="less","Variance ratio is less than 1","Variance ratio is greater than 1"))
          if(CI){
            Fresult <- matrix(c("Comparison groups",paste0(a[1,i],' VS ',a[2,i]),"F statistic",round(result$statistic,digits),"degree of freedom(numerator)",result$parameter[1],"degree of freedom(denominator)",result$parameter[2],"P-value",round(result$p.value,digits),"Ratio Estimate",round(result$estimate,digits),paste0(confi.level*100,"% CI"),paste0("(",round(result$conf.int[1],digits),",",round(result$conf.int[2],digits),")")),byrow=T,ncol=2)
          } else {
            Fresult <- matrix(c("Comparison groups",paste0(a[1,i],' VS ',a[2,i]),"F statistic",round(result$statistic,digits),"degree of freedom(numerator)",result$parameter[1],"degree of freedom(denominator)",result$parameter[2],"P-value",round(result$p.value,digits),"Ratio Estimate",round(result$estimate,digits)),byrow=T,ncol=2)
          }
          cap <- "Ratio Estimate : Estimate of variance ratio between groups"
          R2HTML::HTML(Fresult,file=stdout(),caption=cap, innerBorder = 1,align="left") # 분산 비 추정값
        }
        
      }
    }
    #분석종료시간
    R2HTML::HTMLhr(file=stdout())
    R2HTML::HTML(paste("Analysis is finished at ",Sys.time(),". REx : Test for Homogeneity of Variances",sep=""),file=stdout())
    R2HTML::HTMLhr(file=stdout())   
  })
  
  return(html.output)       
}
