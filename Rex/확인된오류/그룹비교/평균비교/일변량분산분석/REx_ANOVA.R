### 일변량분산분석

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
		qual_var <- unique(unlist(strsplit(vars,"[\\+\\:]")))	# 수정

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
		
		# explanatory variable type, 수정
		if(is.null(qual_var)) {
			error_msg	<- c(error_msg, '\a Error : At least 1 independent variable should be selected. Analysis has been stopped.') 
		} else {
			qll <- c()
			for(ql in qual_var) {
				if (length(intersect(c("character", "factor"), class(dataset[[ql]]))) < 1) {
					qll <- c(qll,ql)
				}
				dataset[,ql] <- as.factor(dataset[,ql]) ;
			}
			if(length(qll)>0)	warn_msg	<- paste0("\a Warning : The variable '", paste(qll,collapse=", "), "' seems not to be categorical: the variable is coerced to the factor, and the fit may be unstable bacause of too many parameters.")
		}
		
		if (posthoc) { 
			 if (!any(grepl(posthoc_var, qual_var))) error_msg	<- c(error_msg, '\a Error : The variable requested for post-hoc analysis is not in the main analysis variable.') 
		
		} ;
		
		formula <- as.formula(paste0(res_var, ' ~ ', paste(c(vars),collapse=' + '))) ;
		if (missing(dataset)) dataset <- environment(formula) ;

		newdat <- dataset[,c(res_var,qual_var),drop=F]
		newdat <- newdat[complete.cases(newdat),,drop=F]
		newdat_2 <- model.matrix(formula, data=newdat)
		newdat_2 <- cbind(newdat_2, newdat[,res_var]) ;
		if (rankMatrix(newdat_2) < ncol(newdat_2)) error_msg 	<- c(error_msg, "\a Error : Linear dependency between columns of the design matrix (including the intercept) detected. Please check the values of indepedent/response variables.")
		
		if(length(error_msg) > 0) {
			R2HTML::HTML(R2HTML::as.title("Warnings"),HR=2,file=stdout())	# 수정
			for (msg in error_msg) {
				R2HTML::HTML(msg, file=stdout()) ;
			}
		} else {
		
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
			data_str	<- matrix(c("Number of observations", "Number of total variables", "Number of used variables", nrow(dataset),ncol(dataset), summ_var), ncol=2) ;#수정
			
			R2HTML::HTML(R2HTML::as.title("Data Structure"), HR=2, file=stdout()) ;
			R2HTML::HTML(data_str, file=stdout(), row.names=FALSE, col.names=TRUE, innerBorder = 1, align="left") ;

			# Model formula
			fit_print	<- matrix(c(variables[2], 
						paste0(variables[2], " = ", variables[3], " + error")
						)) ;
			fit_print	<- cbind(c("Response variable", "Model"), fit_print) ;
			colnames(fit_print)	<- c("Category", "Value") ;
			
			R2HTML::HTML(R2HTML::as.title("Analysis Description"), HR=2, file=stdout()) ;
			R2HTML::HTML(fit_print, file=stdout(), row.names=FALSE, col.names=TRUE, innerBorder = 1, align="left") ;
			if (exists("warn_msg")) R2HTML::HTML(warn_msg, file=stdout()) ;

			R2HTML::HTML(R2HTML::as.title("Results of Analysis of Variance (ANOVA)"), HR=2, file=stdout()) ;
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

			R2HTML::HTML(R2HTML::as.title("Model Fitness Measurements"), HR=3, file=stdout()) ;
			R2HTML::HTML(fitness_print, file=stdout(), na="", innerBorder = 1, align="left") ;
		
			# ANOVA table
			R2HTML::HTML(R2HTML::as.title("ANOVA Table"), HR=3, file=stdout()) ;
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
				
				R2HTML::HTML(R2HTML::as.title(paste0("Post-hoc Analysis for variable '", posthoc_var,"'")), HR=3, file=stdout()) ;
			
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
					
					CD <- post_tmp$statistics$CriticalDifference	# Critical difference, 수정
					if(!is.null(CD)){
						diff_crit	<- format(round(post_tmp$statistics$CriticalDifference, num_digit), nsmall=num_digit) ;
						notice_msg	<- paste0("\a Notice : groups with the difference larger than ", diff_crit, " is considered to be significantly different. (Critical difference is provided only when the sample sizes are same for each group.)") ;
					}
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
				
				post_table	<- cbind(posthoc_method,rownames(post_table), post_table) ;
				names(post_table)[1:2]	<- c("Method","Comparison");
				
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

			if(exists("ANOVA_output")){
				temp.0 <- ANOVA_output
				temp <- as.data.frame(matrix(NA,nrow(dataset),ncol(temp.0)))
				colnames(temp) <- colnames(temp.0)
				row.dataset <- rownames(dataset)
				row.output <- rownames(temp.0)
				for(i in row.dataset) if(i%in%row.output) temp[row.dataset==i,] <- temp.0[row.output==i,]
				ANOVA_output <- temp
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
