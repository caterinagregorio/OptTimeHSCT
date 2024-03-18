define_params <- function(coefs,uncertainty='none',B=0,knots,hsct_aml='no') {
  coefhsct <- as.matrix(1)
  colnames(coefhsct) <- "Tstart"
  
  if(uncertainty=='PSA'){
    coefhsct <- as.matrix(rep(1,B))
    colnames(coefhsct) <- "Tstart"
  }
  
  
  transmod_params <- params_surv_list(
    
    # 1. MDS:AML        (1:2) 
    params_surv(
      coefs = coefs[[1]],
      dist = "survspline",
      aux = list(
        knots = knots[[1]],
        scale = "log_cumhazard",
        timescale = "log"
      )
    ),
    
    # 2. MDS:HSCT       (1:3) deterministic, defined by the treatment policy
    params_surv(coefs = list(est = coefhsct),
                dist = "fixed"),
    
    
    # 3. MDS:DEATH      (1:5) 
    params_surv(
      coefs = coefs[[2]],
      dist = "survspline",
      aux = list(
        knots = knots[[2]],
        scale = "log_cumhazard",
        timescale = "log"
      )
    ),
    # # 4. AML:HSCT       (2:3) deterministic, defined by the treatment policy
    # params_surv(coefs = list(est = coefhsctaml),
    #             dist = "fixed"),
    
    
    # 4. AML:DEATH      (2:5)
    params_surv(
      coefs = coefs[[3]],
      dist = "survspline",
      aux = list(
        knots = knots[[3]],
        scale = "log_cumhazard",
        timescale = "log"
      )
    ),
    # 5. HSCT:RELAPSE   (3:4)
    params_surv(
      coefs = coefs[[4]],
      dist = "survspline",
      aux = list(
        knots = knots[[4]],
        scale = "log_cumhazard",
        timescale = "log"
      )
    ),
    # 6. HSCT:DEATH     (3:5)
    params_surv(
      coefs = coefs[[5]],
      dist = "survspline",
      aux = list(
        knots = knots[[5]],
        scale = "log_cumhazard",
        timescale = "log"
      )
    ),
    # 7. RELAPSE:DEATH  (4:5)
    params_surv(
      coefs = coefs[[6]],
      dist = "survspline",
      aux = list(
        knots = knots[[6]],
        scale = "log_cumhazard",
        timescale = "log"
      )
    )
    
    
    
  )
  
  return(transmod_params)
}




fix_par <- function(index,fit){
  
  if(index==1){
    pars <- t(as.matrix(fit$coefficients[grepl("gamma0",names(fit$coefficients))|!grepl("gamma",names(fit$coefficients))  ]))
  }else{
    pars <- t(as.matrix(fit$coefficients[grepl(paste0("gamma",index-1),names(fit$coefficients))]))
  }
  
  if(ncol(pars)==1)colnames(pars) <- "intercept"
  if(ncol(pars)!=1){
    if(index==1){
      colnames(pars) <- c(
        "intercept",
        names(fit$coefficients)[!grepl("gamma",names(fit$coefficients))])
    }else{
      colnames(pars) <-c(
        "intercept",
        gsub('[()]',"",gsub(paste0("gamma",index-1),"",names(fit$coefficients)[grepl(paste0("gamma",index-1),names(fit$coefficients))][-1]))) 
    }
  }
  
  return(pars)
}



fix_par_boot <- function(index,coefboot){
  
  if(index==1){
    pars <- as.matrix(coefboot[,grepl("gamma0",colnames(coefboot))|!grepl("gamma",colnames(coefboot))  ])
  }else{
    pars <- as.matrix(coefboot[,grepl(paste0("gamma",index-1),colnames(coefboot))])
  }
  
  if(ncol(pars)==1)colnames(pars) <- "intercept"
  if(ncol(pars)!=1){
    if(index==1){
      colnames(pars) <- c(
        "intercept",
        colnames(coefboot)[!grepl("gamma",colnames(coefboot))])
    }else{
      col <- c("intercept",
      gsub('[()]',"",gsub(paste0("gamma",index-1),"",colnames(coefboot)[grepl(paste0("gamma",index-1),colnames(coefboot))][-1]))) 
      
      colnames(pars) <-col
    }
  }
  
  return(pars)
}

coef_to_list <- function(fit){
  
  
  
  coeflist <- lapply(fit$basepars, fix_par,fit=fit)
  names(coeflist) <- paste0("gamma",0:(length(coeflist)-1))
  return(coeflist)
}


coef_to_list_boot <- function(fit,fitorig){
  

  coeflist <- lapply(fitorig$basepars, fix_par_boot,coefboot=fit)
  names(coeflist) <- paste0("gamma",0:(length(coeflist)-1))
  return(coeflist)
}


microsim <- function(thsct,
                      fits,
                      uncertainty = "none",
                      f,
                      nsim,
                      thor,
                      B,
                      prob=F,
                      covs,
                      hsct_aml='no') {
  if(prob){
    transmodprob <- list()
  }
  
  rmst <- list()
  
  
  
  dat <-
    as.data.frame(model.matrix(
      as.formula(paste0("~",f)),
      data = covs
    ))
  
  colnames(dat)[1] <-"intercept"
  
  
  tmat <-matrix(c(NA,1,2,NA,3,
                  NA,NA,NA,NA,4,
                  NA, NA,NA, 5,6,
                  NA, NA, NA,NA,7,
                  NA,NA,NA,NA,NA),
                byrow = T,
                ncol = 5)
  
  
  rownames(tmat) <- colnames(tmat) <-  c("dnh","aml", "hsct","relapse", "dead")
  
  
  
  
  
  if (uncertainty == "none") {
    coefs <- list()
    # 1. MDS:AML        (1:2) 
    coefs[[1]]<-  coef_to_list(fits[[1]])
    # 3. MDS:DEATH      (1:5) 
    coefs[[2]]<-  coef_to_list(fits[[2]])
    # 4. AML:DEATH      (2:5)
    coefs[[3]]<-  coef_to_list(fits[[3]])
    # 5. HSCT:RELAPSE   (3:4)
    coefs[[4]]<-  coef_to_list(fits[[4]])
    # 6. HSCT:DEATH     (3:5)
    coefs[[5]]<-  coef_to_list(fits[[5]])
    # 7. RELAPSE:DEATH  (4:5)
    coefs[[6]]<-  coef_to_list(fits[[6]])
    
    
    
  }
  
  
  
  if (uncertainty == "PSA") {
    
    fboot <- lapply(1:6,function(x) normboot.flexsurvreg(fits[[x]],B=B,raw=T))
    
    
    coefs <- list()
    # 1. MDS:AML        (1:2) 
    coefs[[1]]<-  coef_to_list_boot(fboot[[1]],fits[[1]])
    # 3. MDS:DEATH      (1:5) 
    coefs[[2]]<-  coef_to_list_boot(fboot[[2]],fits[[2]])
    # 4. AML:DEATH      (2:5)
    coefs[[3]]<-  coef_to_list_boot(fboot[[3]],fits[[3]])
    # 5. HSCT:RELAPSE   (3:4)
    coefs[[4]]<-  coef_to_list_boot(fboot[[4]],fits[[4]])
    # 6. HSCT:DEATH     (3:5)
    coefs[[5]]<-  coef_to_list_boot(fboot[[5]],fits[[5]])
    # 7. RELAPSE:DEATH  (4:5)
    coefs[[6]]<-  coef_to_list_boot(fboot[[6]],fits[[6]])
    
    
    
  }
  
  
  knots <- lapply(fits,function(x)x$knots)
  
  transmod_params <- define_params(
    coefs=coefs,
    knots=knots,
  uncertainty,
  B,
  hsct_aml = hsct_aml)
  
  res <- data.table::data.table()
  
  if(nsim>200){
    
    nsim2 <-round(nsim/10,0) 
    
    for (pat in 1:10) {
      print(pat)
      # dat_rep <- dat[rep(which(dat$patient_profile == pat), each = nsim),]
      # Patients
      patients <- dat %>% group_by_all() %>% 
        expand_grid(id=1:nsim2) %>% distinct()
      
      patients$patient_id=1:nrow(patients)
      
      print(nrow(patients))
      # Treatment strategies
      strategies <- data.frame(Tstart = seq(0, 30*48,by=60),
                               strategy_id=1:25)
      age <- patients %>% distinct(patient_id,.keep_all = T) %$%age
      print(length(age))
      # Input data
      hesim_dat <- hesim_data(strategies = strategies,
                              patients = patients)
      
      transmod_data <-hesim::expand(hesim_dat, by = c("strategies", "patients"))
      
      # transmod_data <- bind_cols(transmod_data, dat_rep)
      
      
      
      IndivCtstmTrans_obj <-
        create_IndivCtstmTrans(
          transmod_params,
          input_data = transmod_data,
          trans_mat = tmat,
          clock = "reset",
          start_age = age*365# in days
        )
      
      transmodsim <-
        IndivCtstmTrans_obj$sim_disease(max_t = thor,
                                        max_age=100*365,progress=1)
      
      gc()
      patients <- patients %>% select(patient_id,grp_id)
      transmodsim %<>%select(-grp_id) %>% left_join(patients) %>% select(-patient_id)
      res <- rbind(transmodsim,res)
      
    }}else{
    
      patients <- dat %>% group_by_all() %>% 
      expand_grid(id=1:nsim) %>% distinct()
    
    patients$patient_id=1:nrow(patients)
    
    print(nrow(patients))
    # Treatment strategies
    strategies <- data.frame(Tstart = seq(0, 30*48,by=60),
                             strategy_id=1:25)
    age <- patients %>% distinct(patient_id,.keep_all = T) %$%age
    print(length(age))
    # Input data
    hesim_dat <- hesim_data(strategies = strategies,
                            patients = patients)
    
    transmod_data <-hesim::expand(hesim_dat, by = c("strategies", "patients"))
    
    # transmod_data <- bind_cols(transmod_data, dat_rep)
    
    
    
    IndivCtstmTrans_obj <-
      create_IndivCtstmTrans(
        transmod_params,
        input_data = transmod_data,
        trans_mat = tmat,
        clock = "reset",
        start_age = age*365# in days
      )
    
    transmodsim <-
      IndivCtstmTrans_obj$sim_disease(max_t = thor,
                                      max_age=100*365,progress=1)
    
    gc()
    patients <- patients %>% select(patient_id,grp_id)
    transmodsim %<>%select(-grp_id) %>% left_join(patients) 
    res <- transmodsim
    
  }
 
    suppressMessages(
      if(uncertainty!='PSA'){
        rmst <- res %>%
          group_by(strategy_id, grp_id) %>%
          filter(is.finite(time_stop)) %>% 
          group_by(strategy_id, grp_id,patient_id) %>%
          # c("dnh","aml", "hsct","relapse", "dead")
          mutate(weight=case_when(from==2~0.85,
                                  from>=3~0.90,
                                  TRUE~1)) %>% 
          summarise(util=sum((time_stop-time_start)*weight),
                    util2=sum(time_stop-time_start)) %>% 
          group_by(strategy_id, grp_id) %>%
          summarise(
            rmst = mean(util)/30,
            rmst2 = mean(util2)/30)
      }else{
        rmst <- res %>%
          group_by(strategy_id, grp_id,sample) %>%
          filter(is.finite(time_stop)) %>% 
          group_by(strategy_id, grp_id,patient_id,sample) %>%
          mutate(weight=ifelse(from==2|from==4,0.5,1)) %>% 
          summarise(util=sum((time_stop-time_start)*weight),
                    util2=sum(time_stop-time_start)) %>% 
          group_by(strategy_id, grp_id,sample) %>%
          summarise(
            rmst = mean(util)/30,
            rmst2 = mean(util2)/30)
      }
    )
    
    
    
    
    if(prob){
      transmodprob<-IndivCtstmTrans_obj$sim_stateprobs(disprog = transmodsim,
                                                               t = seq(0, thor, by =
                                                                         0.5))
      
    
    }
    
    
  if(prob){
    rmst <- as.data.frame(do.call('rbind', rmst))
    transmodprob <- as.data.frame(do.call('rbind', transmodprob))
    return(list(prob = transmodprob,rmst=rmst))
  }
  
  else{
    # rmst <- as.data.frame(do.call('rbind', rmst))
    return(rmst)
  }
}


microsim_paths <- function(thsct,
                      fits,
                      uncertainty = "none",
                      f,
                      nsim,
                      thor,
                      B,
                      covs,
                      hsct_aml='no') {
  
  
  rmst <- list()
  
  
  
  dat <-
    as.data.frame(model.matrix(
      as.formula(paste0("~",f)),
      data = covs
    ))
  
  colnames(dat)[1] <-"intercept"
  
  
  tmat <-matrix(c(NA,1,2,NA,3,
                  NA,NA,NA,NA,4,
                  NA, NA,NA, 5,6,
                  NA, NA, NA,NA,7,
                  NA,NA,NA,NA,NA),
                byrow = T,
                ncol = 5)
  
  
  rownames(tmat) <- colnames(tmat) <-  c("dnh","aml", "hsct","relapse", "dead")
  
  
  
  
  
  if (uncertainty == "none") {
    coefs <- list()
    # 1. MDS:AML        (1:2) 
    coefs[[1]]<-  coef_to_list(fits[[1]])
    # 3. MDS:DEATH      (1:5) 
    coefs[[2]]<-  coef_to_list(fits[[2]])
    # 4. AML:DEATH      (2:5)
    coefs[[3]]<-  coef_to_list(fits[[3]])
    # 5. HSCT:RELAPSE   (3:4)
    coefs[[4]]<-  coef_to_list(fits[[4]])
    # 6. HSCT:DEATH     (3:5)
    coefs[[5]]<-  coef_to_list(fits[[5]])
    # 7. RELAPSE:DEATH  (4:5)
    coefs[[6]]<-  coef_to_list(fits[[6]])
    
    
    
  }
  
  
  
  if (uncertainty == "PSA") {
    
    fboot <- lapply(1:6,function(x) normboot.flexsurvreg(fits[[x]],B=B,raw=T))
    
    
    coefs <- list()
    # 1. MDS:AML        (1:2) 
    coefs[[1]]<-  coef_to_list_boot(fboot[[1]],fits[[1]])
    # 3. MDS:DEATH      (1:5) 
    coefs[[2]]<-  coef_to_list_boot(fboot[[2]],fits[[2]])
    # 4. AML:DEATH      (2:5)
    coefs[[3]]<-  coef_to_list_boot(fboot[[3]],fits[[3]])
    # 5. HSCT:RELAPSE   (3:4)
    coefs[[4]]<-  coef_to_list_boot(fboot[[4]],fits[[4]])
    # 6. HSCT:DEATH     (3:5)
    coefs[[5]]<-  coef_to_list_boot(fboot[[5]],fits[[5]])
    # 7. RELAPSE:DEATH  (4:5)
    coefs[[6]]<-  coef_to_list_boot(fboot[[6]],fits[[6]])
    
    
    
  }
  
  
  knots <- lapply(fits,function(x)x$knots)
  
  transmod_params <- define_params(
    coefs=coefs,
    knots=knots,
    uncertainty,
    B,
    hsct_aml = hsct_aml)
  
  
      # dat_rep <- dat[rep(which(dat$patient_profile == pat), each = nsim),]
      # Patients
      patients <- dat %>% group_by_all() %>% 
        expand_grid(id=1:nsim) %>% distinct()
      
      patients$patient_id=1:nrow(patients)
      
      print(nrow(patients))
      # Treatment strategies
      strategies <- data.frame(Tstart = seq(0, 30*48,by=60),
                               strategy_id=1:25)
      age <- patients %>% distinct(patient_id,.keep_all = T) %$%age
      print(length(age))
      # Input data
      hesim_dat <- hesim_data(strategies = strategies,
                              patients = patients)
      
      transmod_data <-hesim::expand(hesim_dat, by = c("strategies", "patients"))
      
      # transmod_data <- bind_cols(transmod_data, dat_rep)
      
      
      
      IndivCtstmTrans_obj <-
        create_IndivCtstmTrans(
          transmod_params,
          input_data = transmod_data,
          trans_mat = tmat,
          clock = "reset",
          start_age = age*365# in days
        )
      
      transmodsim <-IndivCtstmTrans_obj$sim_disease(max_t = thor,
                                        max_age=100*365,progress=1)
      
      gc()
      patients <- patients %>% select(patient_id,grp_id)
      transmodsim %<>%select(-grp_id) %>% left_join(patients) 
     
      
    
  
  
    return(transmodsim)
  
}



microsim_probs <- function(thsct,
                           fits,
                           uncertainty = "none",
                           f,
                           nsim,
                           thor,
                           B,
                           covs,
                           hsct_aml='no') {
  
  
  rmst <- list()
  
  
  
  dat <-
    as.data.frame(model.matrix(
      as.formula(paste0("~",f)),
      data = covs
    ))
  
  colnames(dat)[1] <-"intercept"
  
  
  tmat <-matrix(c(NA,1,2,NA,3,
                  NA,NA,NA,NA,4,
                  NA, NA,NA, 5,6,
                  NA, NA, NA,NA,7,
                  NA,NA,NA,NA,NA),
                byrow = T,
                ncol = 5)
  
  
  rownames(tmat) <- colnames(tmat) <-  c("dnh","aml", "hsct","relapse", "dead")
  
  
  
  
  
  if (uncertainty == "none") {
    coefs <- list()
    # 1. MDS:AML        (1:2) 
    coefs[[1]]<-  coef_to_list(fits[[1]])
    # 3. MDS:DEATH      (1:5) 
    coefs[[2]]<-  coef_to_list(fits[[2]])
    # 4. AML:DEATH      (2:5)
    coefs[[3]]<-  coef_to_list(fits[[3]])
    # 5. HSCT:RELAPSE   (3:4)
    coefs[[4]]<-  coef_to_list(fits[[4]])
    # 6. HSCT:DEATH     (3:5)
    coefs[[5]]<-  coef_to_list(fits[[5]])
    # 7. RELAPSE:DEATH  (4:5)
    coefs[[6]]<-  coef_to_list(fits[[6]])
    
    
    
  }
  
  
  
  if (uncertainty == "PSA") {
    
    fboot <- lapply(1:6,function(x) normboot.flexsurvreg(fits[[x]],B=B,raw=T))
    
    
    coefs <- list()
    # 1. MDS:AML        (1:2) 
    coefs[[1]]<-  coef_to_list_boot(fboot[[1]],fits[[1]])
    # 3. MDS:DEATH      (1:5) 
    coefs[[2]]<-  coef_to_list_boot(fboot[[2]],fits[[2]])
    # 4. AML:DEATH      (2:5)
    coefs[[3]]<-  coef_to_list_boot(fboot[[3]],fits[[3]])
    # 5. HSCT:RELAPSE   (3:4)
    coefs[[4]]<-  coef_to_list_boot(fboot[[4]],fits[[4]])
    # 6. HSCT:DEATH     (3:5)
    coefs[[5]]<-  coef_to_list_boot(fboot[[5]],fits[[5]])
    # 7. RELAPSE:DEATH  (4:5)
    coefs[[6]]<-  coef_to_list_boot(fboot[[6]],fits[[6]])
    
    
    
  }
  
  
  knots <- lapply(fits,function(x)x$knots)
  
  transmod_params <- define_params(
    coefs=coefs,
    knots=knots,
    uncertainty,
    B,
    hsct_aml = hsct_aml)
  
  
  # dat_rep <- dat[rep(which(dat$patient_profile == pat), each = nsim),]
  # Patients
  patients <- dat %>% group_by_all() %>% 
    expand_grid(id=1:nsim) %>% distinct()
  
  patients$patient_id=1:nrow(patients)
  
  print(nrow(patients))
  # Treatment strategies
  strategies <- data.frame(Tstart = seq(0, 30*48,by=60),
                           strategy_id=1:25)
  age <- patients %>% distinct(patient_id,.keep_all = T) %$%age
  print(length(age))
  # Input data
  hesim_dat <- hesim_data(strategies = strategies,
                          patients = patients)
  
  transmod_data <-hesim::expand(hesim_dat, by = c("strategies", "patients"))
  
  # transmod_data <- bind_cols(transmod_data, dat_rep)
  
  
  
  IndivCtstmTrans_obj <-
    create_IndivCtstmTrans(
      transmod_params,
      input_data = transmod_data,
      trans_mat = tmat,
      clock = "reset",
      start_age = age*365# in days
    )
  
  transmodsim <-IndivCtstmTrans_obj$sim_disease(max_t = thor,
                                                max_age=100*365,progress=1)
  
  probs <- IndivCtstmTrans_obj$sim_stateprobs(disprog = transmodsim,
                                              t = seq(0, thor, by =
                                                        60))
  # gc()
  # patients <- patients %>% select(patient_id,grp_id)
  # probs %<>%select(-grp_id) %>% left_join(patients) 
  # 
  
  
  
  
  return(probs)
  
}

