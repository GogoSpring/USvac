#' Add Two Numbers
#'
#' This function takes two numbers and returns their sum.
#'
#' @param x A numeric, the first number to add.
#' @param y A numeric, the second number to add.
#' @return A numeric, the sum of x and y.
#' @export
model_fitting <- function(mod,dat) {
  
  extdata<-dat
  State_name<-c("Alabama",              "Alaska",               "Arizona",              "Arkansas",             "California",           "Colorado",            
                "Connecticut",          "Delaware",             "District of Columbia", "Florida",              "Georgia",              "Hawaii",              
                "Idaho",                "Illinois",             "Indiana",              "Iowa",                 "Kansas",               "Kentucky",            
                "Louisiana",            "Maine",                "Maryland",             "Massachusetts",        "Michigan",             "Minnesota",          
                "Mississippi",          "Missouri",             "Montana",              "Nebraska",             "Nevada",               "New Hampshire",       
                "New Jersey",           "New Mexico",           "New York",             "North Carolina",       "North Dakota",         "Ohio",                
                "Oklahoma",             "Oregon",               "Pennsylvania",         "Rhode Island",         "South Carolina",       "South Dakota",        
                "Tennessee",            "Texas",                "Utah",                 "Vermont",              "Virginia",             "Washington",          
                "West Virginia",        "Wisconsin",            "Wyoming" )

  extdata<-extdata[which(extdata$Adjust_Doses_Vax_Pct_7>43.232 & extdata$Adjust_Doses_Vax_Pct_7<73.343),] # remove 5% outliers
  
  lagknots <- logknots(40,2)
  varknots1<-as.vector(quantile(extdata$Adjust_Doses_Vax_Pct_7,c(0.1,0.75,0.9),na.rm=T))
  cb.vac <- crossbasis(extdata$Adjust_Doses_Vax_Pct_7, group=extdata$State,
                       lag=(40), argvar=list(fun="bs",df=5,knots=varknots1), arglag=list(knots=lagknots)) 
  
  lagknots <- logknots(14,1)
  varknots1<-as.vector(quantile(extdata$log_Cases_MA,c(0.1,0.75,0.9),na.rm=T))
  cb.case <- crossbasis(extdata$log_Cases_MA, group=extdata$State,
                        lag=(14), argvar=list(fun="bs",df=5,knots=varknots1), arglag=list(knots=lagknots)) 
  
  # mutation
  all_mutation_df<-data.frame(Date    =NA,
                              State   =NA,
                              onlyBA1 =NA,
                              bothBA12=NA,
                              onlyBA2 =NA)
  
  
  onlyBA1<-c(2,7,12)
  bothBA12<-c(1,3:6,8:11,13:15)
  onlyBA2<-c(16:19)
  
  for (i in 1:51) {
    State<-State_name[i]
    load_file<-paste("D:/PROJECT/US 2/data/sub_mut_freq_output_MA/",State,"_alg_sub_mut_freq.csv",sep = "",collapse="") # I=42 Date different
    state_mutation<-read.csv(load_file)[,-1];ncol(state_mutation)
    colnames(state_mutation)<-c("G339D", "S371L", "S373P", "S375F", "K417N", "N440K", "G446S",
                                "S477N", "T478K", "E484A", "Q493R", "G496S", "Q498R", "N501Y",
                                "Y505H", "S371F", "T376A", "R408S", "D405N", "Date")
    
    if (i==42) {
      state_mutation$Date<-as.Date(state_mutation$Date,"%m/%d/%Y")
    } else state_mutation$Date<-as.Date(state_mutation$Date)
    
    state_mutation1<-state_mutation[which(state_mutation$Date>as.Date("2021-12-10") & state_mutation$Date<as.Date("2022-03-23")),]
    state_mutation2<-state_mutation1[,onlyBA1] 
    state_mutation2<-apply(state_mutation2,1,mean,na.rm=T)  
    
    state_mutation3<-state_mutation1[,bothBA12] 
    state_mutation3<-apply(state_mutation3,1,mean,na.rm=T)  
    
    state_mutation4<-state_mutation1[,onlyBA2] 
    state_mutation4<-apply(state_mutation4,1,mean,na.rm=T)  
    
    state_df<-data.frame(Date    =state_mutation1$Date,
                         State   =rep(State,102),
                         onlyBA1 =as.vector(state_mutation2),
                         bothBA12=as.vector(state_mutation3),
                         onlyBA2 =as.vector(state_mutation4))
    
    all_mutation_df<-rbind(all_mutation_df,state_df)
  }
  all_mutation_df$Date<-as.Date(all_mutation_df$Date,origin="1970-01-01")
  extdata$Date<-as.Date(extdata$Date,"%Y/%m/%d")
  
  df_mutation<-merge(extdata,all_mutation_df[-1,],by=c("Date","State"),all=T)
  df_mutation<-df_mutation[order(df_mutation$State,df_mutation$Date),]
  
  if (mod=="Basic model") {
    
    
    fit <- suppressWarnings(
      lme4::glmer(Daily_hospital_Admission_pct1 ~ (1 | State) + cb.case + Adjust_Doses_Vax_Pct_7 +
                         inpatient_beds   + Humidity_14MA + TEMP_new_14MA  + as.factor(weekday) + as.factor(holiday),  
                       data = extdata, 
                       weights = New_Cases_MA,
                       family = "binomial",
                       control = glmerControl(optimizer = "bobyqa",
                                              optCtrl = list(maxfun = 2e5)),nAGQ = 0)
    )
    a<-summary(fit)
    cat("Basic model fitting coefficient results:\n")
    print(a$coefficients) #model coefficient
    
    cat("Effects of factors on COVID-19 CHR in the Basic model (OR values and 95% CI):\n")
    print(paste("Vaccination coverage (25%) OR: ", round(exp(25*(-3.990e-02)),3), " (", round(exp(25*(-3.990e-02-1.96*7.638e-04)),3), ",", round(exp(25*(-3.990e-02+1.96*7.638e-04)),3), ")", sep="")) # transform to OR based on the model coefficient
    print(paste("Temperature (5 degree) OR: ", round(exp(5*( 1.683e-02)),3), " (", round(exp(5*(1.683e-02-1.96*4.316e-04)),3), ",", round(exp(5*(1.683e-02+1.96*4.316e-04)),3), ")", sep=""))
    print(paste("Relative humidity (5%) OR: ", round(exp(5*(3.475e-03)),3), " (", round(exp(5*(3.475e-03-1.96*2.075e-04)),3), ",", round(exp(5*(3.475e-03+1.96*2.075e-04)),3), ")", sep=""))
    print(paste("Hospital beds (1000) OR: ", round(exp(1000*(2.838e-06)),3), " (", round(exp(1000*(2.838e-06-1.96*5.595e-07)),3), ",", round(exp(1000*(2.838e-06+1.96*5.595e-07)),3), ")", sep=""))
    print(paste("Weekend (Yes vs No) OR: ", round(exp(1*(1.154e-03)),3), " (", round(exp(1*(1.154e-03-1.96*2.157e-03)),3), ",", round(exp(1*(1.154e-03+1.96*2.157e-03)),3), ")", sep=""))
    print(paste("Holiday (Yes vs No) OR: ", round(exp(1*(2.944e-02)),3), " (", round(exp(1*(2.944e-02-1.96*4.559e-03)),3), ",", round(exp(1*(2.944e-02+1.96*4.559e-03)),3), ")", sep=""))
    
  } else if (mod=="Intermediate model") {
    
    end_date<-as.Date("2022-02-01",origin="1970-01-01")
    df_mutation0<-df_mutation[which(df_mutation$Date<end_date),]
    
    lagknots <- logknots(40,1)
    varknots1<-as.vector(quantile(df_mutation0$Adjust_Doses_Vax_Pct_7,c(0.1,0.75,0.9),na.rm=T))
    cb.vac <- crossbasis(df_mutation0$Adjust_Doses_Vax_Pct_7, group=df_mutation0$State,
                         lag=(40), argvar=list(fun="bs",df=5,knots=varknots1), arglag=list(knots=lagknots)) 

    fit <- suppressWarnings(
      glm(onlyBA1 ~ log_Cases_MA+ cb.vac +ns(Date,2),  
               data = df_mutation0,
               weights = New_Cases_MA,
               family=binomial(link="logit"))
    )
    a<-summary(fit)
    
    cat("Intermediate model (BA.1/BA.1.1-associated mutations) fitting coefficient results:\n")
    print(a$coefficients)
    
    pred_vac <- dlnm::crosspred(cb.vac, fit,by=0.1,cen=48,cumul = T)

    ## matrix
    matRRfit<-pred_vac$matRRfit
    matRRfit<-matRRfit[8:268,]
    
    #the matrix that need be plot
    z0 <- as.matrix(matRRfit[,2:ncol(matRRfit)]);max(z0)
    BA1_OR_log2 <- log2(z0)
    #write.csv(BA1_OR_log2,"D:/PROJECT/US/output/Result/BA1_OR_log2.csv")
    
    BA1_OR_log2<-as.matrix(BA1_OR_log2)
    
    lable        <-c(45,50,55,60,65,70)
    resclae_0to1 <-(lable-44)/(70-44)
    x_axis      <-resclae_0to1
    x_axis_lable<-lable
    y_axis      <-seq(0, 1,length.out = 5)
    y_axis_lable<-round(seq(7, 47,length.out = 5),0)
    z_min<-(-1);z_max<-2
    levels<-seq(z_min, z_max, length.out = 31) 
    
    # color
    pal <- c("#8fa4bc","#a5b6c9","#bcc8d7","#d2dbe4","#e9edf2","#ffffff",
             "#f3d1d2", "#e8a3a4", "#dc7477", "#d14649", "#c5181c")
    col1 <- colorRampPalette(pal[1:6])
    col2 <- colorRampPalette(pal[6:11])
    cols <- c(col1(sum(levels < 0)), col2(sum(levels > 0)))
    
    filled.contour(z = BA1_OR_log2,
                   levels= levels,
                   las = 1,
                   plot.title = title(main ="Intermediate model \n(BA.1/BA.1.1-associated mutations, \nOR using log2 scale)",ylab = "Lag (day)" , xlab = "Vaccine coverage (%)", cex.main = 1),
                   plot.axes = list(axis(2,y_axis,y_axis_lable),  
                                    axis(1,x_axis,x_axis_lable)), 
                   col = cols
    )
    
    ####-----------------------------------------------
    end_date<-as.Date("2022-02-01",origin="1970-01-01")
    df_mutation2<-df_mutation[which(df_mutation$Date<end_date),]
    
    lagknots <- logknots(40,1)
    varknots1<-as.vector(quantile(df_mutation2$Adjust_Doses_Vax_Pct_7,c(0.1,0.75,0.9),na.rm=T))
    cb.vac <- crossbasis(df_mutation2$Adjust_Doses_Vax_Pct_7, group=df_mutation2$State,
                         lag=(40), argvar=list(fun="bs",df=4,knots=varknots1), arglag=list(knots=lagknots)) 
    
    fit <- suppressWarnings(
      glm(bothBA12 ~ log_Cases_MA+ cb.vac+ns(Date,2)  ,  #  
               data = df_mutation2,
               weights = New_Cases_MA,
               family=binomial(link="logit"))
    )
    a<-summary(fit)
    
    cat("Intermediate model (shared mutations between BA.1/BA.1.1 and BA.2) fitting coefficient results:\n")
    print(a$coefficients)
    
    pred_vac <- dlnm::crosspred(cb.vac, fit,by=0.1,cen=48,cumul = T)

    ## PLOT matrix
    matRRfit<-pred_vac$matRRfit
    matRRfit<-matRRfit[8:268,]
    
    #the matrix that need be plot
    z0 <- as.matrix(matRRfit[,2:ncol(matRRfit)]);max(z0)
    BA1BA2_OR_log2 <- log2(z0)
    #write.csv(BA1BA2_OR_log2,"D:/PROJECT/US/output/Result/BA1BA2_OR_log2.csv")
    
    BA1BA2_OR_log2<-as.matrix(BA1BA2_OR_log2)
    
    lable        <-c(45,50,55,60,65,70)
    resclae_0to1 <-(lable-44)/(70-44)
    x_axis      <-resclae_0to1
    x_axis_lable<-lable
    y_axis      <-seq(0, 1,length.out = 5)
    y_axis_lable<-round(seq(7, 47,length.out = 5),0)
    
    z_min<-(-1);z_max<-2
    levels<-seq(z_min, z_max, length.out = 31) #70
    
    # color
    cols <- c(col1(sum(levels < 0)), col2(sum(levels > 0)))
    filled.contour(z = BA1BA2_OR_log2,
                   levels= levels,
                   las = 1,
                   plot.title = title(main ="Intermediate model \n(shared mutations between BA.1/BA.1.1 and BA.2, \nOR using log2 scale)",ylab = "Lag (day)" , xlab = "Vaccine coverage (%)", cex.main = 1),
                   plot.axes = list(axis(2,y_axis,y_axis_lable),  
                                    axis(1,x_axis,x_axis_lable)), 
                   col = cols
    )

    
    ####--------------------------------------------------
    df_mutation1<-df_mutation[which(df_mutation$Date>as.Date("2021-11-15")),] 
    
    lagknots <- logknots(40,1)
    varknots1<-as.vector(quantile(df_mutation1$Adjust_Doses_Vax_Pct_7,c(0.1,0.75,0.9),na.rm=T))
    cb.vac <- crossbasis(df_mutation1$Adjust_Doses_Vax_Pct_7, group=df_mutation1$State,
                         lag=(40), argvar=list(fun="bs",df=2,knots=varknots1), arglag=list(knots=lagknots)) 
    
    fit <- suppressWarnings(
      glm(onlyBA2 ~ log_Cases_MA+ cb.vac,  #   
               data = df_mutation1,
               weights = New_Cases_MA,
               family=binomial(link="logit"))
    )
    a<-summary(fit)
    
    cat("Intermediate model (BA.2-associated mutations) fitting coefficient results:\n")
    print(a$coefficients)
    
    pred_vac <- dlnm::crosspred(cb.vac, fit,by=0.1,cen=48,cumul = T)#model.link="log",
    
    ## PLOT matrix
    matRRfit<-pred_vac$matRRfit
    matRRfit<-matRRfit[8:268,]
    
    #the matrix that need be plot
    z0 <- as.matrix(matRRfit[,2:ncol(matRRfit)]);max(z0)
    BA2_OR_log2 <- log2(z0)
    #write.csv(BA2_OR_log2,"D:/PROJECT/US/output/Result/BA2_OR_log2.csv")
    
    BA2_OR_log2<-as.matrix(BA2_OR_log2)
    
    lable        <-c(45,50,55,60,65,70)
    resclae_0to1 <-(lable-44)/(70-44)
    x_axis      <-resclae_0to1
    x_axis_lable<-lable
    y_axis      <-seq(0, 1,length.out = 5)
    y_axis_lable<-round(seq(7, 47,length.out = 5),0)
    
    BA2_OR_log2[BA2_OR_log2 < (-4)] <- (-4)
    z_min<-(-4);z_max<-4
    levels<-seq(z_min, z_max, length.out = 35) 
    
    # color
    cols <- c(col1(sum(levels < 0)), col2(sum(levels > 0)))
    
    filled.contour(z = BA2_OR_log2,
                   levels= levels,
                   las = 1,
                   plot.title = title(main ="Intermediate model \n(BA.2-associated mutations, \nOR using log2 scale)",ylab = "Lag (day)" , xlab = "Vaccine coverage (%)", cex.main = 1),
                   plot.axes = list(axis(2,y_axis,y_axis_lable),  
                                    axis(1,x_axis,x_axis_lable)), 
                   col = cols
    )
    #dev.off()

    
  }else if (mod=="Mediation model") {
    
    ##### model2: mutation have effect on hospitalization rate
    lagknots <- logknots(14,1)
    varknots1<-as.vector(quantile(df_mutation$log_Cases_MA,c(0.1,0.75,0.9),na.rm=T))
    cb.case <- crossbasis(df_mutation$log_Cases_MA, group=df_mutation$State,
                          lag=(14), argvar=list(fun="bs",df=5,knots=varknots1), arglag=list(knots=lagknots))
    
    fit <- suppressWarnings(
      lme4::glmer(Daily_hospital_Admission_pct1 ~ (1 | State) + cb.case + Adjust_Doses_Vax_Pct_7 + onlyBA1 + onlyBA2 + bothBA12 +
                         inpatient_beds  + Humidity_14MA + TEMP_new_14MA  + as.factor(weekday) + as.factor(holiday),  
                       data = df_mutation, #
                       weights = New_Cases_MA,
                       family = "binomial",
                       control = glmerControl(optimizer = "bobyqa",
                                              optCtrl = list(maxfun = 2e5)),nAGQ = 0)
    )
    a<-summary(fit)
    
    cat("Mediation model fitting coefficient results:\n")
    
    cat("Effects of factors on COVID-19 CHR in the Mediation model (OR values and 95% CI):\n")
    print(paste("Vaccination coverage (25%) OR: ", round(exp(25*(-6.357e-03)),3), " (", round(exp(25*(-6.357e-03-1.96*2.527e-03)),3), ",", round(exp(25*(-6.357e-03+1.96*2.527e-03)),3), ")", sep="")) # transform to OR based on the coefficient
    print(paste("Mean proportion of BA.1/BA.1.1-associated mutations (20%) OR: ", round(exp(0.2*(-6.982e-01)),3), " (", round(exp(0.2*(-6.982e-01-1.96*9.042e-02)),3), ",", round(exp(0.2*(-6.982e-01+1.96*9.042e-02)),3), ")", sep=""))
    print(paste("Mean proportion of BA.2-associated mutations (20%) OR: ", round(exp(0.1*(-1.413e+00)),3), " (", round(exp(0.1*(-1.413e+00-1.96*9.362e-02)),3), ",", round(exp(0.1*(-1.413e+00+1.96*9.362e-02)),3), ")", sep=""))
    print(paste("Mean proportion of shared mutations between BA.1/BA.1.1 and BA.2 (20%) OR: ", round(exp(0.1*(-5.653e-02)),3), " (", round(exp(0.1*(-5.653e-02-1.96*1.002e-01)),3), ",", round(exp(0.1*(-5.653e-02+1.96*1.002e-01)),3), ")", sep=""))
    print(paste("Temperature (5 degree) OR: ", round(exp(5*(1.612e-02)),3), " (", round(exp(5*(1.612e-02-1.96*6.372e-04)),3), ",", round(exp(5*(1.612e-02+1.96*6.372e-04)),3), ")", sep=""))
    print(paste("Relative humidity (5%) OR: ", round(exp(5*(2.563e-03)),3), " (", round(exp(5*(2.563e-03-1.96*2.481e-04)),3), ",", round(exp(5*(2.563e-03+1.96*2.481e-04)),3), ")", sep=""))
    print(paste("Hospital beds (1000) OR: ", round(exp(1000*(2.500e-06)),3), " (", round(exp(1000*(2.500e-06-1.96*6.449e-07)),3), ",", round(exp(1000*(2.500e-06+1.96*6.449e-07)),3), ")", sep=""))
    print(paste("Weekend (Yes vs No) OR: ", round(exp(1*(2.433e-03)),3), " (", round(exp(1*(2.433e-03-1.96*2.315e-03)),3), ",", round(exp(1*(2.433e-03+1.96*2.315e-03)),3), ")", sep=""))
    print(paste("Holiday (Yes vs No) OR: ", round(exp(1*(2.031e-02)),3), " (", round(exp(1*(2.031e-02-1.96*4.880e-03)),3), ",", round(exp(1*(2.031e-02+1.96*4.880e-03)),3), ")", sep=""))
    
  }else {
    
    print("No such model")
    
  } 
  
}
