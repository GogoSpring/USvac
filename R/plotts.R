
#' Add Two Numbers
#'
#' This function takes two numbers and returns their sum.
#'
#' @param x A numeric, the first number to add.
#' @param y A numeric, the second number to add.
#' @return A numeric, the sum of x and y.
#' @export
plot.ts <- function(dat) {

  cases_df1<-dat
  ##########################
  # daily cases
  ##########################
  columun<-"New_Cases_MA_per1000000"
  G1<-IQR_columun_df(cases_df1[cases_df1$State %in% Group1,],columun)
  G2<-IQR_columun_df(cases_df1[cases_df1$State %in% Group2,],columun)
  colnames(G1)<-c("Date","columunG1","columunG1_low","columunG1_up")
  colnames(G2)<-c("Date","columunG2","columunG2_low","columunG2_up")
  All_df<-merge(G1,G2)
  
  
  #---PLOT
  time = 1:nrow(All_df)
  
  Date<- All_df$Date
  actual_case<-All_df$columunG1
  
  plot(time,actual_case,pch=16,xaxt = "n",col="white",yaxt = "n", ylim = c(0,1000),
       xlab= "", ylab="Daily reported cases per million",main=columun,cex.main=1.2,cex.lab=1.2) 
  
  polygon(c(time, rev(time)), c(All_df$columunG1_low[time], rev(All_df$columunG1_up[time])), density = NA,
          col = rgb(244,220,216, 190, maxColorValue=255), border = NA)  
  
  polygon(c(time, rev(time)), c(All_df$columunG2_low[time], rev(All_df$columunG2_up[time])), density = NA,
          col = rgb(204,228,240, 150, maxColorValue=255), border = NA)  
  lines(All_df$columunG1, lwd = 1, col ="#C54F35")  
  lines(All_df$columunG2, lwd = 1, col ="#0473B6") 
  
  legend("top", legend=c("Group1: Low vaccine coverage", "Group2: High vaccine coverage"),
         col=c("#C54F35","#0473B6"), lty=1:1, cex=0.8)
  
  axis(1,at =seq(1,nrow(All_df),10),labels = Date[seq(1,nrow(All_df),10)],las = 2,cex.axis=1.1)
  axis(2,las = 1,cex.axis = 1.2,cex.axis=1.1) 
  
  
  ####################################
  # cumulative cases —— box plot
  ####################################
  df_cum<-aggregate(New_Cases~State,data=df1,sum,na.action = na.omit)
  df_cum1<-merge(df_cum,Pop,by="State")
  df_cum1$cases_per1000000<-(df_cum1$New_Cases/df_cum1$Pop)*1000000
  df_cum1$group[which(df_cum1$State %in% Group1)]<-"Low"
  df_cum1$group[which(df_cum1$State %in% Group2)]<-"High"
  
  df_cum1$group <- as.factor(df_cum1$group)
  df_cum1$group <- factor(df_cum1$group,levels = c("Low", "High"))

  df_cum2<-df_cum1
  
  g1<-df_cum2$cases_per100000[which(df_cum2$group=="Low")]
  g2<-df_cum2$cases_per100000[which(df_cum2$group=="High")]
  t.test(g1,g2)  #p<0.05 有差異
  
  #--- Boxplot1  for cumulative cases
  boxplot(df_cum2$cases_per1000000 ~ df_cum2$group, col=c("#F4DCD8","#CCE4F0"),
          xlab = "Vaccine coverage level",
          ylab="Cumulative cases per million",
          outline = FALSE)
  # Add df_omicron_cum2 points
  mylevels <- c("Low", "High")
  levelProportions <- as.vector(table(as.vector(df_cum2$group))/nrow(df_cum2))
  for(i in 1:length(mylevels)){
    #i<-1
    thislevel <- mylevels[i]
    thisvalues <- df_cum2[df_cum2$group==thislevel, "cases_per1000000"]
    # take the x-axis indices and add a jitter, proportional to the N in each level
    myjitter <- jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/2)
    points(myjitter, thisvalues, pch=20, col=rgb(0,0,0,.9)) 
  }
  
}
