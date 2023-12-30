#' Add Two Numbers
#'
#' This function takes two numbers and returns their sum.
#'
#' @param x A numeric, the first number to add.
#' @param y A numeric, the second number to add.
#' @return A numeric, the sum of x and y.
#' @export
mobility <- function(dat) {
  
  Mob_df<-dat
  head(Mob_df)
  Mob_df$Date<-as.Date(Mob_df$Date)
  # group comparison (each state have a value of mobility)-------------------------------
  # t-test
  Mob_mean<-aggregate(transit_stations_percent_change_from_baseline~State,data=Mob_df,mean,na.action = na.omit)
  Mob_mean$group[which(Mob_mean$State %in% Group1)]<-"Low"
  Mob_mean$group[which(Mob_mean$State %in% Group2)]<-"High"
  
  Mob_mean$group <- as.factor(Mob_mean$group)
  Mob_mean$group <- factor(Mob_mean$group,levels = c("Low", "High"))
  
  g1<-Mob_mean$transit_stations_percent_change_from_baseline[which(Mob_mean$group=="Low")]
  g2<-Mob_mean$transit_stations_percent_change_from_baseline[which(Mob_mean$group=="High")]
  t.test(g1,g2) 
  
  boxplot(Mob_mean$transit_stations_percent_change_from_baseline ~ Mob_mean$group, col=c("#F4DCD8","#CCE4F0"),
          xlab = "Group",
          ylab="Transit stations percent change from baseline (%)",
          outline = FALSE)
  # Add df_omicron_cum2 points
  mylevels <- c("Low", "High")
  levelProportions <- as.vector(table(as.vector(Mob_mean$group))/nrow(Mob_mean))
  for(i in 1:length(mylevels)){
    #i<-1
    thislevel <- mylevels[i]
    thisvalues <- Mob_mean[Mob_mean$group==thislevel, "transit_stations_percent_change_from_baseline"]
    # take the x-axis indices and add a jitter, proportional to the N in each level
    myjitter <- jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/2)
    points(myjitter, thisvalues, pch=20, col=rgb(0,0,0,.9)) 
  }
  
  #######################################################################################
  # 3. ## plot IQR 
  #######################################################################################
  head(Mob_df)
  Mob_df$Date<-as.Date(Mob_df$Date)
  columun<-"transit_stations_percent_change_from_baseline"
  G1<-IQR_columun_df(Mob_df[Mob_df$State %in% Group1,],columun)
  G2<-IQR_columun_df(Mob_df[Mob_df$State %in% Group2,],columun)
  colnames(G1)<-c("Date","columunG1","columunG1_low","columunG1_up")
  colnames(G2)<-c("Date","columunG2","columunG2_low","columunG2_up")
  All_df<-merge(G1,G2)
  min(All_df[,2:7]);max(All_df[,2:7])
  
  #---PLOT
  time = 1:nrow(All_df)
  
  Date<- All_df$Date
  actual_case<-All_df$columunG1
  
  plot(time,actual_case,pch=16,xaxt = "n",col="white",ylim = c(-55,60), yaxt = "n", 
       xlab= "", ylab="Daily mobility index(%)",main=,cex.main=1.2,cex.lab=1.2) 
  
  polygon(c(time, rev(time)), c(All_df$columunG1_low[time], rev(All_df$columunG1_up[time])), density = NA,
          col = rgb(244,220,216, 190, maxColorValue=255), border = NA)  
  
  polygon(c(time, rev(time)), c(All_df$columunG2_low[time], rev(All_df$columunG2_up[time])), density = NA,
          col = rgb(204,228,240, 150, maxColorValue=255), border = NA)  
  lines(All_df$columunG1, lwd = 1, col ="#C54F35")  
  lines(All_df$columunG2, lwd = 1, col ="#0473B6") 
  
  legend("bottom", legend=c("Group1: Low vaccine coverage", "Group2: High vaccine coverage"),
         col=c("#C54F35","#0473B6"), lty=1:1, cex=0.8)
  
  axis(1,at =seq(1,nrow(All_df),5),labels = Date[seq(1,nrow(All_df),5)],las = 2,cex.axis=1.1)
  axis(2,las = 1,cex.axis = 1.2,cex.axis=1.1) 
  
}
