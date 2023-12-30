#' Add Two Numbers
#'
#' This function takes two numbers and returns their sum.
#'
#' @param x A numeric, the first number to add.
#' @param y A numeric, the second number to add.
#' @return A numeric, the sum of x and y.
#' @export
IQR.mutation <- function (df,mutation) {
  #df<-state_mutation_all
  #mutation<-"G339D"
  IQR_mutation_df<-df[,c(mutation,"Date")]
  colnames(IQR_mutation_df)<-c("mutation","Date")
  median_df<-aggregate(mutation~Date,IQR_mutation_df,median)
  colnames(median_df)<-c("Date","mutation")
  
  median_df_low<-aggregate(mutation~Date,IQR_mutation_df,quantile,probs = 0.25)# probs = 0.25 probs = 0.025
  colnames(median_df_low)<-c("Date","mutation_low")
  
  median_df_up<-aggregate(mutation~Date,IQR_mutation_df,quantile,probs = 0.75) #probs = 0.75 probs = 0.975
  colnames(median_df_up)<-c("Date","mutation_up")
  
  IQR_mutation_df1<-merge(median_df,median_df_low)
  IQR_mutation_df2<-merge(IQR_mutation_df1,median_df_up)
  
  return(IQR_mutation_df2)
}