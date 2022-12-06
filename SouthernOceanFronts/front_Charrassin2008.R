
front_Charrassin2008 = function(press_prof, temp_prof){
  prof_array = data.frame(cbind(press_prof, temp_prof))
  prof_array = prof_array[complete.cases(prof_array), ]
  colnames(prof_array) = c("press","temp")
  if(max(press_prof,na.rm = T) < 510){return("max press < 500 dbar")}else{
    
    
    
    ## 1
    # T - P
    t1_lb = 0.82
    t1_mb = 1.28
    t1_ub = 2.20
    t1_pub = 510
    t1_plb = 490

    mean_temp = mean(prof_array$temp[which(prof_array$press >= t1_plb & prof_array$press <= t1_pub)], na.rm = T)
    
    is3 = mean_temp > 2.20
    is2 = mean_temp < 2.20 & mean_temp > 1.28
    is1 = mean_temp < 1.28 & mean_temp > 0.82
    is0 = mean_temp <0.82
    
    #combine result of classification tests
    result = c(is0, is1,is2,is3)

    if(sum(result) == 1){return(c(0,1,2,3)[which(result == 1)])}
  }
}  