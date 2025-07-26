
front_class_Sokolov2009_theta = function(press_prof, temp_prof, sal_prof,lat,lon){
  prof_array = data.frame(cbind(press_prof, swTheta(salinity = sal_prof, temperature = temp_prof, pressure = press_prof,latitude = lat, longitude = lon)))
  prof_array = prof_array[complete.cases(prof_array), ]
  colnames(prof_array) = c("press","temp")
  if(max(press_prof,na.rm = T) < 500){return("max press < 500 dbar")}else{
    
    t_min  = min(prof_array$temp, na.rm = T)
    t_max = max(prof_array$temp, na.rm = T)
    ## 0
    # Tmin <= -1.24 C
    t0_ub = -0.95
    
    fun_0 = any(t_min < t0_ub)
    is0 = fun_0
    
    ## 1
    # Tmin 
    t1_lb = -0.95
    t1_ub = -0.56
    
    # Tmax
    t1_ub_max = 1.59
    
    tmin_1 = t_min < t1_ub & t_min > t1_lb
    tmax_1 = t_max < t1_ub_max
    fun_1 = tmin_1 == 1 & tmax_1 == 1
    is1 = fun_1
    
    ## 2
    # Tmin 
    t2_lb = -0.56
    t2_ub = 0
    
    # Tmax
    t2_lb_max = 1.59
    t2_ub_max = 1.93
    
    
    tmin_2 = any(t_min < t2_ub & t_min > t2_lb)
    tmax_2 = any(t_max < t2_ub_max & t_max > t2_lb_max)
    fun_2 = tmin_2 == 1 & tmax_2 == 1
    is2 = fun_2
    
    ## 3
    # Tmin 
    t3_lb = 0
    t3_ub = 0.98
    
    # Tmax
    t3_lb_max = 1.93
    t3_ub_max = 2.11
    
    tmin_3 = any(t_min < t3_ub & t_min > t3_lb)
    tmax_3 = any(t_max < t3_ub_max & t_max > t3_lb_max)
    fun_3 = tmin_3 == 1 & tmax_3 == 1
    is3 = fun_3
    
    ## 4
    # Tmin 
    t4_lb = 0.98
    t4_ub = 1.15
    
    # Tmax
    t4_lb_max = 2.11
    t4_ub_max = 2.25
    
    tmin_4 = any(t_min < t4_ub & t_min > t4_lb)
    tmax_4 = any(t_max < t4_ub_max & t_max > t4_lb_max)
    fun_4 = tmin_4 == 1 & tmax_4 == 1
    is4 = fun_4
    
    
    #combine result of classification tests
    result = c(is0, is1, is2, is3, is4)
    
    isunclass = sum(result) > 1
    isnoclass = sum(result) == 0
    if(isunclass){
      if(sum(result) == 2 & abs(c(0:4)[which(result == 1)[1]] - c(0:4)[which(result == 1)[2]] ) == 1){
        return(sample(c("0","1","2","3","4")[which(result == 1)],1))
      }else{return("unclassified")}
      }
    if(isnoclass){return("no class")}
    if(sum(result) == 1){return(c("0","1","2","3","4")[which(result == 1)])}
    #return(paste(c("0","1","2","3","4")[which(result == 1)],sep = ","))
  }
 }  