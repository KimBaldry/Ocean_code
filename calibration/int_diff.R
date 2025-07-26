int_diff = function(x1,y1){
  df = data.frame(x1, y1)
  idx = complete.cases(x1,y1)
  df = df[idx,]
  f1 = approxfun(df$x1, df$y1, rule = 2)     # piecewise linear function
  f2 = function(x){abs(f1(x))}   # take the positive value
  integrate(f2, min(x1,na.rm = T), max(x1,na.rm = T), subdivisions = 500)$value
}

