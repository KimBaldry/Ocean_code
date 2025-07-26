# calculate dark offset for optics
library(birk)
optic_dark = function(depth, fluor, bp){
  depth = depth[is.finite(fluor)]
  fluor = fluor[is.finite(fluor)]
  m_depth = max(depth,na.rm = T)
  odx = which(depth < 500)
  idxbp = which.closest(depth, bp):odx[length(odx)]
  
  ifelse(m_depth >= 500 & length(which(depth[is.finite(fluor)] > 500)) > 30,
         min(lm(runmed(fluor[-odx],5, na.action = "na.omit") ~ depth[-odx], na.action=na.exclude)$fitted.values,na.rm = T),
         min(lm(runmed(fluor[idxbp],5, na.action = "na.omit") ~ depth[idxbp], na.action=na.exclude)$fitted.values,na.rm = T)
  )
  
}
