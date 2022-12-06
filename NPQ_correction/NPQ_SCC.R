

NPQ_SCC = function(depth, fluor, bbp, mld,fwl, sfm){
  # This NPQ correction is new, adjusted from Sackmann et al. (2008)
  #
  #
  zeu_15 = Zeu(depth, fluor,0.15)
  if(is.na(fwl)){fwl = F}
  if(fwl){int = NPQ_X12(depth,fluor,zeu_15)}else{
    if(mld < sfm & mld < zeu_15){
      if(is.na(bbp)){
        int = NPQ_X12(depth,fluor,mld)
      }else{int = NPQ_S08(depth,fluor,bbp,mld)}}
    if(mld < sfm & mld > zeu_15){
      if(is.na(bbp)){
        int = NPQ_X12(depth,fluor,zeu_15)
      }else{int = NPQ_S08(depth,fluor,bbp,zeu_15)}}
    if(mld > sfm & mld < zeu_15){int = NPQ_X12(depth,fluor,mld)}
    if(mld > sfm & mld > zeu_15){int = NPQ_X12(depth,fluor,zeu_15)}
  }
   return(list("NPQdepth" = int$NPQdepth, "corr_fluor" = int$corr_fluor))
}

