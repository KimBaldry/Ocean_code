

NPQ_SCC = function(depth, fluor, bbp, mld,fwl, sfm){
  # This NPQ correction is new, adjusted from Sackmann et al. (2008)
  #
  #
  zeu_20 = Zeu(depth, fluor,0.2)
  if(is.na(fwl)){fwl = F}
  if(fwl){int = NPQ_X12(depth,fluor,zeu_20)}else{
    if(zeu_20 < sfm & mld < zeu_20){
      if(all(is.na(bbp))){
        int = NPQ_X12(depth,fluor,mld)
      }else{int = NPQ_S08(depth,fluor,bbp,mld)}}
    if(zeu_20 < sfm & mld > zeu_20){
      if(all(is.na(bbp))){
        int = NPQ_X12(depth,fluor,zeu_20)
      }else{int = NPQ_S08(depth,fluor,bbp,zeu_20)}}
    if(zeu_20 > sfm & mld < zeu_20){int = NPQ_P18(depth,fluor,mld)}
    if(zeu_20 > sfm & mld > zeu_20){
      if(all(is.na(bbp))){
        int = NPQ_X12(depth,fluor,zeu_20)
      }else{int = NPQ_S08(depth,fluor,bbp,mld)}}
  }
  return(list("NPQdepth" = int$NPQdepth, "corr_fluor" = int$corr_fluor))
}

