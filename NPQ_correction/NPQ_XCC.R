NPQ_XCC = function(depth, fluor, bbp, mld, fwl,sfm){
  # This NPQ correction is new, adjusted from Xing et al. (2012)
  #
  #
  zeu_15 = Zeu(depth, fluor,0.15)
  if(is.na(fwl)){fwl = F}
  if(fwl){int = NPQ_X12(depth,fluor,zeu_15)}else{
    if(mld < sfm & mld < zeu_15){int = NPQ_X12(depth,fluor,mld)}
    if(mld < sfm & mld > zeu_15){int = NPQ_X12(depth,fluor,zeu_15)}
    if(mld > sfm & mld < zeu_15){int = NPQ_X12(depth,fluor,mld)}
    if(mld > sfm & mld > zeu_15){int = NPQ_X12(depth,fluor,zeu_15)}
  }
  return(list("NPQdepth" = int$NPQdepth, "corr_fluor" = int$corr_fluor))
}

