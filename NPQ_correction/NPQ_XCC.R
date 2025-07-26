NPQ_XCC = function(depth, fluor, bbp, mld, fwl,sfm){
  # This NPQ correction is new, adjusted from Xing et al. (2012)
  #
  #
  zeu_20 = Zeu(depth, fluor,0.2)
  if(is.na(fwl)){fwl = F}
  if(fwl){int = NPQ_X12(depth,fluor,zeu_20)}else{
    if(zeu_20 < sfm & mld < zeu_20){int = NPQ_X12(depth,fluor,mld)}
    if(zeu_20 < sfm & mld > zeu_20){int = NPQ_X12(depth,fluor,zeu_20)}
    if(zeu_20 > sfm & mld < zeu_20){int = NPQ_P18(depth,fluor,mld)}
    if(zeu_20 > sfm & mld > zeu_20){int = NPQ_X12(depth,fluor,zeu_20)}
  }
  return(list("NPQdepth" = int$NPQdepth, "corr_fluor" = int$corr_fluor))
}
