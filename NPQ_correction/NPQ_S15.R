NPQ_S15 <- function(depth, fluor, bbp, ed){
  # this function calculates the NPQ_depth and a NPQ corrected fluoresence profile based on Thomalla et al (2018).
  # They interperated the method of Swart et al. (2015) to use euphotic depth rather than mld.  
  # Swart actually uses a temperature MLD criterion threshold of 0.2 degC
  # adapted from Sackmann et al. (2008)
  # 
  # ref 1:Thomalla, S. J., Moutier, W., Ryan‐Keogh, T. J., Gregor, L., & Schütt, J. (2018). 
  # An optimized method for correcting fluorescence quenching using optical backscattering on autonomous platforms.
  # Limnology and Oceanography: Methods, 16(2), 132-144.
  # 
  # ref 2: Swart, S., Thomalla, S. J., & Monteiro, P. M. S. (2015). 
  # The seasonal cycle of mixed layer dynamics and phytoplankton biomass in the Sub-Antarctic Zone: 
  # A high-resolution glider experiment. 
  # Journal of Marine Systems, 147, 103-115.
  #
  # ref3:  Sackmann, B., Perry, M., and Eriksen, C.: 
  #       Seaglider observations of variability in daytime fluorescence quenching of chlorophyll-a in 
  #       Northeastern Pacific coastal waters, 
  #       Biogeosciences Discussions, 5, 2839-2865, 2008.
  #
  #
  # Thomalla et al 2018 used Ed was defined using in-situ PAR profiles
  # this method originally includes a visual classification into
  # 1. quenched above MLD
  # 2. quenched below MLD, but uniform ratio exissts below
  # 3. quenched below MLD but no uniform region of ratio exists -> these profiles cannot be corrected
  #
  #
  # Note that units for NPQ depth are inherrited from pres. So, m for Depth and dbar for pressure measurements
  #
  df = data.frame(depth, fluor, bbp)
  df = df[complete.cases(df),]
  r = df$fluor/df$bbp # ratio used by Sackmann et al. (2008)
  max_r = max(r[df$depth < min(ed,df$depth[which.max(df$fluor)], na.rm = T)], na.rm = T) # identify maximum r in upper MLD
  NPQ_depth = df$depth[rev(which(r == max_r))[1]]
  corr_fluor = fluor 
  if(length(which(depth <= NPQ_depth)) < 1){corr_fluor = fluor
  NPQ_depth = 0}else{
    fluor_NPQ = max_r*(bbp[depth <= NPQ_depth])
    corr_fluor[depth <= NPQ_depth] = fluor_NPQ # correct above NPQ_depth
  }
  return(list("NPQdepth" = NPQ_depth, "corr_fluor" = corr_fluor))
 }