
chl_dat$CTDFLUOR_10m = unlist(lapply(1:nrow(chl_dat), FUN = function(z){mean_5m(prof_data$DEPTH,prof_data$CTDFLUOR_adj, depth_out = chl_dat$DEPTH[z], window = 10)}))
chl_dat$CTDFLUOR_5m = unlist(lapply(1:nrow(chl_dat), FUN = function(z){mean_5m(prof_data$DEPTH,prof_data$CTDFLUOR_adj, depth_out = chl_dat$DEPTH[z], window = 2)}))
