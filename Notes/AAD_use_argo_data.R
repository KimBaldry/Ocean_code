## on the server
f <- raadfiles::argo_files()

x <- tidync::tidync(f$fullname[1])
x %>% hyper_array(select_var = c("PRES", "TEMP")) 


## in dev
stars::read_ncdf(f$fullname[1])
