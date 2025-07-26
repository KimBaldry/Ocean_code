## This function will create the supplementry plot for front classification from the water colum structure wide format summary table
#
#
# table - the wide format summary table created by

sup_front = function(input_table){
  wide_table = input_table %>% dplyr::filter(is.finite(FZ_sokolov), FZ_Charrassin != "max press < 500 dbar")
  
  t_charrassin <- as.data.frame(prop.table(table(as.factor(wide_table$FZ_sokolov),as.factor(wide_table$FZ_Charrassin)),margin = 1))
  colnames(t_charrassin)[1:2] = c("FZ_sokolov", "FZ_Charrassin")
  
   p_charrassin = ggplot(t_charrassin, aes(FZ_sokolov, FZ_Charrassin, fill= Freq)) + 
    geom_tile()+ scale_fill_viridis(limits = c(0, 1)) + theme(legend.position = "none")
  
  t_sorarc =  as.data.frame(prop.table(table(as.factor(wide_table$FZ_sokolov),as.factor(wide_table$FZ_sorarc)),margin = 1))
  colnames(t_sorarc)[1:2] = c("FZ_sokolov", "FZ_sorarc")
  
  p_sorarc = ggplot(t_sorarc, aes(FZ_sokolov, FZ_sorarc, fill= Freq)) + 
    geom_tile()+ scale_fill_viridis(limits = c(0, 1)) + theme(legend.position = "none")
  
  t_SokolovA = as.data.frame(prop.table(table(as.factor(wide_table$FZ_sokolov),as.factor(wide_table$FZ_SokolovArgo)),margin = 1))
  colnames(t_SokolovA)[1:2] = c("FZ_sokolov", "FZ_SokolovArgo")
    
  p_FZ_sokolov = ggplot(t_SokolovA, aes(FZ_sokolov, FZ_SokolovArgo, fill= Freq)) + 
    geom_tile() + scale_fill_viridis(limits = c(0, 1)) + theme(legend.position = "none")
  
  t_class <- as.data.frame(prop.table(table(as.factor(wide_table$FZ_sokolov),as.factor(wide_table$FZ_insitu)),margin = 1))
  colnames(t_class)[1:2] = c("FZ_sokolov", "FZ_new")
  
  p_class = ggplot(t_class, aes(FZ_sokolov, FZ_new, fill= Freq)) + 
    geom_tile()+ scale_fill_viridis(limits = c(0, 1)) + theme(legend.position = "none")
  
  
  plist = list()
  
  plot.lab = ggplotGrob(p_charrassin + ggtitle("Charrassin et al. (2008)"))
  plot.lab$layout$l[plot.lab$layout$name == "title"] <- 1
  plist[[1]] = plot.lab       
  
  plot.lab = ggplotGrob(p_sorarc + ggtitle("SOARC"))
  plot.lab$layout$l[plot.lab$layout$name == "title"] <- 1
  plist[[2]] = plot.lab  
  
  plot.lab = ggplotGrob(p_FZ_sokolov + ggtitle("Sokolov et al. (2009) from Argo"))
  plot.lab$layout$l[plot.lab$layout$name == "title"] <- 1
  plist[[3]] = plot.lab  
  
  plot.lab = ggplotGrob(p_class + ggtitle("New based on Sokolov et al. (2009) from Argo"))
  plot.lab$layout$l[plot.lab$layout$name == "title"] <- 1
  plist[[4]] = plot.lab  
  
  # guide
  dummy = ggplot(data = t_SokolovA, aes(x=FZ_sokolov,y=FZ_SokolovArgo)) +geom_raster(aes(fill= Freq)) + 
    scale_fill_viridis(name = "Agreement", limits = range(0,1), breaks = c(0,0.25,0.5,0.75,1),labels = c(0,25,50,75,100)) +
    guides(fill = guide_colourbar(title.position = "top",title.hjust = .5,label.position = "right"))+
    theme(legend.position = "right",legend.key.width = unit(3,"line"))
  
  g = ggplotGrob(dummy)$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  plist[[5]] = legend
  lay = cbind(matrix(c(rep(1,3),rep(2,3),rep(3,3), rep(4,3)),nrow = 1, byrow = T),c(5))
  grid.arrange(grobs = plist,layout_matrix = lay)
  
  
}
