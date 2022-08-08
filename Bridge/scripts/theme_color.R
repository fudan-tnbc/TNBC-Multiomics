mytheme<-theme(
  plot.background = element_rect(colour = "white"),
  axis.title.y = element_text(size=16),
  axis.title.x = element_text(size=32),
  axis.text.y = element_text(colour= "black",size=16),
  axis.text.x = element_text(colour= "black",size=16),
  title=element_text(colour= "black",size=16),
  panel.background = element_rect(fill="white"),
  panel.grid = element_blank(),
  strip.text=element_text(size=16))


mytheme2<-theme(
  plot.background = element_rect(colour = "white"),
  axis.title.y = element_text(size=16),
  axis.title.x = element_text(size=16),
  axis.text.y = element_text(colour= "black",size=16),
  axis.text.x = element_text(colour= "black",size=16),
  title=element_text(colour= "black",size=16),
  panel.background = element_rect(fill="white"),
#  panel.grid = element_blank(),
  strip.text=element_text(size=16))
  
elev_palette = c("#e64b35","#f39b7f","#f8d714","#00a087","#91d1c2","#4dbbd5","#3c5488",
                 "#7e6148","#8491b4","#b09c85","#859900","#383838")

tri_palette = c("#915359","#9dcaac","#7799c1")
