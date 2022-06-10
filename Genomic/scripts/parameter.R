
#################set color set theme
color.PT = '#248067' #HAIWANGLV
color.TT = '#c04851' #YUHONG
mypal =pal_npg("nrc", alpha = 1)(10)
## [1] "#E64B35B2" "#4DBBD5B2""#00A087B2" "#3C5488B2" "#F39B7FB2""#8491B4B2"

## [7] "#91D1C2B2" "#DC0000B2" "#7E6148B2"

colors.type.TNBC <- c(color.PT,color.TT)
names(colors.type.TNBC ) <- c('normal', 'tumor')

mytheme<-theme(
  plot.background = element_rect(colour = "white"),
  plot.title = element_text(hjust=0.5,size=12,face='bold'),
  axis.title = element_text(size=12),
  axis.text = element_text(colour= "black",size=11),
  strip.text=element_text(size=11),
  legend.title = element_text(size=10),
  panel.grid = element_blank())

usample <- c('normal', 'tumor')
comb<- combn(usample,2)

