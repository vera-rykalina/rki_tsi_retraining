library(readxl)
library (dplyr)
library (ggplot2)
library (wesanderson)
library(tidyr)
library(ggsci)
library(gg.gap)
library(RColorBrewer)
library(ggpubr)
library(lubridate)

# images  562 x 519 (460 x 434)
# boxplot 766 x 443

names(wes_palettes)

din_table <- read_excel("2020_08_07_biobank_din.xlsx",sheet = "din_labs_r")
head(din_table) 
nrow(din_table)
din_table <- mutate(din_table, dif_abs=abs(Difference))

din_table <- din_table %>%
  mutate(din_levels=ifelse(dif_abs==0,"no difference", ifelse(dif_abs<=0.2 & dif_abs>0, "less or equal to 0.2", ifelse(dif_abs>0.2 & dif_abs<=0.6, "less or equal to 0.6", ifelse(dif_abs>0.6 & dif_abs<=1.0, "less or equal to 1.0", "greater than 1.0")))))

# used for the paper
pal_sc <- wes_palette(5, name = "Cavalcanti1", type = "continuous")
din_scatter1 <- din_table %>%
  select(DIN_BBP, DIN_CNR, Difference, din_levels) %>%
  ggplot() +
  geom_point(aes(y = DIN_BBP, x = DIN_CNR, color=factor(din_levels, levels=c("no difference", "less or equal to 0.2","less or equal to 0.6", "less or equal to 1.0","greater than 1.0")))) +
  labs(y='DIN (Biobanque de Picardie)', x='DIN (CNRGH)') +
  scale_y_continuous(limits = c(5, 10), breaks = seq(5,10,1)) +
  scale_x_continuous(limits = c(5, 10), breaks = seq(5,10,1)) +
  scale_color_manual(name ="DIN difference", values=c("#FFB90F","#5D92B1","#778899", "#104E8B","#6B4226")) +
  theme_minimal() +
  theme(text=element_text(size=14, color="grey15"), panel.background = element_rect(fill = "white", colour = "grey80"), axis.title=element_text(size=12, color="grey15"),legend.text=element_text(size=10), legend.background =element_rect(color="white"))+
  theme(legend.position=c(0.77 , 0.22), legend.key=element_rect(color="grey15"))
din_scatter1
ggsave(filename="~/Desktop/Figure_6.png", width= 14, height = 14, units = "cm", dpi = 400)

pal_sc <- wes_palette(2, name = "Cavalcanti1", type = "continuous")
din_scatter2 <- din_table %>%
  select(DIN_BBP, DIN_CNR, dif_abs) %>%
  ggplot() +
  geom_point(aes(y = DIN_BBP, x = DIN_CNR, color=dif_abs)) +
  labs(y='DIN (Biobanque de Picardie)', x='DIN (Centre National de Recherche en Génomique Humaine)', color="DIN difference") +
  scale_y_continuous(limits = c(5, 10), breaks = seq(5,10,1)) +
  scale_x_continuous(limits = c(5, 10), breaks = seq(5,10,1)) +
  #scale_color_gradient(name ="DIN difference", values=pal_sc) +
  theme_minimal() +
  scale_color_gradientn(colors=pal_sc)+
  theme(text=element_text(size=12, color="grey15"), panel.background = element_rect(fill = "white", colour = "grey80"), axis.title=element_text(size=12, color="grey15"), legend.text=element_text(size=11)) +
  theme(legend.position=c(0.8, 0.2), legend.key=element_rect(color="grey15"))
din_scatter2

pal_sc <- wes_palette(n=2, name = "Cavalcanti1", type = "continuous")
din_scatter3 <- din_table %>%
  select(DIN_BBP, DIN_CNR, dif_abs) %>%
  ggplot() +
  geom_point(aes(y = DIN_BBP, x = DIN_CNR, color=dif_abs)) +
  labs(y='DIN (Biobanque de Picardie)', x='DIN (Centre National de Recherche en Génomique Humaine)') +
  scale_y_continuous(limits = c(5, 10), breaks = seq(5,10,1)) +
  scale_x_continuous(limits = c(5, 10), breaks = seq(5,10,1)) +
  scale_color_gradient(name="DIN difference", low="#C4961A" , high= "#293352", breaks=c(0.5,1.0,1.5, 2.0),labels=c("0.5", "1.0","1.5","2.0"),limits=c(0,2.5)) +
  theme_minimal() +
  theme(text=element_text(size=12, color="grey15"), panel.background = element_rect(fill = "white", colour = "grey80"), axis.title=element_text(size=12, color="grey15"), legend.text=element_text(size=11)) +
  theme(legend.position=c(0.8, 0.2), legend.key=element_rect(color="grey15"))
din_scatter3



# fig 4a (1 to 10 scale)
din_scatter5 <- din_table %>%
  select(DIN_BBP, DIN_CNR, dif_abs) %>%
  ggplot(aes(y = DIN_BBP, x = DIN_CNR)) +
  #geom_abline(slope=1, size=0.6, color="#3A5894") + #linetype="dashed"
  geom_abline(slope=1, size=1, color="#3A5894") + #linetype="dashed"
  geom_point(color="#4682B4", fill="white", shape=21) +
  #labs(y='DIN (Biobanque de Picardie)', x='DIN (Centre National de Recherche en Génomique Humaine)', color="DIN difference") +
  labs(y='DIN (Biobanque de Picardie)', x='DIN (CNRGH)', color="DIN difference") +
  scale_y_continuous(limits = c(1, 10), breaks = seq(1,10,1)) +
  scale_x_continuous(limits = c(1, 10), breaks = seq(1,10,1)) +
  #scale_color_gradient(name ="DIN difference", values=pal_sc) +
  theme_minimal() +
  theme(text=element_text(size=12, color="grey15"), panel.background = element_rect(fill = "white", colour = "grey80"), axis.title=element_text(size=12, color="grey15"), legend.text=element_text(size=11)) +
  theme(legend.position=c(0.8, 0.2), legend.key=element_rect(color="grey15"))
din_scatter5
ggsave(filename="~/Desktop/Figure_6_1.png", width= 14, height = 14, units = "cm", dpi = 400)

# fig 4b (5 to 10 scale)
din_scatter6 <- din_table %>%
  select(DIN_BBP, DIN_CNR, dif_abs) %>%
  ggplot(aes(y = DIN_BBP, x = DIN_CNR)) +
  geom_abline(slope=1, size=0.5, color="darkblue") + #linetype="dashed"
  #geom_point(color="steelblue", fill="white", shape=21) +
  geom_point(color="steelblue", fill="white", shape=21) +
  #labs(y='DIN (Biobanque de Picardie)', x="DIN (Centre National de Recherche en Génomique Humaine)", color="DIN difference") +
  labs(y='DIN (BBP)', x="DIN (CNRGH)", color="DIN difference") +
  scale_y_continuous(limits = c(5, 10), breaks = seq(5,10,1)) +
  scale_x_continuous(limits = c(5, 10), breaks = seq(5,10,1)) +
  #scale_color_gradient(name ="DIN difference", values=pal_sc) +
  theme_minimal() +
  theme(text=element_text(size=12, color="grey15"), panel.background = element_rect(fill = "white", colour = "grey80"), axis.title=element_text(size=12, color="grey15"), legend.text=element_text(size=11)) +
  theme(legend.position=c(0.8, 0.2), legend.key=element_rect(color="grey15"))
din_scatter6




#
######################## Regression ##########################
din_scatter7 <- din_table %>%
  select(DIN_BBP, DIN_CNR, dif_abs) %>%
  ggplot(aes(y = DIN_BBP, x = DIN_CNR)) +
  geom_smooth(method='lm', formula= y~x, se=F, color="darkblue")+
  stat_cor(label.y = 10)+ #this means at 10th unit in the y axis, the r squared and p value will be shown
  stat_regline_equation(label.y = 9.5)+ #this means at 9.5th unit regresion line equation will be shown
  geom_abline(slope=1, size=1, color="red") + #linetype="dashed"
  geom_point(color="steelblue") +
  labs(y='DIN (Biobanque de Picardie)', x='DIN (Centre National de Recherche en Génomique Humaine)', color="DIN difference") +
  scale_y_continuous(limits = c(6, 10), breaks = seq(6,10,1)) +
  scale_x_continuous(limits = c(6, 10), breaks = seq(6,10,1)) +
  #geom_abline(lm(DIN_BBP ~ DIN_CNR), size=1, color="darkblue") + 
  theme_minimal() +
  theme(text=element_text(size=12, color="grey15"), panel.background = element_rect(fill = "white", colour = "grey80"), axis.title=element_text(size=12, color="grey15"), legend.text=element_text(size=11)) +
  theme(legend.position=c(0.8, 0.2), legend.key=element_rect(color="grey15"))
din_scatter7



raw_table <- read_excel("2020_08_07_biobank_din.xlsx",sheet = "raw_data")
head(raw_table) 
tail(raw_table)
colnames(raw_table)
nrow(raw_table)
range(raw_table$prel_dte)


less7 <- raw_table %>%
  select(DIN_BBP1) %>%
  filter (DIN_BBP1 <7)
  

timeAll<-raw_table %>%
  select(DIN_BBP1, delai_pc_hr) %>%
  drop_na()
describe(timeAll)
summary(timeAll)


time6h<-raw_table %>%
  select(DIN_BBP1, delai_pc_hr) %>%
  filter (delai_pc_hr<=6 | delai_pc_hr==0) %>%
  drop_na()

describe(time6h)
summary(time6h)

time12h<-raw_table %>%
  select(DIN_BBP1, delai_pc_hr) %>%
  filter (delai_pc_hr>6 & delai_pc_hr<12) %>%
  drop_na()

describe(time12h)
summary(time12h)

timeLonger12h<-raw_table %>%
  select(DIN_BBP1, delai_pc_hr) %>%
  filter (delai_pc_hr>=12) %>%
  drop_na()


mean(timeLonger12h$DIN_BBP1)
round(mean(timeLonger12h$DIN_BBP1), 1)
median(timeLonger12h$DIN_BBP1)
round(median(timeLonger12h$DIN_BBP1), 1)


describe(timeLonger12h)
summary(timeLonger12h)


crb <- raw_table %>%
  select(code_CRB, DIN_BBP1) %>%
  drop_na() %>%
  group_by(code_CRB) %>%
  count()

sum(crb$n)

tail(crb_boxplot)
head(crb_boxplot)
str(crb_boxplot)


pal <- wes_palette(13, name = "Darjeeling1", type = "continuous")
crb_boxplot1 <- raw_table %>%
  select(code_CRB, DIN_BBP1) %>%
  group_by(code_CRB) %>%
  ggplot(aes(x = factor(code_CRB, levels = c("CRB1","CRB2","CRB3","CRB4","CRB5","CRB6","CRB7","CRB8","CRB9","CRB10","CRB11","CRB12","CRB13")), y = DIN_BBP1)) +
  stat_boxplot(geom ='errorbar', aes(color=factor(code_CRB, levels = c("CRB1","CRB2","CRB3","CRB4","CRB5","CRB6","CRB7","CRB8","CRB9","CRB10","CRB11","CRB12","CRB13")))) +
  geom_boxplot(aes(color=factor(code_CRB, levels = c("CRB1","CRB2","CRB3","CRB4","CRB5","CRB6","CRB7","CRB8","CRB9","CRB10","CRB11","CRB12","CRB13")))) +
  labs(x='Anonymized ID CRB', y='DIN')+
  scale_y_continuous(limits = c(5, 10), breaks = seq(5,10,1))+
  scale_color_manual(label=c(CRB1="444",CRB2="237",CRB3="177",CRB4="70",CRB5="80",CRB6="84",CRB7="154",CRB8="266",CRB9="213",CRB10="99",CRB11="147",CRB12="236",CRB13="429"), name ="Sample size", values = pal) +
  theme_minimal() +
  theme(text=element_text(size=12, color="grey15"), panel.background = element_rect(fill = "white", colour = "grey80"), axis.title=element_text(size=12, color="grey15"), legend.text=element_text(size=11))+
  theme(legend.position="right")
crb_boxplot1


colourCount = length(unique(raw_table$code_CRB))
coul <- brewer.pal(9, "YIGnBu")[2:9]
coul <- colorRampPalette(coul)(13)

crb_boxplot2 <- raw_table %>%
  select(code_CRB, DIN_BBP1) %>%
  group_by(code_CRB) %>%
  ggplot(aes(x = factor(code_CRB, levels = c("CRB1","CRB2","CRB3","CRB4","CRB5","CRB6","CRB7","CRB8","CRB9","CRB10","CRB11","CRB12","CRB13")), y = DIN_BBP1), linetype="solid") +
  #geom_hline(yintercept = 7, size=1, color="#3A5894")+
  stat_boxplot(geom ='errorbar', aes(color=factor(code_CRB, levels = c("CRB1","CRB2","CRB3","CRB4","CRB5","CRB6","CRB7","CRB8","CRB9","CRB10","CRB11","CRB12","CRB13")))) +
  geom_boxplot(aes(color=factor(code_CRB, levels = c("CRB1","CRB2","CRB3","CRB4","CRB5","CRB6","CRB7","CRB8","CRB9","CRB10","CRB11","CRB12","CRB13"))), outlier.shape = 21) +
  labs(x='Anonymized CRB ID', y='DIN')+
  scale_y_continuous(limits = c(5, 10), breaks = seq(5,10,1))+
  scale_color_manual(label=c(CRB1="444",CRB2="237",CRB3="177",CRB4="70",CRB5="80",CRB6="84",CRB7="154",CRB8="266",CRB9="213",CRB10="99",CRB11="147",CRB12="236",CRB13="429"), name ="Number of samples", values = c("#0099CC","#0D4F8B", "#53868B","#FFCC11", "#236B8E", "#33A1DE", "#5D92B1", "#6E7B8B", "#33A1C9", "#4A708B", "#87CEFF", "#008B8B", "#8E8E38")) +
  theme_minimal() +
  theme(text=element_text(size=12, color="grey15"), panel.background = element_rect(fill = "white", colour = "grey80"), axis.title=element_text(size=12, color="grey15"), legend.text=element_text(size=10))+
  theme(legend.position=c(0.38, 0.144), legend.box = "horizontal", legend.direction = "horizontal",legend.background =element_rect(color="white"))+
  guides(color=guide_legend(nrow=1, title.position = "top", label.position = "bottom")) +
  #theme (panel.grid.minor.x =element_blank(), panel.grid.minor.y =element_blank(), panel.grid.major.y =element_blank(), panel.grid.major.x =element_blank())
  theme (panel.grid.minor.x =element_blank(), panel.grid.minor.y =element_blank())
crb_boxplot2




colnames(raw_table)
pal_bc <- wes_palette(1, name = "Darjeeling1", type = "continuous")
bc_storage1 <- raw_table %>%
  select(delai_ec, DIN_BBP1) %>%
  ggplot(aes(x= delai_ec, y = DIN_BBP1)) +
  geom_smooth(method='lm', formula= y~x, se=F, color="darkblue")+
  stat_cor(label.y = 10)+ #this means at 10th unit in the y axis, the r squared and p value will be shown
  stat_regline_equation(label.y = 9.5)+ #this means at 9.5th unit regresion line equation will be shown
  #geom_point(color="darkgreen", fill="white", shape=1, alpha=.5, stroke=1.5) +
  geom_hline(yintercept = 7, size=2, color="orange") + #linetype="dashed"
  geom_point(color="darkgreen", fill="white", shape=1) +
  labs(x='Buffy coat storage time (months)', y='DIN')+
  scale_y_continuous(limits = c(1, 10), breaks = seq(1,10,1))+
  scale_x_continuous(limits = c(15, 60), breaks = seq(15,60,5))+
  theme_minimal() +
  theme(text=element_text(size=12, color="grey15"),
        panel.background = element_rect(fill = "white", colour = "grey80"),
        legend.position = "none")
bc_storage1


bc_storage2 <- raw_table %>%
  select(delai_ec, DIN_BBP1) %>%
  ggplot(aes(x= delai_ec, y = DIN_BBP1)) +
  geom_smooth(method='lm', formula= y~x, se=F, color="darkblue")+
  stat_cor(label.y = 3)+ #this means at 10th unit in the y axis, the r squared and p value will be shown
  stat_regline_equation(label.y = 2.5)+ #this means at 9.5th unit regresion line equation will be shown
  #geom_point(color="darkgreen", fill="white", shape=1, alpha=.5, stroke=1.5) +
  geom_hline(yintercept = 7, size=1, color="#00688B") + #linetype="dashed"
  geom_point(color="#00B2EE", fill="white", shape=1) +
  labs(x='Buffy coat storage time (months)', y='DIN')+
  scale_y_continuous(limits = c(5, 10), breaks = seq(5,10,1))+
  scale_x_continuous(limits = c(15, 60), breaks = seq(15,60,5))+
  theme_minimal() +
  theme(text=element_text(size=12, color="grey15"),
        panel.background = element_rect(fill = "white", colour = "grey80"),
        legend.position = "none")
bc_storage2

# fig1 
colnames(raw_table)
bc_storage3 <- raw_table %>%
  select(delai_ec, DIN_BBP1) %>%
  ggplot(aes(x= delai_ec, y = DIN_BBP1)) +
  geom_smooth(method='lm', formula= y~x, se=F, color="black", size=0.6)+
  stat_regline_equation(label.y = 3)+ #this means at 9.5th unit regresion line equation will be shown
  stat_cor(label.y = 2)+ #this means at 10th unit in the y axis, the r squared and p value will be shown
  geom_hline(yintercept = 7, size=1, color="#3A5894") + #linetype="dashed" ,"#3A5894"
  #geom_point(color="darkgreen", fill="white", shape=1, alpha=.5, stroke=1.5) +
  geom_point(color="#4682B4",fill="white", shape=21) +
  labs(x='Buffy coat storage before gDNA extraction (months)', y='DIN')+
  scale_y_continuous(limits = c(1, 10), breaks = seq(1,10,1))+
  scale_x_continuous(limits = c(15, 60), breaks = seq(15,60,5))+
  theme_minimal() +
  theme(text=element_text(size=12, color="grey15"),
        panel.background = element_rect(fill = "white", colour = "grey80"),
        legend.position=c(0.8, 0.2))
bc_storage3

# fig1a
  select(delai_ec, DIN_BBP1) %>%
bc_storage3a <- raw_table %>%
  ggplot(aes(x= delai_ec, y = DIN_BBP1)) +
  geom_hline(yintercept = 7, size=0.6, color="#3A5894") + #linetype="dashed"
  #geom_point(color="darkgreen", fill="white", shape=1, alpha=.5, stroke=1.5) +
  #geom_point(color="#4682B4",fill="white", shape=21) +
  geom_point(color="#4682B4", alpha=0.5, shape=1) +
  labs(x='Buffy coat storage before gDNA extraction (months)', y='DIN')+
  scale_y_continuous(limits = c(1, 10), breaks = seq(1,10,1))+
  scale_x_continuous(limits = c(15, 60), breaks = seq(15,60,5))+
  theme_minimal() +
  theme(text=element_text(size=12, color="grey15"),
        panel.background = element_rect(fill = "white", colour = "grey80"),
        legend.position=c(0.8, 0.2))
bc_storage3a
###########################################################
# DIN < 7
less_than7 <- raw_table %>%
  select(delai_ec, DIN_BBP1,delai_pc_hr) %>%
  filter(DIN_BBP1<=7.5) %>%
  ggplot(aes(x= delai_ec, y = DIN_BBP1, size=delai_pc_hr)) +
  geom_point(color="steelblue", shape=21, fill="white") +
  labs(x='Buffy coat storage time (months)', y='DIN', size="Blood storage (hours)")+
  scale_y_continuous(limits = c(5.5,7.5), breaks = seq(5.5,8,0.5))+
  scale_x_continuous(limits = c(25, 55), breaks = seq(25,55,5))+
  theme_minimal() +
  theme(text=element_text(size=12, color="grey15"),
        panel.background = element_rect(fill = "white", colour = "grey80"),
        legend.position=c(0.8, 0.2))
less_than7

####################################################
timeVStime <- raw_table %>%
  select(delai_ec,delai_pc_hr, DIN_BBP1) %>%
  filter(DIN_BBP1<=7.5) %>%
  ggplot(aes(x= delai_ec, y =delai_pc_hr, size=DIN_BBP1)) +
  geom_point(color="steelblue", shape=21, fill="white") +
  labs(x='Buffy coat storage time (months)', y='Time between blood collection and buffy coat freezing (hours)', size="DIN")+
  scale_y_continuous(limits = c(0,25), breaks = seq(0,25, 5))+
  scale_x_continuous(limits = c(15, 65), breaks = seq(15,65,5))+
  theme_minimal() +
  theme(text=element_text(size=12, color="grey15"),
        panel.background = element_rect(fill = "white", colour = "grey80"),
        legend.position=c(0.8, 0.8))
timeVStime



#segments list cantains more than one number vectors
gg.gap(plot=bc_storage3,segments=c(5,6),tick_width = c(5,1),
       rel_heights=c(0.3,0.005,1),
       ylim=c(1,10)) 


colnames(raw_table)
blood_time1 <- raw_table %>%
  select(delai_pc_hr, DIN_BBP1) %>%
  ggplot(aes(x= delai_pc_hr, y = DIN_BBP1)) +
  #geom_point(color="darkgreen", fill="white", shape=1, alpha=.5, stroke=1.5) +
  geom_vline(xintercept = 6, size=1, color="#3A5894") + #linetype="dashed"
  geom_point(color="#4682B4",fill="white", shape=21) +
  labs(x='Blood storage before buffy coat generation (hours)', y='DIN')+
  scale_y_continuous(limits = c(1, 10), breaks = seq(1,10,1))+
  scale_x_continuous(limits = c(0, 36), breaks = seq(0,36,6))+
  theme_minimal() +
  theme(text=element_text(size=12, color="grey15"),
        panel.background = element_rect(fill = "white", colour = "grey80"),
        legend.position=c(0.8, 0.2))
blood_time1


blood_time2 <- raw_table %>%
  select(delai_pc_hr, DIN_BBP1) %>%
  ggplot(aes(x= delai_pc_hr, y = DIN_BBP1)) +
  #geom_point(color="darkgreen", fill="white", shape=1, alpha=.5, stroke=1.5) +
  geom_vline(xintercept = 6, size=1, color="#3A5894") + #linetype="dashed"
  #geom_point(color="#4682B4",fill="white", shape=21) +
  geom_point(color="#4682B4", alpha=0.5, shape=1) +
  labs(x='Blood storage before buffy coat generation (hours)', y='DIN')+
  scale_y_continuous(limits = c(1, 10), breaks = seq(1,10,1))+
  scale_x_continuous(limits = c(0, 36), breaks = seq(0,36,6))+
  theme_minimal() +
  theme(text=element_text(size=12, color="grey15"),
        panel.background = element_rect(fill = "white", colour = "grey80"),
        legend.position=c(0.8, 0.2))
blood_time2



####################################################################
sop_table <- read_excel("Desktop/2020_08_07_biobank_din.xlsx",sheet = "sop_h")
head(sop_table)
dim(sop_table)


sop6 <- sop_table %>%  # 1820
  filter(time_h <=6) %>%
  count()
v6_sop <- (1820/2636)*100 # 69%

sop_more6 <- sop_table %>%  # 816
  filter(time_h >6) %>%
  count()
v12_sop <- (816/2636)*100 # 31% (30.955999 %)

sop12 <- sop_table %>%  # 742
  filter(time_h >6 &time_h <=12) %>%
  count()
v12_sop <- (742/2636)*100 # 28%

sopGreater12 <- sop_table %>%  # 742
  filter(time_h > 12) %>%
  count()
greater12 <- (74/2636)*100 # 69%

dim(din_table)
same_din <- din_table %>%  # 32 have identical DIN
  filter(Difference == 0) %>%
  count()
v1 <- (32/394)*100 # 8.12%
v1

din_0.5 <- din_table %>%  # 297 have DIN diff >= 0.5 (265)
  filter(abs(Difference) <= 0.5 & abs(Difference)>0) %>%
  count()
v2 <- (265/394)*100 # 67.26%
v2

din_0.2 <- din_table %>%  # 94 have DIN diff <= 0.2 but >0 
  filter(abs(Difference) <= 0.2 & abs(Difference)>0) %>%
  count()
v4 <- (94/394)*100 # 23.85%
v4

din_more_0.2 <- din_table %>%  #268  have DIN diff >0.2 
  filter(abs(Difference) >= 0.2)  %>%
  count()
v5 <- (268/394)*100 # 68.0203%
v5

din_0.6 <- din_table %>%  #197  have DIN diff >0.2 but <= 0.6 
  filter(abs(Difference) <= 0.6 & abs(Difference) > 0.2)  %>%
  count()
v6 <- (197/394)*100 # 50%
v6

din_less_1 <- din_table %>%  #61  have DIN diff >0.6 and <=1
  filter(abs(Difference) > 0.6 & abs(Difference) <=1)  %>%
  count()
v7 <- (61/394)*100 # 15.48223%
v7



din_less_than_1 <- din_table %>%  # 87 have DIN diff < 1
  filter(abs(Difference) <1) %>%
  count()
vv <- (377/394)*100 # 
vv 

din_more_than_1 <- din_table %>%  # 10 have DIN diff > 1
  filter(abs(Difference) > 1) %>%
  count()
v8 <- (10/394)*100 # 2.538071
v8
head(din_more_than_1, 10)

concHigher100 <-raw_table %>% # 1769 samples
  select(conc_ADN, DIN_BBP1) %>%
  filter(conc_ADN>100) %>%
  count()

################################################################
colnames(raw_table)
str(raw_table)
raw.table.date <- raw_table

raw.table.date$congel_dte <-lubridate::dmy(raw.table.date$congel_dte)
raw.table.date$date_extract1 <-lubridate::dmy(raw.table.date$date_extract1)
raw.table.date$date_extract2 <-lubridate::dmy(raw.table.date$date_extract2)
str(raw.table.date)



raw.table.date <- raw.table.date %>% 
  mutate(Days=date_extract1-congel_dte)
colnames(raw.table.date)
select(raw.table.date, prelevement_id, Days)
str(raw.table.date)
raw.table.date$Days <- as.integer(raw.table.date$Days)

raw.table.date <- raw.table.date %>% 
  mutate(Years=Days/365)
select(raw.table.date,prelevement_id, Days, Years)
round(raw.table.date$Years, digits = 3)

# Update

bc.storage.years <- raw.table.date %>%
  select(Years, DIN_BBP1) %>%
  ggplot(aes(x= Years, y = DIN_BBP1)) +
  geom_smooth(method='lm', formula= y~x, se=F, color="black", size=0.6)+
  stat_regline_equation(label.y = 3)+ #this means at 9.5th unit regresion line equation will be shown
  stat_cor(label.y = 2)+ #this means at 10th unit in the y axis, the r squared and p value will be shown
  geom_hline(yintercept = 7, size=0.7, color="#3A5894") + #linetype="dashed" ,"#3A5894"
  #geom_point(color="darkgreen", fill="white", shape=1, alpha=.5, stroke=1.5) +
  geom_point(color="#4682B4",fill="white", shape=21) +
  labs(x='Buffy coat storage before gDNA extraction (years)', y='DIN')+
  scale_y_continuous(limits = c(1, 10), breaks = seq(1,10,1))+
  scale_x_continuous(limits = c(1, 5), breaks = seq(1,5,1))+
  theme_minimal() +
  theme(text=element_text(size=12, color="grey15"),
        panel.background = element_rect(fill = "white", colour = "grey80"),
        legend.position=c(0.8, 0.2))
bc.storage.years          
ggsave(filename="~/Desktop/Figure_3.png", width= 12, height = 11, units = "cm", dpi = 400)

max(raw_table$delai_pc_hr)

blood.time.thres <- raw_table %>%
  select(delai_pc_hr, DIN_BBP1) %>%
  ggplot(aes(x= delai_pc_hr, y = DIN_BBP1)) +
  #geom_point(color="darkgreen", fill="white", shape=1, alpha=.5, stroke=1.5) +
  geom_vline(xintercept = 6, size=0.7, color="#3A5894") + #linetype="dashed"
  geom_hline(yintercept = 7, size=0.7, color="#3A5894") + #linetype="dashed"
  geom_point(color="#4682B4",fill="white", shape=21) +
  labs(x='Blood storage before buffy coat generation (hours)', y='DIN')+
  scale_y_continuous(limits = c(1, 10), breaks = seq(1,10,1))+
  scale_x_continuous(limits = c(0, 36), breaks = seq(0,36,6))+
  theme_minimal() +
  theme(text=element_text(size=12, color="grey15"),
        panel.background = element_rect(fill = "white", colour = "grey80"),
        legend.position=c(0.8, 0.2))
blood.time.thres
ggsave(filename="~/Desktop/Figure_4.png", width= 12, height = 11, units = "cm", dpi = 400)

#######################################
crb <- raw_table %>%
  select(code_CRB, DIN_BBP1) %>%
  drop_na() %>%
  group_by(code_CRB) %>%
  count()

crb <- crb[c(1, 6:13, 2:5),]  
crb

crb.boxplot.size <- raw_table %>%
  select(code_CRB, DIN_BBP1) %>%
  drop_na() %>%
  group_by(code_CRB) %>%
  ggplot(aes(x = factor(code_CRB, levels = c("CRB1","CRB2","CRB3","CRB4","CRB5","CRB6","CRB7","CRB8","CRB9","CRB10","CRB11","CRB12","CRB13")), y = DIN_BBP1), linetype="solid") +
  #geom_hline(yintercept = 7, size=1, color="#3A5894")+
  stat_boxplot(geom ='errorbar', aes(color=factor(code_CRB, levels = c("CRB1","CRB2","CRB3","CRB4","CRB5","CRB6","CRB7","CRB8","CRB9","CRB10","CRB11","CRB12","CRB13")))) +
  geom_boxplot(aes(color=factor(code_CRB, levels = c("CRB1","CRB2","CRB3","CRB4","CRB5","CRB6","CRB7","CRB8","CRB9","CRB10","CRB11","CRB12","CRB13"))), outlier.shape = 21) +
  labs(x='Anonymized CRB ID', y='DIN')+
  scale_y_continuous(limits = c(1, 10), breaks = seq(1,10,1))+
  scale_color_manual(label=c(CRB1="444",CRB2="237",CRB3="177",CRB4="70",CRB5="80",CRB6="84",CRB7="154",CRB8="266",CRB9="213",CRB10="99",CRB11="147",CRB12="236",CRB13="429"), name ="Number of samples", values = c("#0099CC","#0D4F8B", "#53868B","#FFCC11", "#236B8E", "#33A1DE", "#5D92B1", "#6E7B8B", "#33A1C9", "#4A708B", "#87CEFF", "#008B8B", "#8E8E38")) +
  geom_text(data=crb, aes(code_CRB, Inf, label = n), vjust = 1.5, color=c("#0099CC","#0D4F8B", "#53868B","#FFCC11", "#236B8E", "#33A1DE", "#5D92B1", "#6E7B8B", "#33A1C9", "#4A708B", "#87CEFF", "#008B8B", "#8E8E38")) +
  theme_minimal() +
  theme(text=element_text(size=12, color="grey15"), panel.background = element_rect(fill = "white", colour = "grey80"), axis.title=element_text(size=12, color="grey15"), legend.text=element_text(size=10))+
  theme(legend.position = "none") +
  theme (panel.grid.minor.x =element_blank(), panel.grid.minor.y =element_blank())
crb.boxplot.size



crb <- raw_table %>%
  select(code_CRB, DIN_BBP1) %>%
  drop_na() %>%
  group_by(code_CRB) %>%
  count()

crb <- crb[c(1, 6:13, 2:5),]  
crb

crb.boxplot.size <- raw_table %>%
  #select(code_CRB, DIN_BBP1) %>%
  #drop_na() %>%
  #group_by(code_CRB) %>%
  left_join(crb) %>%
  mutate(code_CRBn = paste0(code_CRB, "\n", "(n=", n,')')) %>%
  ggplot(aes(x=factor(code_CRBn, levels = c("CRB1\n(n=444)","CRB2\n(n=237)","CRB3\n(n=177)",
                                            "CRB4\n(n=70)","CRB5\n(n=80)","CRB6\n(n=84)",
                                            "CRB7\n(n=154)","CRB8\n(n=266)","CRB9\n(n=213)",
                                            "CRB10\n(n=99)","CRB11\n(n=147)","CRB12\n(n=236)",
                                            "CRB13\n(n=429)")), y = DIN_BBP1),linetype="solid") +
  geom_hline(yintercept = 7, size=0.7, color="#3A5894")+
  stat_boxplot(geom ='errorbar', aes(color=factor(code_CRB, levels = c("CRB1","CRB2","CRB3",
                                                                       "CRB4","CRB5","CRB6",
                                                                       "CRB7","CRB8","CRB9",
                                                                       "CRB10","CRB11","CRB12"
                                                                       ,"CRB13"))), lwd=0.75) + #lwd=1
  geom_boxplot(aes(color=factor(code_CRB, levels = c("CRB1","CRB2","CRB3","CRB4",
                                                     "CRB5","CRB6","CRB7","CRB8",
                                                     "CRB9","CRB10","CRB11",
                                                     "CRB12","CRB13"))), outlier.shape = 21, 
               lwd=0.65, fatten=1.1)+ #lwd=0.8, fatten=0.9
  labs(x='Anonymized CRB ID', y='DIN')+
  scale_y_continuous(limits = c(1, 10), breaks = seq(1,10,1))+
  scale_color_manual(label=c(CRB1="444",CRB2="237",CRB3="177",CRB4="70",
                             CRB5="80",CRB6="84",CRB7="154",CRB8="266",
                             CRB9="213",CRB10="99",CRB11="147",CRB12="236",
                             CRB13="429"),
                     name ="Number of samples", 
                     values = c("#0099CC", "#8E8E38", "#53868B",
                                "#008B8B", "#236B8E", "#33A1DE", 
                                "#5D92B1", "#6E7B8B", "#33A1C9", 
                                "#4A708B", "#87CEFF", "#FFCC11", "#0D4F8B")) +
  theme_minimal() +
  theme(text=element_text(size=12, color="grey15"), panel.background = element_rect(fill = "white", colour = "grey80"), axis.title=element_text(size=12, color="grey15"), legend.text=element_text(size=10))+
  theme(legend.position = "none") +
  theme (panel.grid.minor.x =element_blank(), panel.grid.minor.y =element_blank())
crb.boxplot.size
ggsave(filename="~/Desktop/Figure_5_larger_2.png", width= 22, height = 11, units = "cm", dpi = 400)

