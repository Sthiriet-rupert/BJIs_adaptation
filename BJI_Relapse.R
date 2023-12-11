
#Cytotox + Hemolysis
library(ggpubr)

dataP<-read.table("Cytotox.txt",header = T,sep = "\t")

p1<-ggplot(dataP, aes(x=Strain, y=Point, fill=Strain)) + stat_summary(geom="bar",fun = mean,colour="black",  position='dodge', width=0.7)+ theme_linedraw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ stat_summary(geom="errorbar",fun.data = "mean_sd",aes(width=0.5))+ geom_point(position=position_dodge(width=.7)) +scale_y_continuous(expand = c(0,0), limits = c(0,1.5)) + xlab("")+ylab("Cytotoxicity (relative fluorescence)")+ scale_fill_manual(values=c("firebrick4","dodgerblue4","seagreen4","#D3D3D3","indianred3","dodgerblue2","seagreen3"), name="Strain")+ theme(text = element_text(size=15),legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1))+ geom_hline(yintercept=1, linetype="dashed", color = "black") +facet_grid(.~Group,space = "free",scales = "free")

dataH<-read.table("Hemolysis_Mean.txt",header = T,sep = "\t")
dataPH<-read.table("Hemolysis_Point.txt",header = T,sep = "\t")
order=c("1h30","6h","24h")

p2<-ggplot(dataH, aes(x=Strain, y=Mean, fill=factor(Time, levels = order))) + geom_bar(colour="black", stat="identity", position='dodge', width=0.7)+theme_classic()+ geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.2, position=position_dodge(0.7))+ geom_point(data=dataPH, aes(x=Strain, y=Point, fill=factor(Time, levels = order)), position=position_dodge(width=.7)) +scale_y_continuous(expand = c(0,0), limits = c(0,1.5)) + xlab("")+ylab(expression(atop("Hemolysis capacity relative to", paste("Infection strain (a.u.)"))))+ scale_fill_manual(values=c("white",  "lightgrey",  "darkgrey"), name="")+ theme(text = element_text(size=15),strip.text = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.margin=unit(c(0.2,0,1,0.5), "cm"))+ geom_hline(yintercept=1, linetype="dashed", color = "black")

ggarrange(p1,p2,labels=c("A","B"),nrow = 1)
ggsave("Cytotox_Hemolysis.tiff",width=7,height = 5,device = "tiff",dpi = 700)



#Osteoblasts
dataP<-read.table("osteoblasts_points.txt",header = T,sep = "\t")

ggplot(dataP, aes(x=Sample, y=cfu, fill=Sample)) + stat_summary(geom="bar",fun = mean,colour="black",  position='dodge', width=0.7)+ theme_linedraw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ stat_summary(geom="errorbar",fun.data = "mean_sd",aes(width=0.5))+ geom_point(position=position_dodge(width=.7)) +scale_y_continuous(trans = "log",expand = c(0,0), limits = c(NA,15000000),breaks = c(0, 1e+02, 1e+04, 1e+06, 1e+08)) + xlab("")+ylab(expression(paste("CFU / ", 10^6, " osteoblast cells")))+ scale_fill_manual(values=c("firebrick4","dodgerblue4","seagreen4","#D3D3D3","indianred3","dodgerblue2","seagreen3","white"), name="Strain")+ theme(text = element_text(size=15),legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1))+ geom_hline(yintercept=1, linetype="dashed", color = "black") +facet_grid(.~Group,space = "free",scales = "free")

ggsave("Osteoblasts.tiff",width=6,height = 5,device = "tiff",dpi = 700)


#Biofilm formation
##I1-R1

order=c("I1","R1","MG1655","TG1")
dataCFU<-read.table("I1-R1_Biof_mean.txt",header = T,sep = "\t")
dataCFUp<-read.table("I1-R1_Biof_points.txt",header = T,sep = "\t")

ggplot(dataCFU, aes(x=Group, y=Mean,fill=factor(Sample, levels = order))) + geom_bar(colour="black", stat="identity", position=position_dodge(), width=0.7)+theme_classic()+ geom_errorbar(aes(ymin=pmax(0,Mean-SD), ymax=Mean+SD), width=.2, position=position_dodge(0.7))+ geom_point(data=dataCFUp, aes(x=Group, y=cfu, fill=factor(Sample, levels = order)), position=position_dodge(width=.7)) + ylab("CFU / mL")+ xlab("") + scale_fill_manual(values=c("firebrick4", "indianred3", "#D3D3D3", "#FFFFFF"), name="Strain")+ theme(text = element_text(size=15),axis.text.x=element_blank(),axis.ticks.x=element_blank())+facet_grid(.~Group,space = "free",scales = "free") +scale_y_continuous(trans = "log",expand=c(0,0), breaks = c(0, 1e+02, 1e+04, 1e+06, 1e+08), limits = c(NA, 1e+10))
ggsave("I1-R1_CFUlog.tiff",width=7,height = 5,device = "tiff",dpi = 700)


##I2-R2

dataCFU<-read.table("I2-R2_Biof_mean.txt",header = T,sep = "\t")
dataCFUp<-read.table("I2-R2_Biof_points.txt",header = T,sep = "\t")
order=c("I2","R2","MG1655","TG1")

ggplot(dataCFU, aes(x=Group, y=Mean,fill=factor(Sample, levels = order))) + geom_bar(colour="black", stat="identity", position=position_dodge(), width=0.7)+theme_classic()+ geom_errorbar(aes(ymin=pmax(0,Mean-SD), ymax=Mean+SD), width=.2, position=position_dodge(0.7))+ geom_point(data=dataCFUp, aes(x=Group, y=cfu, fill=factor(Sample, levels = order)), position=position_dodge(width=.7)) + ylab("CFU / mL")+ xlab("") + scale_fill_manual(values=c("dodgerblue4", "dodgerblue2", "#D3D3D3", "#FFFFFF"), name="Strain")+ theme(text = element_text(size=15),axis.text.x=element_blank(),axis.ticks.x=element_blank())+facet_grid(.~Group,space = "free",scales = "free") +scale_y_continuous(trans = "log",expand=c(0,0), breaks = c(0, 1e+02, 1e+04, 1e+06, 1e+08), limits = c(NA, 1e+10))
ggsave("I2-R2_CFUlog.tiff",width=7,height = 5,device = "tiff",dpi = 700)


##I3-R3

dataCFU<-read.table("I3-R3_Biof_mean.txt",header = T,sep = "\t")
dataCFUp<-read.table("I3-R3_Biof_points.txt",header = T,sep = "\t")
order=c("I3","R3","MG1655","TG1")

ggplot(dataCFU, aes(x=Group, y=Mean,fill=factor(Sample, levels = order))) + geom_bar(colour="black", stat="identity", position=position_dodge(), width=0.7)+theme_classic()+ geom_errorbar(aes(ymin=pmax(0,Mean-SD), ymax=Mean+SD), width=.2, position=position_dodge(0.7))+ geom_point(data=dataCFUp, aes(x=Group, y=cfu, fill=factor(Sample, levels = order)), position=position_dodge(width=.7)) + ylab("CFU / mL")+ xlab("") + scale_fill_manual(values=c("seagreen4", "seagreen3", "#D3D3D3", "#FFFFFF"), name="Strain")+ theme(text = element_text(size=15),axis.text.x=element_blank(),axis.ticks.x=element_blank())+facet_grid(.~Group,space = "free",scales = "free") +scale_y_continuous(trans = "log",expand=c(0,0), breaks = c(0, 1e+02, 1e+04, 1e+06, 1e+08), limits = c(NA, 1e+10))
ggsave("I3-R3_CFUlog.tiff",width=7,height = 5,device = "tiff",dpi = 700)


#Serum survival - All
dataP<-read.table("Serum_All_Points.txt",header = T,sep = "\t")

ggplot(dataP, aes(x=Strain, y=Value, fill=Strain)) + stat_summary(geom="bar",fun = mean,colour="black",  position='dodge', width=0.7)+ theme_linedraw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ stat_summary(geom="errorbar",fun.data = "mean_sd",aes(width=0.5))+ geom_point(position=position_dodge(width=.7)) +scale_y_continuous(expand = c(0,0), limits = c(0,100)) + xlab("")+ylab("Normalized Serum Survival (%)")+ scale_fill_manual(values=c("white","firebrick4","dodgerblue4","seagreen4","#D3D3D3","indianred3","dodgerblue2","seagreen3"), name="Strain")+ theme(text = element_text(size=15),legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1))+facet_grid(.~Group,space = "free",scales = "free")
ggsave("Serum_survival_All.tiff",width=7,height = 5,device = "tiff",dpi = 700)


#Serum survival - mutants
data<-read.table("SerumKilling_Mean.txt",header = T,sep="\t")
datPoint<-read.table("SerumKilling_Points.txt",header = T,sep = "\t")

ggplot(data, aes(x=Strain, y=Mean, fill=Strain)) + geom_bar(colour="black", stat="identity", position='dodge', width=0.7)+theme_classic()+ geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.2, position=position_dodge(0.7))+ geom_point(data=datPoint, aes(x=Strain, y=Point, fill=Strain), position=position_dodge(width=.7)) +scale_y_continuous(expand = c(0,0), limits = c(0,1.2)) + ylab(expression(atop("Serum survival relative to", paste("R2 strain (a.u.)"))))+ scale_fill_manual(values=c("dodgerblue4", "white", "grey", "dodgerblue2"), name="Strain")+ theme(text = element_text(size=20),legend.position ="none",axis.text.x = element_text(angle = 45, hjust = 1))+ geom_hline(yintercept=1, linetype="dashed", color = "black")+xlab("")+ scale_x_discrete(limits=c("R2","I2","I2∆afaBCDE","I2∆yadA-like"))
ggsave("Serum_mutants.tiff",width=6,height = 5,device = "tiff",dpi = 700)



#Zebrafish I1-R1
order=c("I1","R1","MG1655","CFT073")

dataCFU<-read.table("CFU_I1_R1_Mean.txt",header = T,sep = "\t")
dataCFUp<-read.table("CFU_I1_R1_Points.txt",header = T,sep = "\t")
orderhpi=c("0 hpi","6 hpi","24 hpi")
ggplot(dataCFU, aes(x=Time, y=Mean,fill=factor(Strain, levels = order))) + geom_bar(colour="black", stat="identity", position=position_dodge(), width=0.7)+theme_classic()+ geom_errorbar(aes(ymin=pmax(0,Mean-SD), ymax=Mean+SD), width=.2, position=position_dodge(0.7))+ geom_point(data=dataCFUp, aes(x=Time, y=CFU, fill=factor(Strain, levels = order)), position=position_dodge(width=.7)) + ylab("CFU / zebrafish")+ xlab("") + scale_fill_manual(values=c("firebrick4", "indianred3", "#D3D3D3", "azure4"), name="Strain")+ theme(text = element_text(size=15),axis.text.x=element_blank(),axis.ticks.x=element_blank())+facet_grid(.~factor(Time,levels=orderhpi),space = "free",scales = "free") +scale_y_continuous(trans = "log",expand=c(0,0), breaks = c(0, 1e+01, 1e+03, 1e+05, 1e+07), limits = c(NA, 1e+08))

ggsave("Zebra_CFU_I1-R1.tiff",width=7,height = 5,device = "tiff",dpi = 700)



#Zebrafish I2-R2 #1
order=c("I2","R2","MG1655","CFT073")

dataCFU<-read.table("CFU_I2_R2_Mean.txt",header = T,sep = "\t")
dataCFUp<-read.table("CFU_I2_R2_Points.txt",header = T,sep = "\t")
orderhpi=c("0 hpi","6 hpi","24 hpi")
ggplot(dataCFU, aes(x=Time, y=Mean,fill=factor(Strain, levels = order))) + geom_bar(colour="black", stat="identity", position=position_dodge(), width=0.7)+theme_classic()+ geom_errorbar(aes(ymin=pmax(0,Mean-SD), ymax=Mean+SD), width=.2, position=position_dodge(0.7))+ geom_point(data=dataCFUp, aes(x=Time, y=CFU, fill=factor(Strain, levels = order)), position=position_dodge(width=.7)) + ylab("CFU / zebrafish")+ xlab("") + scale_fill_manual(values=c("dodgerblue4", "dodgerblue2", "#D3D3D3", "azure4"), name="Strain")+ theme(text = element_text(size=15),axis.text.x=element_blank(),axis.ticks.x=element_blank())+facet_grid(.~factor(Time,levels=orderhpi),space = "free",scales = "free") +scale_y_continuous(trans = "log",expand=c(0,0), breaks = c(0, 1e+01, 1e+03, 1e+05, 1e+07), limits = c(NA, 1e+08))

ggsave("Zebra_CFU_I2-R2.tiff",width=7,height = 5,device = "tiff",dpi = 700)


#Zebrafish I2-R2 #2
order=c("I2","R2","MG1655","CFT073")

dataCFU<-read.table("CFU_I2_R2_manip2_Mean.txt",header = T,sep = "\t")
dataCFUp<-read.table("CFU_I2_R2_manip2_Points.txt",header = T,sep = "\t")
orderhpi=c("0 hpi","6 hpi","24 hpi")
ggplot(dataCFU, aes(x=Time, y=Mean,fill=factor(Strain, levels = order))) + geom_bar(colour="black", stat="identity", position=position_dodge(), width=0.7)+theme_classic()+ geom_errorbar(aes(ymin=pmax(0,Mean-SD), ymax=Mean+SD), width=.2, position=position_dodge(0.7))+ geom_point(data=dataCFUp, aes(x=Time, y=CFU, fill=factor(Strain, levels = order)), position=position_dodge(width=.7)) + ylab("CFU / zebrafish")+ xlab("") + scale_fill_manual(values=c("dodgerblue4", "dodgerblue2", "#D3D3D3", "azure4"), name="Strain")+ theme(text = element_text(size=15),axis.text.x=element_blank(),axis.ticks.x=element_blank())+facet_grid(.~factor(Time,levels=orderhpi),space = "free",scales = "free") +scale_y_continuous(trans = "log",expand=c(0,0), breaks = c(0, 1e+01, 1e+03, 1e+05, 1e+07), limits = c(NA, 1e+08))

ggsave("Zebra_CFU_I2-R2_2.tiff",width=7,height = 5,device = "tiff",dpi = 700)


#Zebrafish I3-R3
order=c("I3","R3","MG1655","CFT073")

dataCFU<-read.table("CFU_I3_R3_Mean.txt",header = T,sep = "\t")
dataCFUp<-read.table("CFU_I3_R3_Points.txt",header = T,sep = "\t")
orderhpi=c("0 hpi","6 hpi","24 hpi")
ggplot(dataCFU, aes(x=Time, y=Mean,fill=factor(Strain, levels = order))) + geom_bar(colour="black", stat="identity", position=position_dodge(), width=0.7)+theme_classic()+ geom_errorbar(aes(ymin=pmax(0,Mean-SD), ymax=Mean+SD), width=.2, position=position_dodge(0.7))+ geom_point(data=dataCFUp, aes(x=Time, y=CFU, fill=factor(Strain, levels = order)), position=position_dodge(width=.7)) + ylab("CFU / zebrafish")+ xlab("") + scale_fill_manual(values=c("seagreen4", "seagreen3", "#D3D3D3", "azure4"), name="Strain")+ theme(text = element_text(size=15),axis.text.x=element_blank(),axis.ticks.x=element_blank())+facet_grid(.~factor(Time,levels=orderhpi),space = "free",scales = "free") +scale_y_continuous(trans = "log",expand=c(0,0), breaks = c(0, 1e+01, 1e+03, 1e+05, 1e+07), limits = c(NA, 1e+08))

ggsave("Zebra_CFU_I3-R3.tiff",width=7,height = 5,device = "tiff",dpi = 700)


#Competitions Plk Biof
data<-read.table("Competition_Plk_Biof_Mean.txt",header = T,sep="\t")
datPoint<-read.table("Competition_Plk_Biof_Point.txt",header = T,sep = "\t")

ggplot(data, aes(x=Condition, y=Mean, fill=Competition)) + geom_bar(colour="black", stat="identity", position='dodge', width=0.7)+theme_classic()+ geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.2, position=position_dodge(0.7))+ geom_point(data=datPoint, aes(x=cond, y=fit, fill=strain), position=position_dodge(width=.7)) +scale_y_continuous(expand = c(0,0), limits = c(0,1.64)) + ylab("Fitness coefficient of relapse strain")+ scale_fill_manual(values=c("indianred3","dodgerblue2","seagreen3"), name="Competition")+ xlab("")+theme(text = element_text(size=20))+ geom_hline(yintercept=1, linetype="dashed", color = "black")
ggsave("Competitions_Plk-Biof.tiff",width=5,height = 7,device = "tiff",dpi = 700)


#Competitions I2-R2
data<-read.table("Competitions_Sucres_Mean.txt",header = T,sep="\t")
datPoint<-read.table("Competitions_Sucres_Points.txt",header = T,sep = "\t")

ggplot(data, aes(x=Condition, y=Mean, fill=Condition)) + geom_bar(colour="black", stat="identity", position='dodge', width=0.7)+theme_classic()+ geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.2, position=position_dodge(0.7))+ geom_point(data=datPoint, aes(x=Condition, y=Point, fill=Condition), position=position_dodge(width=.7)) +scale_y_continuous(expand = c(0,0), limits = c(0,2.05)) + ylab("Fitness coefficient of R2 relapse strain")+ scale_fill_manual(values=c("dodgerblue2","dodgerblue2","dodgerblue2","dodgerblue2","dodgerblue2"), name="Competition")+ theme(text = element_text(size=20),legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1))+ geom_hline(yintercept=1, linetype="dashed", color = "black") + scale_x_discrete(limits=c("LB", "LB + BIP","Glucose","Maltose","Lactate"))+xlab("")
ggsave("Competitions_R2.tiff",width=5,height = 7,device = "tiff",dpi = 700)



#Effect k7
dataP<-read.table("Effect_k7.txt",header = T,sep = "\t")

ggplot(dataP, aes(x=Competition, y=Fitness, fill=Competition)) + stat_summary(geom="bar",fun = mean,colour="black",  position='dodge', width=0.7)+ theme_linedraw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ stat_summary(geom="errorbar",fun.data = "mean_sd",aes(width=0.5))+ geom_point(position=position_dodge(width=.7)) +scale_y_continuous(expand = c(0,0), limits = c(0,1.5)) + xlab("")+ylab("Relative Fitness")+ scale_fill_manual(values=c("firebrick4", "seagreen4", "indianred3", "seagreen3"), name="Competition")+ theme(text = element_text(size=15),legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1))+ geom_hline(yintercept=1, linetype="dashed", color = "black")+scale_x_discrete(limits=c("I1-GFP vs I1-RFP","R1-GFP vs R1-RFP","I3-GFP vs I3-RFP","R3-GFP vs R3-RFP"))
ggsave("Effect-k7.tiff",width=5,height = 5,device = "tiff",dpi = 700)

#Agregation I1-R1
data<-read.table("Agregation_I1-R1.txt",header = T,sep="\t")
ggplot(data, aes(x=Time, y=OD, color=Sample)) +geom_line(size=1) + geom_point(aes(shape=Sample, size=0.1)) + theme_classic() + geom_errorbar(aes(ymin=OD-SD, ymax=OD+SD), width=.6, position=position_dodge(0)) + ylab(expression("%Initial OD"[600])) + xlab("Time (hours)") + theme(text = element_text(size=20)) +scale_x_continuous(expand = c(0,0)) +scale_color_manual(name="Sample",values=c("black", "darkgrey","firebrick4", "indianred3"),limits=c("∆flu","PCL-flu","I1","R1"))+guides(size=FALSE)+scale_shape_manual(name="Sample",labels=c("∆flu","PCL-flu","I1","R1"),values=c(15,16,17,18))

ggsave("Agregation_I1.tiff",width=7,height = 5,device = "tiff",dpi = 700)


#Agregation I2-R2
data<-read.table("Agregation_I2-R2.txt",header = T,sep="\t")
ggplot(data, aes(x=Time, y=OD, color=Sample)) +geom_line(size=1) + geom_point(aes(shape=Sample, size=0.1)) + theme_classic() + geom_errorbar(aes(ymin=OD-SD, ymax=OD+SD), width=.6, position=position_dodge(0)) + ylab(expression("%Initial OD"[600])) + xlab("Time (hours)") + theme(text = element_text(size=20)) +scale_x_continuous(expand = c(0,0)) +scale_color_manual(name="Sample",values=c("black", "darkgrey","dodgerblue4", "dodgerblue2"),limits=c("∆flu","PCL-flu","I2","R2"))+guides(size=FALSE)+scale_shape_manual(name="Sample",labels=c("∆flu","PCL-flu","I2","R2"),values=c(15,16,17,18))
ggsave("Agregation_I2-R2.tiff",width=7,height = 5,device = "tiff",dpi = 700)


#Agregation I3-R3
data<-read.table("Agregation_I3-R3.txt",header = T,sep="\t")
ggplot(data, aes(x=Time, y=OD, color=Sample)) +geom_line(size=1) + geom_point(aes(shape=Sample, size=0.1)) + theme_classic() + geom_errorbar(aes(ymin=OD-SD, ymax=OD+SD), width=.6, position=position_dodge(0)) + ylab(expression("%Initial OD"[600])) + xlab("Time (hours)") + theme(text = element_text(size=20)) +scale_x_continuous(expand = c(0,0)) +scale_color_manual(name="Sample",values=c("black", "darkgrey","seagreen4", "seagreen3"),limits=c("∆flu","PCL-flu","I3","R3"))+guides(size=FALSE)+scale_shape_manual(name="Sample",labels=c("∆flu","PCL-flu","I3","R3"),values=c(15,16,17,18))

ggsave("Agregation_I3-R3.tiff",width=7,height = 5,device = "tiff",dpi = 700)


#Growth curves
#LB
data<-read.table("GC_LB_I2-R2.txt",header = T,sep="\t")
ggplot(data, aes(x=Time, y=OD600, color=Strain)) +geom_line(size=1) + geom_point(aes(shape=Strain, size=0.1)) + theme_classic() + geom_errorbar(aes(ymin=OD600-SD, ymax=OD600+SD), width=.6, position=position_dodge(0)) + ylab(expression("OD"[600])) + xlab("Time (min)") + theme(plot.title = element_text(size=20,hjust = 0.5),text = element_text(size=20)) +scale_x_continuous(expand = c(0,0)) +scale_color_manual(name="Sample",values=c("dodgerblue4", "dodgerblue2"),limits=c("I2","R2"))+guides(size=FALSE)+scale_shape_manual(name="Sample",labels=c("I2","R2"),values=c(17,18))+scale_y_log10()+ggtitle("LB broth")

ggsave("GC_LB_I2-R2.tiff",width=7,height = 5,device = "tiff",dpi = 700)


#LB + BIP
data<-read.table("GC_BIP_I2-R2.txt",header = T,sep="\t")
ggplot(data, aes(x=Time, y=OD600, color=Strain)) +geom_line(size=1) + geom_point(aes(shape=Strain, size=0.1)) + theme_classic() + geom_errorbar(aes(ymin=OD600-SD, ymax=OD600+SD), width=.6, position=position_dodge(0)) + ylab(expression("OD"[600])) + xlab("Time (min)") + theme(plot.title = element_text(size=20,hjust = 0.5),text = element_text(size=20)) +scale_x_continuous(expand = c(0,0)) +scale_color_manual(name="Sample",values=c("dodgerblue4", "dodgerblue2"),limits=c("I2","R2"))+guides(size=FALSE)+scale_shape_manual(name="Sample",labels=c("I2","R2"),values=c(17,18))+scale_y_log10()+ggtitle("LB broth + 500µM 2-2'BIP")

ggsave("GC_BIP_I2-R2.tiff",width=7,height = 5,device = "tiff",dpi = 700)

#Glu
data<-read.table("GC_Glu_I2-R2.txt",header = T,sep="\t")
ggplot(data, aes(x=Time, y=OD600, color=Strain)) +geom_line(size=1) + geom_point(aes(shape=Strain, size=0.1)) + theme_classic() + geom_errorbar(aes(ymin=OD600-SD, ymax=OD600+SD), width=.6, position=position_dodge(0)) + ylab(expression("OD"[600])) + xlab("Time (min)") + theme(plot.title = element_text(size=20,hjust = 0.5),text = element_text(size=20)) +scale_x_continuous(expand = c(0,0)) +scale_color_manual(name="Sample",values=c("dodgerblue4", "dodgerblue2"),limits=c("I2","R2"))+guides(size=FALSE)+scale_shape_manual(name="Sample",labels=c("I2","R2"),values=c(17,18))+scale_y_log10()+ggtitle("M63B1 + 0.4% Glucose")

ggsave("GC_Glu_I2-R2.tiff",width=7,height = 5,device = "tiff",dpi = 700)

#Mal
data<-read.table("GC_Mal_I2-R2.txt",header = T,sep="\t")
ggplot(data, aes(x=Time, y=OD600, color=Strain)) +geom_line(size=1) + geom_point(aes(shape=Strain, size=0.1)) + theme_classic() + geom_errorbar(aes(ymin=OD600-SD, ymax=OD600+SD), width=.6, position=position_dodge(0)) + ylab(expression("OD"[600])) + xlab("Time (min)") + theme(plot.title = element_text(size=20,hjust = 0.5),text = element_text(size=20)) +scale_x_continuous(expand = c(0,0)) +scale_color_manual(name="Sample",values=c("dodgerblue4", "dodgerblue2"),limits=c("I2","R2"))+guides(size=FALSE)+scale_shape_manual(name="Sample",labels=c("I2","R2"),values=c(17,18))+scale_y_log10()+ggtitle("M63B1 + 0.4% Maltose")

ggsave("GC_Mal_I2-R2.tiff",width=7,height = 5,device = "tiff",dpi = 700)

#Lactate
data<-read.table("GC_Lactate_I2-R2.txt",header = T,sep="\t")
ggplot(data, aes(x=Time, y=OD600, color=Strain)) +geom_line(size=1) + geom_point(aes(shape=Strain, size=0.1)) + theme_classic() + geom_errorbar(aes(ymin=OD600-SD, ymax=OD600+SD), width=.6, position=position_dodge(0)) + ylab(expression("OD"[600])) + xlab("Time (min)") + theme(plot.title = element_text(size=20,hjust = 0.5),text = element_text(size=20)) +scale_x_continuous(expand = c(0,0)) +scale_color_manual(name="Sample",values=c("dodgerblue4", "dodgerblue2"),limits=c("I2","R2"))+guides(size=FALSE)+scale_shape_manual(name="Sample",labels=c("I2","R2"),values=c(17,18))+scale_y_log10()+ggtitle("M63B1 + 0.4% Lactate")

ggsave("GC_Lactate_I2-R2.tiff",width=7,height = 5,device = "tiff",dpi = 700)

#Pyr
data<-read.table("GC_Pyr_I1-R1.txt",header = T,sep="\t")
ggplot(data, aes(x=Time, y=OD600, color=Strain)) +geom_line(size=1) + geom_point(aes(shape=Strain, size=0.1)) + theme_classic() + geom_errorbar(aes(ymin=OD600-SD, ymax=OD600+SD), width=.6, position=position_dodge(0)) + ylab(expression("OD"[600])) + xlab("Time (min)") + theme(plot.title = element_text(size=20,hjust = 0.5),text = element_text(size=20)) +scale_x_continuous(expand = c(0,0)) +scale_color_manual(name="Sample",values=c("firebrick4", "indianred3"),limits=c("I1","R1"))+guides(size=FALSE)+scale_shape_manual(name="Sample",labels=c("I1","R1"),values=c(17,18))+scale_y_log10()+ggtitle("M63B1 + 0.4% Pyruvate")

ggsave("GC_Pyr_I1-R1.tiff",width=7,height = 5,device = "tiff",dpi = 700)


