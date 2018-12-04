data.del<-data[-which(data$plot=="L2-19" & data$TimeStep=="0")]
data.del<-data.del[-which(data.del$plot=="L2-23" & data.del$TimeStep=="0")]
data.del<-data.del[-which(data.del$plot=="L2-24" & data.del$TimeStep=="0")]
data.del<-data.del[-which(data.del$plot=="L2-26" & data.del$TimeStep=="0")]

library(vegan)
names(data.del)
community<-data.del[,6:16]
id.tags<-data.del[,c(1:4,17:18)]

community.mds<-metaMDS(community,distance="bray",trymax=50,k=3)
### check stress plot
stressplot(community.mds)


t0<-c(which(data.del$TimeStep==0))
t1<-c(which(data.del$TimeStep==1))
t4<-c(which(data.del$TimeStep==4))

community0.mds3<-metaMDS(community[t0,],k=3,trymax=100)
stressplot(community0.mds3)

treat<-id.tags$removal[t0]
colvec <- c("green4", "red2", "mediumblue","orange")
pchs<-c(2,8,16)
gr.removal<-factor(data$removal)
gr.time<-factor(data$TimeStep)

ordiplot(community0.mds3,type="n",main="Initial")
ordiellipse(community0.mds3, display = "sites", groups=treat, col = colvec, kind="sd",lty=2)
ordiellipse(community0.mds3, display = "sites", groups=treat, col = colvec, kind="se")
orditorp(community0.mds3,display="species",col="red",air=0.01)
legend("topright", legend = levels(gr.removal), bty = "n", col = colvec, pch=c(20),title="Removal")


MDS1<-community.mds$points[,1]
MDS2<-community.mds$points[,2]
MDS3<-community.mds$points[,3]

NMDS <- data.frame(MDS1 = MDS1, MDS2 = MDS2, MDS3 = MDS3, 
                   Plot = data.del$plot, Site = data.del$site, TimeStep = data.del$TimeStep,
                   Removal = data.del$removal, Pred = data.del$pred)


library(ggplot2)

# Site groupings
ggplot(NMDS[t0,], aes(x=MDS1, y=MDS2, col=Site)) +
  geom_point() +
  stat_ellipse() +
  theme_bw() +
  labs(title = "Initial (pre-treatment)")

ggplot(NMDS[t0,], aes(x=MDS1, y=MDS2, col=Removal)) +
  geom_point() +
  stat_ellipse() +
  theme_bw() +
  labs(title = "Initial (pre-treatment), restricted")

ggplot(NMDS[t0,], aes(x=MDS1, y=MDS2, col=Pred)) +
  geom_point() +
  stat_ellipse() +
  theme_bw() +
  labs(title = "Initial (pre-treatment), restricted")

centroids <- aggregate(cbind(MDS1,MDS2,MDS3)~Removal*Pred*TimeStep,NMDS,mean)

ggplot(NMDS[t0,], aes(x=MDS1, y=MDS2, col=Removal,pch=Pred)) +
  scale_color_manual(values = c("#7CAE00", "#F8766D", "#00BFC4", "#C77CFF"))+
  scale_linetype_manual(values = c(1,3,2))+
  stat_ellipse(aes(linetype=Pred)) +
  theme_bw() +
  labs(title = "Initial (pre-treatment), restricted")+
  geom_point(data=centroids[which(centroids$TimeStep==0),],size=2)

t<-t0 # choose from c(t0, t1, t4)

adonis.t0<-adonis(vegdist(data.del[t,c(6:16)],method="bray")
                    ~pred*removal, strata=data.del$site[t], data=data.del[t], permutations=9999)

adonis.t0
