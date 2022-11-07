
#color pallette

pal9 <- c("#920000","#db6d00","#24ff24","#009292","#004949","#b6dbff","#006ddb","#b66dff","#490092")
pal4 <- c("#920000","#db6d00","#009292","#490092")
pal2a<- c("#004949","#009292")
pal2b <- c("#6db6ff","#490092")

pal <- c("#000000","#004949","#009292","#ff6db6","#ffb6db",
         "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
         "#920000","#924900","#db6d00","#24ff24","#ffff6d")

pie(rep(1,2), col=pal2b)


# Day 1
MunidaD1
Surf_MunidaD1 <- subset_samples(MunidaD1, !Location=="9")


### How do things separate out?
###################################################################################################
## nMDS

# OTUs

Mun_D1_OTUs <- ordinate(Surf_MunidaD1, "NMDS", "jaccard", binary=TRUE, trymax = 500)
Mun_D1_OTU_plot = plot_ordination(Surf_MunidaD1, Mun_D1_OTUs, type="samples", color="Samp",shape="WaterType")
Mun_D1_OTU_plot + geom_point(size=5) + ggtitle("nMDS of surface samples - Day 1 - OTUs")+theme(plot.title = element_text(size=30))+ scale_color_manual(values=pal9)+theme_bw()
#+stat_ellipse(type = "norm",linetype = 5)+ facet_wrap(Day~.)

# genus

MunidaD1_genus = tax_glom(Surf_MunidaD1, "Genus")
phylaGFilter = c("", "g__")
MunidaD1_genus = subset_taxa(MunidaD1_genus, !Genus %in% phylaGFilter)

Mun_D1_G <- ordinate(MunidaD1_genus, "NMDS", "jaccard", binary=TRUE, trymax = 500)
Mun_D1_G_plot = plot_ordination(MunidaD1_genus, Mun_D1_G, type="samples", color="Samp",shape="WaterType")
Mun_D1_G_plot + geom_point(size=5) + ggtitle("nMDS of surface samples - Day 1 - Genus")+theme(plot.title = element_text(size=30))+ scale_color_manual(values=pal9)+theme_bw()
#+stat_ellipse(type = "norm",linetype = 5)+ facet_wrap(Day~.)

##################################################################
## Permanova for Day 1

library(vegan)

# Run a PERMANOVA on species but with Jaccard distances for OTUs
MunidaSS_df <- data.frame(sample_data(Surf_MunidaD1))
MunidaSS_jac <- phyloseq::distance(Surf_MunidaD1, method = "jaccard", binary=T)
MunidaSS_adonis_J <- adonis2(MunidaSS_jac ~ WaterType, data = MunidaSS_df)

MunidaSS_adonis_J   # significant @ F(3,39)=10.526, p<0.001

#and permDisp
groups <- MunidaSS_df[["WaterType"]]
class(groups)
groups
MunidaSS_PermDisp <- betadisper(MunidaSS_jac, groups)

anova(MunidaSS_PermDisp)      # significant @ F(3,39)=37.525, p<0.00000001

# plot
PD_SurfSamp_J<-plot(MunidaSS_PermDisp,main="PERMDISP for Day 1 surface-level OTUs", col=pal4)
### This is the jaccard, species level taxonomy, with depth 9 sample graph


boxplot(MunidaSS_PermDisp)
MunidaSS_PermDisp.HSD <- TukeyHSD(MunidaSS_PermDisp)
plot(MunidaSS_PermDisp.HSD)


# Run a PERMANOVA on species but with Jaccard distances for genus
MunidaSS_g_df <- data.frame(sample_data(MunidaD1_genus))
MunidaSS_g_jac <- phyloseq::distance(MunidaD1_genus, method = "jaccard", binary=T)
MunidaSS_g_adonis_J <- adonis2(MunidaSS_g_jac ~ WaterType, data = MunidaSS_g_df)

MunidaSS_g_adonis_J   # significant @ F(3,39)=10.62, p<0.001

#and permDisp
groups <- MunidaSS_g_df[["WaterType"]]
class(groups)
groups
MunidaSS_g_PermDisp <- betadisper(MunidaSS_g_jac, groups)

anova(MunidaSS_g_PermDisp)      # significant @ F(4,79)=3.14, p<0.019

# plot

SM_D1_PCoA = ordinate(Surf_MunidaD1, "PCoA", "jaccard", binary=TRUE)
p3<-plot_ordination(Surf_MunidaD1, SM_D1_PCoA, color="Samp", shape="WaterType")+ geom_point(size=5)+ scale_color_manual(values=pal9)+theme_bw(base_size=25)+scale_x_reverse()

SM_D1_PCoA_G <- ordinate(MunidaD1_genus, "PCoA", "jaccard", binary=TRUE)
p4<-plot_ordination(Surf_MunidaD1, SM_D1_PCoA_G, color="Samp", shape="WaterType")+ geom_point(size=5)+ scale_color_manual(values=pal9)+theme_bw(base_size=25)

p3
p4

PD_SurfSamp_g_J<-plot(MunidaSS_g_PermDisp,main="PERMDISP for Day 1 surface-level Genus", col=pal4)   ### This is the jaccard, species level taxonomy, with depth 9 sample graph



boxplot(MunidaSS_sp_PermDisp)
permutest(MunidaSS_sp_PermDisp)
TukeyHSD(MunidaSS_sp_PermDisp)
MunidaSS_sp_PermDisp.HSD <- TukeyHSD(MunidaSS_sp_PermDisp)
plot(MunidaSS_sp_PermDisp.HSD)


## plots only for day 1

library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

p1<-Mun_D1_OTU_plot + geom_point(size=5) +theme(plot.title = element_text(size=30))+ scale_color_manual(values=pal9)+theme_bw(base_size=25)
p2<-Mun_D1_G_plot + geom_point(size=5)+theme(plot.title = element_text(size=30))+ scale_color_manual(values=pal9)+theme_bw(base_size=25)

SM_D1_PCoA = ordinate(Surf_MunidaD1, "PCoA", "jaccard", binary=TRUE)
p3<-plot_ordination(Surf_MunidaD1, SM_D1_PCoA, color="Samp", shape="WaterType")+ geom_point(size=5)+ scale_color_manual(values=pal9)+theme_bw(base_size=25)+scale_x_reverse()

SM_D1_PCoA_G <- ordinate(MunidaD1_genus, "PCoA", "jaccard", binary=TRUE)
p4<-plot_ordination(Surf_MunidaD1, SM_D1_PCoA_G, color="Samp", shape="WaterType")+ geom_point(size=5)+ scale_color_manual(values=pal9)+theme_bw(base_size=25)

p3
p4


grid_arrange_shared_legend_D1 <-
  function(...,
           ncol = 2,
           nrow = 2,
           position = c("right")) {

    plots <- list(p1,p2,p3,p4)
    position <- match.arg(position)
    g <-
      ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x)
      x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x)
      x + theme(legend.position = "none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)

    combined <- switch(
      position,
      "bottom" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 2,
        heights = unit.c(unit(1, "npc") - lheight, lheight)
      ),
      "right" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 2,
        widths = unit.c(unit(1, "npc") - lwidth, lwidth)
      )
    )

    grid.newpage()
    grid.draw(combined)

    # return gtable invisibly
    invisible(combined)

  }

grid_arrange_shared_legend_D1(p1, p2, p3, p4)



## PCoA only of Day 1



grid_arrange_shared_legend_D2 <-
  function(...,
           ncol = 2,
           nrow = 1,
           position = c("right")) {
    
    plots <- list(p3,p4)
    position <- match.arg(position)
    g <-
      ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x)
      x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x)
      x + theme(legend.position = "none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)
    
    combined <- switch(
      position,
      "bottom" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 2,
        heights = unit.c(unit(1, "npc") - lheight, lheight)
      ),
      "right" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 2,
        widths = unit.c(unit(1, "npc") - lwidth, lwidth)
      )
    )
    
    grid.newpage()
    grid.draw(combined)
    
    # return gtable invisibly
    invisible(combined)
    
  }

grid_arrange_shared_legend_D2(p3, p4)

## other permdisp plots

par(mfrow=c(1,2))

p3<-PD_SurfSamp_J<-plot(MunidaSS_PermDisp,main="Day1", col=pal4)
p4<-PD_SurfSamp_g_J<-plot(MunidaSS_g_PermDisp,col=pal4)

par(mfrow=c(1,1))



###################################################################################################
# Day 2
MunidaD2
Surf_MunidaD2 <- subset_samples(MunidaD2, !Location=="9")
### How do things separate out?

## nMDS

# OTUs

Mun_D2_OTUs <- ordinate(Surf_MunidaD2, "NMDS", "jaccard", binary=TRUE, trymax = 1500)
Mun_D2_OTU_plot = plot_ordination(Surf_MunidaD2, Mun_D2_OTUs, type="samples", color="Samp",shape="WaterType")
Mun_D2_OTU_plot + geom_point(size=5) + ggtitle("nMDS of surface samples - Day 2 - OTUs")+theme(plot.title = element_text(size=30))+ scale_color_manual(values=pal9)+theme_bw()
#+stat_ellipse(type = "norm",linetype = 5)+ facet_wrap(Day~.)

# genus

MunidaD2_genus = tax_glom(Surf_MunidaD2, "Genus")
phylaGFilter = c("", "g__")
MunidaD2_genus = subset_taxa(MunidaD2_genus, !Genus %in% phylaGFilter)

Mun_D2_G <- ordinate(MunidaD2_genus, "NMDS", "jaccard", binary=TRUE, trymax = 500)
Mun_D2_G_plot = plot_ordination(MunidaD2_genus, Mun_D2_G, type="samples", color="Samp",shape="WaterType")
Mun_D2_G_plot + geom_point(size=5) + ggtitle("nMDS of surface samples - Day 2 - Genus")+theme(plot.title = element_text(size=30))+ scale_color_manual(values=pal9)+theme_bw()
#+stat_ellipse(type = "norm",linetype = 5)+ facet_wrap(Day~.)

##################################################################
## Permanova for Day 2

## Permanova

# Run a PERMANOVA on species but with Jaccard distances for OTUs
Munida_D2_df <- data.frame(sample_data(Surf_MunidaD2))
Munida_D2_jac <- phyloseq::distance(Surf_MunidaD2, method = "jaccard", binary=T)
Munida_D2_adonis_J <- adonis2(Munida_D2_jac ~ WaterType, data = Munida_D2_df)

Munida_D2_adonis_J   # significant @ F(3,39)=10.545, p<0.001

#and permDisp
groups <- Munida_D2_df[["WaterType"]]
class(groups)
groups
Munida_D2_PermDisp <- betadisper(Munida_D2_jac, groups)

anova(Munida_D2_PermDisp)      # significant @ F(3,39)=51.211, p<0.00000001

# plot
D2_SurfSamp_J<-plot(Munida_D2_PermDisp,main="PERMDISP for Day 2 surface-level OTUs")   ### This is the jaccard, species level taxonomy, with depth 9 sample graph

boxplot(Munida_D2_PermDisp)
Munida_D2_PermDisp.HSD <- TukeyHSD(Munida_D2_PermDisp)
plot(Munida_D2_PermDisp.HSD)


# Run a PERMANOVA on species but with Jaccard distances for genus
Munida_D2_g_df <- data.frame(sample_data(MunidaD2_genus))
Munida_D2_g_jac <- phyloseq::distance(MunidaD2_genus, method = "jaccard", binary=T)
Munida_D2_g_adonis_J <- adonis2(Munida_D2_g_jac ~ WaterType, data = Munida_D2_g_df)

Munida_D2_g_adonis_J   # significant @ F(3,39)=12.88, p<0.001

#and permDisp
groups <- Munida_D2_g_df[["WaterType"]]
class(groups)
groups
Munida_D2_g_PermDisp <- betadisper(Munida_D2_g_jac, groups)

anova(Munida_D2_g_PermDisp)      # significant @ F(4,79)=4.7244, p<0.007015

# plot
D2_SurfSamp_g_J<-plot(Munida_D2_g_PermDisp,main="PERMDISP for Day 2 surface-level Genus")   ### This is the jaccard, species level taxonomy, with depth 9 sample graph
boxplot(Munida_D2_sp_PermDisp)
permutest(Munida_D2_sp_PermDisp)
TukeyHSD(Munida_D2_sp_PermDisp)
Munida_D2_sp_PermDisp.HSD <- TukeyHSD(Munida_D2_sp_PermDisp)
plot(MunidaSS_sp_PermDisp.HSD)



## plots only ##

p5<-Mun_D2_OTU_plot + geom_point(size=5)+theme(plot.title = element_text(size=30))+ scale_color_manual(values=pal9)+theme_bw(base_size=25)
p6<-Mun_D2_G_plot + geom_point(size=5) +theme(plot.title = element_text(size=30))+ scale_color_manual(values=pal9)+theme_bw(base_size=25)


SM_D2_PCoA = ordinate(Surf_MunidaD2, "PCoA", "jaccard", binary=TRUE)
p7<-plot_ordination(Surf_MunidaD2, SM_D2_PCoA, color="Samp", shape="WaterType")+ geom_point(size=5)+ scale_color_manual(values=pal9)+theme_bw(base_size=25)

SM_D2_PCoA_G <- ordinate(MunidaD2_genus, "PCoA", "jaccard", binary=TRUE)
p8<-plot_ordination(Surf_MunidaD2, SM_D2_PCoA_G, color="Samp", shape="WaterType")+ geom_point(size=5)+ scale_color_manual(values=pal9)+theme_bw(base_size=25)+scale_y_reverse()

p7
p8

## all

grid_arrange_shared_legend_D3 <-
  function(...,
           ncol = 2,
           nrow = 2,
           position = c("right")) {

    plots <- list(p5,p6, p7, p8)
    position <- match.arg(position)
    g <-
      ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x)
      x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x)
      x + theme(legend.position = "none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)

    combined <- switch(
      position,
      "bottom" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 2,
        heights = unit.c(unit(1, "npc") - lheight, lheight)
      ),
      "right" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 2,
        widths = unit.c(unit(1, "npc") - lwidth, lwidth)
      )
    )

    grid.newpage()
    grid.draw(combined)

    # return gtable invisibly
    invisible(combined)

  }

grid_arrange_shared_legend_D3(p5, p6, p7, p8)



## only PCoA


grid_arrange_shared_legend_D4 <-
  function(...,
           ncol = 2,
           nrow = 1,
           position = c("right")) {
    
    plots <- list(p7, p8)
    position <- match.arg(position)
    g <-
      ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x)
      x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x)
      x + theme(legend.position = "none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)
    
    combined <- switch(
      position,
      "bottom" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 2,
        heights = unit.c(unit(1, "npc") - lheight, lheight)
      ),
      "right" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 2,
        widths = unit.c(unit(1, "npc") - lwidth, lwidth)
      )
    )
    
    grid.newpage()
    grid.draw(combined)
    
    # return gtable invisibly
    invisible(combined)
    
  }

grid_arrange_shared_legend_D4(p7, p8)













par(mfrow=c(1,2))

p7<-plot(Munida_D2_PermDisp,main="PERMDISP for Day 2 surface-level OTUs",col=pal4, font = 6)
p8<-plot(Munida_D2_g_PermDisp,main="PERMDISP for Day 2 surface-level Genus",col=pal4)

par(mfrow=c(1,1))

p7

ordu = ordinate(Surf_MunidaD2, "PCoA", "jaccard", binary=TRUE)
plot_ordination(Surf_MunidaD2, ordu, color="Samp", shape="WaterType")+ geom_point(size=5)+ scale_color_manual(values=pal9)+theme_bw(base_size=25)










###################################################################################################
#################################   Differences between days     #################################
###################################################################################################

Munida
Surface_Munida <- subset_samples(Munida, !Location=="9")

MunidaBothDays_df <- data.frame(sample_data(Surface_Munida))
MunidaBothDays_jac <- phyloseq::distance(Surface_Munida, method = "jaccard", binary=T)
MunidaBothDays_adonis_day_J <- adonis2(MunidaBothDays_jac ~ Day, data = MunidaBothDays_df)

MunidaBothDays_adonis_day_J ## Significant @ F(1,79)=10.402; p < 0.001

#and permDisp
groups <- MunidaBothDays_df[["Day"]]
class(groups)
groups
MunidaBothDays_PermDisp_day <- betadisper(MunidaBothDays_jac, groups)

anova(MunidaBothDays_PermDisp_day) ## non-siginficant, difference is only a grouping effect of day
## insignificant @ F(1,79)=0.0099, p=0.9212

MunidaBothDays_PermDisp_plot_day<-plot(MunidaBothDays_PermDisp_day,main="PERMDISP by day",col=pal2b) # very different


SM_pcoa = ordinate(Surface_Munida, "PCoA", "jaccard", binary=TRUE)
plot_ordination(Surface_Munida, SM_pcoa, color="Day", shape="WaterType")+ geom_point(size=5)+ scale_color_manual(values=pal9)+theme_bw(base_size=25)







###################################################################################################
#################################   surface vs depth      #################################
###################################################################################################

# surface v depth - all otus
SvD_Munida<-subset_samples(Munida, Distance=="8")


Mun_SvD_OTUs <- ordinate(SvD_Munida, "NMDS", "jaccard", binary=TRUE, trymax = 5000)
Mun_SvD_OTU_plot = plot_ordination(SvD_Munida, Mun_SvD_OTUs, type="samples", color="Day",shape="WaterType")
Mun_SvD_OTU_plot + geom_point(size=5) + ggtitle("nMDS of Surface v Depth - OTUs")+theme(plot.title = element_text(size=30))
#+stat_ellipse(type = "norm",linetype = 5)+ facet_wrap(Day~.)

Munida_SvD_genus = tax_glom(SvD_Munida, "Genus")
phylaGFilter = c("", "g__")
#Munida_SvD_genus = subset_taxa(Munida_SvD_genus, !Genus %in% phylaGFilter)

Munida_SvD_genus_ord <- ordinate(Munida_SvD_genus, "NMDS", "jaccard", binary=TRUE, trymax = 500)
Munida_SvD_genus_ord_plot = plot_ordination(Munida_SvD_genus, Munida_SvD_genus_ord, type="samples", color="WaterType",shape="Day")
Munida_SvD_genus_ord_plot + geom_point(size=5) + ggtitle("nMDS of Surface v Depth - Genus")+theme(plot.title = element_text(size=30))

# SS for surface samples
# Run a PERMANOVA on species but with Jaccard distances
Munida_SubAnt_df <- data.frame(sample_data(SvD_Munida))
Munida_SubAnt_jac <- phyloseq::distance(SvD_Munida, method = "jaccard", binary=T)
Munida_SubAnt_adonis_J <- adonis2(Munida_SubAnt_jac ~ Samp, data = Munida_SubAnt_df)

Munida_SubAnt_adonis_J   # significant @ F(1,19)=10.022, p<0.001

#and permDisp
groups <- Munida_SubAnt_df[["Samp"]]
class(groups)
groups
Munida_SubAnt_PermDisp <- betadisper(Munida_SubAnt_jac, groups)

anova(Munida_SubAnt_PermDisp)      # significant @ F(1,19)=35.704, p<0.00001185

# plot
Munida_SubAnt_J<-plot(Munida_SubAnt_PermDisp,main="PERMDISP surface v depth OTUs")   ### This is the jaccard, species level taxonomy, with depth 9 sample graph
boxplot(MunidaSS_PermDisp)
MunidaSS_PermDisp.HSD <- TukeyHSD(MunidaSS_PermDisp)
plot(MunidaSS_PermDisp.HSD)


#
# Run a PERMANOVA on species but with Jaccard distances
Munida_SubAnt_g_df <- data.frame(sample_data(Munida_SvD_genus))
Munida_SubAnt_g_jac <- phyloseq::distance(Munida_SvD_genus, method = "jaccard", binary=T)
Munida_SubAnt_adonis_g_J <- adonis2(Munida_SubAnt_g_jac ~ WaterType, data = Munida_SubAnt_g_df)

Munida_SubAnt_adonis_g_J   # significant @ F(1,19)=12.69, p<0.001

#and permDisp
groups <- Munida_SubAnt_g_df[["WaterType"]]
class(groups)
groups
Munida_SubAnt_PermDisp_g <- betadisper(Munida_SubAnt_g_jac, groups)

anova(Munida_SubAnt_PermDisp_g)      # significant @ F(1,19)=35.704, p<0.00001185

# plot
Munida_SubAnt_J<-plot(Munida_SubAnt_PermDisp_g,main="PERMDISP surface v depth Genus")   ### This is the jaccard, species level taxonomy, with depth 9 sample graph
boxplot(MunidaSS_PermDisp)
MunidaSS_PermDisp.HSD <- TukeyHSD(MunidaSS_PermDisp)
plot(MunidaSS_PermDisp.HSD)




##plot depth


Munida_SvD_genus_ord_plot + geom_point(size=5) + ggtitle("nMDS of Surface v Depth - Genus")+theme(plot.title = element_text(size=30))+ scale_color_manual(values=pal2b)+theme_bw(base_size=25)



par(mfrow=c(1,2))


Munida_SubAnt_J<-plot(Munida_SubAnt_PermDisp,main="PERMDISP surface v depth OTUs", col=pal2b)
Munida_SubAnt_J<-plot(Munida_SubAnt_PermDisp_g,main="PERMDISP surface v depth Genus", col=pal2b)

par(mfrow=c(1,1))



SvD_Munida<-subset_samples(Munida, Distance=="8")
SvD_Munida_D1<-subset_samples(MunidaD1, Distance=="8")
SvD_Munida_D2<-subset_samples(MunidaD2, Distance=="8")

SvD_pcoa_D1 = ordinate(SvD_Munida_D1, "PCoA", "jaccard", binary=TRUE)
pd1<-plot_ordination(SvD_Munida_D1, SvD_pcoa_D1 , color="WaterType", shape="WaterType")+ geom_point(size=5)+ scale_color_manual(values=pal2b)+theme_bw(base_size=25)
SvD_pcoa_D2 = ordinate(SvD_Munida_D2, "PCoA", "jaccard", binary=TRUE)
pd2<-plot_ordination(SvD_Munida_D2, SvD_pcoa_D2 , color="WaterType", shape="WaterType")+ geom_point(size=5)+ scale_color_manual(values=pal2b)+theme_bw(base_size=25)


par(mfrow=c(1,2))
pd1
pd2
par(mfrow=c(1,1))


## Graphing PCoA of samples - Surface vs Depth

SvD_Munida<-subset_samples(Munida, Distance=="8")
SvD_pcoa = ordinate(SvD_Munida, "PCoA", "jaccard", binary=TRUE)

## both days together, surface v depth, OTUs

pd2<-plot_ordination(SvD_Munida, SvD_pcoa, color="Day", shape="WaterType")+ geom_point(size=5)+ scale_color_manual(values=pal2b)+theme_bw(base_size=25)
pd2

# split out by day

pd<-plot_ordination(SvD_Munida, SvD_pcoa, color="WaterType", shape="WaterType")+ geom_point(size=5)+ scale_color_manual(values=pal2b)+theme_bw(base_size=25)+facet_wrap("Day")
pd

## both days together, Surface v depth, Genus

SvD_pcoa_genus = ordinate(Munida_SvD_genus, "PCoA", "jaccard", binary=TRUE)
pd3<-plot_ordination(Munida_SvD_genus, SvD_pcoa_genus, color="Day", shape="WaterType")+ geom_point(size=5)+ scale_color_manual(values=pal2b)+theme_bw(base_size=25)
pd3






###################################################################################################
#################################   Indicator Species Analyses     #################################
###################################################################################################

phylaSFilter = c("", "s__")
Munida_species = subset_taxa(MunidaD1, !Species %in% phylaSFilter)
tax_table(Munida_species)
OTU_Tab_IndVal<-otu_table(Munida_species)

Munida_otu<-Munida_species@otu_table  # get the otu table out of the phyloseq object
Munida_otu_T<-t(Munida_otu)            # transpose (rotate) the otu table
Munida_otu_T[1:4,]                  # double check
metadata<-read.csv("Fixedmetadata2.csv")  # read in metadata
metadata2<-subset(metadata, JulianDate!=0)
metadata2[85:90,]                  # double check

library(tibble)                 # change row to an actual column and name it
Munida_otu_T_2<-data.frame(Munida_otu_T)
Munida_otu_T_3<-rownames_to_column(Munida_otu_T_2, var="SampleID")

# Join based on sample ID
Munida_otu_T_4 <- left_join(metadata2, Munida_otu_T_3, by = c("SampleID" = "SampleID"))
Munida_otu_T_4<-subset(Munida_otu_T_4, JulianDate!=0) # get rid of negatives
Mun_D1_IndicSpe<-na.omit(Munida_otu_T_4)
# Load indicspecies package
library(indicspecies); packageVersion('indicspecies')



# set indicipieces factors
abund = Mun_D1_IndicSpe[,20:112] #select just the ASV table part of the file
wat = Mun_D1_IndicSpe$WaterType
#IndicatorSpecies = multipatt(abund, wat, func = "IndVal.g", control = how(nperm=9999))
#summary(IndicatorSpecies, indvalcomp = TRUE)

# without group combinations
IndicatorSpecies = multipatt(abund, wat, duleg = TRUE, func = "IndVal.g", control = how(nperm=9999))
summary(IndicatorSpecies, indvalcomp = TRUE)



### Correction for a p-value
library(data.table)

#extract table of stats
indisp.sign<-as.data.table(IndicatorSpecies$sign, keep.rownames=TRUE)
#add adjusted p-value
indisp.sign[ ,p.value.b:=p.adjust(p.value, method="bonferroni")]
#now can select only the indicators with adjusted significant p-values
indisp.sign[p.value.b<=0.05, ]


#### Day 2

phylaSFilter = c("", "s__")
Munida_species = subset_taxa(MunidaD2, !Species %in% phylaSFilter)
tax_table(Munida_species)
OTU_Tab_IndVal<-otu_table(Munida_species)

Munida_otu<-Munida_species@otu_table  # get the otu table out of the phyloseq object
Munida_otu_T<-t(Munida_otu)            # transpose (rotate) the otu table
Munida_otu_T[1:4,]                  # double check
metadata<-read.csv("Fixedmetadata2.csv")  # read in metadata
metadata2<-subset(metadata, JulianDate!=0)
metadata2[85:90,]                  # double check

library(tibble)                 # change row to an actual column and name it
Munida_otu_T_2<-data.frame(Munida_otu_T)
Munida_otu_T_3<-rownames_to_column(Munida_otu_T_2, var="SampleID")

# Join based on sample ID
Munida_otu_T_4 <- left_join(metadata2, Munida_otu_T_3, by = c("SampleID" = "SampleID"))
Munida_otu_T_4<-subset(Munida_otu_T_4, JulianDate!=0) # get rid of negatives
Mun_D2_IndicSpe<-na.omit(Munida_otu_T_4)
# Load indicspecies package
library(indicspecies); packageVersion('indicspecies')

# without group combinations
abund2 = Mun_D2_IndicSpe[,20:112] #select just the ASV table part of the file
wat2 = Mun_D2_IndicSpe$WaterType
IndicatorSpecies2 = multipatt(abund2, wat2, duleg = TRUE, func = "IndVal.g", control = how(nperm=9999))
summary(IndicatorSpecies2, indvalcomp = TRUE)

#extract table of stats
indisp.sign<-as.data.table(IndicatorSpecies2$sign, keep.rownames=TRUE)
#add adjusted p-value
indisp.sign[ ,p.value.b:=p.adjust(p.value, method="bonferroni")]
#now can select only the indicators with adjusted significant p-values
indisp.sign[p.value.b<=0.05, ]