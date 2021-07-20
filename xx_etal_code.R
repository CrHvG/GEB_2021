# XX et al. Climate mediates the relationship between plant diversity and forest structure across the United States. GEB submission
  
# misc. functions ------
  
  my_theme <- function(base_size = 12, base_family = "",mycolor_background="black",plot_background="white",mycolor_elements="white",mycolor_text="white",font="Times") {
    theme_grey(base_size = base_size, base_family = base_family) %+replace%
      theme(axis.line = element_blank(),  
        axis.text.x = element_text(family=font,size = base_size*0.6, color = mycolor_text, lineheight = 0.9,angle=25, hjust = 1),  
        axis.text.y = element_text(family=font,size = base_size*0.8, color = mycolor_text, lineheight = 0.9),  
        axis.ticks = element_line(color = mycolor_elements, size  =  0.2),  
        axis.title.x = element_text(family=font,size = base_size, color = mycolor_text, margin = ggplot2::margin(0, 10, 0, 0)),  
        axis.title.y = element_text(family=font,size = base_size, color = mycolor_text, angle = 90, margin = ggplot2::margin(0, 10, 0, 0)),  
        axis.ticks.length = unit(0.3, "lines"),   
        legend.background = element_rect(color = plot_background, fill = plot_background),  
        legend.key = element_rect(color = plot_background,  fill = plot_background),  
        legend.key.size = unit(1.2, "lines"),  
        legend.key.height = NULL,  
        legend.key.width = NULL,      
        legend.text = element_text(family=font,size = base_size*0.8, color = mycolor_text),  
        legend.title = element_text(family=font,size = base_size*0.8, face = "bold", hjust = 0, color = mycolor_text),  
        legend.position = "right",  
        legend.text.align = NULL,  
        legend.title.align = NULL,  
        legend.direction = "vertical",  
        legend.box = NULL, 
        panel.background = element_rect(fill = mycolor_background, color  =  NA),  
        panel.border = element_rect(fill = NA, color = mycolor_elements),  
        panel.grid.major = element_line(color = mycolor_background),  
        panel.grid.minor = element_blank(),  
        panel.spacing = unit(0.5, "lines"),   
        plot.background = element_rect(color = plot_background, fill = plot_background),  
        plot.title = element_text(size = base_size*1.2, color = mycolor_text,face="bold",family=font),  
        plot.margin = unit(rep(1, 4), "lines")
      )
  }
  
  
# Fig. 3 Figure 3. Vascular plant and tree diversity covaries with forest structural conditions among NEON plots and sites -----

  # load data
    library(MASS)
    full_data<-read.csv("/... xx_etal_data.csv")
    colors<-read.csv("/... xx_etal_site_colors.csv")

  # 3a analysis -----------

    div_vars<-c("plant_SR_400","plant_H_400","tree_SR_400","tree_H_400")  
    strct_vars<-c("CHM_nni_3_2512m2","pc_skew_DBH_2512m2","z_ent1_ent2_2512m2","z_1m_mode_mean_2512m2","pc_ent1_01_2512m2","z_cc_1_3_2512m2")
    
    plot_list<-list()  
    for (strct_var in strct_vars){
      for (div_var in div_vars){
        tmp<-full_data[,c("siteID",div_var,strct_var,"latitude","longitude")]
        tmp<-tmp[complete.cases(tmp),]
        names(tmp)<-c("siteID","div","strct","latitude","longitude")
        tmp<-tmp[tmp$siteID%in%names(table(tmp$siteID)[table(tmp$siteID)>=9]),]
        site_means<-aggregate(.~siteID,data=tmp,FUN=mean, na.rm = TRUE)
        
        colors$siteID<-as.character(colors$siteID)
        colors$color<-as.character(colors$color)
        bubble_colors<-colors[which(colors$siteID%in%tmp$siteID),"color"]
        label_colors<-rep("black",length(bubble_colors))
        label_colors[which(site_means$siteID%in%c("BART","BONA","DEJU","HEAL","NIWO","RMNP","STEI","UNDE"))]<-"white"
  
        ifelse(div_var=="tree_SR_400",mymax<-21,mymax<-max(tmp$div[!is.na(tmp$div)]))
        
        # run GLMM models
          if(div_var=="tree_SR_400"){
            glm_re_spatial_freq<- glmmPQL(div ~ strct,random=~1|siteID,data=tmp,family=poisson,correlation=corExp(form=~longitude+latitude,nugget = TRUE))
          }
          if(div_var=="plant_SR_400"){
            m0 <- glm.nb(div ~ strct + siteID , data=tmp)
            theta<-m0$theta
            glm_re_spatial_freq<- glmmPQL(div ~ strct,random=~1|siteID,data=tmp,family=negative.binomial(theta),correlation=corExp(form=~longitude+latitude,nugget = TRUE))
          }
          #cor_r<- round(r.squaredGLMM(glm_re_spatial_freq)[2,2],2)
        
          if(div_var=="plant_H_400" | div_var=="tree_H_400"){
            glm_re_spatial_freq<- glmmPQL(div ~ strct,random=~1|siteID,data=tmp,family=gaussian,correlation=corExp(form=~longitude+latitude,nugget = TRUE))
            cor_r<- round(r.squaredGLMM(glm_re_spatial_freq)[1,2],2)
          }
        
        # org results  
          results<-as.data.frame(summary(glm_re_spatial_freq)$tTable)
          results$lowC<-results[,1]-results[,2]*1.96
          results$hiC<-results[,1]+results[,2]*1.96
          results<- results[2,c(1,6,7)]
          
          if(all(results<0)){cor_r<-cor_r} else if (all(results>0)) {cor_r<-cor_r} else {cor_r<-0}
          pframe <- with(tmp,expand.grid(strct=seq(min(strct),max(strct),length=nrow(tmp))))
          pframe$div <- predict(glm_re_spatial_freq, pframe, level = 0, type="response")
          pframe$siteID<-tmp$siteID
        
      # create raw charts
        
        # plot elements   
          axis_text_size<-12
          axis_title_size<-15
          plot_point<-.5
          site_point<-1
          r_box_size<-5
          reg_line_size<-1.7
          line_cols<- mycols<-c(brewer.pal(n = 11, name = "PiYG")[10],brewer.pal(n = 11, name = "BrBG")[2])
          lr_margins<-0.1
          ud_margins<-0.1
          mytit_h<- -.5
          mytit_v<- .7
          
        tmp_plot <- ggplot(tmp,aes(x = strct, y = div,colour = siteID)) +
          geom_point(size=plot_point) +    
          scale_colour_manual(values=bubble_colors) +
          annotate("label",x=min(tmp$strct), y = max(tmp$div), colour="black",fill = "white",family="Times",vjust=1,hjust=0,label = paste0("R^2 == ",cor_r),parse=T,size=r_box_size) +  
          my_theme(mycolor_background="grey95",mycolor_elements="grey60",mycolor_text="black") +
          theme(legend.position = "none",
                axis.text.x = element_text(angle=25,size=axis_text_size,family="Helvetica",hjust=.7,vjust=.7),
                axis.text.y = element_text(size=axis_title_size,family="Helvetica"),
                axis.title.x = element_text(size=axis_title_size,vjust=0,family="Helvetica"),
                plot.margin = unit(c(ud_margins,lr_margins,ud_margins,lr_margins), "cm"),
                axis.title.y = element_text(size=axis_title_size,family="Helvetica")) + #hjust=0,margin = ggplot2::margin(b = -20)),
          ylab("test") + xlab(strct_var) 

        if (cor_r==0){
          plot_list[[plot_num]]<-tmp_plot
        } else {
          plot_list[[plot_num]] <-tmp_plot + geom_line(data=pframe,aes(x=strct, y=div),col="black",size=2)
        }
        plot_num<-plot_num+1
      }
    }

  # 3b final charts ----------- 

    quartz()
    plot_list
   

# Figure 4. Climate’s effects on continental patterns of diversity and forest structure conditions across the United States ###########  

  # load data
    library(MASS)
    full_data<-read.csv("... xx_etal_data.csv")

  # Fig 4a analysis -------
    div_vars<-c("plant_SR_400","plant_H_400","tree_SR_400","tree_H_400")  
    strct_vars<-c("CHM_nni_3_2512m2","pc_skew_DBH_2512m2","z_ent1_ent2_2512m2","z_1m_mode_mean_2512m2","pc_ent1_01_2512m2","z_cc_1_3_2512m2")
    clim_vars<-c("MAT","AET","MTCM","alpha","TSeas","PSeas")
    beta_mat<-as.data.frame(matrix(NA,length(clim_vars)*length(div_vars),6))
    colnames(beta_mat)<-c("clim_var","div_var","beta_025","beta_5","beta_975","sig_test")
    beta_mat$clim_var<-rep(clim_vars,length(div_vars))
    beta_mat$div_var<-rep(div_vars,each=length(clim_vars))
    div_var<-"tree_H_400";clim_var<-"MTCM";sd<-1;index<-1;subset<-9 #@@

    for (div_var in div_vars){
      for (clim_var in clim_vars){
        
        tmp<-full_data[,c("plotID","siteID",div_var,clim_var,"longitude","latitude")]
        names(tmp)<-c("plotID","siteID","div","clim","longitude","latitude")
        tmp2<-tmp[complete.cases(tmp),]
        
        # glm + RE + spatial
          tmp2<-tmp[complete.cases(tmp),]
          tmp3<-tmp2[tmp2$siteID%in%names(table(tmp2$siteID)[table(tmp2$siteID)>=9]),]
          tmp3[,colnames(tmp3)%in%c("clim")]<-scale(tmp3[,colnames(tmp3)%in%c("clim")])
          
          if(div_var=="tree_SR_400"){
            glm_re_spatial_freq<- glmmPQL(div ~ clim,random=~1|siteID,data=tmp3,family=poisson,correlation=corExp(form=~longitude+latitude,nugget = TRUE))
            
          }
          if(div_var=="plant_SR_400"){
            m0 <- glm.nb(div ~ clim + siteID , data=tmp3)
            theta<-m0$theta
            glm_re_spatial_freq<- glmmPQL(div ~ clim,random=~1|siteID,data=tmp3,family=negative.binomial(theta),correlation=corExp(form=~longitude+latitude,nugget = TRUE))
          }
          if(div_var=="tree_H_400" | div_var=="plant_H_400"){
            tmp3[,colnames(tmp3)%in%c("div","clim")]<-scale(tmp3[,colnames(tmp3)%in%c("div","clim")])
            glm_re_spatial_freq<- glmmPQL(div ~ clim,random=~1|siteID,data=tmp3,family=gaussian,correlation=corExp(form=~longitude+latitude,nugget = TRUE))
          }
          results<-as.data.frame(summary(glm_re_spatial_freq)$tTable)
          results$lowC<-results[,1]-results[,2]*1.96
          results$hiC<-results[,1]+results[,2]*1.96
          results<- results[2,c(1,6,7)]
          beta_mat[beta_mat$clim_var==clim_var & beta_mat$div_var==div_var,3:5]<-round(results[c(2,1,3)],2)
        
          if (results[2]<0 & results[3]>0){
            beta_mat[beta_mat$clim_var==clim_var & beta_mat$div_var==div_var,6]<-"non"
          } else {beta_mat[beta_mat$clim_var==clim_var & beta_mat$div_var==div_var,6]<-div_var}
          
      }
    }
  
    beta_mat$clim_var<-as.factor(beta_mat$clim_var)
    beta_mat$div_var<-as.factor(beta_mat$div_var)
    
    beta_mat$div_var<-revalue(beta_mat$div_var, c("plant_SR_400"="SR[plant]","tree_SR_400"="SR[tree]","plant_H_400"="H[plant]","tree_H_400"="H[tree]"))
    beta_mat$div_var<-factor(beta_mat$div_var, levels = c("SR[plant]","H[plant]","SR[tree]","H[tree]"))
    
    beta_mat$clim_var<-revalue(beta_mat$clim_var, c("MAT"="T[mean]","AET"="AET","alpha"="aridity","MTCM"="coldness","TSeas"="T[Seas]","PSeas"="P[Seas]"))
    beta_mat$clim_var<-factor(beta_mat$clim_var, levels = rev(c("coldness","aridity","T[mean]","AET", "T[Seas]", "P[Seas]")))
    
    beta_mat$sig_test<-revalue(beta_mat$sig_test, c("plant_SR_400"="SR[plant]","tree_SR_400"="SR[tree]","plant_H_400"="H[plant]","tree_H_400"="H[tree]"))
    beta_mat$sig_test<-factor(beta_mat$sig_test, levels = c("non","SR[plant]","H[plant]","SR[tree]","H[tree]"))


  # Fig 4a chart -------

    # chart elements
      mycols<-c(brewer.pal(n = 9, name = "Greens")[c(9,5)],brewer.pal(n = 9, name = "BrBG")[c(1,3)])
      size_tmp<-10
      axis.text<-9
      myfills<-c("grey98",mycols)
      annotate_pos=min(beta_mat$beta_025)*1.35
      label_cols<-c(brewer.pal(n = 9, name = "Reds")[8],brewer.pal(n = 9, name = "Greens")[8],brewer.pal(n = 9, name = "Blues")[8])

    fig4a_plot <-  ggplot(beta_mat, aes(x = as.factor(clim_var),y=beta_5,colour=factor(div_var),shape=factor(div_var),fill=factor(sig_test))) +
      geom_hline(yintercept = 0, linetype = "longdash")  +  
      geom_errorbar(aes(ymin=beta_025, ymax=beta_975), size=.6,width=0,position=position_dodge(width=-0.6)) +
      geom_point(position=position_dodge(width=-0.6),size=1.5) +
      guides(fill=FALSE,shape=FALSE,colour = guide_legend(override.aes = list(shape =21:24,fill=mycols))) +    
      scale_shape_manual(values = c(21:24)) +
      scale_colour_manual(values=mycols,labels = parse(text = levels(beta_mat$div_var))) +
      scale_fill_manual(values = c(myfills)) +
      my_theme(mycolor_background="grey98",mycolor_elements="grey60",mycolor_text="black",plot_background="white") +
      theme(axis.text.y=element_text(angle=0,family="Helvetica",size=size_tmp,hjust=1,vjust=0.5),
            axis.text.x=element_text(angle=15,family="Helvetica",size=size_tmp-1),
            axis.title.x=element_text(angle=0,family="Helvetica",size=size_tmp+2),
            axis.title.y=element_text(family="Helvetica"),
            plot.title = element_text(size=size_tmp,vjust=2),
            plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
            legend.title=element_blank(),
            legend.spacing.x = unit(0.1, 'cm'),
            legend.spacing.y = unit(0, 'cm'),
            legend.position=c(1,1),
            legend.background = element_rect(colour = 'transparent', fill = 'transparent'),
            legend.box.background = element_rect(colour = "grey60"),
            legend.key = element_rect(colour = 'transparent', fill = 'transparent'),
            legend.key.size = unit(0.5, "cm"),
            legend.key.height = unit(0, "cm"),
            legend.text.align = 0,
            legend.justification=c(1,1),
            legend.text=element_text(size=axis.text-2,family="Helvetica")) +
      xlab("") +  ylab(expression(beta[climate])) +
      geom_vline(xintercept=c(2.5,4.5), color = "grey80",size=.3) +
      annotate(geom = "text", x = c(1.5,3.5,5.5), y =annotate_pos-0.06, label = c("seasonality","energy","stress"),angle=90,hjust=.5,vjust=0.5, size = (size_tmp-6),parse=T,family="Helvetica") +  
      annotate("segment", x = 0.7, xend = 2.3, y = annotate_pos, yend = annotate_pos,color="black") +
      annotate("segment", x = 2.7, xend = 4.3, y = annotate_pos, yend = annotate_pos,color="black") +
      annotate("segment", x = 4.7, xend = 6.3, y = annotate_pos, yend = annotate_pos,color="black") +
      coord_flip(xlim=c(0.5,6.5),ylim =  c(min(beta_mat$beta_025)-0.001,max(beta_mat$beta_975)+0.001), expand = FALSE, clip = "off") + 
      scale_x_discrete("", labels = parse(text = levels(beta_mat$clim_var))) +
      annotate(geom = "text", x = 6.4, y =min(beta_mat$beta_025)*.95, label = "(a)",hjust=0,vjust=1,size = (size_tmp-6),family="Helvetica",fontface="bold") 

  # Fig 4b analysis --------
  
    beta_mat4b<-as.data.frame(matrix(NA,length(clim_vars)*length(strct_vars),6))
    colnames(beta_mat4b)<-c("clim_var","strct_var","beta_025","beta_5","beta_975","sig_test")
    beta_mat4b$clim_var<-rep(clim_vars,length(strct_vars))
    beta_mat4b$strct_var<-rep(strct_vars,each=length(clim_vars))
    
    strct_var<-"z_cc_3_2512m2";clim_var<-"alpha";sd<-2;index<-1
    for (strct_var in strct_vars){
      for (clim_var in clim_vars){
    
        tmp<-full_data[,c("plotID","siteID",strct_var,clim_var,"longitude","latitude")]
        names(tmp)<-c("plotID","siteID","strct","clim","longitude","latitude")
        tmp2<-tmp[complete.cases(tmp),]
    
        tmp2[,colnames(tmp2)%in%c("strct","clim")]<-scale(tmp2[,colnames(tmp2)%in%c("strct","clim")])
        tmp3<-tmp2[tmp2$siteID%in%names(table(tmp2$siteID)[table(tmp2$siteID)>=9]),]
    
        skip_to_next <- FALSE
        tryCatch(re_model_freq_SAC <- lme(strct ~ clim,random=~1|siteID,data=tmp3,correlation=corExp(form=~longitude+latitude,nugget = TRUE)), error = function(e) {skip_to_next <<- TRUE})
        if(skip_to_next) { next }
        results<-as.data.frame(summary(re_model_freq_SAC)$tTable)
        results$lowC<-results[,1]-results[,2]*1.96
        results$hiC<-results[,1]+results[,2]*1.96
        results<- results[2,c(1,6,7)]
        print(clim_var)
        beta_mat4b[beta_mat4b$strct_var==strct_var & beta_mat4b$clim_var==clim_var,3:5]<-round(results[c(2,1,3)],3)
        if (results[2]<=0 & results[3]>=0){
          beta_mat4b[beta_mat4b$clim_var==clim_var & beta_mat4b$strct_var==strct_var,6]<-"non"
        } else {beta_mat4b[beta_mat4b$clim_var==clim_var & beta_mat4b$strct_var==strct_var,6]<-strct_var}
      }
    }
  
    beta_mat4b$clim_var<-as.factor(beta_mat4b$clim_var)
    beta_mat4b$strct_var<-as.factor(beta_mat4b$strct_var)
    
    beta_mat4b$clim_var<-revalue(beta_mat4b$clim_var, c("MAT"="T[mean]","AET"="AET","alpha"="aridity","MTCM"="coldness","TSeas"="T[Seas]","PSeas"="P[Seas]"))
    beta_mat4b$clim_var<-factor(beta_mat4b$clim_var, levels = rev(c("coldness","aridity","T[mean]","AET", "T[Seas]", "P[Seas]")))
    
    beta_mat4b$strct_var<-revalue(beta_mat4b$strct_var, c("CHM_nni_3_2512m2"="C[clumping]","pc_skew_DBH_2512m2"="H[skew]","z_ent1_ent2_2512m2"="C[hetero]",
                                                          "z_1m_mode_mean_2512m2"="H[mean]","pc_ent1_01_2512m2"="H[hetero]","z_cc_1_3_2512m2"="C[cover]"))
    beta_mat4b$strct_var<-factor(beta_mat4b$strct_var, levels = c("C[cover]","H[mean]","C[hetero]","H[hetero]","C[clumping]","H[skew]"))
    
    beta_mat4b$sig_test<-revalue(beta_mat4b$sig_test, c("CHM_nni_3_2512m2"="C[clumping]","pc_skew_DBH_2512m2"="H[skew]","z_ent1_ent2_2512m2"="C[hetero]",
                                                        "z_1m_mode_mean_2512m2"="H[mean]","pc_ent1_01_2512m2"="H[hetero]","z_cc_1_3_2512m2"="C[cover]"))
    beta_mat4b$sig_test<-factor(beta_mat4b$sig_test, levels = c("non","C[cover]","H[mean]","C[hetero]","H[hetero]","C[clumping]","H[skew]"))

  # Fig 4b chart -------

    # chart elements
      mycols<-viridis_pal()(11)[c(1,3,5,7,9,11)]
      size_tmp<-10
      myfills<-c("grey98",mycols)
      annotate_pos=min(beta_mat4b$beta_025)*1.4

    fig4b_plot <- ggplot(beta_mat4b, aes(x = as.factor(clim_var),y=beta_5,colour=factor(strct_var),shape=factor(strct_var),fill=factor(sig_test))) +
      geom_hline(yintercept = 0, linetype = "longdash")  +  
      geom_errorbar(aes(ymin=beta_025, ymax=beta_975), width=0,position=position_dodge(width=-0.7),size=.6) +
      geom_point(position=position_dodge(width=-0.7),size=1.2) +
      guides(fill=FALSE,shape=FALSE,colour = guide_legend(override.aes = list(shape =c(21:25,21),colour=mycols,fill=mycols))) +    
      scale_shape_manual(values = c(21:25,21)) +
      scale_colour_manual(values=mycols,labels = parse(text = levels(beta_mat4b$strct_var))) +
      scale_fill_manual(values = myfills) +
      my_theme(mycolor_background="grey98",mycolor_elements="grey60",mycolor_text="black",plot_background="white") +
      theme(axis.text.y=element_blank(),
            axis.text.x=element_text(angle=15,family="Helvetica",size=size_tmp-1),
            axis.title.x=element_text(angle=0,family="Helvetica",size=size_tmp+2),
            axis.title.y=element_text(family="Helvetica"),
            plot.title = element_text(size=size_tmp+1,vjust=2),
            plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
            legend.title=element_blank(),
            legend.spacing.x = unit(0, 'cm'),
            legend.position=c(1,1),
            legend.justification=c(1,1),
            legend.background = element_blank(),
            legend.box.background = element_rect(colour = "grey60"),
            legend.spacing.y = unit(-.1, "mm"),
            legend.key = element_rect(colour = 'transparent', fill = 'transparent'),
            legend.key.size = unit(0.5, "cm"),
            legend.key.height = unit(0, "cm"),
            legend.text.align = 0,
            legend.text=element_text(size=size_tmp-2,family="Helvetica")) +
      xlab("") +  ylab(expression(beta[climate])) +
      geom_vline(xintercept=c(2.5,4.5), color = "grey80",size=.3) +
      coord_flip(xlim=c(0.5,6.5),ylim =  c(min(beta_mat4b$beta_025)-0.01,max(beta_mat4b$beta_975)+0.01), expand = FALSE, clip = "off") + 
      scale_x_discrete("", labels = parse(text = levels(beta_mat4b$clim_var)))+
      annotate(geom = "text", x = 6.4, y =min(beta_mat4b$beta_025)*.95, label = "(b)",hjust=0,vjust=1,size = (size_tmp-6),family="Helvetica",fontface="bold") 

  # Fig 4 final subplots 
    quartz()
    fig4a_plot
    fig4b_plot


# Figure 5. Individual and joint effects of climate and forest structure variables on plant and tree species richness ###########  
    
  # load data
    library(MASS)
    full_data<-read.csv("... xx_etal_data.csv")  # @@
    panel_vars<-sort(c("tree_SR_400","plant_SR_400"));panel_label<-"diversity"
    clim_vars<-c("MAT","AET","MTCM","alpha","TSeas","PSeas");y_label<-"climate"
    treatment_vars<-c("CHM_nni_3_2512m2","pc_skew_DBH_2512m2","z_ent1_ent2_2512m2","z_1m_mode_mean_2512m2","pc_ent1_01_2512m2","z_cc_1_3_2512m2")
    treatment_label<-"structure"           

  # Fig 5 analysis ------

    # model output matrix
      interactions_mat<-as.data.frame(matrix(NA,(length(clim_vars)*length(treatment_vars)),2+length(panel_vars)*4))
      colnames(interactions_mat)<-c("clim_var","treatment_var",sort(paste(panel_vars,rep(c("_025","_5","_975","_sig"),each=length(panel_vars)),sep="")))
      interactions_mat$treatment_var<-rep(treatment_vars,length(clim_vars))
      interactions_mat$clim_var<-rep(clim_vars,each=length(treatment_vars))
      pred1_mat<-pred2_mat<-interactions_mat

    # run all models  
      treatment_var<-"z_ent1_ent2_2512m2";panel_var<-"plant_SR_400";clim_var<-"PSeas";sd<-2;site_index<-1; subset<-9;siteID<-"epa_merge";model_type<-3 #@@
      progress.bar <- create_progress_bar("text")
      progress.bar$init(length(treatment_vars)*length(panel_vars)*length(clim_vars))
  
      for (treatment_var in treatment_vars){
        for (panel_var in panel_vars){
          for (clim_var in clim_vars){
            
            tmp<-full_data[,c("siteID","plotID",treatment_var,panel_var,clim_var,"longitude","latitude")]
            names(tmp)[grep("bird|mam|plant|tree", names(tmp))]<-"tmp_div"
            
            if ((names(tmp)[-c(1:2)][names(tmp)[-c(1:2)]!="tmp_div"]%in%treatment_vars)[1]==T){pred1<-"treatment"}
            if ((names(tmp)[-c(1:2)][names(tmp)[-c(1:2)]!="tmp_div"]%in%treatment_vars)[2]==T){pred2<-"treatment"}
            if ((names(tmp)[-c(1:2)][names(tmp)[-c(1:2)]!="tmp_div"]%in%clim_vars)[1]==T){pred1<-"y"}
            if ((names(tmp)[-c(1:2)][names(tmp)[-c(1:2)]!="tmp_div"]%in%clim_vars)[2]==T){pred2<-"y"}
            if ((names(tmp)[-c(1:2)][names(tmp)[-c(1:2)]!="tmp_div"]%in%panel_vars)[1]==T){pred1<-"panel"}
            if ((names(tmp)[-c(1:2)][names(tmp)[-c(1:2)]!="tmp_div"]%in%panel_vars)[2]==T){pred2<-"panel"}
            
            names(tmp)<-c("siteID","plotID","pred1","tmp_div","pred2","longitude","latitude")
            tmp2<-tmp[complete.cases(tmp),]
            tmp2[,colnames(tmp2)%in%c("pred1","pred2")]<-scale(tmp2[,colnames(tmp2)%in%c("pred1","pred2")])
            tmp3<-tmp2[tmp2$siteID%in%names(table(tmp2$siteID)[table(tmp2$siteID)>=9]),]
            
            if(panel_var=="tree_SR_400"){
              glm_re_spatial_freq<- glmmPQL(tmp_div ~ pred1*pred2,random=~1|siteID,data=tmp3,family=poisson,correlation=corExp(form=~longitude+latitude,nugget = TRUE))
            }
            if(panel_var=="plant_SR_400"){
              m0 <- glm.nb(tmp_div ~ pred1*pred2 + siteID , data=tmp3)
              theta<-m0$theta
              glm_re_spatial_freq<- glmmPQL(tmp_div ~ pred1*pred2,random=~1|siteID,data=tmp3,family=negative.binomial(theta),correlation=corExp(form=~longitude+latitude,nugget = TRUE))
            }
            
            results<-as.data.frame(summary(glm_re_spatial_freq)$tTable)
            results$lowC<-results[,1]-results[,2]*1.96
            results$hiC<-results[,1]+results[,2]*1.96
            results<- results[2:4,c(1,6,7)]
      
            pred1_mat[pred1_mat$treatment_var==treatment_var & pred1_mat$clim_var==clim_var,colnames(pred1_mat) %ilike% panel_var & ! colnames(pred1_mat) %ilike% "sig"]<-round(c(results[1,2],results[1,1],results[1,3]),5)
            pred2_mat[pred2_mat$treatment_var==treatment_var & pred2_mat$clim_var==clim_var,colnames(pred2_mat) %ilike% panel_var & ! colnames(pred2_mat) %ilike% "sig"]<-round(c(results[2,2],results[2,1],results[2,3]) ,5)
            interactions_mat[interactions_mat$treatment_var==treatment_var & interactions_mat$clim_var==clim_var,colnames(interactions_mat) %ilike% panel_var & ! colnames(interactions_mat) %ilike% "sig"]<-round(c(results[3,2],results[3,1],results[3,3]) ,5)
            
            if (results[1,2]<0 & results[1,3]>0){
              pred1_mat[pred1_mat$treatment_var==treatment_var & pred1_mat$clim_var==clim_var,colnames(pred1_mat) %ilike% panel_var & colnames(pred1_mat) %ilike% "sig"]<-"non"
            } else {pred1_mat[pred1_mat$treatment_var==treatment_var & pred1_mat$clim_var==clim_var,colnames(pred1_mat) %ilike% panel_var & colnames(pred1_mat) %ilike% "sig"]<-clim_var}
            
            if (results[2,2]<0 & results[2,3]>0){
              pred2_mat[pred1_mat$treatment_var==treatment_var & pred2_mat$clim_var==clim_var,colnames(pred2_mat) %ilike% panel_var & colnames(pred2_mat) %ilike% "sig"]<-"non"
            } else {pred2_mat[pred1_mat$treatment_var==treatment_var & pred2_mat$clim_var==clim_var,colnames(pred2_mat) %ilike% panel_var & colnames(pred2_mat) %ilike% "sig"]<-clim_var}
            
            if (results[3,2]<0 & results[3,3]>0){
              interactions_mat[pred1_mat$treatment_var==treatment_var & interactions_mat$clim_var==clim_var,colnames(interactions_mat) %ilike% panel_var & colnames(interactions_mat) %ilike% "sig"]<-"non"
            } else {interactions_mat[interactions_mat$treatment_var==treatment_var & interactions_mat$clim_var==clim_var,colnames(interactions_mat) %ilike% panel_var & colnames(interactions_mat) %ilike% "sig"]<-clim_var}
            site_index<-site_index+1
            progress.bar$step()
          }
        }
      }
      
      pred1_mat[names(pred1_mat)[which(names(pred1_mat)%ilike%"sig")]] <- lapply(pred1_mat[names(pred1_mat)[which(names(pred1_mat)%ilike%"sig")]], function(x){factor(x,levels = c("non",clim_vars))})
      pred2_mat[names(pred2_mat)[which(names(pred2_mat)%ilike%"sig")]] <- lapply(pred2_mat[names(pred2_mat)[which(names(pred2_mat)%ilike%"sig")]], function(x){factor(x,levels = c("non",clim_vars ))})
      interactions_mat[names(interactions_mat)[which(names(interactions_mat)%ilike%"sig")]] <- lapply(interactions_mat[names(interactions_mat)[which(names(interactions_mat)%ilike%"sig")]], function(x){factor(x,levels = c("non",clim_vars))})
      
      if(pred1=="treatment"){pred1_label<-treatment_label};if(pred2=="treatment"){pred2_label<-treatment_label}
      if(pred1=="y"){pred1_label<-y_label};if(pred2=="y"){pred2_label<-y_label}
      if(pred1=="panel"){pred1_label<-panel_label};if(pred2=="panel"){pred2_label<-panel_label}
      
      mat_list<-list(pred1_mat,pred2_mat,interactions_mat)
      names(mat_list)<-c(pred1_label,pred2_label,"interaction")
      
  # Fig 5 org results ------

    mat_list[[1]]<-mat_list[[1]][which(mat_list[[1]]$treatment_var%in%treatment_vars),]
    mat_list[[2]]<-mat_list[[2]][which(mat_list[[2]]$treatment_var%in%treatment_vars),]
    mat_list[[3]]<-mat_list[[3]][which(mat_list[[3]]$treatment_var%in%treatment_vars),]
    
    i<-1
    for (i in 1:length(mat_list)){
      mat_list[[i]]$treatment_var<-revalue( mat_list[[i]]$treatment_var, c("CHM_nni_3_2512m2"="C[clumping]","pc_skew_DBH_2512m2"="H[skew]",
                                                                           "z_ent1_ent2_2512m2"="C[hetero]",
                                                                           "z_1m_mode_mean_2512m2"="H[mean]","pc_ent1_01_2512m2"="H[hetero]","z_cc_1_3_2512m2"="C[cover]"))
      
      mat_list[[i]]$treatment_var<-factor( mat_list[[i]]$treatment_var, levels = rev(c( "C[cover]","H[mean]","C[hetero]","H[hetero]","C[clumping]","H[skew]")))
      mat_list[[i]]$clim_var<-revalue( mat_list[[i]]$clim_var, c("MAT"="T[mean]","AET"="AET","alpha"="aridity","MTCM"="coldness","PSeas"="P[seas]","TSeas"="T[seas]"))
      mat_list[[i]]$clim_var<-factor( mat_list[[i]]$clim_var, levels = c("coldness","aridity","T[mean]","AET","T[seas]","P[seas]"))
      
      mat_list[[i]][,6]<-factor(mat_list[[i]][,6],levels=c("non",rev(c("MTCM","alpha","MAT","AET","TSeas","PSeas"))))
      
      mat_list[[i]][,10]<-factor(mat_list[[i]][,10],levels=c("non",rev(c("MTCM","alpha","MAT","AET","TSeas","PSeas"))))
    }

    # Fig 5 charts ------
    
      # chart elements
        bars<-length(treatment_vars)
        mycols<-c(brewer.pal(n = 9, name = "Reds")[c(5,8)],brewer.pal(n = 9, name = "Greens")[c(5,8)],brewer.pal(n = 9, name = "Blues")[c(5,8)])
        myfills<-c("grey98",rev(mycols))
        if (bars<6){my_shapes<-21:(20+bars)}
        if (bars>5){my_shapes<-c(21:25,21:(bars+15))}
        plot_list<-list()
        large_plots<-list();index<-1
        panel_titles<-c("tree species richness","plant species richness")
        fig<-8
        axis_text<-9
        annotate_label_size<-axis_text-6
        erroBar<-.5
        point_szie<-1.2
        axis_title_size<-axis_text+1
        legend_text_size<-axis_text-2
        title_size<-axis_text+1
        panel_num<-2;term<-3;index<-1
      
      # run charts
        
        for (panel_num in 1:length(panel_vars)){
          if (panel_num==1){offset<-.085} else {offset<-.22}
          for (term in 1:3){
            
            myxmin<-min(mat_list[[term]][(panel_num*4)-1])
            myxmax<-max(mat_list[[term]][(panel_num*4)+1])
            label_offset<-(myxmax-myxmin)*0.02
            
            if (panel_num==1 & term==1){label_x<-"(a)"}
            if (panel_num==1 & term==2){label_x<-"(b)"}
            if (panel_num==1 & term==3){label_x<-"(c)"}
            if (panel_num==2 & term==1){label_x<-"(d)"}
            if (panel_num==2 & term==2){label_x<-"(e)"}
            if (panel_num==2 & term==3){label_x<-"(f)"}
            
            plot_list[[index]] <-  ggplot(mat_list[[term]], aes_(x = as.name(colnames(mat_list[[term]])[panel_num*4]),y=~as.factor(treatment_var),
                                                                 colour=~factor(clim_var),shape=~factor(clim_var),fill=as.name(colnames(mat_list[[term]])[(panel_num*4)+2]))) +
              geom_vline(xintercept = 0, linetype = "longdash",colour="grey80") +
              geom_errorbar(aes_(xmin=as.name(colnames(mat_list[[term]])[(panel_num*4)-1]), 
                                 xmax=as.name(colnames(mat_list[[term]])[(panel_num*4)+1])), width=0,position=position_dodge(width=-0.7),size=erroBar)  +
              geom_point(position=position_dodge(width=-0.7),size=point_szie) +
              guides(fill=FALSE,shape=FALSE,colour = guide_legend(override.aes = list(shape = my_shapes,fill=mycols))) +    
              scale_shape_manual(values = c(my_shapes)) +
              scale_colour_manual(values=mycols,labels = parse(text = levels(mat_list[[term]]$clim_var)))  +
              scale_fill_manual(values = c(myfills[c(which(c("non",rev(clim_vars))%in%unique(mat_list[[term]][,(panel_num*4)+2])))])) +
              my_theme(mycolor_background="grey98",mycolor_elements="grey60",mycolor_text="black",plot_background="white") +
              geom_hline(yintercept = c(2.5,4.5), colour="grey80") +
              theme(legend.position = "none",
                    axis.text.x = element_text(hjust = .5,angle=15,family="Helvetica",size=axis_text-2),
                    axis.text.y = element_text(family="Helvetica"),
                    axis.title.x=element_text(family="Helvetica"),
                    axis.title.y=element_text(family="Helvetica")) + xlab("") + ylab("") +
              coord_cartesian(xlim= c(myxmin-label_offset,myxmax+label_offset),ylim =c(.5,(length(clim_vars)+.5)) , expand = FALSE, clip = "off") + ggtitle("") +
              scale_y_discrete("", labels = parse(text = levels(mat_list[[term]]$treatment_var))) +
              annotate(geom = "text", x = myxmin, y =6.4, label = label_x,hjust=0,vjust=1,size = annotate_label_size,family="Helvetica",fontface="bold") 
            
            if (term==1){
              plot_list[[index]]<- plot_list[[index]] +
                annotate(geom = "segment", x =myxmin-offset, xend = myxmin-offset , y = .6, yend = 2.4,colour="black") +
                annotate(geom = "segment", x = myxmin-offset , xend = myxmin-offset, y = 2.6, yend = 4.4,colour="black") +
                annotate(geom = "segment", x = myxmin-offset, xend = myxmin-offset, y = 4.6, yend = 6.4,colour="black") +
                annotate(geom = "text", x = myxmin-offset-offset*0.4, y =c(1.5,3.5,5.5),
                         label = c("spatial","structural","canopy"),angle=90,hjust=0.5,vjust=0.5,parse=T,
                         size = annotate_label_size,colour="black",family="Helvetica") +
                annotate(geom = "text", x = myxmin-offset-offset*0.15, y =c(1.5,3.5,5.5),
                         label = c("configuration","heterogeneity","dimensions"),angle=90,hjust=0.5,vjust=0.5,parse=T,
                         size = annotate_label_size,colour="black",family="Helvetica")
            }
            index<-index+1
          }
        }

  # final charts       
    quartz()
    plot_list
    