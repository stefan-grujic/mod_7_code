
# Produce tSNE plots for pre and post scoary COGs

rootloc <- readline(prompt="Enter Project Directory: ")
outdir <- readline(prompt="Enter Output Directory: ")

# load packages: ----------------------------------------------------------------------------
library(pacman)
p_load(corrplot)
p_load(tidyverse)
p_load(Hmisc)
p_load(stats)
p_load(plotly)
p_load(htmlwidgets)
p_load(Rtsne)
library(data.table)

# read file: --------------------------------------------------------------------------------
cogs <- read.csv(paste(rootloc,'/COG/expression_table/cog_table.csv',sep=''))
meta <- read.csv(paste(rootloc,'/COG/expression_table/cog_multinomial_metadata.csv',sep=''))

scoary_ae <- read.csv(paste(rootloc,'/SCOARY/COG/aerobic_cogs_scoary.csv',sep=''))
scoary_an <- read.csv(paste(rootloc,'/SCOARY/COG/anaerobic_cogs_scoary.csv',sep=''))
scoary_ep <- read.csv(paste(rootloc,'/SCOARY/COG/epibiotic_cogs_scoary.csv',sep=''))
scoary_fa <- read.csv(paste(rootloc,'/SCOARY/COG/facultative_cogs_scoary.csv',sep=''))
scoary_gr <- read.csv(paste(rootloc,'/SCOARY/COG/group attack_cogs_scoary.csv',sep=''))
scoary_np <- read.csv(paste(rootloc,'/SCOARY/COG/non predator_cogs_scoary.csv',sep=''))
scoary_pr <- read.csv(paste(rootloc,'/SCOARY/COG/predator_cogs_scoary.csv',sep=''))


# convert to df: ----------------------------------------------------------------------------
cogs_df <- as.data.frame(cogs)
meta_df <- as.data.frame(meta)

cogs_ae <- select(cogs_df, name, colnames(scoary_ae))
cogs_an <- select(cogs_df, name, colnames(scoary_an))
cogs_ep <- select(cogs_df, name, colnames(scoary_ep))
cogs_fa <- select(cogs_df, name, colnames(scoary_fa))
cogs_gr <- select(cogs_df, name, colnames(scoary_gr))
cogs_np <- select(cogs_df, name, colnames(scoary_np))
cogs_pr <- select(cogs_df, name, colnames(scoary_pr))

# tSNE: ---------------------------------------------------------------------------------------
do_the_tsne = function(cog_data, saveloc, tol, title){
  
  labels <- cog_data$name
  cog_data$name <- as.factor(cog_data$name)
  
  tsne = Rtsne(cog_data[,-1], dims = 2, perplexity=30, verbose=TRUE, max_iter = 5000, check_duplicates = F)
  tsne_scores = data.frame(meta_df, tsne$Y)
  
  if (tol == 'oxygen_tolerance'){
    col_list = c('lightseagreen', 'purple', 'royalblue')
    legend_list = c('Anaerobic', 'Facultative', 'Aerobic')
    tol_id = tsne_scores$oxygen_tolerance
  }
  else{
    col_list = c('gray44', 'green', 'deepskyblue','red')
    legend_list = c('Non-Predator', 'Group-Attack', 'Epibiotic', 'Endobiotic')
    tol_id = tsne_scores$predation_strategy
  }
    
  ggplot(tsne_scores, aes(x=X1, y=X2,
                         col=as.factor(tol_id),
                         size=2))+
    geom_point()+
    ggtitle(title)+
    theme_bw()+ theme(plot.title = element_text(hjust = 0.5))+
    scale_color_manual(labels = legend_list, values = col_list)+
    guides(size=F, col=guide_legend(title=NULL), shape=guide_legend(title=NULL))
    #xlab(paste("PC1:", toString(format(round(prop_variance[1], 2), nsmall = 2)), sep=" "))+
    #ylab(paste("PC2:", toString(format(round(prop_variance[2], 2), nsmall = 2)), sep=" "))+
  
  if (saveloc != 0){
    ggsave(saveloc)
    
  }
}

do_the_tsne(cogs_ae,
            paste(outdir,"/aerobic.png",sep=""),
            'oxygen_tolerance',
            'Scoary Aerobic')

do_the_tsne(cogs_an,
            paste(outdir,"/anaerobic.png",sep=""),
            'oxygen_tolerance',
            'Scoary Anaerobic')

do_the_tsne(cogs_fa,
            paste(outdir,"/facultative.png",sep=""),
            'oxygen_tolerance',
            'Scoary Facultative')

do_the_tsne(cogs_gr,
            paste(outdir,"/group_attack.png",sep=""),
            'predation_strategy',
            'Scoary Group-Attack')

do_the_tsne(cogs_ep,
            paste(outdir,"/epibiotic.png",sep=""),
            'predation_strategy',
            'Scoary Epibiotic')

do_the_tsne(cogs_pr,
            paste(outdir,"/predator.png",sep=""),
            'predation_strategy',
            'Scoary Predator')

do_the_tsne(cogs_np,
            paste(outdir,"/non_predator.png",sep=""),
            'predation_strategy',
            'Scoary Non-Predator')

do_the_tsne(cogs_df,
            paste(outdir,"/pre_scoary_oxy.png",sep=""),
            'oxygen_tolerance',
            'Pre-Scoary Oxygen Tolerance')

do_the_tsne(cogs_df,
            paste(outdir,"/pre_scoary_pre.png",sep=""),
            'predation_strategy',
            'Pre-Scoary Predation Strategy')

