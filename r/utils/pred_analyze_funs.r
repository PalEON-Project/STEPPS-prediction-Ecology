bw <- readOGR('data/map/bw/bw_albers.shp')
bw_fort <- fortify(bw)
bw_buff = 108000
bw_lims <- list(xlims=c(min(bw_fort$long)-bw_buff, max(bw_fort$long)+bw_buff)/1e6, 
                ylims=c(min(bw_fort$lat)-bw_buff, max(bw_fort$lat)+64000)/1e6)


get_r_quants <- function(rIts){
  
  N = dim(rIts)[1]
  K = dim(rIts)[2]
  T = dim(rIts)[3]
  
  r_int = apply(rIts, c(1,2,3), function(x) quantile(x, probs=c(0.1, 0.5, 0.9)))
  r_sd = apply(rIts, c(1,2,3), function(x) sd(x))
  
  r_quants = array(NA, c(4, N*T, K))
  for (k in 1:K){
    for (n in 1:N){
      for (t in 1:T){
        r_quants[1:3,(n-1)*T+t,k] = r_int[,n,k,t]
        r_quants[4,(n-1)*T+t,k] = r_sd[n,k,t]
      }
    }
  }
  return(r_quants)
}


# compute posterior differences and then screen for statistical and ecological significance
# rIts is an N, K, T, niter array with the proportional composition estimates
# t_diffs: vector of times to be included in the comparison
post_diffs <- function(rIts, t_diffs){
  
  N=dim(rIts)[1]
  K=dim(rIts)[2]
  niter=dim(rIts)[4]
  n_diffs = length(t_diffs)
  
  diffs = array(NA, c(N, K, n_diffs, n_diffs, niter))
  for (i in 1:n_diffs){
    for (j in 1:n_diffs){
      diffs[,,i,j,] = rIts[,,t_diffs[i],] - rIts[,,t_diffs[j],] 
    }
  }
  
  
  return(diffs)
}  

# compute posterior differences and then screen for statistical and ecological significance
# rIts is an N, K, T, niter array with the proportional composition estimates
# t_diffs: vector of times to be included in the comparison
post_diffs_bw <- function(rIts, t_diffs){
  
  N=dim(rIts)[1]
  niter=dim(rIts)[3]
  n_diffs = length(t_diffs)
  
  diffs = array(NA, c(N, n_diffs, n_diffs, niter))
  for (i in 1:n_diffs){
    for (j in 1:n_diffs){
      diffs[,i,j,] = rIts[,t_diffs[i],] - rIts[,t_diffs[j],] 
    }
  }
  
  
  return(diffs)
}  

# identify the estimates that are above eco_cut*100 % relative composition
# modele estimates non-zero proportions everywhere for all taxa, must screen to retain only those that are meaningful
# rIts is an N, K, T, niter array with the proportional composition estimates
# t_diffs: vector of times to be included in the comparison
# eco_cut is the cutoff to determine ecological significan (ex: eco_cut = 0.03 means that we ignore any taxa who are less than 3%)
eco_sig <- function(rIts, t_diffs, eco_cut){
  
  N=dim(rIts)[1]
  K = dim(rIts)[2]
  niter=dim(rIts)[4]
  n_diffs = length(t_diffs)
  
  es  = array(NA, c(N, K, n_diffs, n_diffs))
  for (i in 1:n_diffs){
    for (j in 1:n_diffs){
      
      es_i = apply(rIts[,,t_diffs[i],], c(1,2), function(x) mean(x) > eco_cut)
      es_j = apply(rIts[,,t_diffs[j],], c(1,2), function(x) mean(x) > eco_cut)
      
      es[,,i,j] = es_i | es_j
    }
  }
  
  return(es)
}  

# identify the estimates that are above eco_cut*100 % relative composition
# modele estimates non-zero proportions everywhere for all taxa, must screen to retain only those that are meaningful
# rIts is an N, K, T, niter array with the proportional composition estimates
# t_diffs: vector of times to be included in the comparison
# eco_cut is the cutoff to determine ecological significan (ex: eco_cut = 0.03 means that we ignore any taxa who are less than 3%)
eco_sig_bw <- function(rIts_bw, t_diffs, eco_cut){
  
  N=dim(rIts_bw)[1]
  niter=dim(rIts_bw)[4]
  n_diffs = length(t_diffs)
  
  es  = array(NA, c(N, n_diffs, n_diffs))
  for (i in 1:n_diffs){
    for (j in 1:n_diffs){
      
      es_i = apply(rIts_bw[,t_diffs[i],], c(1), function(x) mean(x) > eco_cut)
      es_j = apply(rIts_bw[,t_diffs[j],], c(1), function(x) mean(x) > eco_cut)
      
      es[,i,j] = es_i | es_j
    }
  }
  
  return(es)
}  


diff_q <- function(x, prob){
  n_pos = length(which(x>0))
  n_neg = length(which(x<0))
  n_x = length(x)
  if (n_pos/n_x > prob){
    return(n_pos/n_x)
  } else if (n_neg/n_x > prob){
    return(-n_neg/n_x) 
  } else {
    return(NA)
  }
}

diff_mag <- function(x, prob){
  n_pos = length(which(x>0))
  n_neg = length(which(x<0))
  n_x = length(x)
  if (n_pos/n_x > prob){
    return(median(x[which(x>0)]))
  } else if (n_neg/n_x > prob){
    return(median(x[which(x<0)])) 
  } else {
    return(NA)
  }
}

# # plot_post_change(post_prob_es, t_diffs, centers_veg, N, taxa, rescale, type='prob', subDir)
# plot_post_change <- function(post_change, t_diffs, ages, centers_veg, N, taxa, taxa_sub, limits, rescale, type, subDir){
#   
#   K = length(taxa)
#   if (type=='bw'){
#     taxa='BW'
#     K=1
#   }
#   
#   pp_melted = data.frame(change=numeric(0), x=integer(0), y=integer(0), t1=numeric(0), t2=numeric(0), taxon=character())
#   for (k in 1:K){
#     for (t2 in 2:length(t_diffs)){
#       for (t1 in (t2-1):1){
#         pp_melted  = rbind(pp_melted , data.frame(prob = as.vector(post_change[,k,t1,t2]), 
#                                                   x     = centers_veg[,1]*rescale, 
#                                                   y     = centers_veg[,2]*rescale, 
#                                                   # t1 = rep(t_diffs[t1], N),
#                                                   # t2 = rep(t_diffs[t2], N),
#                                                   t1 = rep(ages[t_diffs[t1]]*100, N),
#                                                   t2 = rep(ages[t_diffs[t2]]*100, N),
#                                                   taxon = rep(taxa[k], N)))
#       }
#     }
#   }
#   
#   if (type %in% c('prob', 'bw')){
#     values=c(-1, -0.9, -0.899, 0.899, 0.9, 1)
#     lims = c(-1.0,1.0)
#   } else if (type=='mag') {  
#     bound = ceiling(10*max(abs(c(max(post_change, na.rm=TRUE), min(post_change, na.rm=TRUE)))))/10
#     values=c(-bound, -0.3, -0.001, 0.001, 0.3, bound)
#     lims = c(-bound,bound)
#   }
#   
#   print('here')
#   
#   # pdf(file=paste0(subDir, '/post_diff_', type, '.pdf'))
#   
#   for (taxon in taxa_sub){
#     pdf(file=paste0(subDir, '/post_diff_', taxon, '_', type, '.pdf'))
#     
#     pp_taxon = pp_melted[pp_melted$taxon==taxon,]
#     
#     p <- ggplot(data=pp_taxon) + geom_tile(aes(x=x, y=y, fill=prob)) 
#     p <- p + scale_fill_gradientn(colours=c("blue", "lightskyblue", "white", "pink", "red"),  
#                                   values=values, limits=lims, na.value="white",
#                                   rescaler = function(x, ...) x, oob = identity, name='Probability') 
#     p <- p + coord_fixed()
#     p <- add_map_albers(p, map_data=us.fort, limits)
#     p <- p + facet_grid(t2~t1)
#     # p <- p + ggtitle(taxon)
#     # p <- theme_clean(p)
#     p <- theme_clean(p) + theme_bw() + theme(strip.text = element_text(size = 14),
#                                              strip.background = element_rect(colour = 'grey'),
#                                              axis.title.y = element_blank(),
#                                              axis.text = element_blank(),
#                                              axis.title = element_blank(),
#                                              axis.ticks = element_blank())
#     print(p)
#     
#     dev.off()
#   }
#   
#   # dev.off()
#   
# }

# plot_post_change(post_prob_es, t_diffs, centers_veg, N, taxa, rescale, type='prob', subDir)
plot_class_change <- function(class_change, t_diffs, ages, centers_veg, N, taxa, taxa_sub, cuts, limits, rescale, subDir, suff=''){
  
  K = length(taxa)
  
  cc_melted = data.frame(change=numeric(0), x=integer(0), y=integer(0), t1=numeric(0), t2=numeric(0), taxon=character())
  for (k in 1:K){
    for (t2 in 2:length(t_diffs)){
      for (t1 in (t2-1):1){
        cc_melted  = rbind(cc_melted , data.frame(change = as.vector(class_change[,k,t1,t2]), 
                                                  x     = centers_veg[,1]*rescale, 
                                                  y     = centers_veg[,2]*rescale, 
                                                  # t1 = rep(t_diffs[t1], N),
                                                  # t2 = rep(t_diffs[t2], N),
                                                  t1 = rep(ages[t_diffs[t1]]*100, N),
                                                  t2 = rep(ages[t_diffs[t2]]*100, N),
                                                  taxon = rep(taxa[k], N)))
      }
    }
  }
  
  cc_melted$change = factor(cc_melted$change, levels=c('no_change', 'intermediate', 'change'))#unique(cc_melted$change))
  
  if (nchar(suff>=1)){suff=paste0('_', suff)}
  pdf(file=paste0(subDir, '/class_change', suff, '.pdf'))
  
  for (taxon in taxa_sub){
    
    cc_taxon = cc_melted[cc_melted$taxon==taxon,]
    # cc_taxon$change = factor(cc_taxon$change, levels=unique(cc_melted$change))
    
    p <- ggplot(data=cc_taxon) + geom_tile(aes(x=x, y=y, fill=change)) 
    p <- p + scale_fill_manual(values=c('blue', 'orange' , 'red'), labels=levels(cc_melted$change), drop=FALSE)
    # p <- p + scale_fill_manual(drop=FALSE)
    
    p <- p + coord_fixed()
    p <- add_map_albers(p, map_data=us.fort, limits)
    p <- p + facet_grid(t2~t1)
    p <- p + ggtitle(paste0(taxon))
    p <- theme_clean(p) + theme_bw() + theme(strip.text = element_text(size = 14),
                                             strip.background = element_rect(colour = 'grey'),
                                             axis.title.y = element_blank(),
                                             axis.text = element_blank(),
                                             axis.title = element_blank(),
                                             axis.ticks = element_blank())
    print(p)
    
  }
  
  dev.off()
  
}

# plot_post_change(post_prob_es, t_diffs, ages, centers_veg, N, taxa, taxa_sub=taxa, limits, rescale, type='prob', subDir=subDir)
plot_post_change <- function(post_change, t_diffs, ages, centers_veg, N, taxa, taxa_sub, limits, rescale, type, subDir, suff=''){
  
  K = length(taxa)
  if (type=='bw'){
    taxa='BW'
    K=1
  }
  
  pp_melted = data.frame(change=numeric(0), x=integer(0), y=integer(0), t1=numeric(0), t2=numeric(0), taxon=character())
  for (k in 1:K){
    for (t2 in 2:length(t_diffs)){
      for (t1 in (t2-1):1){
        pp_melted  = rbind(pp_melted , data.frame(prob = as.vector(post_change[,k,t1,t2]), 
                                                  x     = centers_veg[,1]*rescale, 
                                                  y     = centers_veg[,2]*rescale, 
                                                  # t1 = rep(t_diffs[t1], N),
                                                  # t2 = rep(t_diffs[t2], N),
                                                  t1 = rep(ages[t_diffs[t1]]*100, N),
                                                  t2 = rep(ages[t_diffs[t2]]*100, N),
                                                  taxon = rep(taxa[k], N)))
      }
    }
  }
  
  if (type %in% c('prob', 'bw')){
    # values=c(-1, -0.9, -0.899, 0.899, 0.9, 1)
    values=c(-1, -0.9, -0.5, 0.5, 0.9, 1)
    lims = c(-1.0,1.0)
    cols=brewer.pal(11, "RdBu")
    # cols=rev(brewer.pal(11, "RdBu"))
    # cols=rev(terrain.colors(24))
    # cols=c("blue", "lightskyblue", "white", "pink", "red")
    legend_label = 'Probability'
  } else if (type=='mag') {  
    bound = ceiling(10*max(abs(c(max(post_change, na.rm=TRUE), min(post_change, na.rm=TRUE)))))/10
    values=c(-bound, -0.25, -0.001, 0.001, 0.25, bound)
    lims = c(-bound,bound)
    # cols=brewer.pal(6, "RdBu")
    cols=brewer.pal(11, "RdBu")
    # cols=c("blue", "lightskyblue", "white", "pink", "red")
    legend_label = 'Magnitude'
  } else if (type=='sd') {  
    bound = ceiling(10*max(abs(c(max(post_change, na.rm=TRUE), min(post_change, na.rm=TRUE)))))/10
    values=c(0, 0.01, 0.05, 0.10, 0.15, bound)
    lims = c(0,bound)
    cols=brewer.pal(6, "Reds")
    legend_label = 'SD'
  } else if (type=='cov'){
    bound = ceiling(10*max(abs(c(max(post_change, na.rm=TRUE), min(post_change, na.rm=TRUE)))))/10
    values=c(0, 0.5, 1, 1.5, 3, bound)
    lims = c(0,bound)
    cols=brewer.pal(6, "YlOrBr")
    legend_label = 'COV'
  }
  
  # pdf(file=paste0(subDir, '/post_diff_', type, '.pdf'))
  
  for (taxon in taxa_sub){
    pdf(file=paste0(subDir, '/post_diff_', taxon, '_', type, suff, '.pdf'))
    
    pp_taxon = pp_melted[pp_melted$taxon==taxon,]
    pp_taxon = pp_taxon[which(pp_taxon$t2>pp_taxon$t1),]
    pp_taxon$t1 = factor(pp_taxon$t1, levels=sort(unique(pp_taxon$t1), decreasing=TRUE))
    
    p <- ggplot(data=pp_taxon) + geom_tile(aes(x=x, y=y, fill=prob)) 
    p <- p + scale_fill_gradientn(colours=cols,  
                                  values=values, limits=lims, na.value="white",
                                  rescaler = function(x, ...) x, name=legend_label)#, oob = identity) 
    p <- p + coord_fixed()
    p <- add_map_albers(p, map_data=us.fort, limits)
    p <- p + facet_grid(t2~t1, switch='x')
    # p <- p + ggtitle(taxon)
    # p <- theme_clean(p)
    p <- theme_clean(p) + theme_bw() + theme(strip.text = element_text(size = 14),
                                             strip.background = element_rect(colour = 'grey'),
                                             axis.title.y = element_blank(),
                                             axis.text = element_blank(),
                                             axis.title = element_blank(),
                                             axis.ticks = element_blank(),
                                             panel.grid.major = element_blank(), 
                                             panel.grid.minor = element_blank())
    topright=c(0.9,0.9)
    # topright=c(0.25,0.08)
    p <- p + ggsn::scalebar(dist=100, 
                            transform=FALSE, 
                            dist_unit="km",
                            location="topright",
                            height=0.03,
                            st.dist=0.03,
                            st.size=2,
                            border.size=0.4,
                            anchor=c(x=(1-topright[1])*limits$xlims[1]*1000000+topright[1]*limits$xlims[2]*1000000,
                                     y=(1-topright[2])*limits$ylims[1]*1000000+topright[2]*limits$ylims[2]*1000000),
                            x.min=limits$xlims[1]*1000000,
                            x.max=limits$xlims[2]*1000000,
                            y.min=limits$ylims[1]*1000000,
                            y.max=limits$ylims[2]*1000000)
    print(p)
    
    dev.off()
  }
  
  # dev.off()
  
}

plot_post_change_meanmag<- function(post_change, t_diffs, ages, centers_veg, N, taxa, taxa_sub, limits, rescale, type, subDir, suff=''){
  
  K = length(taxa)
  if (type=='bw'){
    taxa='BW'
    K=1
  }
  
  pp_melted = data.frame(change=numeric(0), x=integer(0), y=integer(0), t1=numeric(0), t2=numeric(0), taxon=character())
  for (k in 1:K){
    for (t2 in 2:length(t_diffs)){
      for (t1 in (t2-1):1){
        pp_melted  = rbind(pp_melted , data.frame(prob = as.vector(post_change[,k,t1,t2]), 
                                                  x     = centers_veg[,1]*rescale, 
                                                  y     = centers_veg[,2]*rescale, 
                                                  # t1 = rep(t_diffs[t1], N),
                                                  # t2 = rep(t_diffs[t2], N),
                                                  t1 = rep(ages[t_diffs[t1]]*100, N),
                                                  t2 = rep(ages[t_diffs[t2]]*100, N),
                                                  taxon = rep(taxa[k], N)))
      }
    }
  }
  
  # values=c(-1, -0.9, -0.899, 0.899, 0.9, 1)
  values=c(-1, -0.9, -0.5, 0.5, 0.9, 1)
  lims = c(-1.0,1.0)
  # cols=brewer.pal(11, "RdBu")
  cols=rev(brewer.pal(11, "RdBu"))
  # cols=rev(terrain.colors(24))
  # cols=c("blue", "lightskyblue", "white", "pink", "red")
  legend_label = 'Probability'
  
 
  # pdf(file=paste0(subDir, '/post_diff_', type, '.pdf'))
  
  for (taxon in taxa_sub){
    pdf(file=paste0(subDir, '/post_diff_', taxon, '_', type, suff, '.pdf'))
    
    pp_taxon = pp_melted[pp_melted$taxon==taxon,]
    pp_taxon = pp_taxon[which(pp_taxon$t2>pp_taxon$t1),]
    pp_taxon$t1 = factor(pp_taxon$t1, levels=sort(unique(pp_taxon$t1), decreasing=TRUE))
    
    p <- ggplot(data=pp_taxon) + geom_tile(aes(x=x, y=y, fill=prob)) 
    p <- p + scale_fill_gradientn(colours=cols,  
                                  values=values, limits=lims, na.value="white",
                                  rescaler = function(x, ...) x, name=legend_label)#, oob = identity) 
    p <- p + coord_fixed()
    p <- add_map_albers(p, map_data=us.fort, limits)
    p <- p + facet_grid(t2~t1, switch='x')
    # p <- p + ggtitle(taxon)
    # p <- theme_clean(p)
    p <- theme_clean(p) + theme_bw() + theme(strip.text = element_text(size = 14),
                                             strip.background = element_rect(colour = 'grey'),
                                             axis.title.y = element_blank(),
                                             axis.text = element_blank(),
                                             axis.title = element_blank(),
                                             axis.ticks = element_blank())
    print(p)
    
    dev.off()
  }
  
  # dev.off()
  
}


# plot_post_change(post_prob_es, t_diffs, centers_veg, N, taxa, rescale, type='prob', subDir)
plot_post_change_cov <- function(post_change, pol, t_diffs, ages, centers_veg, N, taxa, taxa_sub, limits, rescale, subDir){
  
  K = length(taxa)
  
  pp_melted = data.frame(change=numeric(0), x=integer(0), y=integer(0), t1=numeric(0), t2=numeric(0), taxon=character())
  for (k in 1:K){
    for (t2 in 2:length(t_diffs)){
      for (t1 in (t2-1):1){
        pp_melted  = rbind(pp_melted , data.frame(prob = as.vector(post_change[,k,t1,t2]), 
                                                  x     = centers_veg[,1]*rescale, 
                                                  y     = centers_veg[,2]*rescale, 
                                                  # t1 = rep(t_diffs[t1], N),
                                                  # t2 = rep(t_diffs[t2], N),
                                                  t1 = rep(ages[t_diffs[t1]]*100, N),
                                                  t2 = rep(ages[t_diffs[t2]]*100, N),
                                                  taxon = rep(taxa[k], N)))
      }
    }
  }
  
  bound = ceiling(10*max(abs(c(max(post_change, na.rm=TRUE), min(post_change, na.rm=TRUE)))))/10
  values=c(0, 0.5, 1, 1.5, 3, bound)
  lims = c(0,bound)
  cols=brewer.pal(6, "YlOrBr")
  legend_label = 'COV'
  
  # pdf(file=paste0(subDir, '/post_diff_', type, '.pdf'))
  
  for (taxon in taxa_sub){
    pdf(file=paste0(subDir, '/post_diff_cov_', taxon, '.pdf'))
    
    pp_taxon = pp_melted[pp_melted$taxon==taxon,]
    
    p <- ggplot(data=pp_taxon) + geom_tile(aes(x=x, y=y, fill=prob)) 
    p <- p + scale_fill_gradientn(colours=cols,  
                                  values=values, limits=lims, na.value="white",
                                  rescaler = function(x, ...) x, name=legend_label) #oob = identity, 
    p <- p + geom_point(data=pol, aes(x=x, y=y), alpha=0.4, colour="grey35", shape=1)
    p <- p + coord_fixed()
    p <- add_map_albers(p, map_data=us.fort, limits)
    p <- p + facet_grid(t2~t1)
    p <- p + ggtitle(taxon)
    # p <- theme_clean(p)
    p <- theme_clean(p) + theme_bw() + theme(strip.text = element_text(size = 14),
                                             strip.background = element_rect(colour = 'grey'),
                                             axis.title.y = element_blank(),
                                             axis.text = element_blank(),
                                             axis.title = element_blank(),
                                             axis.ticks = element_blank())
    print(p)
    
    dev.off()
  }
  
  # dev.off()
  
}


# plot_post_change(post_prob_es, t_diffs, centers_veg, N, taxa, rescale, type='prob', subDir)
plot_post_change_group <- function(post_change, t_diffs, ages, centers_veg, N, 
                                   taxa, rescale, limits, type, group, subDir){
  
  K = 1
  
  pp_melted = data.frame(change=numeric(0), x=integer(0), y=integer(0), t1=numeric(0), t2=numeric(0), taxon=character())
  for (t2 in 2:length(t_diffs)){
    for (t1 in (t2-1):1){
      pp_melted  = rbind(pp_melted , data.frame(prob = as.vector(post_change[,t1,t2]), 
                                                x     = centers_veg[,1]*rescale, 
                                                y     = centers_veg[,2]*rescale, 
                                                t1 = rep(ages[t_diffs[t1]], N)*100,
                                                t2 = rep(ages[t_diffs[t2]], N)*100,
                                                taxon = rep(group, N)))
    }
  }
  
  if (type %in% c('prob')){
    values=c(-1, -0.9, -0.899, 0.899, 0.9, 1)
    lims = c(-1.0,1.0)
    cols=brewer.pal(11, "RdBu")
  } else if (type=='mag') {  
    bound = ceiling(10*max(abs(c(max(post_change, na.rm=TRUE), min(post_change, na.rm=TRUE)))))/10
    values=c(-bound, -0.3, -0.001, 0.001, 0.3, bound)
    lims = c(-bound,bound)
  }
  
  print('here')
  
  # pdf(file=paste0(subDir, '/post_diff_', type, '.pdf'))
  
  pdf(file=paste0(subDir, '/post_diff_', type, '_', group,'.pdf'))
  
  pp_taxon = pp_melted
  pp_taxon$t1 = factor(pp_taxon$t1, levels=sort(unique(pp_taxon$t1), decreasing=TRUE))
  
  p <- ggplot(data=pp_taxon) + geom_tile(aes(x=x, y=y, fill=prob)) 
  p <- p + scale_fill_gradientn(colours=cols,  
                                values=values, limits=lims, na.value="white",
                                rescaler = function(x, ...) x)#, oob = identity, name='Probability') 
  p <- p + coord_fixed()
  p <- add_map_albers(p, map_data=us.fort, limits)
  p <- p + geom_path(data=bw, aes(x=long, y=lat),  colour='grey', alpha=0.8)
  p <- p + facet_grid(t2~t1, switch='x')
  p <- theme_clean(p) + theme_bw() + theme(strip.text = element_text(size = 14),
                                           strip.background = element_rect(colour = 'grey'),
                                           axis.title.y = element_blank(),
                                           axis.text = element_blank(),
                                           axis.title = element_blank(),
                                           axis.ticks = element_blank(),
                                           panel.grid.major = element_blank(), 
                                           panel.grid.minor = element_blank())
  topright=c(0.9,0.9)
  # topright=c(0.25,0.08)
  p <- p + ggsn::scalebar(dist=100, 
                          transform=FALSE, 
                          dist_unit="km",
                          location="topright",
                          height=0.03,
                          st.dist=0.03,
                          st.size=2,
                          border.size=0.4,
                          anchor=c(x=(1-topright[1])*limits$xlims[1]*1000000+topright[1]*limits$xlims[2]*1000000,
                                   y=(1-topright[2])*limits$ylims[1]*1000000+topright[2]*limits$ylims[2]*1000000),
                          x.min=limits$xlims[1]*1000000,
                          x.max=limits$xlims[2]*1000000,
                          y.min=limits$ylims[1]*1000000,
                          y.max=limits$ylims[2]*1000000)
  print(p)
  dev.off()
  
}

# distance funs
squared_chord <- function(x, y){
  sum((sqrt(x) - sqrt(y))^2)
}

canberra <- function(x, y){
  sum(abs(x - y)/(abs(x) + abs(y)))
}

hellinger <- function(x, y){
  1/sqrt(2)*sqrt(sum((sqrt(x) - sqrt(y))^2))
}

bray_curtis <- function(x, y){
  sum(abs(x - y))/sum(x + y)
}



get_diags <- function(centers_veg){
  
  # if (n_diags==1){y_start=0.862}else{y_start=0.862+0.024}
  
  # if (n_diags>5){stop("n_diags must be less than 4!")}
  
  p_start = data.frame(rbind(c(0,15),
                             c(0,12),
                             c(0, 9),
                             c(0, 6),
                             c(0, 3),
                             c(0, 0), 
                             c(3, 0),
                             c(6, 0),
                             c(9, 0),
                             c(12,0),
                             c(15,0)))
  
  coords = data.frame(x=numeric(0), y=numeric(0), diag=numeric(0))
  for (i in 1:nrow(p_start)){
    
    p1 = c(-0.059, 0.862) + p_start[i,]*0.024 + c(0.024,0.024)
    # p1 = c(-0.059, 0.862) + c(0.024,0.024)
    coord_old = p1
    in_domain = TRUE
    coords_diag = p1
    while (in_domain) {
      coord_new = coord_old + c(0.024,0.024)
      # coord_new = coord_old + c(0.024,2*0.024)
      if (rdist(matrix(coord_new, nrow=1), matrix(c(0.205, 1.126), nrow=1)) < 1e-10) {
        in_domain == TRUE
      } else {
        in_domain = any(rdist(matrix(coord_new, nrow=1), centers_veg) < 1e-10)
      }
      if (in_domain) {
        coords_diag = rbind(coords_diag, coord_new)
        coord_old = coord_new
      }
    }
    
    coords = rbind(coords, data.frame(x=coords_diag[,1], y=coords_diag[,2], diag=rep(i)))
  }
  
  return(coords)
}

diag_transect <- function(dat, coords_diag, rescale_coords=1e6, rescale_dat=1){
  dat = dat[with(dat, order(x)), ]
  d_diag = rdist(coords_diag[,1:2]*rescale_coords, matrix(cbind(dat$x, dat$y)*rescale_dat, ncol=2))
  idx_keep = apply(d_diag, 2, function(x) any(x < 1e-8))
  idx_tran = apply(d_diag, 2, function(x) {if(any(x < 1e-8)){ which(x < 1e-8)}else{NA}})
  idx_tran = unlist(idx_tran)
  
  dat_tran = dat[idx_keep, ]
  
  dat_tran = data.frame(dat_tran, diag=coords_diag[idx_tran[idx_keep],3])
  
  return(dat_tran)
}

diff_change <- function(x, cut_lo, cut_hi, prob){
  n_hi  = length(which(abs(x)>cut_hi))
  n_lo  = length(which(abs(x)<cut_lo))
  n_mid = length(which((cut_lo<abs(x)) & (abs(x)<cut_hi)))
  
  n_x = length(x)
  if (n_hi/n_x > prob){
    return('change')
  } else if (n_lo/n_x > prob){
    return('no_change') 
  } else if (n_mid/n_x > prob) {
    return('intermediate')
  } else {
    return(NA)
  }
  
}

diff_change_above <- function(x, cut_change, prob){
  n_change  = length(which(abs(x)>cut_change))
  
  n_x = length(x)
  if (n_change/n_x > prob){
    return('change')
  } else if (n_change/n_x < prob){
    return('no_change')
  } else {
    return(NA)
  }
}


diff_change_below <- function(x, cut_change, prob){
  n_change  = length(which(abs(x)<cut_change))
  
  n_x = length(x)
  if (n_change/n_x > prob){
    return('no_change')
  } else if (n_change/n_x < prob){
    return('uncertain')
  } else {
    return(NA)
  }
}