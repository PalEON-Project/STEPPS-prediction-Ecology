library(gridExtra)
library(maptools)

us.shp <- readOGR('data/map/us/us_alb.shp')
us.shp@data$id <- rownames(us.shp@data)
us.fort <- fortify(us.shp, region='id') 

bw <- readOGR('data/map/bw/bw_albers.shp')

create_figure_path <- function(subDir){
  mainDir <- getwd()
  
  figure_path = file.path(mainDir, subDir)
  
  if (file.exists(subDir)){
    print('Folder already exists: contents may get overwritten!')
  } else {
    dir.create(figure_path)
  }
  
  print(paste('All figures will be saved in: ', figure_path, sep=''))
  
#    return(figure_path)
}

# trace plots
# fit is a stanfit object
# can pass true values as a list if known
trace_pars <- function(post, N_knots, T, N_pars, taxa, suff=suff, save_plots=TRUE, fpath=subDir){
  
#   post = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
#   post = post[250:500,,]
  n    = dim(post)[3]
  
  #don't plot knots here
  idx_pars = c(seq(1,N_pars),n)
  
  labels = colnames(post[,1,])[idx_pars]
  
  labels[2:(length(labels)-1)] =  paste(labels[2:(length(labels)-1)], taxa[1:(K-1)], sep=' ')

#   avg = summary(fit)$summary[,"mean"]
  
  par(mfrow=c(4,2))
  if (save_plots){
    pdf(paste(fpath, "/trace_pars.pdf", sep=""), width=8, height=12)
    par(mfrow=c(5,2))
  }
  
  
  for (i in idx_pars){
    plot(post[,1,i], type="l", ylab=labels[i], xlab="iter")
#     plot(post[,i], type="l", ylab=labels[i], xlab="iter")
    abline(h=mean(post[,1,i]), col="blue")
    abline(h=quantile(post[,1,i],probs=0.025), col='blue', lty=2)
    abline(h=quantile(post[,1,i],probs=0.975), col='blue', lty=2)
    #     if (dim(post)[2] >= 2){
    #       lines(post[,2,i], col="blue")
    #     }
  }
  
  if (save_plots){
    dev.off()
  }

}

# trace plots
# fit is a stanfit object
# can pass true values as a list if known
trace_plot_mut <- function(post, N_knots, T, N_pars, taxa, mean_type='other', suff=suff, save_plots=TRUE, fpath=subDir){
  
  n    = dim(post)[3]
  
  #idx_pars = seq(3+N_pars+1, 3+N_pars+W*T)
  #idx_pars = seq(N_pars+1, N_pars+W*T)
  idx_pars = seq(N_pars+1, N_pars+W*(T-1))
  
  if (mean_type=='MRF'){
    idx_pars = seq(3, 3+W*T)
  }
#   
#   labels = colnames(post[,1,])[idx_pars]
  
#   labels[4:(length(labels)-1)] =  paste(labels[4:(length(labels)-1)], taxa[1:(K-1)], sep=' ')
  
#   avg = summary(fit)$summary[,"mean"]
  
  
  if (save_plots){
    pdf(paste(fpath, "/trace_mut.pdf", sep=""), width=8, height=12)
    par(mfrow=c(5,2))
  }
  
  for (k in 1:W){
    par(mfrow=c(5,2))
    print(k)
    idx_taxon = seq(k, W*(T-1), by=W)
    #idx_taxon = seq(k, W*T, by=W)
    idx = idx_pars[idx_taxon]
    labels = colnames(post[,1,])[idx]
    draws = post[,1,idx]
    
    for (i in 1:length(idx)){
      plot(draws[,i], type="l", ylab=labels[i], xlab="iter")
      #     plot(post[,i], type="l", ylab=labels[i], xlab="iter")
      abline(h=mean(draws[,i]), col="blue")
      abline(h=quantile(draws[,i],probs=0.025), col='blue', lty=2)
      abline(h=quantile(draws[,i],probs=0.975), col='blue', lty=2)
      #     if (dim(post)[2] >= 2){
      #       lines(post[,2,i], col="blue")
      #     }
    }
  }

  if (save_plots){
    dev.off()
  }
  
}


knot_idx <- function(w, n, t,K){
  K + 1 + (n-1)*T*W + (t-1)*W + w-1
}

# trace plots of knots
# each knot on a different panels, multiple times per pane
trace_plot_knots <- function(fit, N_knots, T, K, N_pars, suff=suff, save_plots=TRUE, fpath=subDir){
  
  post = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
  labels = colnames(post[,1,])
  niter = dim(post)[1]
  
#   t = 1
#   n = 2
#   w = 1
  
  t = seq(1,T)
  
  
#   labels[6 + (n-1)*T*W + (t-1)*W + w-1]
#   knot = post[,1,6 + (n-1)*T*W + (t-1)*W + w-1]
  
  par(mfrow=c(3,2))
  
  if (save_plots){
    pdf(paste(fpath, "/trace_knots_", suff, ".pdf", sep=""), width=8, height=4)
    par(mfrow=c(1,1))
  }
  
 
  for (w in 1:W){
    for (n in 1:N_knots){
      
      post_knots = post[,1,knot_idx(w,n,t,N_pars)]
  
      ylo = min(post_knots)
      yhi = max(post_knots)
  
      plot(c(0,0), type="n", ylab=paste('Taxon', w, '; Knot', n ,sep=' '), xlab="iter", ylim=c(ylo, yhi), xlim=c(0,niter))
      #plot(c(0,0), type="n", ylab=paste('Knot', n ,sep=' '), xlab="iter", ylim=c(0.0379, 0.038), xlim=c(0,niter))
      cols=c("red","blue", "black")
      for (i in 1:T){
        lines(post_knots[,i], col=cols[i])
      }
    }
  }
  
  if (save_plots){
    dev.off()
  }

}

# trace plots
# fit is a stanfit object
# can pass true values as a list if known
trace_plot_process <- function(r, suff=suff, save_plots, fpath=subDir){
  
  K = dim(r)[2] # n taxa
  NT = dim(r)[1]
  skip = ceiling(NT/15)
  cells = seq(1, NT, by=skip)
  
  cols = rep('black', K)#rainbow(K)
  
  par(mfrow=c(4,2))
  if (save_plots){
    if (nchar(suff)>0) suff1 = paste0('_', suff)
    pdf(paste(fpath, "/trace", suff1, ".pdf", sep=""), width=8, height=4)
    par(mfrow=c(1,1))
  }
  
  if (suff=='r'){
    minr = 0
    maxr = 1
  } else {
    minr = min(r)
    maxr = max(r)
  }
  
  for (i in 1:length(cells)){
    for (k in 1:K){
      plot(r[cells[i],k,], type="l", col=cols[k], ylab=paste(suff, cells[i], 'taxon', k, sep=' '), 
           xlab="iter", ylim=c(minr,maxr))
    }
  }

  
  if (save_plots){
    dev.off()
  }
}

get_limits <- function(centers){
  xlo = min(centers[,1])
  xhi = max(centers[,1])

  ylo = min(centers[,2])
  yhi = max(centers[,2])

  return(list(xlims=c(xlo,xhi),ylims=c(ylo, yhi)))
}  

plot_pred_maps <- function(r_mean, centers, taxa, t, N, K, T, thresh, limits, type, suff, save_plots, fpath=subDir){
  
  rescale=1000000
  
  if (is.null(taxa)){taxa=seq(1,K)}
  if (type=='prop'){bar_title='Proportions'}else {bar_title='Values'}
  
  prop_dat = data.frame(props=numeric(0), x=integer(0), y=integer(0),time=integer(0), taxon=character())
  for (k in 1:K){
    prop_dat = rbind(prop_dat, data.frame(props = r_mean[,k], 
                                          x     = rep(centers[,1], each=T)*rescale, 
                                          y     = rep(centers[,2], each=T)*rescale, 
                                          time  = rep(t,times=N), 
                                          taxon = rep(taxa[k], N*T)))
  }
  
  if (!is.na(thresh)){
    prop_dat$props[which(prop_dat$props > thresh)] = thresh
  }
    
  p <- ggplot() + geom_tile(data=prop_dat, aes(x=x, y=y, fill=props)) + 
    scale_fill_gradientn(colours=tim.colors(), name=bar_title) + coord_fixed() +
    scale_x_continuous(limits$xlims) + scale_y_continuous(limits$ylims)
  p <- add_map_albers(p, map_data=us.fort, limits)
  p <- p + facet_grid(time~taxon)
  p <- theme_clean(p) + theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
  print(p)
  Sys.sleep(2)
  if (save_plots){
    ggsave(file=paste(fpath, '/veg_maps_', suff, '.pdf', sep=''), scale=1, width=12)
#     dev.off()
  }
  return(p)
}


# plot_pred_maps_select(r_mean, centers_veg, taxa=taxa, ages, N, K, T, thresh=0.5, limits, type='prop', suff=suff,  save_plots=save_plots)

plot_pred_maps_select <- function(r_mean, centers, taxa, ages, N, K, T, thresh, limits, type, suff, save_plots, fpath=subDir){
  
  if (is.null(taxa)){taxa=seq(1,K)}
  if (!is.null(taxa)){
    taxa[which(taxa == 'OTHER.HARDWOOD')] = 'OTHER\nHW'
    taxa[which(taxa == 'OTHER.CONIFER')] = 'OTHER\nCON'
    taxa[which(taxa == 'TAMARACK')] = 'LARCH'
    
#     taxa = as.vector(sapply(taxa, simpleCap))
  }
  
  if (type=='prop'){bar_title='Proportions'}else {bar_title='Values'}
  
  # find a better way to do this later
  if (length(ages) > 2){
    idx.keep  = c(1,2,length(ages)/2,T)
    ages.keep = ages[idx.keep]
    T.keep    = length(ages.keep)
  } else if (length(ages) == 1){
    idx.keep = 1
    ages.keep = ages
    T.keep = length(ages.keep)
  }
  
  idx_r_keep = vector(length=0)
  for (i in 1:T.keep){
    idx_orig   = seq(idx.keep[i], N*T, by=T)
    idx_r_keep = c(idx_r_keep, idx_orig)
  }
  
  idx_r_keep = sort(idx_r_keep)
  r_mean_keep = r_mean[idx_r_keep,]
  
  prop_dat = data.frame(props=numeric(0), x=integer(0), y=integer(0),time=integer(0), taxon=character())
  for (k in 1:K){
    prop_dat = rbind(prop_dat, data.frame(props = r_mean_keep[,k], 
                                          x     = rep(centers[,1]*1000000, each=T.keep), 
                                          y     = rep(centers[,2]*1000000, each=T.keep), 
                                          time  = rep(ages.keep, times=N), 
                                          taxon = rep(taxa[k], N*T.keep)))
  }
  
  if (!is.na(thresh)){
    prop_dat$props[which(prop_dat$props > thresh)] = thresh
  }

  p <- ggplot() + geom_tile(data=prop_dat, aes(x=x, y=y, fill=props)) + 
    scale_fill_gradientn(colours=tim.colors(), name=bar_title) + coord_fixed() #+
    #scale_x_continuous(limits$xlims*1000) + scale_y_continuous(limits$ylims*1000)
  p <- add_map_albers(p, map_data=us.fort, limits)
  p <- p + facet_grid(time~taxon)
  p <- theme_clean(p) + theme(strip.text.y = element_text(size = rel(1.3)), strip.text.x = element_text(size = rel(1.0)))
#   print(p)
  Sys.sleep(2)
  if (save_plots){
    fname = paste0(fpath, '/veg_maps_', suff, '.pdf')
    ggsave(file=fname, scale=1, width=12, height=12)
    sys_str = paste("pdfcrop", fname, fname, sep=' ')
    system(sys_str)
    #     dev.off()
  }

  return(p)
}

# plot_pred_maps_binned_select(r_mean, centers_veg, breaks, taxa, ages, N, K, T, limits, suff=suff_figs, save_plots, fpath=subDir)

plot_pred_maps_binned_select <- function(r_mean, centers, breaks, taxa, taxa_sub, ages, N, K, T, limits, suff, save_plots, fpath=subDir){
  
  if (is.null(taxa)){taxa=seq(1,K)}
  if (!is.null(taxa)){
    taxa[which(taxa == 'OTHER.HARDWOOD')] = 'OTHER\nHW'
    taxa[which(taxa == 'OTHER.CONIFER')] = 'OTHER\nCON'
    taxa[which(taxa == 'TAMARACK')] = 'LARCH'
    
    #     taxa = as.vector(sapply(taxa, simpleCap))
  }
  
  if (!is.null(taxa_sub)){
    taxa_sub[which(taxa_sub == 'OTHER.HARDWOOD')] = 'OTHER\nHW'
    taxa_sub[which(taxa_sub == 'OTHER.CONIFER')] = 'OTHER\nCON'
    taxa_sub[which(taxa_sub == 'TAMARACK')] = 'LARCH'
    
    #     taxa = as.vector(sapply(taxa, simpleCap))
  }
  
  if (length(ages) > 2){
    idx.keep  = c(1,2,T/4,length(ages)/2,3*T/4, T)
    ages.keep = ages[idx.keep]
    T.keep    = length(ages.keep)
  } else if (length(ages) == 1){
    idx.keep  = 1
    ages.keep = ages
    T.keep    = length(ages.keep)
  }
  
  idx_r_keep = vector(length=0)
  for (i in 1:T.keep){
    idx_orig   = seq(idx.keep[i], N*T, by=T)
    idx_r_keep = c(idx_r_keep, idx_orig)
  }
  
  idx_r_keep = sort(idx_r_keep)
  r_mean_keep = r_mean[idx_r_keep,]

  r_mean_binned = matrix(0, nrow=nrow(r_mean_keep), ncol=ncol(r_mean_keep))
  colnames(r_mean_binned) <- colnames(r_mean_keep)
  
  for (i in 1:ncol(r_mean)){
    r_mean_binned[,i] = cut(r_mean_keep[,i], breaks, include.lowest=TRUE, labels=FALSE)
  }
  
  breaklabels = apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1, 
                      function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })
  
  prop_dat = data.frame(props=numeric(0), x=integer(0), y=integer(0),time=integer(0), taxon=character())
  for (k in 1:K){
    prop_dat = rbind(prop_dat, data.frame(props = r_mean_binned[,k], 
                                          x     = rep(centers[,1]*1000000, each=T.keep), 
                                          y     = rep(centers[,2]*1000000, each=T.keep), 
                                          time  = rep(ages.keep*100, times=N), 
                                          taxon = rep(taxa[k], N*T.keep)))
  }
  
  if (length(taxa_sub)>0){  
    prop_dat = prop_dat[prop_dat$taxon %in% taxa_sub, ]
  }
  
#   prop_dat$taxon[prop_dat$taxon == 'OTHER.HARDWOOD'] = 'OTHER HW'
#   prop_dat$taxon[prop_dat$taxon == 'OTHER.CONIFER']  = 'OTHER CON'
#   prop_dat$taxon[prop_dat$taxon == 'TAMARACK']  = 'LARCH'
    
  p <- ggplot() + geom_tile(data=prop_dat, aes(x=x, y=y, fill=factor(props))) + 
    scale_fill_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportion') + 
    coord_fixed() #+ scale_x_continuous(limits$xlims*1000) + scale_y_continuous(limits$ylims*1000)
  p <- add_map_albers(p, map_data=us.fort, limits)
  p <- p + facet_grid(time~taxon)
#   p <- theme_clean(p) + theme(strip.text.x = element_blank(), strip.text.y = element_text(size = rel(1.5)))
  # p <- theme_clean(p) + theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.2)))
  #theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
  p <- theme_clean(p) + theme_tufte() + theme(strip.text = element_text(size = 12),
                                           strip.background = element_rect(colour = 'grey'),
                                           axis.title.y = element_blank(),
                                           axis.text = element_blank(),
                                           axis.title = element_blank(),
                                           axis.ticks = element_blank())
#   print(p)
  Sys.sleep(2)
  if (save_plots){
    fname = paste0(fpath, '/veg_maps_props_binned_', suff, '.pdf')	
    ggsave(file=fname, scale=1, width=12, height=12)
    sys_str = paste("pdfcrop", fname, fname, sep=' ')
    system(sys_str)
  }
  return(p)
}

# # plot_pred_maps_binned_gif
# p_binned <- plot_pred_maps_binned_gif(r_mean2, centers_veg, breaks, taxa, taxa_sub=taxa, ages, N, K, T, limits, 
#                                       suff=run$suff_figs, save_plots=TRUE, fpath=subDir)
plot_pred_maps_binned_gif <- function(r_mean, centers, breaks, taxa, taxa_sub, ages, N, K, T, limits, suff, save_plots, fpath=subDir){
  
  if (is.null(taxa)){taxa=seq(1,K)}
  ages = ages * 100

  r_mean_binned = matrix(0, nrow=nrow(r_mean), ncol=ncol(r_mean))
  colnames(r_mean_binned) <- colnames(r_mean)
  
  for (i in 1:ncol(r_mean)){
    r_mean_binned[,i] = cut(r_mean[,i], breaks, include.lowest=TRUE, labels=FALSE)
  }
  
  breaklabels = apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1, 
                      function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })
  
  prop_dat = data.frame(props=numeric(0), x=integer(0), y=integer(0),time=integer(0), taxon=character())
  for (k in 1:K){
    prop_dat = rbind(prop_dat, data.frame(props = r_mean_binned[,k], 
                                          x     = rep(centers[,1]*1000000, each=T), 
                                          y     = rep(centers[,2]*1000000, each=T), 
                                          time  = rep(ages, times=N), 
                                          taxon = rep(taxa[k], N*T)))
  }
  
  prop_dat$props = as.factor(prop_dat$props)
  
  if (length(taxa_sub)>0){  
    prop_dat = prop_dat[prop_dat$taxon %in% taxa_sub, ]
  }
  
  dir.create(file.path(subDir, 'gif'), showWarnings = FALSE)
  
  for (taxon in taxa_sub){
    taxonnm = taxon
    if (taxon == 'OTHER.HARDWOOD'){taxonnm = 'OTHER_HARDWOOD'}
    if (taxon == 'OTHER.CONIFER'){taxonnm = 'OTHER_CONIFER'}
    dir.create(file.path(subDir, 'gif', taxonnm), showWarnings = FALSE)
  # taxon = 'ASH'
    for (t in ages){
      dat_sub = prop_dat[which((prop_dat$taxon == taxon) & (prop_dat$time == t)),]
      
      print(taxon)
      
      if (nchar(t) <= max(nchar(ages))){tnm = paste0(paste(rep(0, max(nchar(ages))-nchar(t)), collapse="") , t)} else {tnm = t}
      png(paste0(file.path(subDir, 'gif', taxonnm), '/veg_props_', taxonnm, '_', tnm, '.png'), width=960, height=960)
      p <- ggplot() + geom_tile(data=dat_sub, aes(x=x, y=y, fill=props)) + 
        scale_fill_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportion', drop=FALSE) + 
        coord_fixed() #+ scale_x_continuous(limits$xlims*1000) + scale_y_continuous(limits$ylims*1000)
      p <- add_map_albers(p, map_data=us.fort, limits)
      # p <- p + ggtitle(paste0(taxon, ' @ ', t*100-50, ' - ', t*100-50, ' YBP'))
      # p <- p + labs(title = paste0(taxon, ' @ ', t*100-50, ' - ', t*100+50, ' YBP'))
      p <- p + labs(title = paste0(taxon, ' @ ', t, ' YBP'))
      #   p <- theme_clean(p) + theme(strip.text.x = element_blank(), strip.text.y = element_text(size = rel(1.5)))
      p <- theme_clean(p) + theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.2))) +
           theme(plot.title = element_text(size=22))
      #theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
      print(p)
      dev.off()
    }
    
    #"ls *.png | sort -r"
    fnames = sort(list.files(file.path(getwd(), subDir, 'gif', taxonnm), '*.png'), decreasing=TRUE)
    fnames_abs = paste(file.path(getwd(), subDir, 'gif', taxonnm, fnames), collapse=' ')
    sys_call = paste0("convert -delay 80 ", fnames_abs, " ", file.path(getwd(), subDir, 'gif'), "/", taxonnm, ".gif")
    system(sys_call)
    
  }
  
  return(p)
}

plot_pred_maps_binned_gifun <- function(r_mean, r_sd, centers, breaks, taxa, taxa_sub, ages, N, K, T, limits, suff, save_plots, fpath=subDir){
  
  if (is.null(taxa)){taxa=seq(1,K)}
  
  r_mean_binned = matrix(0, nrow=nrow(r_mean), ncol=ncol(r_mean))
  colnames(r_mean_binned) <- colnames(r_mean)
  
  for (i in 1:ncol(r_mean)){
    r_mean_binned[,i] = cut(r_mean[,i], breaks, include.lowest=TRUE, labels=FALSE)
  }
  
  # classIntervals(r_sd, 10, stype="quantile")
  breaks_sd = c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25) 
  
  r_sd_binned = matrix(0, nrow=nrow(r_sd), ncol=ncol(r_sd))
  colnames(r_sd_binned) <- colnames(r_sd)
  for (i in 1:ncol(r_sd)){
    r_sd_binned[,i] = cut(r_sd[,i], breaks_sd, include.lowest=TRUE, labels=FALSE)
  }
  
  
  
  breaklabels = apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1, 
                      function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })
  
  prop_dat = data.frame(props=numeric(0), sd=numeric(0), x=integer(0), y=integer(0),time=integer(0), taxon=character())
  for (k in 1:K){
    prop_dat = rbind(prop_dat, data.frame(props = r_mean_binned[,k], 
                                          sd    = r_sd_binned[,k],
                                          x     = rep(centers[,1]*1000000, each=T), 
                                          y     = rep(centers[,2]*1000000, each=T), 
                                          time  = rep(ages, times=N), 
                                          taxon = rep(taxa[k], N*T)))
  }
  
  prop_dat$props = as.factor(prop_dat$props)
  prop_dat$sd    = as.factor(prop_dat$sd)
  
  if (length(taxa_sub)>0){  
    prop_dat = prop_dat[prop_dat$taxon %in% taxa_sub, ]
  }
  
  dir.create(file.path(fpath, 'gifsd'), showWarnings = FALSE)
  
  for (taxon in taxa_sub){
    taxonnm = taxon
    if (taxon == 'OTHER.HARDWOOD'){taxonnm = 'OTHER_HARDWOOD'}
    if (taxon == 'OTHER.CONIFER'){taxonnm = 'OTHER_CONIFER'}
    dir.create(file.path(fpath, 'gifsd', taxonnm), showWarnings = FALSE)
    # taxon = 'ASH'
    for (t in ages){
      dat_sub = prop_dat[which((prop_dat$taxon == taxon) & (prop_dat$time == t)),]
      
      print(taxon)
      
      if (nchar(t) == 3){tnm = paste0(0, t)} else {tnm = t}
      png(paste0(file.path(subDir, 'gifsd', taxonnm), '/veg_props_', taxonnm, '_', tnm, '.png'), width=960, height=960)
      p <- ggplot() + geom_tile(data=dat_sub, aes(x=x, y=y, fill=props)) + 
        geom_point(data=dat_sub, aes(x=x, y=y, shape=sd, size=sd)) + 
        scale_fill_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportion', drop=FALSE) + 
        scale_shape_manual(values=rep(1,5), drop=FALSE, name='sd') + 
        scale_size_manual(values=seq(1,5), name='sd', labels=seq(1,5), drop=FALSE) + 
        # scale_shape_manual(values=c(17, 19, 1, 95, 20), drop=FALSE) + 
        coord_fixed() #+ scale_x_continuous(limits$xlims*1000) + scale_y_continuous(limits$ylims*1000)
      p <- add_map_albers(p, map_data=us.fort, limits)
      # p <- p + ggtitle(paste0(taxon, ' @ ', t*100-50, ' - ', t*100-50, ' YBP'))
      p <- p + labs(title = paste0(taxon, ' @ ', t*100-50, ' - ', t*100+50, ' YBP'))
      #   p <- theme_clean(p) + theme(strip.text.x = element_blank(), strip.text.y = element_text(size = rel(1.5)))
      p <- theme_clean(p) + theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.2))) +
        theme(plot.title = element_text(size=22))
      #theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
      print(p)
      dev.off()
    }
    
    #"ls *.png | sort -r"
    fnames = sort(list.files(file.path(getwd(), subDir, 'gifsd', taxonnm), '*.png'), decreasing=TRUE)
    fnames_abs = paste(file.path(getwd(), subDir, 'gifsd', taxonnm, fnames), collapse=' ')
    sys_call = paste0("convert -delay 100 ", fnames_abs, " ", file.path(getwd(), subDir, 'gifsd'), "/", taxonnm, ".gif")
    system(sys_call)
    
  }
  
  return(p)
}

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), tolower(substring(s, 2)),
        sep="", collapse=" ")
}

plot_data_maps <- function(y, centers, taxa, N, K, thresh, limits, suff, save_plots, fpath=subDir){
  
  if (is.null(taxa)){taxa=seq(1,K)}
  
  props_data = t(apply(y, 1, function(x) x/sum(x)))
  #colnames(props_data) = taxa
  
  prop_dat = data.frame(props=numeric(0), x=integer(0), y=integer(0), taxon=character())
  for (k in 1:K){
    prop_dat = rbind(prop_dat, data.frame(props = props_data[,k], 
                                          x     = centers[,1], 
                                          y     = centers[,2], 
                                          taxon = rep(taxa[k], N)))
  }
  
  if (!is.na(thresh)){
    prop_dat$props[which(prop_dat$props > thresh)] = thresh
  }
    
  p <- ggplot() + geom_tile(data=prop_dat, aes(x=x*1000, y=y*1000, fill=props)) + 
    scale_fill_gradientn(colours=tim.colors(), guide='none') + coord_fixed() + 
    scale_x_continuous(limits$xlims*1000) + scale_y_continuous(limits$ylims*1000)
  p <- add_map_albers(p, map_data=us.fort, limits)
  p <- p + facet_grid(~taxon)
  p <- theme_clean(p) + theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
  print(p)
  Sys.sleep(2)
  if (save_plots){
    fname = paste0(fpath, '/veg_maps_data_', suff, '.pdf')
    ggsave(file=fname, scale=1, width=12)
    sys_str = paste("pdfcrop", fname, fname, sep=' ')
    system(sys_str)
#     dev.off()
  }
   return(p)
}

# plot_data_maps_bw(y_bw, centers_pls, taxa='bw', N_pls, K=1, thresh=0.6, limits, suff='bw', save_plots, fpath=subDir)
plot_data_maps_bw <- function(props_data, centers, taxa, N, K, thresh, limits, suff, save_plots, fpath=subDir){
  
  if (is.null(taxa)){taxa=seq(1,K)}

  prop_dat = data.frame(props = props_data, 
                        x     = centers[,1], 
                        y     = centers[,2], 
                        taxon = rep(taxa[k], N))
  
  if (!is.na(thresh)){
    prop_dat$props[which(prop_dat$props > thresh)] = thresh
  }
  
  p <- ggplot() + geom_tile(data=prop_dat, aes(x=x*1000, y=y*1000, fill=props)) + 
    scale_fill_gradientn(colours=tim.colors(), guide='none') + coord_fixed() + 
    scale_x_continuous(limits$xlims*1000) + scale_y_continuous(limits$ylims*1000)
  p <- add_map_albers(p, map_data=us.fort, limits)
  p <- p + facet_grid(~taxon)
  p <- theme_clean(p) + theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
  print(p)
  Sys.sleep(2)
  if (save_plots){
    fname = paste0(fpath, '/veg_maps_data_', suff, '.pdf')
    ggsave(file=fname, scale=1, width=12)
    sys_str = paste("pdfcrop", fname, fname, sep=' ')
    system(sys_str)
    #     dev.off()
  }
  return(p)
}

# plot_both_maps <- function(r_mean, y, centers, taxa, t, N, K, T, thresh, suff, save_plots, fpath=subDir){
#   
#   if (is.null(taxa)){taxa=seq(1,K)}
#   
#   props_data = t(apply(y, 1, function(x) x/sum(x)))
#   #colnames(props_data) = taxa
#   
#   prop_dat = data.frame(props=numeric(0), x=integer(0), y=integer(0), taxon=character())
#   for (k in 1:K){
#     prop_dat = rbind(prop_dat, data.frame(props = props_data[,k], 
#                                           x     = centers[,1], 
#                                           y     = centers[,2], 
#                                           taxon = rep(taxa[k], N)))
#   }
#   
#   if (!is.na(thresh)){
#     prop_dat$props[which(prop_dat$props > thresh)] = thresh
#   }
#   
#   p1 <- ggplot() + geom_raster(data=prop_dat, aes(x=x, y=y, fill=props)) + 
#     scale_fill_gradientn(colours=tim.colors()) + coord_fixed()
#   p1 <- p1 + facet_grid(~taxon)
#   p1 <- theme_clean(p1) + theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
# #   print(p)
#   
#   
#   prop_preds = data.frame(props=numeric(0), x=integer(0), y=integer(0),time=integer(0), taxon=character())
#   for (k in 1:K){
#     prop_preds = rbind(prop_preds, data.frame(props = r_mean[,k], 
#                                           x     = rep(centers[,1], each=T), 
#                                           y     = rep(centers[,2], each=T), 
#                                           time  = rep(t,times=N), 
#                                           taxon = rep(taxa[k], N*T)))
#   }
#   
#   if (!is.na(thresh)){
#     prop_preds$props[which(prop_preds$props > thresh)] = thresh
#   }
#   
#   p2 <- ggplot() + geom_raster(data=prop_preds, aes(x=x, y=y, fill=props)) + 
#     scale_fill_gradientn(colours=tim.colors()) + coord_fixed()
#   p2 <- p2 + facet_grid(time~taxon)
#   p2 <- theme_clean(p2) + theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
# #   print(p2)
#   
#   grid.arrange(p1, p2, ncol=1)
#   
#   
# #   Sys.sleep(1)
# #   if (save_plots){
# #     ggsave(file=paste('figures/pred_model/veg_maps_', suff, '.pdf', sep=''), scale=1)
# #   }
#   #   print(p)
# }


plot_pred_maps_binned <- function(r_mean, centers, breaks, taxa, taxa_sub, ages, N, K, T, limits, suff, save_plots, fpath=subDir){
  
  rescale=1000000
  
  if (is.null(taxa)){taxa=seq(1,K)}
  
  r_mean_binned = matrix(0, nrow=nrow(r_mean), ncol=ncol(r_mean))
  colnames(r_mean_binned) <- colnames(r_mean)
  
  for (i in 1:ncol(r_mean)){
    r_mean_binned[,i] = cut(r_mean[,i], breaks, include.lowest=TRUE, labels=FALSE)
  }
  
  breaklabels = apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1, 
                      function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })
  
  prop_dat = data.frame(props=numeric(0), x=integer(0), y=integer(0),time=integer(0), taxon=character())
  for (k in 1:K){
    prop_dat = rbind(prop_dat, data.frame(props = r_mean_binned[,k], 
                                          x     = rep(centers[,1], each=T)*rescale, 
                                          y     = rep(centers[,2], each=T)*rescale, 
                                          time  = rep(ages,times=N), 
                                          taxon = rep(taxa[k], N*T)))
  }
  
  prop_dat = prop_dat[which(prop_dat$taxon %in% taxa_sub),]
  
  p <- ggplot() + geom_tile(data=prop_dat, aes(x=x, y=y, fill=factor(props))) + 
    scale_fill_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportions') + 
    coord_fixed() + scale_x_continuous(limits$xlims) + scale_y_continuous(limits$ylims)
  p <- add_map_albers(p, map_data=us.fort, limits)
  p <- p + facet_grid(time~taxon)
  p <- theme_clean(p) + theme(strip.text.x = element_blank(), strip.text.y = element_text(size = rel(1.5)))
    #theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
  print(p)
  Sys.sleep(2)
  if (save_plots){
    ggsave(file=paste(fpath, '/veg_maps_binned_', suff, '.pdf', sep=''), scale=1)
#     dev.off()
  }
  return(p)
}

theme_clean <- function(plot_obj){
  plot_obj <- plot_obj + theme(axis.ticks = element_blank(), 
                               axis.text.y = element_blank(), 
                               axis.text.x = element_blank(),
                               axis.title.x = element_blank(),
                               axis.title.y = element_blank(),
                               plot.background = element_rect(fill = "transparent",colour = NA))
  
  return(plot_obj)
}

# plot_core_locations_select <- function(y, centers_pol, centers_pls, ages, idx.keep, limits, suff, fpath)
plot_core_locations <- function(y, centers_pol, centers_pls, ages, limits, fpath){
  
  centers = centers_pol*1000
  
  K = ncol(y)
  N_cores = nrow(centers_pol)
  
  buffer = 0#100000
  
  # xlo = min(centers[,1]) - buffer
  # xhi = max(centers[,1]) + buffer
  # ylo = min(centers[,2]) - buffer
  # yhi = max(centers[,2]) + buffer
  
  xlo = min(centers_pls[,1])*1000 - buffer
  xhi = max(centers_pls[,1])*1000 + buffer
  ylo = min(centers_pls[,2])*1000 - buffer
  yhi = max(centers_pls[,2])*1000 + buffer
  
  core_dat = data.frame(x=integer(0), y=integer(0), age=character(), core_sites=integer(0))
  core_missing = data.frame(x=integer(0), y=integer(0), age=character(), core_sites=integer(0))
  for (i in 1:length(ages)){
    print(i)
    print((N_cores*(i-1) + 1))
    
    y_sub    = y[(N_cores*(i-1) + 1):(N_cores*i),]
    idx_data = which(rowSums(y_sub) != 0) 
    idx_missing = which(rowSums(y_sub) == 0) 
    
    print(idx_data)
    
    core_dat = rbind(core_dat, data.frame(x     = centers_pol[idx_data,1]*1000, 
                                          y     = centers_pol[idx_data,2]*1000, 
                                          age   = rep(ages[i], length(idx_data))
    ))
    
    core_missing = rbind(core_missing, data.frame(x     = centers_pol[idx_missing,1]*1000, 
                                          y     = centers_pol[idx_missing,2]*1000, 
                                          age   = rep(ages[i], length(idx_missing))
    ))
  }
  
#   col = rep(NA, nrow(core_dat))
#   for (i in 1:nrow(core_dat)) {
#     idx = which( (core_dat[,1] == core_dat[i,1]) & (core_dat[,2] == core_dat[i,2]))
#     core_ages = core_dat[idx, 'age']
#     print(length(idx))
#     if (length(idx) > 1){
#       age_diff = (max(core_ages) - min(core_ages))
#       if (age_diff >= 10) {
#         col[i] = 4
#       } else if (length(idx) <= 4 ) {
#         col[i] = 1
#       } else {
#         col[i] = 2
#       }
#     } else {
#       col[i] = 1
#     }
# #     if (length(idx) >= 10) {
# #       col[i] = 4
# #     } else if (length(idx) <= 4 ) {
# #       col[i] = 1
# #     } else {
# #       col[i] = 2
# #     }
#   }
  
  core_dat = cbind(core_dat, col)
  
  p <- ggplot() + geom_point(data=core_dat, aes(x=x, y=y), colour=orange, shape=19) + 
    coord_fixed() #+ scale_x_continuous(limits$xlims) + scale_y_continuous(limits$ylims)
  p <- add_map_albers(p, map_data=us.fort, limits)#, xlims=c(xlo,xhi), ylims=c(ylo, yhi))
  p <- p + facet_grid(age~.)
  p <- p + theme(strip.text.x = element_blank(),
                 strip.text.y = element_blank())
  p <- p + theme(strip.background = element_blank())
  p <- theme_clean(p) #+ theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
  
  
#   p <- ggplot() + geom_point(data=core_missing, aes(x=x, y=y), colour='#FF6600', shape=19) + 
#     coord_fixed() #+ scale_x_continuous(limits$xlims) + scale_y_continuous(limits$ylims)
#   p <- add_map_albers(p, map_data=us.fort, limits)#, xlims=c(xlo,xhi), ylims=c(ylo, yhi))
#   p <- p + facet_grid(age~.)
#   p <- p + theme(strip.text.x = element_blank(),
#                  strip.text.y = element_blank())
#   p <- p + theme(strip.background = element_blank())
#   p <- theme_clean(p)
#   
  Sys.sleep(2)
  print(p)
  
  if (save_plots){
    ggsave(file=paste(fpath, '/core_locs_', suff, '.pdf', sep=''), scale=1)
    #     dev.off()
  }
  
  return(p)
}

plot_core_locations_sj <- function(y, centers_pol, centers_sj, centers_pls, ages, limits, fpath, save_plots=TRUE){
  
  #centers = centers_pol*1000
  #centers_sj = centers_sj*1000
  rescale = 1000000
  centers_pol = rescale*centers_pol
  
  K = ncol(y)
  N_cores = nrow(centers_pol)
  
  buffer = 0#100000
  
  xlo = min(centers_pls[,1])*rescale - buffer
  xhi = max(centers_pls[,1])*rescale + buffer
  ylo = min(centers_pls[,2])*rescale - buffer
  yhi = max(centers_pls[,2])*rescale + buffer
  
  core_dat = data.frame(x=integer(0), y=integer(0), age=character(), core_sites=integer(0))
  core_missing = data.frame(x=integer(0), y=integer(0), age=character(), core_sites=integer(0))
  for (i in 1:length(ages)){
    print(i)
    print((N_cores*(i-1) + 1))
    
    y_sub    = y[(N_cores*(i-1) + 1):(N_cores*i),]
    idx_data = which(rowSums(y_sub) != 0) 
    idx_missing = which(rowSums(y_sub) == 0) 
    
    print(idx_data)
    
    core_dat = rbind(core_dat, data.frame(x     = centers_pol[idx_data,1]*rescale, 
                                          y     = centers_pol[idx_data,2]*rescale, 
                                          age   = rep(ages[i], length(idx_data))
    ))
    
    core_missing = rbind(core_missing, data.frame(x     = centers_pol[idx_missing,1]*rescale, 
                                                  y     = centers_pol[idx_missing,2]*rescale, 
                                                  age   = rep(ages[i], length(idx_missing))
    ))
  }
  
  # core_dat = cbind(core_dat)
  
  p <- ggplot() + geom_point(data=core_dat, aes(x=x, y=y), colour='orange', shape=19) + 
    coord_fixed() #+ scale_x_continuous(limits$xlims) + scale_y_continuous(limits$ylims)
  p <- p + geom_point(data=centers_sj, aes(x=x, y=y), colour='blue')
  p <- add_map_albers(p, map_data=us.fort, limits)#, xlims=c(xlo,xhi), ylims=c(ylo, yhi))
  p <- p + facet_grid(age~.)
  p <- p + theme(strip.text.x = element_blank(),
                 strip.text.y = element_blank())
  p <- p + theme(strip.background = element_blank())
  p <- theme_clean(p) #+ theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))

  Sys.sleep(2)
  print(p)
  
  if (save_plots){
    ggsave(file=paste(fpath, '/core_locs_sj.pdf', sep=''), scale=1)
  }
  
  return(p)
}


# plot_core_locations_select <- function(y, centers_pol, centers_pls, ages, idx.keep, limits, suff, fpath)
plot_core_intervals <- function(y, centers_pol, centers_pls, ages, limits, fpath){
  
  centers = centers_pol*1000
  
  K = ncol(y)
  N_cores = nrow(centers_pol)
  
  buffer = 0#100000
  
  # xlo = min(centers[,1]) - buffer
  # xhi = max(centers[,1]) + buffer
  # ylo = min(centers[,2]) - buffer
  # yhi = max(centers[,2]) + buffer
  
  xlo = min(centers_pls[,1])*1000 - buffer
  xhi = max(centers_pls[,1])*1000 + buffer
  ylo = min(centers_pls[,2])*1000 - buffer
  yhi = max(centers_pls[,2])*1000 + buffer
  
  core_num = rep(seq(1, N_cores), each=T)
  age      = rep(ages, N_cores)
  
  missing  = rowSums(y) != 0
  core_dat = data.frame(core = missing*1, core_num = core_num, age = age*100)
  core_dat = core_dat[core_dat$core == 1,]
  
  age_max = aggregate(age ~ core_num, core_dat, function(x) max(x))
  age_min = aggregate(age ~ core_num, core_dat, function(x) min(x))
  # tenth = rep(c(rep(0, 11), 1), 12)
  age_seg = data.frame(core_num = age_max[order(age_max$age, decreasing=TRUE),'core_num'], 
                       age_max  = age_max[order(age_max$age, decreasing=TRUE), 'age']+50, 
                       age_min  = age_min[order(age_max$age, decreasing=TRUE),'age']-50, 
                       # tenth,
                       core_lab = seq(1, N_cores)) 
  core_dat$core_lab = age_seg$core_lab[match(core_dat$core_num, age_seg$core_num)]
  age_seg$core_num = factor(age_seg$core_lab, levels=age_seg$core_lab[order(age_seg$age_max, decreasing=TRUE)])
  p <- ggplot(age_seg) + 
    geom_segment(data=age_seg, aes(x=age_min, y=core_lab, xend=age_max, yend=core_lab), colour='grey', size=0.8) + #, colour=factor(tenth))) +
    xlab('Age (YBP)') + ylab('Core number') #+ theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
  p <- p + scale_y_discrete(breaks=seq(0,N_cores, by=20), labels=seq(0,N_cores, by=20)) + #scale_x_reverse(breaks=seq(0, 2000, by=250))
       scale_x_continuous(breaks=seq(0, 2000, by=250)) #+ 
  # p <- theme_clean(p)
  p <- p + geom_point(data=core_dat, aes(x=age, y=factor(core_lab)), shape=3, size=1.2) #+ coord_flip()
#   ggplot(core_dat) + geom_line(data=core_dat, aes(x=age, y=factor(core_num)))
#   
  # core_count = aggregate(. ~ age, core_dat, function(x) length(unique(x)))
  # p <- p + geom_text(data=NULL, aes(x=c(150, 950, 1950), y=c(124,34,31), label=c(124,34,31)))
#   ggplot(core_count) + geom_point(data=core_count, aes(x=age, y=core_num))
  p <- theme_clean(p) + theme_bw() + theme(strip.text = element_text(size = 12),
                                           axis.text = element_text(size=14),
                                           axis.title = element_text(size=14))
  
  Sys.sleep(2)
  print(p)
  
  if (save_plots){
    ggsave(file=paste(fpath, '/core_time_spans.pdf', sep=''), scale=1)
    #     dev.off()
  }
  
  return(p)
}

# plot_core_locations_select <- function(y, centers_pol, centers_pls, ages, idx.keep, limits, suff, fpath)
# plot_core_intervals3(y, meta_pol, centers_pol, centers_pls, ages, limits, fpath=subDir)
plot_core_intervals3 <- function(meta_all, fpath){
  
  N_cores = length(unique(meta_all$stat_id))
  
  age_max = aggregate(age_bacon ~ stat_id, meta_all, function(x) max(x, na.rm=TRUE))
  age_min = aggregate(age_bacon ~ stat_id, meta_all, function(x) min(x, na.rm=TRUE))
  age_seg = data.frame(stat_id = age_max[order(age_max$age_bacon, decreasing=TRUE),'stat_id'], 
                       age_max  = age_max[order(age_max$age_bacon, decreasing=TRUE), 'age_bacon'], 
                       age_min  = age_min[order(age_max$age_bacon, decreasing=TRUE),'age_bacon'], 
                       core_num = seq(1, N_cores)) 
  age_seg$stat_id = factor(age_seg$stat_id, levels=age_seg$stat_id[order(age_seg$age_max, decreasing=TRUE)])
  
  # meta_all = meta_all[match(meta_all$stat_id, age_seg$stat_id),]
  meta_all$stat_id = factor(meta_all$stat_id, levels=levels(age_seg$stat_id))
  meta_all$core_num = age_seg[match(meta_all$stat_id, age_seg$stat_id), 'core_num']
  
  age_seg$age_max[age_seg$age_max>2150] = 2150
  age_seg$age_min[age_seg$age_min<150] = 150
  
  p <- ggplot(age_seg) + 
    geom_vline(xintercept=c(150,500,1000,1500, 2000), colour="lightgrey") +
    geom_segment(data=age_seg, aes(x=age_min, y=core_num, xend=age_max, yend=core_num), colour='grey', size=0.5) + #, colour=factor(tenth))) +
    xlab('YB1950') + ylab('Core') #+ theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
  # p <- theme_clean(p)
  p <- p + geom_point(data=meta_all, aes(x=age_bacon, y=core_num), shape=3, size=1.2) #+ coord_flip()
  p <- p  + 
    scale_x_continuous(breaks=c(150,500,1000,1500,2000), labels=c(150,500,1000,1500,2000), limits=c(150, 2150)) 
  # p <- p + scale_y_discrete(breaks=as.character(seq(0,N_cores, by=20)), labels=as.character(seq(0,N_cores, by=20))) + expand_limits(y=c(-5, 210)) #+ #scale_x_reverse(breaks=seq(0, 2000, by=250))
  
  # p <- theme_clean(p) + theme_bw() + theme(strip.text = element_text(size = 12),
  #                                          axis.text = element_text(size=14),
  #                                          axis.title = element_text(size=14))
  
  # like this one
  p <- p + theme_light() + theme(strip.text = element_text(size = 14),
                                           axis.text = element_text(size=14),
                                           axis.title = element_text(size=14))
  
  Sys.sleep(2)
  print(p)
  
  if (save_plots){
    ggsave(file=paste(fpath, '/core_time_spans_bacon.pdf', sep=''), scale=1)
    #     dev.off()
  }
  
  return(p)
}
plot_core_intervals2 <- function(y, centers_pol, centers_pls, ages, limits, fpath){
  
  centers = centers_pol*1000
  
  K = ncol(y)
  N_cores = nrow(centers_pol)
  
  buffer = 0#100000

  xlo = min(centers_pls[,1])*1000 - buffer
  xhi = max(centers_pls[,1])*1000 + buffer
  ylo = min(centers_pls[,2])*1000 - buffer
  yhi = max(centers_pls[,2])*1000 + buffer
  
  core_num = rep(seq(1, N_cores), each=T)
  age      = rep(ages, N_cores)
  
  missing  = rowSums(y) != 0
  core_dat = data.frame(core = missing*1, core_num = core_num, age = age*100)
  core_dat = core_dat[core_dat$core == 1,]
  
  age_max = aggregate(age ~ core_num, core_dat, function(x) max(x))
  age_min = aggregate(age ~ core_num, core_dat, function(x) min(x))
  # tenth = rep(c(rep(0, 11), 1), 12)
  age_seg = data.frame(core_num = age_max[order(age_max$age, decreasing=TRUE),'core_num'], 
                       age_max  = age_max[order(age_max$age, decreasing=TRUE), 'age']+50, 
                       age_min  = age_min[order(age_max$age, decreasing=TRUE),'age']-50, 
                       # tenth,
                       core_lab = seq(1, N_cores)) 
  core_dat$core_lab = age_seg$core_lab[match(core_dat$core_num, age_seg$core_num)]
  age_seg$core_num = factor(age_seg$core_lab, levels=age_seg$core_lab[order(age_seg$age_max, decreasing=TRUE)])
  p <- ggplot(age_seg) + 
    # geom_segment(data=age_seg, aes(x=age_min, y=core_lab, xend=age_max, yend=core_lab), colour='grey', size=1) + #, colour=factor(tenth))) +
    geom_segment(data=age_seg, aes(x=age_min, y=core_lab, xend=age_max, yend=core_lab), colour='grey', size=1) + #, colour=factor(tenth))) +
    xlab('Age (YBP)') + ylab('Pollen record') #+ arrow()
  p <- p + scale_y_discrete(breaks=seq(0,N_cores, by=20), labels=seq(0,N_cores, by=20))
  # p <- theme_clean(p)
  p <- p + geom_point(data=core_dat, aes(x=age, y=factor(core_lab)), shape=1, size=1.2)
  #   ggplot(core_dat) + geom_line(data=core_dat, aes(x=age, y=factor(core_num)))
  #   
  # core_count = aggregate(. ~ age, core_dat, function(x) length(unique(x)))
  # p <- p + geom_text(data=NULL, aes(x=c(150, 950, 1950), y=c(124,34,31), label=c(124,34,31)))
  #   ggplot(core_count) + geom_point(data=core_count, aes(x=age, y=core_num))
  p <- p + theme(axis.text=element_text(size=12),
                 axis.title=element_text(size=14,face="bold"))
  
  Sys.sleep(2)
  print(p)
  
  # if (save_plots){
  #   ggsave(file=paste(fpath, '/core_time_spans.pdf', sep=''), scale=1)
  # }
  # 
  return(p)
}


# plot_core_locations_select(y, centers_pol, centers_pls, ages, idx.keep, limits, suff=suff_figs, fpat=subDir)

plot_core_locations_select <- function(y, centers_pol, centers_pls, ages, idx.keep, limits, suff, fpath){
 
  rescale = 1000000
  centers = centers_pol*rescale

  K = ncol(y)
  N_cores = nrow(centers_pol)
  
  buffer = 0#100000

  # xlo = min(centers[,1]) - buffer
  # xhi = max(centers[,1]) + buffer
  # ylo = min(centers[,2]) - buffer
  # yhi = max(centers[,2]) + buffer

  xlo = min(centers_pls[,1])*rescale - buffer
  xhi = max(centers_pls[,1])*rescale + buffer
  ylo = min(centers_pls[,2])*rescale - buffer
  yhi = max(centers_pls[,2])*rescale + buffer
  
#   idx.keep  = c(1,2,length(ages)/2,T)
  ages.keep = ages[idx.keep]
  T.keep    = length(ages.keep)
  
  idx_y = vector(length=0)
  for (i in 1:T.keep){
    idx_orig   = seq(idx.keep[i], N_cores*T, by=T)
    idx_y = c(idx_y, idx_orig)
  }
  
  # idx_y = sort(idx_y)
  y_keep = y[idx_y,]
  
  ncores = rep(NA, T.keep)
  core_dat = data.frame(x=integer(0), y=integer(0), age=character(), core_sites=integer(0))
  for (i in 1:length(ages.keep)){
    print(i)
    print((N_cores*(i-1) + 1))
    
    y_sub    = y_keep[(N_cores*(i-1) + 1):(N_cores*i),]
    idx_data = which(rowSums(y_sub) != 0) 
    
    print(idx_data)
    
    core_dat = rbind(core_dat, data.frame(x     = centers_pol[idx_data,1]*rescale, 
                                          y     = centers_pol[idx_data,2]*rescale, 
                                          age   = rep(ages.keep[i], length(idx_data))
                                          ))
    ncores[i] = length(idx_data)
  }
  
  col = rep(NA, nrow(core_dat))
  for (i in 1:nrow(core_dat)) {
    idx = which( (core_dat[,1] == core_dat[i,1]) & (core_dat[,2] == core_dat[i,2]))
    if (length(idx) == 4) {
      col[i] = 4
    } else if (length(idx) == 1) {
      col[i] = 1
    } else {
      col[i] = 2
    }
  }
  
  core_dat = cbind(core_dat, col)
  
  core_dat$age = core_dat$age*100
  
#   p <- ggplot() + geom_point(data=core_dat, aes(x=x, y=y, colour=factor(col)), shape=19) + #, colour='#FF6600', shape=19) + 
#     coord_fixed() #+ scale_x_continuous(limits$xlims) + scale_y_continuous(limits$ylims)
#   p <- add_map_albers(p, map_data=us.fort, limits)#, xlims=c(xlo,xhi), ylims=c(ylo, yhi))
#   p <- p + facet_grid(age~.)
# #   p <- p + theme(strip.text.x = element_blank(),
# #                 strip.text.y = element_blank())
# #   p <- p + theme(strip.background = element_blank())
#   p <- theme_clean(p) + theme(strip.text.y = element_text(size = rel(1.5)), 
#                           strip.text.x = element_text(size = rel(1.5)))
  
  print(head(core_dat))
  
  p <- ggplot() + geom_point(data=core_dat, aes(x=x, y=y), colour='#FF6600', shape=19, size=1.1) + #, colour='#FF6600', shape=19) + 
    coord_fixed() #+ scale_x_continuous(limits$xlims) + scale_y_continuous(limits$ylims)
  p <- add_map_albers(p, map_data=us.fort, limits)#, xlims=c(xlo,xhi), ylims=c(ylo, yhi))
  p <- p + facet_grid(.~age)

  p <- theme_clean(p) + theme_bw() + theme(strip.text = element_text(size = 14),
                                           strip.background = element_rect(colour = 'grey'),
                                           axis.text = element_blank(),
                                           axis.title = element_blank(),
                                           axis.ticks = element_blank(),
                                           panel.grid.major = element_blank(),
                                           panel.grid.minor = element_blank())
  topright=c(0.9,0.9)
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
  
  
  Sys.sleep(2)
  print(p)
  
  if (save_plots){
    fname = paste0(fpath, '/core_locs_', suff, '.pdf')	
    ggsave(file=fname, scale=1)
    sys_str = paste("pdfcrop", fname, fname, sep=' ')
    system(sys_str)
  }
  
  return(p)
}

plot_core_locations_select_sj <- function(y, centers_pol, centers_sj, centers_pls, ages, idx.keep, limits, fpath, save_plots=TRUE){
  
  rescale = 1000000
  #centers = centers_pol*rescale
  # centers_sj = centers_sj
  
  K = ncol(y)
  N_cores = nrow(centers_pol)
  
  buffer = 0#100000

  xlo = min(centers_pls[,1])*rescale - buffer
  xhi = max(centers_pls[,1])*rescale + buffer
  ylo = min(centers_pls[,2])*rescale - buffer
  yhi = max(centers_pls[,2])*rescale + buffer
  
  #   idx.keep  = c(1,2,length(ages)/2,T)
  ages.keep = ages[idx.keep]
  T.keep    = length(ages.keep)
  
  idx_y = vector(length=0)
  for (i in 1:T.keep){
    idx_orig   = seq(idx.keep[i], N_cores*T, by=T)
    idx_y = c(idx_y, idx_orig)
  }
  
  # idx_y = sort(idx_y)
  y_keep = y[idx_y,]
  
  ncores = rep(NA, T.keep)
  core_dat = data.frame(x=integer(0), y=integer(0), age=character(), core_sites=integer(0))
  for (i in 1:length(ages.keep)){
    print(i)
    print((N_cores*(i-1) + 1))
    
    y_sub    = y_keep[(N_cores*(i-1) + 1):(N_cores*i),]
    idx_data = which(rowSums(y_sub) != 0) 
    
    print(idx_data)
    
    core_dat = rbind(core_dat, data.frame(x     = centers_pol[idx_data,1]*rescale, 
                                          y     = centers_pol[idx_data,2]*rescale, 
                                          age   = rep(ages.keep[i], length(idx_data))
    ))
    ncores[i] = length(idx_data)
  }
  
  col = rep(NA, nrow(core_dat))
  for (i in 1:nrow(core_dat)) {
    idx = which( (core_dat[,1] == core_dat[i,1]) & (core_dat[,2] == core_dat[i,2]))
    if (length(idx) == 4) {
      col[i] = 4
    } else if (length(idx) == 1) {
      col[i] = 1
    } else {
      col[i] = 2
    }
  }
  
  core_dat = cbind(core_dat, col)
  
  p <- ggplot() + geom_point(data=core_dat, aes(x=x, y=y), colour='#FF6600', shape=19) + #, colour='#FF6600', shape=19) + 
    coord_fixed()
  p <- p + geom_point(data=centers_sj, aes(x=x*rescale, y=y*rescale), colour='blue') +
  geom_text(data=centers_sj, aes(label=handle,x=x*rescale, y=y*rescale),hjust=0, vjust=0)
  p <- add_map_albers(p, map_data=us.fort, limits)
  p <- p + facet_grid(age~.)
  p <- theme_clean(p) + theme(strip.text.y = element_text(size = rel(1.5)), 
                              strip.text.x = element_text(size = rel(1.5)))
  
  Sys.sleep(2)
  print(p)
  
  if (save_plots){
    fname = paste0(fpath, '/core_locs_sj.pdf')	
    ggsave(file=fname, scale=1)
    sys_str = paste("pdfcrop", fname, fname, sep=' ')
    system(sys_str)
  }
  
  return(p)
}


# 
# plot_core_locations_select <- function(y, centers_pol, centers_pls, ages, idx.keep, limits, suff, fpath){
#   
#   rescale = 1000000
#   centers = centers_pol*rescale
#   
#   K = ncol(y)
#   N_cores = nrow(centers_pol)
#   
#   buffer = 0#100000
#   
#   # xlo = min(centers[,1]) - buffer
#   # xhi = max(centers[,1]) + buffer
#   # ylo = min(centers[,2]) - buffer
#   # yhi = max(centers[,2]) + buffer
#   
#   xlo = min(centers_pls[,1])*rescale - buffer
#   xhi = max(centers_pls[,1])*rescale + buffer
#   ylo = min(centers_pls[,2])*rescale - buffer
#   yhi = max(centers_pls[,2])*rescale + buffer
#   
#   #   idx.keep  = c(1,2,length(ages)/2,T)
#   ages.keep = ages[idx.keep]
#   T.keep    = length(ages.keep)
#   
#   idx_y = vector(length=0)
#   for (i in 1:T.keep){
#     idx_orig   = seq(idx.keep[i], N_cores*T, by=T)
#     idx_y = c(idx_y, idx_orig)
#   }
#   
#   idx_y = sort(idx_y)
#   y_keep = y[idx_y,]
#   
#   core_dat = data.frame(x=integer(0), y=integer(0), age=character(), core_sites=integer(0))
#   for (i in 1:length(ages.keep)){
#     print(i)
#     print((N_cores*(i-1) + 1))
#     
#     y_sub    = y_keep[(N_cores*(i-1) + 1):(N_cores*i),]
#     idx_data = which(rowSums(y_sub) != 0) 
#     
#     print(idx_data)
#     
#     core_dat = rbind(core_dat, data.frame(x     = centers_pol[idx_data,1]*rescale, 
#                                           y     = centers_pol[idx_data,2]*rescale, 
#                                           age   = rep(ages.keep[i], length(idx_data))
#     ))
#   }
#   
#   
#   p <- ggplot() + geom_point(data=core_dat, aes(x=x, y=y), colour='#FF6600', shape=19) + 
#     coord_fixed() #+ scale_x_continuous(limits$xlims) + scale_y_continuous(limits$ylims)
#   p <- add_map_albers(p, map_data=us.fort, limits)#, xlims=c(xlo,xhi), ylims=c(ylo, yhi))
#   p <- p + facet_grid(age~.)
#   #   p <- p + theme(strip.text.x = element_blank(),
#   #                 strip.text.y = element_blank())
#   #   p <- p + theme(strip.background = element_blank())
#   p <- theme_clean(p) + theme(strip.text.y = element_text(size = rel(1.5)), 
#                               strip.text.x = element_text(size = rel(1.5)))
#   
#   Sys.sleep(2)
#   print(p)
#   
#   if (save_plots){
#     fname = paste0(fpath, '/core_locs_', suff, '.pdf')	
#     ggsave(file=fname, scale=1)
#     sys_str = paste("pdfcrop", fname, fname, sep=' ')
#     system(sys_str)
#   }
#   
#   return(p)
# }


add_map_albers <- function(plot_obj, map_data=us.fort, limits){
  p <- plot_obj + geom_path(data=map_data, aes(x=long, y=lat, group=group),  colour='grey55') + 
   scale_x_continuous(limits = limits$xlims*1000000) +
   scale_y_continuous(limits = limits$ylims*1000000) #+ coord_map("albers")
  return(p)
}




poster_fig <- function(y, y_veg, r_mean, centers_veg, centers_pls, centers_polU, taxa, t, N, K, T, thresh, limits, type, suff, save_plots, fpath=subDir){
  
  if (is.null(taxa)){taxa=seq(1,K)}
  
  props_data = t(apply(y_veg, 1, function(x) x/sum(x)))
  
  K = ncol(y)
  N_cores = nrow(centers_polU)
  
  prop_dat = data.frame(props=numeric(0), x=integer(0), y=integer(0), time=character(0), taxon=character())
#   core_dat = data.frame(x=integer(0), y=integer(0), age=character(), core_sites=integer(0))
  for (i in 1:length(ages)){
    print(i)
    print((N_cores*(i-1) + 1))
    
    y_sub    = y[(N_cores*(i-1) + 1):(N_cores*i),]
    idx_data = which(rowSums(y_sub) != 0) 
    
    print(idx_data)
    
    prop_dat = rbind(prop_dat, data.frame(props = rep(0, length(idx_data)),
                                          x     = centers_polU[idx_data,1]*1000, 
                                          y     = centers_polU[idx_data,2]*1000, 
                                          time  = rep(ages[i], length(idx_data)),
                                          taxon = rep('cores', length(idx_data))))

  }
  
  core_dat = data.frame(x=integer(0), y=integer(0), time=character(0), taxon=character())
  #   core_dat = data.frame(x=integer(0), y=integer(0), age=character(), core_sites=integer(0))
  for (i in 1:length(ages)){
    print(i)
    print((N_cores*(i-1) + 1))
    
    y_sub    = y[(N_cores*(i-1) + 1):(N_cores*i),]
    idx_data = which(rowSums(y_sub) != 0) 
    
    print(idx_data)
    
    core_dat = rbind(core_dat, data.frame(
                                          x     = centers_polU[idx_data,1]*1000, 
                                          y     = centers_polU[idx_data,2]*1000, 
                                          time  = rep(ages[i], length(idx_data)),
                                          taxon = rep('cores', length(idx_data))))
    
  }
  #colnames(props_data) = taxa
  
#   prop_dat = data.frame(props=numeric(0), x=integer(0), y=integer(0), time=character(0), taxon=character())
  for (k in 1:K){
    prop_dat = rbind(prop_dat, data.frame(props = props_data[,k], 
                                          x     = centers_pls[,1]*1000, 
                                          y     = centers_pls[,2]*1000,
                                          time  = rep('pls', times = nrow(centers_pls)),
                                          taxon = rep(taxa[k], nrow(centers_pls))))
  }
  

  
#   prop_dat = data.frame(props=numeric(0), x=integer(0), y=integer(0),time=integer(0), taxon=character())
  for (k in 1:K){
    prop_dat = rbind(prop_dat, data.frame(props = r_mean[,k], 
                                          x     = rep(centers_veg[,1]*1000, each=T), 
                                          y     = rep(centers_veg[,2]*1000, each=T), 
                                          time  = rep(as.character(ages),times=N), 
                                          taxon = rep(taxa[k], N*T)))
  }
  
  if (!is.na(thresh)){
    prop_dat$props[which(prop_dat$props > thresh)] = thresh
  }
  
  p <- ggplot() + geom_tile(data=prop_dat, aes(x=x, y=y, fill=props)) + 
    scale_fill_gradientn(colours=tim.colors()) + coord_fixed() +
     scale_x_continuous(limits$xlims*1000) + scale_y_continuous(limits$ylims*1000)
  p <- add_map_albers(p, us.shp, limits)
  p <- p + facet_grid(time~taxon)
  p <- theme_clean(p) + theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
  print(p)
  
  p <- p + geom_point(data = core_dat, aes(x=x, y=y))
  
  print(p)
  #ggsave(file='figures/pred/pred_plot_test.pdf', scale=1)
  Sys.sleep(2)
  if (save_plots){
    ggsave(file=paste(fpath, '/veg_maps_', suff, '.pdf', sep=''), scale=1)
    #     dev.off()
  }
  return(p)
}



# plot(centers_veg[,1]*1000, centers_veg[,2]*1000, col='blue', pch=19)
# points(centers_pls[,1]*1000, centers_pls[,2]*1000)
# 
# us.shp <- readShapeLines('r/data/map_data/us_alb.shp',
#                          proj4string=CRS('+init=epsg:3175'))
# plot(us.shp, add=T, lwd=2)


plot_data_maps <- function(y, centers, taxa, t, N, K, T, thresh, limits, suff, save_plots, fpath=subDir){
  
  rescale=1000000
  
  if (is.null(taxa)){taxa=seq(1,K)}
  
  props_data = t(apply(y, 1, function(x) x/sum(x)))
  #colnames(props_data) = taxa
  
  prop_dat = data.frame(props=numeric(0), x=integer(0), y=integer(0), taxon=character())
  for (k in 1:K){
    prop_dat = rbind(prop_dat, data.frame(props = props_data[,k], 
                                          x     = centers[,1]*rescale, 
                                          y     = centers[,2]*rescale, 
                                          taxon = rep(taxa[k], N)))
  }
  
  if (!is.na(thresh)){
    prop_dat$props[which(prop_dat$props > thresh)] = thresh
  }
  
  prop_dat$type = rep('PLS', nrow(prop_dat))
  
  p <- ggplot() + geom_tile(data=prop_dat, aes(x=x, y=y, fill=props)) + 
    scale_fill_gradientn(colours=tim.colors(), guide='none') + coord_fixed() + 
    scale_x_continuous(limits$xlims) + scale_y_continuous(limits$ylims)
  p <- add_map_albers(p, us.shp, limits)
  p <- p + facet_grid(type~taxon)
  p <- theme_clean(p) #+ theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
  
  p <- p + theme(strip.text.x = element_blank(),
                 strip.text.y = element_blank())
  p <- p + theme(strip.background = element_blank())
  
  print(p)
  Sys.sleep(2)
  if (save_plots){
    ggsave(file=paste(fpath, '/veg_maps_', suff, '.pdf', sep=''), scale=1)
    ggsave(file=paste(fpath, '/veg_maps_', suff, '.eps', sep=''), scale=1)
    #     dev.off()
  }
  return(p)
}

plot_data_maps_binned_bw <- function(r, centers, taxa, N, K, breaks, limits, suff, save_plots, fpath=subDir){
  
  rescale=1e6
  props_data_binned = matrix(0, nrow=length(r), ncol=1)
  colnames(props_data_binned) <- taxa
  

  props_data_binned = cut(r, breaks, include.lowest=TRUE, labels=FALSE)
  
  breaklabels = apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1, 
                      function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })
    
  prop_dat = data.frame(props = props_data_binned, 
                        x     = centers[,1]*rescale, 
                        y     = centers[,2]*rescale, 
                        taxon = rep(taxa[k], N))
  
  p <- ggplot() + geom_tile(data=prop_dat, aes(x=x, y=y, fill=factor(props)), alpha=0.8) + 
    scale_fill_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportion', drop=FALSE) + 
    coord_fixed() + 
    scale_x_continuous(limits$xlims*rescale) + scale_y_continuous(limits$ylims*rescale)
  p <- add_map_albers(p, us.shp, limits)
  p <- p + geom_path(data=bw, aes(x=long, y=lat),  colour='grey31', alpha=0.7,size=1.8)
  p <- p + theme_bw()
  p <- theme_clean(p) + 
    guides(fill = guide_legend(reverse=T)) + 
    theme_bw() + theme(strip.text = element_text(size = 14),
                                           strip.background = element_rect(colour = 'grey'),
                                           axis.title.y = element_blank(),
                                           axis.text = element_blank(),
                                           axis.title = element_blank(),
                                           axis.ticks = element_blank(),
                                           panel.grid.major = element_blank(), 
                                           panel.grid.minor = element_blank())
  topright=c(0.9,0.9)
  #topright=c(0.4,0.08)
  p <- p + ggsn::scalebar(dist=100, 
                          transform=FALSE, 
                          dist_unit="km",
                          location="topright",
                          height=0.03,
                          st.dist=0.03,
                          st.size=5,
                          border.size=0.4,
                          anchor=c(x=(1-topright[1])*limits$xlims[1]*1000000+topright[1]*limits$xlims[2]*1000000,
                                   y=(1-topright[2])*limits$ylims[1]*1000000+topright[2]*limits$ylims[2]*1000000),
                          x.min=limits$xlims[1]*1000000,
                          x.max=limits$xlims[2]*1000000,
                          y.min=limits$ylims[1]*1000000,
                          y.max=limits$ylims[2]*1000000)
  
  print(p)
  Sys.sleep(2)
  if (save_plots){
    fname = paste0(fpath, '/veg_maps_data_binned_', suff, '.pdf')	
    ggsave(file=fname, scale=1, width=12)
    sys_str = paste("pdfcrop", fname, fname, sep=' ')
    system(sys_str)
  }
  return(p)
}

plot_data_maps_binned <- function(y, centers, taxa, N, K, T, breaks, limits, suff, save_plots, fpath=subDir){
  
  rescale=1000000
  
  if (is.null(taxa)){taxa=seq(1,K)}
  
  props_data = t(apply(y, 1, function(x) if (sum(x) > 0) {x/sum(x)} else {x}))
  colnames(props_data) = taxa
  
  props_data_binned = matrix(0, nrow=nrow(props_data), ncol=ncol(props_data))
  colnames(props_data_binned) <- colnames(props_data)
  
  for (i in 1:ncol(props_data)){
    props_data_binned[,i] = cut(props_data[,i], breaks, include.lowest=TRUE, labels=FALSE)
  }
  
  breaklabels = apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1, 
                      function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })
  
  prop_dat = data.frame(props=numeric(0), x=integer(0), y=integer(0), taxon=character())
  for (k in 1:K){
    prop_dat = rbind(prop_dat, data.frame(props = props_data_binned[,k], 
                                          x     = centers[,1]*rescale, 
                                          y     = centers[,2]*rescale, 
                                          taxon = rep(taxa[k], N)))
  }
  
  prop_dat$type = rep('PLS', nrow(prop_dat))
  
  p <- ggplot() + geom_tile(data=prop_dat, aes(x=x, y=y, fill=factor(props))) + 
    scale_fill_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportions') + 
    coord_fixed() + 
    scale_x_continuous(limits$xlims) + scale_y_continuous(limits$ylims)
  p <- add_map_albers(p, us.shp, limits)
  p <- p + facet_grid(type~taxon)
  p <- theme_clean(p) #+ theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
  
  p <- p + theme(strip.text.x = element_blank(),
                 strip.text.y = element_blank())
  p <- p + theme(strip.background = element_blank())
  
  print(p)
  Sys.sleep(2)
  if (save_plots){
    fname = paste0(fpath, '/veg_maps_data_binned_', suff, '.pdf')	
    ggsave(file=fname, scale=1, width=12)
    sys_str = paste("pdfcrop", fname, fname, sep=' ')
    system(sys_str)
  }
  return(p)
}