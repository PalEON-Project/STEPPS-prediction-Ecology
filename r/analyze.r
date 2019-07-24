library(ggplot2)
library(fields)
library(reshape2)
library(RColorBrewer)
library(maptools)
library(plyr)
library(rgdal)

source('r/utils/pred_analyze_funs.r')
source('r/utils/pred_helper_funs.r')
source('r/utils/pred_plot_funs.r')

PL_var = list(suff_run = '120knots_150to2150ybp_PL_umw_3by_v2.3_ar_set1',
              suff_figs = 'v2.3')

run = PL_var
run_num  = '1'
subDir <- paste0('runs/figures')

rescale=1e6
save_plots = TRUE
breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)

#######################################################################################################################
## load the data
#######################################################################################################################

# comp data
load("data/calibration/cal_data_12taxa_mid_comp_ALL_v0.3.rdata")
centers_comp=centers_veg

# pred run data
load(paste0('runs/input.rdata'))

rIts = readRDS(file=paste0('runs/rIts_sub.RDS'))

# work with something smaller until production
rIts = rIts[,,,seq(1, dim(rIts)[4], by=2)]

r_quants = get_r_quants(rIts)
r_mean2  = r_quants[2,,]
r_sd2    = r_quants[4,,]

if (!file.exists(subDir)){
  dir.create(file.path(subDir))
} 

limits   = get_limits(centers_pls)
taxa_sub = taxa[c(2,5,7,10,11)]
#######################################################################################################################
## Figure 1
#######################################################################################################################

load('runs/input.rdata')
pollen_ts = readRDS(file='data/pollen_ts.RDS')

pollen_ts = pollen_ts[pollen_ts$id %in% meta_pol$id,]
pollen_ts$stat_id = meta_pol$stat_id[match(pollen_ts$id, meta_pol$id)]

states_pol = c('minnesota', 'wisconsin', 'michigan:north')
pollen_ts_new = pollen_ts[which(pollen_ts$state %in% states_pol),]
plot_core_intervals3(pollen_ts_new, fpath=subDir)
plot_core_locations_select(y, centers_pol, centers_pls, ages, idx.keep=c(2,5,10,15,20), limits, suff='', fpath=subDir)

#######################################################################################################################
## Figure 2
#######################################################################################################################
t_diffs = c(2, 5, 10, 15, 20)
es      = eco_sig(rIts, t_diffs, eco_cut=0.03)
diffs   = post_diffs(rIts, t_diffs)

post_mag    = apply(diffs, c(1,2,3,4), mean, na.rm=TRUE)
post_mag_es = post_mag
post_mag_es[which(!es)] = NA

post_sd    = apply(diffs, c(1,2,3,4), sd, na.rm=TRUE)
post_sd_es = post_sd
post_sd_es[which(!es)] = NA

post_prob    = apply(diffs, c(1,2,3,4), diff_q, prob=0.85)
post_prob_es = post_prob
post_prob_es[which(!es)] = NA

counts_bool = rowSums(y) > 0 
pol_samples = data.frame(counts_bool=counts_bool, 
                         age=rep(ages, nrow(centers_pol))*100,
                         x = rep(centers_pol$x, each=T)*rescale,
                         y = rep(centers_pol$y, each=T)*rescale)

pol_samples = pol_samples[which(pol_samples$counts_bool==TRUE),]
pol_samples$t1 = pol_samples$age
pol_samples = pol_samples[which(pol_samples$t1 %in% (ages[t_diffs]*100)),]
pol_samples$t2 = ages[t_diffs][match(pol_samples$t1, ages[t_diffs]*100)+1]*100
pol_samples = pol_samples[which(!is.na(pol_samples$t2)),]

cutoff = 80 
cuts   = c(20, 50, 100, 150, 200) 
dists  = rdist(centers_veg*rescale, pol_samples[,c('x', 'y')])/1e3

## LEFT PANEL
## Volcano plots for paper: Beech and Pine
taxa_keep = taxa[c(2, 10)]

mag_sd_all = melt(post_mag_es)
colnames(mag_sd_all) = c('cell', 'taxon', 't1', 't2', 'mag')
mag_sd_all$sd = melt(post_sd_es)$value
mag_sd_all$taxon = taxa[mag_sd_all$taxon]
mag_sd_all$sig = !is.na(melt(post_prob_es)$value)
mag_sd_all$sig = factor(mag_sd_all$sig, levels=c('FALSE', 'TRUE'))
mag_sd_all$mag = abs(mag_sd_all$mag)
mag_sd_subset = mag_sd_all[which(mag_sd_all$taxon %in% taxa_keep),]
mag_sd_subset$t1 = ages[t_diffs[mag_sd_subset$t1]]*100
mag_sd_subset$t2 = ages[t_diffs[mag_sd_subset$t2]]*100
mag_sd_subset$x = centers_veg$x*1e6
mag_sd_subset$y = centers_veg$y*1e6
mag_sd_subset = mag_sd_subset[which((mag_sd_subset$t1==300)&(mag_sd_subset$t2==2100)),]

# get n sites
counts_bool = rowSums(y) > 0 
pol_samples = data.frame(counts_bool=counts_bool, 
                         age=rep(ages, nrow(centers_pol))*100,
                         x = rep(centers_pol$x, each=T)*rescale,
                         y = rep(centers_pol$y, each=T)*rescale,
                         core = rep(seq(1, nrow(centers_pol)), each=T))

pol_samples = pol_samples[which(pol_samples$counts_bool==TRUE),]
# pol_samples = pol_samples[which(pol_samples$age %in% (ages[t_diffs]*100)),]

pol_samples$extent = NA
for (core in 1:N_cores){
  modern=any(c(200, 300, 400, 500) %in% pol_samples[which(pol_samples$core == core),'age'])
  historic=any(c(1900, 2000, 2100) %in% pol_samples[which(pol_samples$core == core),'age'])
  if (modern&historic) {
    pol_samples$extent[which(pol_samples$core == core)] = 'full'
  } else if (any((ages*100) %in% pol_samples[which(pol_samples$core == core),'age'])) {
    pol_samples$extent[which(pol_samples$core == core)] = 'partial'
  } else {
    pol_samples$extent[which(pol_samples$core == core)] = 'none'
  }
}

pol_extent = pol_samples[!duplicated(pol_samples[,'core']),]
pol_extent = pol_extent[which(pol_extent$extent %in% c('full', 'partial')),]
pol_extent$extent <- factor(pol_extent$extent)
# pol_extent$extent <- factor(pol_extent$extent, levels=c('full', 'partial'))


# add number of sites within radius
radius = 50 

dists  = rdist(centers_veg, centers_pol)
nsite_close = apply(dists, 1, function(x) sum(x<radius/1e3))

mag_sd_subset$nsites = nsite_close[mag_sd_subset$cell]
mag_sd_subset$nsites_bool = mag_sd_subset$nsites > 1

p1 <- ggplot(data=mag_sd_subset) + geom_point(data=mag_sd_subset, aes(x=abs(mag), y=sd, colour=sig, shape=20), size=2) + 
  scale_shape_identity() + coord_fixed() + 
  xlab('Mean magnitude proportional change') + ylab('Standard deviation of proportional change') + 
  facet_grid(taxon~.) 
print(p1)
# p1_record <- recordPlot()
ggsave(paste0(subDir, '/figure2_mag_vs_sd.pdf'))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(2)

## PLOT FIGURE 2: LEFT PANEL
p1 <- ggplot(data=mag_sd_subset) + 
  geom_point(data=mag_sd_subset, aes(x=abs(mag), y=sd, shape=nsites_bool, colour=sig), alpha=0.8, size=2.5) + 
  # coord_fixed() + 
  scale_color_manual(name="Sig.", values=cols)+
  scale_shape_manual(values = c(20, 4), name="Data-rich") +
  xlab('Mean magnitude proportional change') + 
  ylab('Standard deviation of proportional change') + 
  facet_grid(taxon~.) + theme_bw()
print(p1)
ggsave(paste0(subDir, '/figure2_mag_vs_sd_nsite.pdf'), width=13, units="cm")

# MIDDLE PANEL
# plot mag 
bound = ceiling(10*max(abs(c(max(mag_sd_subset$mag, na.rm=TRUE), min(mag_sd_subset$mag, na.rm=TRUE)))))/10
values=c(0, 0.05, 0.1, 0.25, bound)
lims = c(0,bound)
cols=brewer.pal(11, "BrBG")
cols=rev(terrain.colors(20))
legend_label = 'Magnitude'

mag_sd_subset_sig = mag_sd_subset[which(mag_sd_subset$sig==TRUE), ]
# mag_sd_subset = factor(mag_sd_subset)

p2 <- ggplot() + geom_tile(data=mag_sd_subset, aes(x=x, y=y, fill=mag), alpha=0.9)
p2 <- p2 + geom_tile(data=mag_sd_subset_sig, aes(x=x, y=y, group=taxon), alpha=0, color='black', size=0.8) 
p2 <- p2 + geom_tile(data=mag_sd_subset_sig, aes(x=x, y=y, group=taxon, fill=mag)) 
p2 <- p2 + scale_fill_gradientn(colours=cols,
                                values=values, 
                                limits=lims,
                                na.value="white",
                                rescaler = function(x, ...) x,
                                #                               oob = identity, 
                                name=legend_label)
p2 <- p2 + geom_point(data=pol_extent, aes(x=x, y=y, shape=extent), colour="#252525", alpha=0.9, size=1.5) +
  scale_shape_manual(values = c(1, 19), name="Coverage")
p2 <- p2 + coord_fixed()
p2 <- add_map_albers(p2, map_data=us.fort, limits)
p2 <- p2 + facet_grid(taxon~.)
p2 <- theme_clean(p2) + theme_bw() + theme(strip.text = element_text(size = 14),
                                           strip.background = element_rect(colour = 'grey'),
                                           axis.title.y = element_blank(),
                                           axis.text = element_blank(),
                                           axis.title = element_blank(),
                                           axis.ticks = element_blank(),
                                           panel.grid.major = element_blank(), 
                                           panel.grid.minor = element_blank())
topright=c(0.9,0.9)
p2 <- p2 + ggsn::scalebar(dist=100, 
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
print(p2)
ggsave(paste0(subDir, '/figure2_mag_map.pdf'))

# RIGHT PANEL
### plot sd
bound = ceiling(10*max(abs(c(max(mag_sd_subset$sd, na.rm=TRUE), min(mag_sd_subset$sd, na.rm=TRUE)))))/10
values=c(0, 0.05, 0.15, 0.2, bound)
labels=format(breaks)
lims = c(0,bound)
cols=brewer.pal(9, "OrRd")
legend_label = 'SD'

p3 <- ggplot() + geom_tile(data=mag_sd_subset, aes(x=x, y=y, fill=sd), alpha=0.9) 
p3 <- p3 + geom_tile(data=mag_sd_subset_sig, aes(x=x, y=y, group=taxon), alpha=0, color='black', size=0.8) 
p3 <- p3 + geom_tile(data=mag_sd_subset_sig, aes(x=x, y=y, group=taxon, fill=sd)) 
p3 <- p3 + scale_fill_gradientn(colours=cols,  
                                values=values, limits=lims, na.value="white",
                                rescaler = function(x, ...) x, oob = identity, name=legend_label) 
p3 <- p3 + geom_point(data=pol_extent, aes(x=x, y=y, shape=extent), colour="#252525", alpha=0.9, size=1.5) +
  scale_shape_manual(values = c(1, 19), name="Coverage")
p3 <- p3 + coord_fixed()
p3 <- add_map_albers(p3, map_data=us.fort, limits)
p3 <- p3 + facet_grid(taxon~.)
p3 <- theme_clean(p3) + theme_bw() + theme(strip.text = element_text(size = 14),
                                           strip.background = element_rect(colour = 'grey'),
                                           axis.title.y = element_blank(),
                                           axis.text = element_blank(),
                                           axis.title = element_blank(),
                                           axis.ticks = element_blank(),
                                           panel.grid.major = element_blank(), 
                                           panel.grid.minor = element_blank())
print(p3)
ggsave(paste0(subDir, '/figure2_sd_map.pdf'))

#######################################################################################################################
## Figure 3
#######################################################################################################################
suff = '-percent-abs'
t_diffs = c(1, 2, 5, 10, 15, 20)

es    = eco_sig(rIts, t_diffs, eco_cut=0.03)
diffs = post_diffs(rIts, t_diffs)

# k_cut = rbind(rep(0.02, K), rep(0.1, K))
k_cut = rbind(rep(0.05, K), rep(0.1, K))
k_cut_low = rbind(rep(0.03, K), rep(0.1, K))

class_change = array(NA, c(N, K, length(t_diffs), length(t_diffs)))
class_no_change = array(NA, c(N, K, length(t_diffs), length(t_diffs)))
for (k in 1: K){
  class_change[,k,,] = apply(diffs[,k,,,], c(1,2,3), diff_change_above, cut_change=k_cut[1,k], prob=0.85)
  class_no_change[,k,,] = apply(diffs[,k,,,], c(1,2,3), diff_change_below, cut_change=k_cut_low[1,k], prob=0.85)
}

## now add the mag of the change
# post_mag  = apply(diffs, c(1,2,3,4), diff_mag, prob=0.5)
post_mag  = apply(diffs, c(1,2,3,4), mean, na.rm=TRUE)
# post_mag  = apply(diffs, c(1,2,3,4), median, na.rm=TRUE)
post_mag[which(!es)] = NA

cc_melted = data.frame(change=numeric(0), no_change=character(0), x=integer(0), y=integer(0), t1=numeric(0), t2=numeric(0), taxon=character(), mag=numeric(0))
for (k in 1:K){
  for (t2 in 2:length(t_diffs)){
    for (t1 in (t2-1):1){
      cc_melted  = rbind(cc_melted , data.frame(change = as.vector(class_change[,k,t1,t2]), 
                                                no_change = as.vector(class_no_change[,k,t1,t2]), 
                                                x     = centers_veg[,1]*rescale, 
                                                y     = centers_veg[,2]*rescale, 
                                                t1 = rep(ages[t_diffs[t1]]*100, N),
                                                t2 = rep(ages[t_diffs[t2]]*100, N),
                                                taxon = rep(taxa[k], N),
                                                mag = post_mag[,k,t1,t2]))
    }
  }
}

write.csv(cc_melted, paste0('runs/class_change.csv'))

## aggregate to grid cell where categories are:
# at least one taxon changed
# no taxa changed
# a taxon might have changed
cc_sig = cc_melted[which(cc_melted$change %in% c('change', 'intermediate')),]
cc_all_sig = ddply(cc_sig, .(x,y,t1,t2), function(df) df[df$mag == max(df$mag, na.rm=TRUE),])
cc_all_sig = cc_all_sig[which(!is.na(cc_all_sig$t1)), ] 
cc_all_sig = cc_all_sig[which(!is.na(cc_all_sig$t2)), ] 
cc_all_sig = cc_all_sig[which((cc_all_sig$t1>200)&(cc_all_sig$t2>200)),]

write.csv(cc_all_sig, paste0('runs/class_change_sig.csv'))


cc_no_sig = ddply(cc_melted, .(x,y,t1,t2), function(df) df[all(df$no_change == 'no_change'),])
cc_no_sig = cc_no_sig[which(!is.na(cc_no_sig$t1)), ] 
cc_no_sig = cc_no_sig[which(!is.na(cc_no_sig$t2)), ] 
cc_no_sig = cc_no_sig[which((cc_no_sig$t1>200)&(cc_no_sig$t2>200)),]

cc_all = ddply(cc_melted, .(x,y,t1,t2), function(df) df[df$mag == max(df$mag, na.rm=TRUE),])
cc_all = cc_all[which(!is.na(cc_all$t1)), ] 
cc_all = cc_all[which(!is.na(cc_all$t2)), ] 
cc_all = cc_all[which(abs(cc_all$mag)>0.05), ] 
cc_all = cc_all[which((cc_all$t1>200)&(cc_all$t2>200)),]

## PLOT RESULTS
n_sig = length(unique(cc_all$taxon))

# plot what change the most (same as above) but also outline regions of sig change and sig no change 
cc_all$t1 = factor(cc_all$t1, levels=sort(unique(cc_all$t1), decreasing=TRUE))
cc_all_sig$t1 = factor(cc_all_sig$t1, levels=sort(unique(cc_all_sig$t1), decreasing=TRUE))
cc_no_sig$t1 = factor(cc_no_sig$t1, levels=sort(unique(cc_no_sig$t1), decreasing=TRUE))

cols1=rev(brewer.pal(n_sig+1, "Paired"))[c(1,2,3,4,5,10,11,8,9,6,7)]
# cols1=rev(brewer.pal(n_sig+1, "Paired"))
cols=tim.colors(n_sig)
cols[1]=cols1[1]
cols=cols1

cols = c('#b15928', '#6a3d9a', '#cab2d6', '#1f78b4', 
         '#ff7f00', '#33a02c', '#b2df8a', '#e31a1c', 
         '#fb9a99')

pdf(file=paste0(subDir, '/class_change_most_outline_sig_mask_no_sig', suff, '.pdf'))
p <- ggplot(data=cc_all) + geom_tile(data=cc_all, aes(x=x, y=y, fill=taxon, color=taxon))#, alpha=0.7) 
p <- add_map_albers(p, map_data=us.fort, limits)
p <- p + geom_tile(data=cc_all_sig, aes(x=x, y=y, group=taxon), alpha=0, color="black", size=0.8)
p <- p + geom_tile(data=cc_all_sig, aes(x=x, y=y, group=taxon, fill=taxon, color=taxon))
p <- p + scale_fill_manual(values=cols)
p <- p + scale_color_manual(values=cols)
p <- p + geom_tile(data=cc_no_sig, aes(x=x, y=y), alpha=0, color="black", size=0.8)
p <- p + geom_tile(data=cc_no_sig, aes(x=x, y=y), fill="white", color="white")
p <- p + coord_fixed()
p <- p + facet_grid(t2~t1, switch='x')
p <- theme_clean(p) + theme_bw() + theme(strip.text = element_text(size = 14),
                                         strip.background = element_rect(colour = 'grey'),
                                         axis.title.y = element_blank(),
                                         axis.text = element_blank(),
                                         axis.title = element_blank(),
                                         axis.ticks = element_blank(),
                                         panel.grid.major = element_blank(), 
                                         panel.grid.minor = element_blank())
# topright=c(0.9,0.9)
topright=c(0.25,0.08)
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

# get proportion of cells in each category
prop_change    = aggregate(change~t1+t2, cc_all_sig, FUN=function(x) length(x))
prop_change$percent_change = prop_change$change/704*100

prop_no_change = aggregate(change~t1+t2, cc_no_sig, FUN=function(x) length(x)/12)
colnames(prop_no_change)[3] = 'no_change'
prop_no_change$percent_no_change = prop_no_change$no_change/704*100

change_percent = merge(prop_change, prop_no_change, by=c('t1', 't2'), all.x=TRUE, all.y=TRUE)
change_percent = change_percent[which(change_percent$t1 != 200),]
write.csv(change_percent, file=paste0('runs/', run$suff_run, '/change_percent.csv'))

change_percent$x = rep(0.65*1e6, nrow(change_percent))
change_percent$y = rep(1.45*1e6, nrow(change_percent))
change_percent$percent_change[which(is.na(change_percent$percent_change))] = 0
change_percent$percent_no_change[which(is.na(change_percent$percent_no_change))] = 0
change_percent$label_change = paste0("US: ", as.character(round(change_percent$percent_change)), " %; S: ", as.character(round(change_percent$percent_no_change)), " %")

# which experienced large change
prop_change    = aggregate(change~t1+t2, cc_all, FUN=function(x) length(x))
prop_change$percent_change = prop_change$change/704*100


pdf(file=paste0(subDir, '/class_change_most_outline_sig_mask_no_sig', suff, '.pdf'))
p <- ggplot(data=cc_all) + geom_tile(data=cc_all, aes(x=x, y=y, fill=taxon, color=taxon))#, alpha=0.7) 
p <- add_map_albers(p, map_data=us.fort, limits)
p <- p + geom_tile(data=cc_all_sig, aes(x=x, y=y, group=taxon), alpha=0, color="black", size=0.8)
p <- p + geom_tile(data=cc_all_sig, aes(x=x, y=y, group=taxon, fill=taxon, color=taxon))
p <- p + scale_fill_manual(values=cols)
p <- p + scale_color_manual(values=cols)
p <- p + geom_tile(data=cc_no_sig, aes(x=x, y=y), alpha=0, color="black", size=0.8)
p <- p + geom_tile(data=cc_no_sig, aes(x=x, y=y), fill="white", color="white")
p <- p + coord_fixed()
p <- p + facet_grid(t2~t1, switch='x')
p <- p + geom_text(data=change_percent, aes(x=x, y=y, label=label_change), size=2.5)
p <- theme_clean(p) + theme_bw() + theme(strip.text = element_text(size = 14),
                                         strip.background = element_rect(colour = 'grey'),
                                         axis.title.y = element_blank(),
                                         axis.text = element_blank(),
                                         axis.title = element_blank(),
                                         axis.ticks = element_blank(),
                                         panel.grid.major = element_blank(), 
                                         panel.grid.minor = element_blank())
# topright=c(0.9,0.9)
topright=c(0.25,0.08)
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


#######################################################################################################################
## Figure 4
#######################################################################################################################

#UPPER PANEL
t_diffs = c(2, 5, 10, 15, 20)

es    = eco_sig(rIts, t_diffs, eco_cut=0.03)
diffs = post_diffs(rIts, t_diffs)

post_prob = apply(diffs, c(1,2,3,4), diff_q, prob=0.85)
post_mag  = apply(diffs, c(1,2,3,4), diff_mag, prob=0.85)

post_prob_es = post_prob
post_prob_es[which(!es)] = NA

post_mag_es = post_mag
post_mag_es[which(!es)] = NA

tsuga_lims = list(xlims=c(3*10^5/rescale, max(centers_veg$x)), ylims=c(min(centers_veg$y), 1200000/rescale))
plot_post_change(post_prob_es, t_diffs, ages, centers_veg, N, taxa, taxa_sub=taxa[5], tsuga_lims, rescale, type='prob', subDir=subDir)

# LOWER PANEL
p1 = c(-0.059, 0.862) + 14*c(0.024, 0) - 2*c(0, 0.024)
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

coords_diag = data.frame(coords_diag)*1e6
colnames(coords_diag) = c('x', 'y')

## ecotone proportions
prop_dat = data.frame(props=numeric(0), x=integer(0), y=integer(0),time=integer(0), taxon=character())
for (k in 1:K){
  prop_dat = rbind(prop_dat, data.frame(props = r_mean2[,k], 
                                        x     = rep(centers_veg[,1]*1000000, each=T), 
                                        y     = rep(centers_veg[,2]*1000000, each=T), 
                                        time  = rep(ages, times=N), 
                                        taxon = rep(taxa[k], N*T)))
}

prop_dat = prop_dat[with(prop_dat, order(x)), ]
d_diag = rdist(coords_diag[,1:2], matrix(cbind(prop_dat$x, prop_dat$y), ncol=2))
idx_keep = apply(d_diag, 2, function(x) any(x < 1e-8))
prop_dat1 = prop_dat[idx_keep, ]
prop_dat2 = prop_dat1[which(prop_dat1$taxon %in% c('HEMLOCK')), ]
# prop_dat3 = prop_dat2[order(prop_dat2$x),]

prop_dat3 = data.frame(props=numeric(0), x=numeric(0), y=numeric(0), time=numeric(0), taxon=character(0))
for (t in 1:T){
  prop_sub = prop_dat2[which(prop_dat2$time == ages[t]),]
  
  yvals = seq(min(prop_sub$y), max(prop_sub$y), length=201)
  
  prop_sub_ss = spline(prop_sub$x, prop_sub$props, n=201)
  
  prop_dat3 = rbind(prop_dat3, data.frame(props=prop_sub_ss$y, 
                                          x=prop_sub_ss$x, 
                                          y=yvals, 
                                          time=rep(ages[t]), 
                                          taxon=rep(prop_sub$taxon[1])))
}


prop_dat2$time = prop_dat2$time * 100
prop_dat3$time = prop_dat3$time * 100
prop_dat2 = prop_dat2[which(prop_dat2$time != 200),]
prop_dat3 = prop_dat3[which(prop_dat3$time != 200),]
prop_dat2$time = factor(prop_dat2$time)
prop_dat3$time = factor(prop_dat3$time)

prop_dat2$dist_origin = as.vector(rdist(matrix(p1, nrow=1), matrix(cbind(prop_dat2$x, prop_dat2$y), ncol=2)/1e6))*1e3
prop_dat3$dist_origin = as.vector(rdist(matrix(p1, nrow=1), matrix(cbind(prop_dat3$x, prop_dat3$y), ncol=2)/1e6))*1e3

colfunc <- colorRampPalette(c('#fee391', '#fec44f', '#fe9929', '#ec7014', '#cc4c02', '#993404', '#662506'))

q <- ggplot() + geom_point(data=prop_dat2, aes(x=dist_origin, y=props, colour=time)) + 
  geom_line(data=prop_dat3, aes(x=dist_origin, y=props, colour=time), size=1) +
  scale_colour_manual(values=colfunc(20), name='YB1950') +
  ylim(c(0,0.65)) + 
  ylab('Proportion') + 
  xlab('Distance along transect (km)') 
q <- q + theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14)) + guides(colour=guide_legend(nrow=10))
print(q)
ggsave(q, file=paste0(subDir,'/hemlock_diag_trans_props.pdf'))


coords_diag$dist_origin = as.vector(rdist(matrix(p1, nrow=1), matrix(cbind(coords_diag$x, coords_diag$y), ncol=2)/1e6))*1e3

p <- ggplot(data=subset(prop_dat, taxon == 'HEMLOCK')) + geom_tile(aes(x=x, y=y, fill=props)) 
p <- p + scale_fill_gradientn(colours=c("white", "grey80", "grey58", "grey22"),  
                              values=c(0,0.05,0.15, 1), na.value="white",
                              rescaler = function(x, ...) x, name='Proportion', guide=FALSE) 
# p <- p + geom_line(data=coords_diag, aes(x=x, y=y, colour=dist_origin), size=2)
p <- p + geom_line(data=coords_diag, aes(x=x, y=y, colour=dist_origin), size=2)
p <- p + scale_colour_gradientn(colors=tim.colors(20), name='Distance along transect') 
p <- p + coord_fixed()
p <- add_map_albers(p, map_data=us.fort, limits)
p <- p + theme_bw()
p <- theme_clean(p)
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  
               panel.border = element_rect(colour = "black", fill=NA, size=0.8))
print(p)
ggsave(p, file=paste0(subDir,'/hemlock_trans_map.pdf'))


#A viewport taking up a fraction of the plot area
vp <- viewport(width = 0.43, height = 0.43, x = 0.7, y = 0.75)

#Just draw the plot twice
pdf(paste0(subDir, "/hemlock_diag_transect_panel.pdf"), width=8, height=6)
print(q)
print(p, vp = vp)
dev.off()

#######################################################################################################################
## Figure 5
#######################################################################################################################

t_diffs = c(2, 5, 10, 15, 20)

es    = eco_sig(rIts, t_diffs, eco_cut=0.03)
diffs = post_diffs(rIts, t_diffs)

post_prob = apply(diffs, c(1,2,3,4), diff_q, prob=0.85)
post_mag  = apply(diffs, c(1,2,3,4), diff_mag, prob=0.85)

post_prob_es = post_prob
post_prob_es[which(!es)] = NA

post_mag_es = post_mag
post_mag_es[which(!es)] = NA

pine_lims = list(xlims=c(-0.05, 0.9), ylims=c(0.670, 1.45))
plot_post_change(post_prob_es, t_diffs, ages, centers_veg, N, taxa, taxa_sub=taxa[10], pine_lims, rescale, type='prob', subDir=subDir)

# LOWER PANELA (eof)
source('r/utils/util.r')

r_mean = apply(rIts, c(1,2,3), mean, na.rm=TRUE)

niter=dim(rIts)[4]
burn=1500
rItsSub = rIts[,,,burn:niter]

myTaxaEOFs=taxaEOFs(r_mean,rItsSub)

# scores
scores_dat = melt(myTaxaEOFs$scoreMean)
colnames(scores_dat) = c('cell', 'EOF', 'time', 'score')
cveg = do.call("rbind", replicate(nrow(scores_dat)/nrow(centers_veg), centers_veg, simplify = FALSE))
scores_dat = cbind(scores_dat, cveg)
scores_dat$time = ages[scores_dat$time]
scores_dat = scores_dat[which(scores_dat$EOF %in% seq(1,5)),]
scores_dat$time = scores_dat$time * 100

coords_diag = get_diags(centers_veg)
cdiag=data.frame(x=coords_diag[,1], y=coords_diag[,2], diag=coords_diag[,3])
cdiag1 = cdiag[cdiag$diag==6, ]

limits2 = limits
limits2$xlims[2] = 0.750
limits2$ylims[1] = 0.690
limits2$ylims[2] = 1.460

scores_dat1 = scores_dat[scores_dat$EOF == 1,]
scores_dat1 = scores_dat1[which(scores_dat1$time != 200),]
s <- ggplot(scores_dat1) + 
  geom_contour(data=scores_dat1, aes(z=score, x=x*rescale, y=y*rescale, group=factor(time),colour=factor(time)), size=1,breaks=0) + 
  # geom_contour(data=scores_dat1, aes(z=score, x=x, y=y, group=factor(time),colour=factor(time)), size=2,breaks=-1) + 
  coord_fixed() + geom_line(data=cdiag1, aes(x=x*rescale, y=y*rescale, group=diag), linetype=2, size=1.3, colour='lightgrey')#+
#geom_dl(data=cdiag, aes(x=x*rescale, y=y*rescale, group=diag, label = diag), method = list("last.points", dl.trans(x=x+0.1, y=y+0.1) ), cex = 0.8)
s <- s + scale_colour_manual(values = tim.colors(19), name='Cal Yr BP')
s <- add_map_albers(s, map_data=us.fort, limits2)
s <- s + theme_bw() + theme(legend.text=element_text(size=14), legend.title=element_text(size=16))
s <- theme_clean(s) 
topright=c(0.95,0.95)
# topright=c(0.25,0.08)
s <- s + ggsn::scalebar(dist=100, 
                        transform=FALSE, 
                        dist_unit="km",
                        location="topright",
                        height=0.03,
                        st.dist=0.03,
                        st.size=3,
                        border.size=0.4,
                        anchor=c(x=(1-topright[1])*limits2$xlims[1]*1000000+topright[1]*limits2$xlims[2]*1000000,
                                 y=(1-topright[2])*limits2$ylims[1]*1000000+topright[2]*limits2$ylims[2]*1000000),
                        x.min=limits2$xlims[1]*1000000,
                        x.max=limits2$xlims[2]*1000000,
                        y.min=limits2$ylims[1]*1000000,
                        y.max=limits2$ylims[2]*1000000)
print(s)
ggsave(file=paste0(subDir, '/ecotone_eof_sign_change_paper.pdf'))

# find where the sign change happens
#find where the sign changes
scores_dat = melt(myTaxaEOFs$scoreMean)
colnames(scores_dat) = c('cell', 'EOF', 'time', 'score')
cveg = do.call("rbind", replicate(nrow(scores_dat)/nrow(centers_veg), centers_veg, simplify = FALSE))
scores_dat = cbind(scores_dat, cveg)
scores_dat$time = ages[scores_dat$time]
scores_dat = scores_dat[which(scores_dat$EOF %in% seq(1,5)),]
scores_dat$time = scores_dat$time * 100

scores_sub = diag_transect(scores_dat, coords_diag, rescale_coords=1e6, rescale_dat=1e6)
scores_sub = scores_sub[scores_sub$EOF == 1,]
# scores_sub = scores_sub[scores_sub$time %in% c(50, seq(150, 1950, by=200)),]
scores_sub = scores_sub[scores_sub$time %in% (ages[c(1,2,5,10,15,20, 25)]*100),]
ts=unique(scores_sub$time)

scores_keep = data.frame(matrix(NA, nrow=0, ncol=4))
colnames(scores_keep)= c('time', 'y', 'score', 'diag')

scores_zero = data.frame(matrix(NA, nrow=0, ncol=4))
colnames(scores_zero)= c('time', 'y', 'score', 'diag')

n_diags = length(unique(coords_diag$diag))
for (j in 1:n_diags){
  d_diag=rdist(coords_diag[coords_diag$diag==j,1:2],matrix(cbind(scores_sub$x, scores_sub$y), ncol=2))
  idx_keep = apply(d_diag, 2, function(x) any(x < 1e-8))
  # id_keep  = apply(d_diag, 2, function(x){if(any(x < 1e-10)){ which(x < 1e-10)}})
  score_sub = scores_sub[idx_keep, ]
  for (t in 1:length(ts)){
    score_add = score_sub[score_sub$time == ts[t],]
    score_add = score_add[order(score_add$y),]
    score_add$d = seq(0, sqrt(2*0.024^2)*(nrow(score_add)-1), by=sqrt(2*0.024^2))
    app = approx(score_add$d, score_add$score, n=500)
    foo = data.frame(time=rep(ts[t]), y=app$x, score=app$y, diag=rep(j))
    
    scores_zero = rbind(scores_zero, foo[which.min(abs(foo$score)),])
    
    scores_keep = rbind(scores_keep, foo)
  }
}

scores_zero_re = ddply(scores_zero, .(diag), summarize, time=time, value=(y-head(y, 1))*1e3)
scores_zero_re1 = scores_zero_re[scores_zero_re$diag==6,] 
scores_zero_re1 = scores_zero_re1[which(scores_zero_re1$time != 200),]
q <- ggplot() + 
  geom_point(data=scores_zero_re1, aes(x=time, y=value), size=3) + 
  scale_x_continuous(breaks=c(200, 500, 1000, 1500, 2000, 1000)) +
  scale_y_continuous(breaks=c(0, 2.5, 5, 7.5, 10, 12.5)) + 
  xlab("Cal Yr BP") + 
  ylab("Distance from present\n position (km)") + 
  theme_bw() +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14), 
        panel.border=element_rect(color="black", fill=NA))
print(q)
ggsave(paste0(subDir, "/ecotone_distance_from_present.pdf"), width=4.5, height=3)

#######################################################################################################################
## Figure 6: Big Woods
#######################################################################################################################

bw_buff = 164000
bw_lims <- list(xlims=c(min(bw_fort$long)-bw_buff, max(bw_fort$long)+bw_buff)/1e6, 
                ylims=c(min(bw_fort$lat)-64000, max(bw_fort$lat)+bw_buff)/1e6)

r_bw = rowSums(r[,c(1,4)])
plot_data_maps_binned_bw(r_bw, centers_comp/1e6, taxa='bw', 
                         nrow(centers_comp), K=1, breaks, bw_lims, 
                         suff='pls-comp_bw', save_plots, fpath=subDir)

rp_bw = rowSums(r_mean2[seq(2,N*T, T),c(1,4)])
plot_data_maps_binned_bw(rp_bw, centers_veg, taxa='pbw', 
                         N, K=1, breaks, bw_lims, 
                         suff='pred-comp_bw', save_plots, fpath=subDir)

# visual checks
r2_bw = rbind(data.frame(props=r_bw, x=centers_comp[,1], y=centers_comp[,2], type='PLS composition'),
              data.frame(props=rp_bw, x=centers_veg[,1]*rescale, y=centers_veg[,2]*rescale, type='STEPPS'))

p <- ggplot() + geom_tile(data=r2_bw, aes(x=x, y=y, fill=factor(props))) + 
  scale_fill_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportions') + 
  coord_fixed() + 
  scale_x_continuous(limits$xlims*rescale) + scale_y_continuous(limits$ylims*rescale)
p <- add_map_albers(p, us.shp, limits)
p <- p + geom_path(data=bw, aes(x=long, y=lat),  colour='grey', alpha=0.8,size=1.8)
p <- theme_clean(p) #+ theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
p <- p + theme(strip.text.x = element_blank(),
               strip.text.y = element_blank())
p <- p + theme(strip.background = element_blank())
print(p)

p <- ggplot() + geom_tile(data=r2_bw, aes(x=x, y=y, fill=props)) + 
  coord_fixed() + 
  scale_x_continuous(limits$xlims*rescale) + scale_y_continuous(limits$ylims*rescale)
p <- add_map_albers(p, us.shp, limits)
p <- p + geom_path(data=bw, aes(x=long, y=lat),  colour='grey', alpha=0.8,size=1.8)
p <- p + facet_grid(.~type)
p <- theme_clean(p) 
p <- p + theme(strip.text.x = element_blank(),
               strip.text.y = element_blank())
p <- p + theme(strip.background = element_blank())
print(p)

# posterior difference plot
t_diffs = c(2, 5, 10, 15, 20)

es    = eco_sig(rIts, t_diffs, eco_cut=0.03)
diffs = post_diffs(rIts, t_diffs)

post_prob = apply(diffs, c(1,2,3,4), diff_q, prob=0.85)
post_mag  = apply(diffs, c(1,2,3,4), diff_mag, prob=0.85)

post_prob_es = post_prob
post_prob_es[which(!es)] = NA

post_mag_es = post_mag
post_mag_es[which(!es)] = NA

plot_post_change(post_prob_es, t_diffs, ages, centers_veg, N, taxa, taxa_sub=taxa, limits, rescale, type='prob', subDir=subDir)
plot_post_change(post_mag_es, t_diffs, ages, centers_veg, N, taxa, taxa_sub=taxa, limits, rescale, type='mag', subDir=subDir)

rIts_bw = apply(rIts[,c(1,4),,], c(1,3,4), sum, na.rm=TRUE)
es_bw = eco_sig_bw(rIts_bw, t_diffs, eco_cut=0.03)
diffs_bw = post_diffs_bw(rIts_bw, t_diffs)

post_prob_bw = apply(diffs_bw, c(1,2,3), diff_q, prob=0.85)
post_mag_bw  = apply(diffs_bw, c(1,2,3), diff_mag, prob=0.85)

post_prob_bw_es = post_prob_bw
# post_prob_bw_es[which(!es_bw)] = NA
# 
post_mag_bw_es = post_mag_bw
# post_mag_bw_es[which(!es_bw)] = NA

plot_post_change_group(post_prob_bw_es, t_diffs, ages, centers_veg, N, taxa[c(1,4)], 
                       rescale, limits=bw_lims, type='prob', group='bw', subDir)
plot_post_change_group(post_mag_bw_es, t_diffs, ages, centers_veg, N, taxa[c(1,4)], 
                       rescale, limits=bw_lims, type='mag', group='bw', subDir)
#######################################################################################################################
## Appendix Figures
#######################################################################################################################
p_binned <- plot_pred_maps_binned_select(r_mean2, centers_veg, breaks, taxa, taxa_sub=taxa_sub, ages, N, K, T, limits, 
                                         suff=run$suff_figs, save_plots=TRUE, fpath=subDir)

# plot_pred_maps_binned_gif
p_binned <- plot_pred_maps_binned_gif(r_mean2, centers_veg, breaks, taxa, taxa_sub=taxa, ages, N, K, T, limits, 
                                      suff=run$suff_figs, save_plots=TRUE, fpath=subDir)

p_binned <- plot_pred_maps_binned_gifun(r_mean2, r_sd2, centers_veg, breaks, taxa, taxa_sub=taxa, ages, N, K, T, limits, 
                                        suff=run$suff_figs, save_plots=TRUE, fpath=subDir)

suff = paste0('sd_', run$suff_figs)
p_binned <- plot_pred_maps_binned_select(r_sd2, centers_veg, breaks, taxa, taxa_sub=taxa, ages, N, K, T, limits, 
                                         suff=suff, save_plots, fpath=subDir)

# posterior difference plots
t_diffs = c(2, 5, 10, 15, 20)

es    = eco_sig(rIts, t_diffs, eco_cut=0.03)
diffs = post_diffs(rIts, t_diffs)

post_prob = apply(diffs, c(1,2,3,4), diff_q, prob=0.85)
post_mag  = apply(diffs, c(1,2,3,4), diff_mag, prob=0.85)

post_prob_es = post_prob
post_prob_es[which(!es)] = NA

post_mag_es = post_mag
post_mag_es[which(!es)] = NA

plot_post_change(post_prob_es, t_diffs, ages, centers_veg, N, taxa, taxa_sub=taxa, limits, rescale, type='prob', subDir=subDir)
plot_post_change(post_mag_es, t_diffs, ages, centers_veg, N, taxa, taxa_sub=taxa, limits, rescale, type='mag', subDir=subDir)






##########################################################################################################################################
## Additional code for figures not included in paper
##########################################################################################################################################
t_diffs = c(2, 5, 10, 15, 20)

es    = eco_sig(rIts, t_diffs, eco_cut=0.03)
diffs = post_diffs(rIts, t_diffs)

post_prob = apply(diffs, c(1,2,3,4), diff_q, prob=0.85)
post_mag  = apply(diffs, c(1,2,3,4), diff_mag, prob=0.85)

post_prob_es = post_prob
post_prob_es[which(!es)] = NA

post_mag_es = post_mag
post_mag_es[which(!es)] = NA

# deciduous
plot_post_change_group(post_prob_es[,c(1,2,3,4,6,7,9),,], t_diffs, ages, centers_veg, N, taxa[c(1,2,3,4,6,7,9)],
                       rescale, limits=limits, type='prob', group='deciduous', subDir=subDir)
# deciduous
plot_post_change_group(post_mag_es[,c(1,2,3,4,6,7,9),,], t_diffs, ages, centers_veg, N, taxa[c(1,2,3,4,6,7,9)],
                       rescale, limits=limits, type='mag', group='deciduous', subDir=subDir)

# conifer
plot_post_change_group(post_prob_es[,c(5,8,10,11,12),,], t_diffs, ages, centers_veg, N, taxa[c(5,8,10,11,12)],
                       rescale, limits=limits, type='prob', group='conifer', subDir)
plot_post_change_group(post_mag_es[,c(5,8,10,11,12),,], t_diffs, ages, centers_veg, N, taxa[c(5,8,10,11,12)],
                       rescale, limits=limits, type='mag', group='conifer', subDir)

#################################################################################################################################################
# hists of changes for 500 and 1000 year intervals
#################################################################################################################################################

eco_sig_consec <- function(rIts, t_diffs, eco_cut){
  
  niter=dim(rIts)[4]
  n_diffs = length(t_diffs)
  
  es  = array(NA, c(N, K, n_diffs-1))
  for (i in 1:(n_diffs-1)){
      
      es_i = apply(rIts[,,t_diffs[i],], c(1,2), function(x) mean(x) > eco_cut)
      es_ip1 = apply(rIts[,,t_diffs[i+1],], c(1,2), function(x) mean(x) > eco_cut)
      
      es[,,i] = es_i | es_ip1
  }
  
  return(es)
}

post_diffs_consec <- function(rIts, t_diffs){
  
  N=dim(rIts)[1]
  K=dim(rIts)[2]
  niter=dim(rIts)[4]
  n_diffs = length(t_diffs)
  
  diffs = array(NA, c(N, K, n_diffs-1, niter))
  for (i in 1:(n_diffs-1)){
      diffs[,,i,] = rIts[,,t_diffs[i],] - rIts[,,t_diffs[i+1],] 
  }
  
  return(diffs)
}

t_diffs = c(2,7,12,17)
es    = eco_sig_consec(rIts, t_diffs, eco_cut=0.03)
diffs = post_diffs_consec(rIts, t_diffs)

post_mag  = apply(diffs, c(1,2,3), mean, na.rm=TRUE)
post_mag_es = post_mag
post_mag_es[which(!es)] = NA

cc_melted = data.frame(x=integer(0), y=integer(0), t1=numeric(0), t2=numeric(0), scale=numeric(0), taxon=character(), mag=numeric(0))
for (k in 1:K){
  for (t1 in 1:(length(t_diffs)-1)){
      cc_melted  = rbind(cc_melted , data.frame(x     = centers_veg[,1]*rescale, 
                                                y     = centers_veg[,2]*rescale, 
                                                t1 = rep(ages[t_diffs[t1]]*100, N),
                                                t2 = rep(ages[t_diffs[t1+1]]*100, N),
                                                scale = rep(500, N),
                                                taxon = rep(taxa[k], N),
                                                mag = post_mag_es[,k,t1]))
  }
}

ggplot() + geom_histogram(data=cc_melted, aes(x=mag)) + geom_vline(xintercept=0, colour="red")+ facet_wrap(~taxon)
ggsave(paste0(subDir, '/mag_histogram.pdf'))

## differences in magnitude for 500 and 1000 years
cc_melted = data.frame(x=integer(0), y=integer(0), t1=numeric(0), t2=numeric(0), scale=numeric(0), taxon=character(), mag=numeric(0))

# for 500
for (i in 1:15){
  t_diffs = c(i, i+5)
  es    = eco_sig_consec(rIts, t_diffs, eco_cut=0.03)
  diffs = post_diffs_consec(rIts, t_diffs)
  
  post_mag  = apply(diffs, c(1,2,3), mean, na.rm=TRUE)
  post_mag_es = post_mag
  post_mag_es[which(!es)] = NA
  
  for (k in 1:K){
    for (t1 in 1:(length(t_diffs)-1)){
      cc_melted  = rbind(cc_melted , data.frame(x     = centers_veg[,1]*rescale, 
                                                y     = centers_veg[,2]*rescale, 
                                                t1 = rep(ages[t_diffs[t1]]*100, N),
                                                t2 = rep(ages[t_diffs[t1+1]]*100, N),
                                                scale = rep(500, N),
                                                taxon = rep(taxa[k], N),
                                                mag = post_mag_es[,k,t1]))
    }
  }
}


# for 1000
for (i in 1:10){
  t_diffs = c(i, i+10)
  es    = eco_sig_consec(rIts, t_diffs, eco_cut=0.03)
  diffs = post_diffs_consec(rIts, t_diffs)
  
  post_mag  = apply(diffs, c(1,2,3), mean, na.rm=TRUE)
  post_mag_es = post_mag
  post_mag_es[which(!es)] = NA
  
  for (k in 1:K){
    for (t1 in 1:(length(t_diffs)-1)){
      cc_melted  = rbind(cc_melted , data.frame(x     = centers_veg[,1]*rescale, 
                                                y     = centers_veg[,2]*rescale, 
                                                t1 = rep(ages[t_diffs[t1]]*100, N),
                                                t2 = rep(ages[t_diffs[t1+1]]*100, N),
                                                scale = rep(1000, N),
                                                taxon = rep(taxa[k], N),
                                                mag = post_mag_es[,k,t1]))
    }
  }
}

xlims=c(-0.05, 0.05)

ggplot() + geom_vline(xintercept=0, colour="red") + 
  geom_histogram(data=cc_melted, aes(x=mag, ..density.., fill=factor(scale)), bins=20) + xlim(xlims) + facet_wrap(~taxon)
ggsave(paste0(subDir, '/mag_histogram.pdf'))

ggplot() + geom_vline(xintercept=0, colour="grey") + 
  geom_freqpoly(data=cc_melted, aes(x=mag, ..density.., color=factor(scale)), bins=20) + 
  xlim(xlims) + facet_wrap(~taxon)
ggsave(paste0(subDir, '/mag_freqpoly.pdf'))

ggplot() + geom_vline(xintercept=0, colour="grey") + 
  geom_density(data=cc_melted, aes(x=mag, ..density.., color=factor(scale))) + xlim(xlims) + facet_wrap(~taxon)
ggsave(paste0(subDir, '/mag_density.pdf'))

quants = aggregate(mag~scale+taxon, cc_melted, function(x) quantile(x, probs=c(0.05,0.5, 0.95)))

write.csv(quants, paste0('runs/', run$suff_run, '/mag_quants_millenial-scale.csv'))

quants_abs = aggregate(mag~scale+taxon, cc_melted, function(x) quantile(abs(x), probs=c(0.05,0.5, 0.95)))

write.csv(quants_abs, paste0('runs/', run$suff_run, '/mag_quants_millenial-scale-abs.csv'))

#################################################################################################################################################
# plot mag of all changes, not just sig
#################################################################################################################################################
t_diffs = c(2, 5, 10, 15, 20)

es    = eco_sig(rIts, t_diffs, eco_cut=0.03)
diffs = post_diffs(rIts, t_diffs)
post_mag  = apply(diffs, c(1,2,3,4), mean, na.rm=TRUE)
post_mag_es = post_mag
post_mag_es[which(!es)] = NA

plot_post_change(post_mag_es, t_diffs, ages, centers_veg, N, taxa, taxa_sub=taxa, limits, rescale, type='mag', subDir=subDir, suff='')

post_sd  = apply(diffs, c(1,2,3,4), sd, na.rm=TRUE)
post_sd_es = post_sd
post_sd_es[which(!es)] = NA

plot_post_change(post_sd_es, t_diffs, ages, centers_veg, N, taxa, taxa_sub=taxa, limits, rescale, type='sd', subDir=subDir, suff='')


##########################################################################################################################################
## maximum probability of change
##########################################################################################################################################
t_diffs = c(1, 2, 5, 10, 15, 20)
n_diffs = length(t_diffs)
niter   = dim(rIts)[4]

es    = eco_sig(rIts, t_diffs, eco_cut=0.03)
diffs = post_diffs(rIts, t_diffs)

post_prob = apply(diffs, c(1,2,3,4), diff_q, prob=0.85)

post_prob_es = post_prob
post_prob_es[which(!es)] = NA

post_prob_max = apply(post_prob_es, c(1,3,4), function(x) if (all(is.na(x))){NA} else {max(x, na.rm=TRUE)})

# start at 2?
prob_max = data.frame(p=numeric(0), x=integer(0), y=integer(0), t1=numeric(0), t2=numeric(0))
for (t2 in 2:n_diffs){
  for (t1 in (t2-1):1){
    prob_max  = rbind(prob_max, data.frame(p = post_prob_max[,t1,t2],
                                           x = centers_veg[,1]*rescale,
                                           y = centers_veg[,2]*rescale,
                                           t1 = rep(t_diffs[t1], N),
                                           t2 = rep(t_diffs[t2], N)))
  }
}

prob_max$t1 = ages[prob_max$t1]*100
prob_max$t2 = ages[prob_max$t2]*100

prob_max$d_bin = rep(NA, nrow(prob_max))
prob_max$d_bin[which(prob_max$d>=0.85)] = 1

values = c(0, 0.5, 0.9, 1)

p <- ggplot(data=prob_max) + geom_tile(aes(x=x, y=y, fill=p)) 
p <- p + scale_fill_gradientn(colours=c("white", "white", "white", "red", "red"),
                              values=values, limits=c(0,1), na.value="white",
                              rescaler = function(x, ...) x, oob = identity)
p <- p + coord_fixed()
p <- add_map_albers(p, map_data=us.fort, limits)
p <- p + facet_grid(t2~t1, switch='x')
p <- theme_clean(p) + theme_bw() + theme(strip.text = element_text(size = 12),
                                         strip.background = element_rect(colour = 'grey'),
                                         axis.ticks = element_blank(), 
                                         axis.text = element_blank(),
                                         axis.title = element_blank(),
                                         legend.text = element_text(size = 12))+ labs(color='Prob') 
print(p)
ggsave(paste0(subDir, '/post_diffs_prob_max.pdf'))
