#' Plot zygosity HMM
#' @description a 3-4 track visualization of an HMM fitted to a series of observed
#' heterozygous counts. The top track will be a simple count of heterozygous counts
#' along the genome. The middle track(s) will be the estimated state based on a
#' fitted HMM model. If the actual state is known, it will also include the visualization
#' for the actual state on top of the estimated states. The bottom track contains
#' line graphs for the posterior probability of each state.
#'
#' @param model Data frame output fitHmm() function with obs, states, and PP
#' @param est_cols 2 element character vector for the estimated [state, state.label] columns
#' @param act_cols If known, 2 element character vector for the actual [state, state.label] columns
#' @param xval String for column ID corresponding to x-axis value (e.g. bin)
#' 
#' @importFrom dplyr %>%
#' @importFrom purrr map
#' @importFrom assertthat assert_that
#' @import ggplot2
#' @import reshape2 
#' @import gridExtra
#' 
#' @return
#' a grid object for visualization
#' @export
plotHMM <- function(model, est_cols, act_cols=NULL, xval='bin'){
  known_actual <- if(!is.null(act_cols)) TRUE else FALSE
  mycols <- c("#d73027",  
              "#7fcdbb", "#41b6c4", "#1d91c0",
              "#225ea8", "#253494", "#081d58") #red, blue1, blue2, blue3...
  # Ensure the standards of est_cols
  .check(est_cols, 'plotHMM.est_cols', m=model)
  if(known_actual) .check(act_cols, 'plotHMM.est_cols', m=model)
  
  # Ensures there are enough colours to go around
  num_uniq_states <- length(unique(model[,est_cols[1]]))
  if(!is.null(act_cols)) num_uniq_states <- max(c(num_uniq_states,
                                                  length(unique(model[,act_cols[1]]))))
  stopifnot(length(mycols) >= num_uniq_states)
  mycols          <- mycols[1:num_uniq_states]
  
  # Identify the chromosome transition points
  idxs <- WadingPool:::.getChrIdx(model)
  
  # Attributing state to state-label
  state_map <- WadingPool:::.genStateMap(model, est_cols, act_cols)
  uniq_map  <- state_map[['map']]
  act_map   <- state_map[['act.map']]
  model     <- state_map[['model']]
  model$bin <- as.integer(gsub("^bin", "", model$bin))
  
  # Plot observed Het cnts
  gobs <- (ggplot(model, aes_string(x = xval, y = 'obs')) + 
             geom_line() +
             geom_vline(xintercept = model[idxs$chrs, xval], size=1, colour='grey') +
             theme_classic() +
             lims(y=c(-3,3)) +
             theme(axis.ticks = element_blank(), 
                   axis.title.y = element_blank()) +
             scale_x_continuous(breaks=model[idxs$mids, xval],
                                labels=names(idxs$mids))) %>% ggplotGrob
  pdf("~/test.pdf")
  grid.arrange(gobs)
  dev.off()
  
  # Barplot of different states (actual or estimated)
  .barplot_hmm <- function(model, est_cols, typelab="Estimated", labels, idxs){
    idxs <- .getChrIdx(model)
    
    (ggplot(model, aes_string(x = xval, y = est_cols[1], fill = est_cols[2], col = est_cols[2])) + 
        geom_bar(stat = "identity", alpha = I(0.7), size=0) + 
        geom_vline(xintercept = model[idxs$chrs, xval], size=1.5, colour='grey') +
        theme_classic() +
        scale_fill_manual(values = mycols, name = "State:", labels = labels) +
        scale_color_manual(values = mycols, name = "State:", labels = labels) +
        theme(axis.ticks = element_blank(), axis.text.y = element_blank()) +
        labs(y = paste0(typelab, " State")) +
        scale_x_continuous(breaks=model[idxs$mids, xval],
                           labels=names(idxs$mids))) %>% ggplotGrob
  }
  gest <- .barplot_hmm(model, est_cols, typelab="Estimated", uniq_map[,est_cols[1]], idxs)
  if(known_actual) gact <- .barplot_hmm(model, act_cols, typelab="Actual", act_map[,act_cols[1]], idxs)
  
  # Line Plot of the posterior probabilities
  melt_model <- melt(model[,c(xval, paste0("S", uniq_map[,est_cols[1]]), 'chrom')], 
                     id = c(xval, 'chrom'))
  gpp <- (ggplot(melt_model, aes_string(x = xval, y = 'value', col = 'variable')) + 
            geom_line() +
            geom_vline(xintercept = model[idxs$chrs, xval], size=1.5, colour='grey') +
            theme_classic() +
            lims(y=c(0,1)) +
            scale_color_manual(values = mycols, name = "State:", labels = uniq_map[,est_cols[1]]) +
            theme(axis.ticks = element_blank(), axis.text.y = element_blank()) + 
            labs(y = "Posterior Prob.") +
            scale_x_continuous(breaks=model[idxs$mids, xval],
                               labels=names(idxs$mids))) %>%
    ggplotGrob()
  
  gobs$widths <- gest$widths
  if(known_actual){
    grid.arrange(gobs, gact, gest, gpp, widths = 1, nrow = 4)
  } else {
    grid.arrange(gobs, gest, gpp, widths = 1, nrow = 3)
  }
}

#' Cn-Zyg plotter
#' @description Plots the zygosity and cn state
#'
#' @param cn GRanges object with a discrete integer cncol column
#' @param zyg GRanges object with a discrete 'state' column of c('LOH', 'H1', etc...)
#' @param tiles Reference GRanges object from tileGenome()
#' @param cncol character: copy number column for cn GRanges object
#' @param ... pass in for plot()
#'
#' @importFrom S4Vectors queryHits
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors mcols
#' @importFrom GenomicRanges findOverlaps
#' @importFrom graphics plot
#' @importFrom graphics rect
#' @importFrom graphics abline
#' @importFrom graphics axis
#' 
#' @return
#' data frame containing the segment-level CN and zygosity state
#' 
#' @export
plotCnZyg <- function(cn, zyg, tiles, cncol='copy.number', ...){
  # cn <- segs[['net-001a']]
  # zyg <- combZygPos(zyg=hmm1$df, zygpos=tiles_ov)
  # cncol <- 'copy.number'
  mycols <- c("darkgrey", "#d73027",  
              "#7fcdbb", "#41b6c4", "#1d91c0",
              "#225ea8", "#253494", "#081d58") #red, blue1, blue2, blue3...
  
  # For plotting, range of CN from 0 to maxCN
  yrange <- c(0, max(mcols(cn)[,cncol]))
  
  # Initialize a blank genomic-positional data.frame for zygosity and CN to map to
  cnzyg <- as.data.frame(matrix(nrow=length(tiles), ncol=2))
  colnames(cnzyg) <- c('zyg', 'cn')
  cnzyg$xpos <- c(1:length(tiles))
  
  # Map zygosity to the reference genomic bins (tiles)
  ov <- findOverlaps(tiles, zyg, select='all')
  if(any(duplicated(queryHits(ov)))) {
    warning("Duplicates in overlap were found, the tiles genomic bins may not be
            granular enough to capture all the data")
    ov <- ov[-which(duplicated(queryHits(ov)))]
  }
  zygstat <- zyg[subjectHits(ov)]
  zygstat <- factor(zygstat$state, levels=c("LOH", paste0("H", c(1:(length(table(zygstat$state))-1)))))
  cnzyg[queryHits(ov),]$zyg <- zygstat
  if(any(is.na(cnzyg$zyg))) cnzyg$zyg[is.na(cnzyg$zyg)] <- 0
  
  # Map CN to the reference genomic bins (tiles)
  ov <- findOverlaps(tiles, cn)
  if(any(duplicated(queryHits(ov)))) {
    warning("Duplicates in overlap were found, the tiles genomic bins may not be
            granular enough to capture all the data")
    ov <- ov[-which(duplicated(queryHits(ov)))]
  }
  cnstat <- mcols(cn[subjectHits(ov)])[,cncol]
  cnzyg[queryHits(ov),]$cn <- cnstat
  
  # Reduce bin level data to segment level cn/zyg states
  cnzyg_simp <- rle(apply(cnzyg[,c('cn', 'zyg')], 1, paste, collapse="_"))
  cnzyg_cum <- cumsum(c(0, cnzyg_simp$lengths))
  cnzyg_simp <- do.call(rbind, strsplit(cnzyg_simp$values, "_"))
  storage.mode(cnzyg_simp) <- 'integer'
  
  # Annotate the reduced seg with chr, start, end, binpos
  cnzyg_red <- as.data.frame(cnzyg_simp)
  colnames(cnzyg_red) <- c('cn', 'zyg')
  cnzyg_red$sbin <- cnzyg[cnzyg_cum[-length(cnzyg_cum)]+1,'xpos']
  cnzyg_red$ebin <- cnzyg[cnzyg_cum[-1],'xpos']
  cnzyg_red$chr <- as.character(seqnames(tiles[cnzyg_red$sbin]))
  cnzyg_red$start <- start(tiles[cnzyg_red$sbin])
  cnzyg_red$end <- end(tiles[cnzyg_red$ebin])
  
  
  # Plot the CN/Zyg data
  chr_brk <- rle(as.character(seqnames(tiles)))
  plot(0, type='n',
       pch=16, ylim=yrange, xlim=c(0, length(tiles)), 
       xlab = 'Chromosomes', ylab="Copy number", xaxt='n',
       cex.lab=1,cex.axis=1, cex=0.7, las=1, ...)
  rect(xleft = cnzyg_cum[-length(cnzyg_cum)], 
       ybottom = (cnzyg_simp[-nrow(cnzyg_simp),1]-0.3), 
       xright = cnzyg_cum[-1], 
       ytop = (cnzyg_simp[-nrow(cnzyg_simp),1]+0.3), 
       col=mycols[cnzyg_simp[,2]+1], border = NA)
  abline(v = cumsum(chr_brk$length))
  axis(side = 1, at = cumsum(chr_brk$length)-(chr_brk$length/2),
       labels = gsub("chr", "", chr_brk$values, ignore.case = T), 
       tick = FALSE, cex.axis=0.7)
  
  return(cnzyg_red)
}


# Identify the chromosome transition points
.getChrIdx <- function(model){
  chr_idx <- cumsum(rle(model$chrom)$lengths)
  mid_chr <- ceiling(chr_idx - diff(c(1,chr_idx)) /2)
  mid_chr <- setNames(mid_chr, gsub("^chr([om]*)?", "", rle(model$chrom)$values))
  list("chrs"=chr_idx, "mids"=mid_chr)
}

