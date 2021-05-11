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
  mycols <- c("#1b9e77", "#7570b3", "#d95f02",
              "#e7298a", "#66a61e", "#e6ab02")
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
  idxs <- .getChrIdx(model)
  
  # Attributing state to state-label
  state_map <- .genStateMap(model, est_cols, act_cols)
  uniq_map  <- state_map[['map']]
  act_map   <- state_map[['act.map']]
  model     <- state_map[['model']]
  model$bin <- as.integer(gsub("^bin", "", model$bin))
  
  # Plot observed Het cnts
  gobs <- (ggplot(model, aes_string(x = xval, y = 'obs')) + 
             geom_line() +
             geom_vline(xintercept = idxs$chrs, size=1.5, colour='grey') +
             theme_classic() +
             lims(y=c(-2.5,2.5)) +
             theme(axis.ticks = element_blank(), 
                   axis.title.y = element_blank()) +
             scale_x_continuous(breaks=idxs$mids,
                                labels=names(idxs$mids))) %>% ggplotGrob
  
  # Barplot of different states (actual or estimated)
  .barplot_hmm <- function(model, est_cols, typelab="Estimated", labels, idxs){
    idxs <- .getChrIdx(model)
    
    (ggplot(model, aes_string(x = xval, y = est_cols[1], fill = est_cols[2], col = est_cols[2])) + 
       geom_bar(stat = "identity", alpha = I(0.7), size=0) + 
       geom_vline(xintercept = idxs$chrs, size=1.5, colour='grey') +
       theme_classic() +
       scale_fill_manual(values = mycols, name = "State:", labels = labels) +
       scale_color_manual(values = mycols, name = "State:", labels = labels) +
       theme(axis.ticks = element_blank(), axis.text.y = element_blank()) +
       labs(y = paste0(typelab, " State")) +
       scale_x_continuous(breaks=idxs$mids,
                          labels=names(idxs$mids))) %>% ggplotGrob
  }
  gest <- .barplot_hmm(model, est_cols, typelab="Estimated", uniq_map[,est_cols[1]], idxs)
  if(known_actual) gact <- .barplot_hmm(model, act_cols, typelab="Actual", act_map[,act_cols[1]], idxs)
  
  # Line Plot of the posterior probabilities
  melt_model <- melt(model[,c(xval, paste0("S", uniq_map[,est_cols[1]]), 'chrom')], 
                     id = c(xval, 'chrom'))
  gpp <- (ggplot(melt_model, aes_string(x = xval, y = 'value', col = 'variable')) + 
            geom_line() +
            geom_vline(xintercept = idxs$chrs, size=1.5, colour='grey') +
            theme_classic() +
            lims(y=c(0,1)) +
            scale_color_manual(values = mycols, name = "State:", labels = uniq_map[,est_cols[1]]) +
            theme(axis.ticks = element_blank(), axis.text.y = element_blank()) + 
            labs(y = "Posterior Prob.") +
            scale_x_continuous(breaks=idxs$mids,
                               labels=names(idxs$mids))) %>%
    ggplotGrob()
  
  gobs$widths <- gest$widths
  if(known_actual){
    grid.arrange(gobs, gact, gest, gpp, widths = 1, nrow = 4)
  } else {
    grid.arrange(gobs, gest, gpp, widths = 1, nrow = 3)
  }
}


# Identify the chromosome transition points
.getChrIdx <- function(model){
  chr_idx <- cumsum(rle(model$chrom)$lengths)
  mid_chr <- ceiling(chr_idx - diff(c(1,chr_idx)) /2)
  mid_chr <- setNames(mid_chr, gsub("^chr([om]*)?", "", rle(model$chrom)$values))
  list("chrs"=chr_idx, "mids"=mid_chr)
}

