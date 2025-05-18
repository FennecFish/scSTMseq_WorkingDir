structure_plot <-
  function (scSTMobj, topics, grouping, # loadings_order = "embed", 
            n = 2000,
            colors, gap = 1,
            embed_method = structure_plot_default_embed_method,
            ggplot_call = structure_plot_ggplot_call, ...) {
    
      # check to see scSTMobj is correct
      if(!class(scSTMobj) == "STM") stop("Input \"scSTMobj\" should be an output from \"scSTMseq\".")
           
      #
      n0 <- scSTMobj$settings$dim$N # number of doc/cell
      k <- scSTMobj$settings$dim$K
      theta <- scSTMobj$theta
      colnames(theta) <- paste0("topic_",1:k)
      rownames(theta) <- scSTMobj$DocName
    # Check and process input argument "topics".
    if (missing(topics))
      topics <- order(colMeans(theta))
    if (!is.character(topics))
      topics <- colnames(theta[,topics,drop = FALSE])

    # Check and process input argument "grouping".
    if (missing(grouping))
      grouping <- factor(rep(1,n0))
    if (!is.factor(grouping))
      grouping <- as.factor(grouping)
    if (length(grouping) != n0)
      stop("Input argument \"grouping\" should be a factor with one entry ",
           "for each row of scSTMobj$theta")

    # Check and process input argument "colors". colors9 is the from
    # colorbrewer2.org (qualitative data, 9-class Set1).
    colors9 <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00",
                 "#ffff33","#a65628","#f781bf","#999999")
    if (missing(colors)) {
      if (k < 10)
        colors <- colors9
      else if (k < 22)
        colors <- kelly()[-1]
      else      
        colors <- glasbey()[-1]
    }

    if (length(colors) < k)
      stop("There must be at least as many colours as topics")
    names(colors) <- colnames(theta)
    colors <- colors[topics]

    # Check and process input arguments "loadings_order" and "n".
    # If necessary, randomly subsample the rows of L.
    if (n < n0) {
        rows <- sample(n0,n)
        theta <- theta[rows,]
        grouping <- grouping[rows,drop = FALSE]
    }
   
    # The ordering of the rows is not provided, so determine an
    # ordering by computing a 1-d embedding of L.
    if (nlevels(grouping) == 1) {
        y <- rnorm(nrow(theta))
        loadings_order <- order(y)
    } else {
        loadings_order <- NULL
        for (group in levels(grouping)) {
            i <- which(grouping == group)
            if (length(i) > 0)
                y <- rnorm(nrow(theta[i, , drop = FALSE]))
            # y <- embed_method(select_loadings(scSTMobj,i),...)
            loadings_order <- c(loadings_order,i[order(y)])
        }
    }
   
    # Prepare the data for plotting and create the structure plot.
    theta <- theta[loadings_order,]
    grouping <- grouping[loadings_order,drop = TRUE]
    if (nlevels(grouping) == 1) {
      dat <- compile_structure_plot_data(theta,topics)
      return(ggplot_call(dat,colors))
    } else {
      out <- compile_grouped_structure_plot_data(theta,topics,grouping,gap)
      return(ggplot_call(out$dat,colors,out$ticks))
    }
  }

#' @rdname structure_plot
#'
#' @importFrom stats rnorm
#' 
#' @export
#' 
# structure_plot_default_embed_method <- function (fit,...) {
#     
#   if (nrow(fit) < 20)
#     return(rnorm(nrow(fit)))
#   else {
#     d <- dim(fit$L)
#     message(sprintf("Running tsne on %s x %s matrix.",d[1],d[2]))
#     return(drop(suppressMessages(tsne_from_topics(fit,dims = 1,...))))
#   }
# }

structure_plot_ggplot_call <- function (dat, colors, ticks = NULL,
                                        font.size = 9)
  ggplot(dat,aes_string(x = "sample",y = "prop",color = "topic",
                        fill = "topic")) +
  geom_col() +
  scale_x_continuous(limits = c(0,max(dat$sample) + 1),breaks = ticks,
                     labels = names(ticks)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = "",y = "topic proportion") +
  theme_cowplot(font.size) +
  theme(axis.line   = element_blank(),
        axis.ticks  = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1))

# This is used by structure_plot to create a data frame suitable for
# plotting with 'ggplot'. Input argument L is the topic proportions
# matrix. Input argument "topics" is the vector of the selected topics
# (that is, selected columns of L). The output is a data frame with
# three columns: "sample", a row of L (numeric); "topic", a topic
# (factor); and "prop", the topic proportion for the given sample
# (numeric).
compile_structure_plot_data <- function (L, topics) {

  n <- nrow(L)
  k <- length(topics)
  dat <- data.frame(sample = rep(1:n,times = k),
                    topic  = rep(topics,each = n),
                    prop   = c(L[,topics]))
  dat$topic <- factor(dat$topic,topics)
  dat$source_Doc <- rownames(L[,topics])
  return(dat)
}

# This is used by structure_plot to create a data frame suitable for
# plotting with 'ggplot' when the data are grouped. Input argument L is
# the topic proportions matrix. Input argument "topics" is the vector
# of selected topics (that is, selected columns of L).  Input argument
# "grouping" is a factor with one entry for each row of L giving the
# assigned group. The rows of L (and, correspondingly, the grouping
# vector) should already ordered by the groups; that is, grouping =
# sort(grouping). Finally, a "gap" is added to the sample indices in
# each group to provide a visual spacing of the groups.
compile_grouped_structure_plot_data <- function (L, topics, grouping,
                                                 gap = 0) {

  ticks <- rep(0,nlevels(grouping))
  names(ticks) <- levels(grouping)
  dat <- NULL
  m <- 0
  for (j in levels(grouping)) {
    i          <- which(grouping == j)
    out        <- compile_structure_plot_data(L[i,,drop = FALSE],topics)
    out$sample <- out$sample + m
    n          <- length(i)
    dat        <- rbind(dat,out)
    ticks[j]   <- m + n/2
    m          <- m + n + gap
  }
  return(list(dat = dat,ticks = ticks))
}
