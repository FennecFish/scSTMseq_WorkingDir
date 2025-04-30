# this script is used to plot fixed effect size and its standard error

plot.fixedEffect <- function(estobj, covariate) {
    
    fixEffect_cov <- lapply(estobj$fixedEffect, function(x) x[covariate, ])
    fixEffect <- do.call(rbind, fixEffect_cov) %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column("Topic")
    
    p <- ggplot(fixEffect, aes(x = Topic, y = Estimate)) +
        geom_point(size = 4) + 
        geom_errorbar(aes(ymin = Estimate - `Std. Error`, ymax = Estimate + `Std. Error`), width = 0.2, size = 1) + 
        theme_minimal() +
        coord_flip() + 
        labs(y = "Estimate", x = "Fixed Effect", 
             title = paste("Fixed Effects Estimates with Standard Errors for", covariate)) +
        theme(text = element_text(size = 16), # Increase text size
              axis.title = element_text(size = 18), # Increase axis title size
              axis.text = element_text(size = 14)) + # Increase axis text size
        geom_text(aes(y = Estimate, label = ifelse(`Pr(>|t|)` < 0.05, 
                                                   paste0("pValue=", sprintf("%.2f", `Pr(>|t|)`), " **"), 
                                                   paste0("pValue=", sprintf("%.2f", `Pr(>|t|)`)))), 
                  vjust = 3, size = 5) 
    return(p)
}