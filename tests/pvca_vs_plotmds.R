# Compare pvca to plotMDS

library(golubEsets)
library(pvca) # original boioconductor package
library(limma)
library(ggplot2)


data(Golub_Merge)
pct_threshold <- 0.6
batch.factors <- c("ALL.AML", "BM.PB", "Source")


message("Testing PVCA with ", length(batch.factors), " factors")




# default - with interactions
start_def <- Sys.time()
pvcaObj <- pvcaMod::pvcaBatchAssess(Golub_Merge, batch.factors, pct_threshold)
elapsed_def  <- Sys.time () - start_def

bp <- barplot(pvcaObj$dat,  xlab = "Effects",
              ylab = "Weighted average proportion variance", ylim= c(0,1.1),
              col = c("blue"), las=2, main="Default PVCA bar chart")
axis(1, at = bp, labels = pvcaObj$label, xlab = "Effects", cex.axis = 0.5, las=2)
values = pvcaObj$dat
new_values = round(values , 3)
text(bp,pvcaObj$dat,labels = new_values, pos=3, cex = 0.8) 

# message("PVCA with interactions took ", elapsed_def)

# run limma::plotMDS
mdsres <- plotMDS(Golub_Merge)

# merge mds results with pData
mdsr_df <- data.frame(x = mdsres$x, y = mdsres$y)
colnames(mdsr_df) <- make.names(paste(mdsres$axislabel, mdsres$dim.plot))

comb_res <- merge(mdsr_df, pData(Golub_Merge), by = "row.names")


# build base scatter plot
p <- ggplot(data = comb_res, 
            aes(x = Leading.logFC.dim.1, y = Leading.logFC.dim.2))

sapply(batch.factors, function(x) {
  print( p + geom_point(aes_string(color = x)) )
})
