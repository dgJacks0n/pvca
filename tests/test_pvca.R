# test script for pvca package mods
# copied verbatim from example in pvca-batch man page

library(golubEsets)
library(pvcaMod) # my mod for regression testing
data(Golub_Merge)
pct_threshold <- 0.6
batch.factors <- c("ALL.AML", "BM.PB", "Source")



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

message("PVCA with interactions took ", elapsed_def)

# modified - without interactions
start_noi <- Sys.time()
pvcaObj_noi <- pvcaBatchAssess (Golub_Merge, batch.factors, pct_threshold,
                            include.interactions = FALSE)
elapsed_noi <- Sys.time() - start_noi

bp_noi <- barplot(pvcaObj_noi$dat,  xlab = "Effects",
              ylab = "Weighted average proportion variance", ylim= c(0,1.1),
              col = c("blue"), las=2, main="PVCA bar chart w/o Interactions")
axis(1, at = bp_noi, labels = pvcaObj_noi$label, xlab = "Effects", cex.axis = 0.5, las=2)
values = pvcaObj_noi$dat
new_values = round(values , 3)
text(bp_noi, pvcaObj_noi$dat, labels = new_values, pos=3, cex = 0.8) 

message("PVCA without interactions took ", elapsed_noi)

