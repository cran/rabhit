## ---- eval=TRUE, message=FALSE, warning=FALSE---------------------------------
library(rabhit)

## ---- eval=TRUE, warning=FALSE------------------------------------------------
# Load example sequence data and example germline database
data(samples_db, HVGERM, HDGERM)
# Selecting a single individual
clip_db <- samples_db[samples_db$SUBJECT=='I5', ]
# Inferring haplotype using J6 as anchor
haplo_db_J6 <- createFullHaplotype(clip_db, toHap_col=c("V_CALL","D_CALL"),
                                hapBy_col="J_CALL", hapBy="IGHJ6", 
                                toHap_GERM=c(HVGERM, HDGERM))

## ---- eval=TRUE, warning=FALSE------------------------------------------------
head(haplo_db_J6,3)

## ---- eval=TRUE, warning=FALSE,fig.width=15,fig.height=10---------------------
# Plotting the haplotype map
plotHaplotype(haplo_db_J6)

## ---- eval=TRUE, warning=FALSE, fig.width=12, fig.height=8--------------------
# Inferring haplotype using D2-21 as anchor
haplo_db_D2_21 <- createFullHaplotype(clip_db, toHap_col="V_CALL",
                        hapBy_col="D_CALL", hapBy="IGHD2-21", toHap_GERM=HVGERM)

## ---- eval=TRUE, warning=FALSE------------------------------------------------
haplo_db_J6[haplo_db_J6$GENE == "IGHD2-21", ]

## ---- eval=TRUE, warning=FALSE------------------------------------------------
# rename the subject
haplo_db_J6$SUBJECT <- 'J6'
haplo_db_D2_21$SUBJECT <- 'D2-21'

# change the anchor gene columns
# For D2-21*01 and J6*03 we will change the column to Anchor_J03_D01
names(haplo_db_J6)[which(names(haplo_db_J6)=='IGHJ6_03')] <- "AnchorJ03D01"
names(haplo_db_D2_21)[which(names(haplo_db_D2_21)=='IGHD2-21_01')] <- "AnchorJ03D01"

# For D2-21*02 and J6*02 we will change the column to Anchor_J02_D02
names(haplo_db_J6)[which(names(haplo_db_J6)=='IGHJ6_02')] <- "AnchorJ02D02"
names(haplo_db_D2_21)[which(names(haplo_db_D2_21)=='IGHD2-21_02')] <- "AnchorJ02D02"

# Subseting the haplo_db_J6 dataset to include only the V genes
haplo_db_J6 <- haplo_db_J6[grep('IGHV',haplo_db_J6$GENE),]

# Combining the datasets rowwise
haplo_comb <- rbind(haplo_db_J6,haplo_db_D2_21)

## ---- eval=TRUE, warning=FALSE, fig.width=14, fig.height=6--------------------
# Plot the haplotype inferred dendogram
hapDendo(haplo_comb)

## ---- cache = T, eval=TRUE, warning=FALSE, fig.width=12, fig.height=12--------
# Removing the individual I5_FR1 with the paritial V coverage sequence.
clip_dbs <- samples_db[samples_db$SUBJECT!='I5_FR2', ]

# Inferred haplotype summary table
haplo_db <- createFullHaplotype(clip_dbs, toHap_col=c("V_CALL","D_CALL"),
                        hapBy_col="J_CALL", hapBy="IGHJ6", toHap_GERM=c(HVGERM, HDGERM))
# Plot the haplotype inferred heatmap
p.list <- hapHeatmap(haplo_db)

# The function return a list with the plot and the optimal width and height
# we can use both parametrs to render the plot to the desired size.
width <- p.list$width
height <- p.list$height


## ---- echo = F, cache = T, eval=TRUE, warning=FALSE, fig.width=width, fig.height=height----
# Plotting the heatmap
p.list$p

## ---- eval=TRUE, warning=FALSE------------------------------------------------
# Infering double choromosome deletions
del_binom_db <- deletionsByBinom(clip_dbs)
head(del_binom_db)

## ---- eval=TRUE, warning=FALSE,fig.height=9,fig.width=15----------------------
# Don't plot IGHJ
del_binom_db <- del_binom_db[grep('IGHJ', del_binom_db$GENE, invert = T),]
# Inferred deletion summary table
plotDeletionsByBinom(del_binom_db) 

## ---- eval=TRUE, warning=FALSE, cache=TRUE------------------------------------
# Selecting a single individual with partial V coverage
clip_db <- samples_db[samples_db$SUBJECT=='I5_FR2', ]

# Detecting non reliable genes
nonReliable_Vgenes <- nonReliableVGenes(clip_db)

## ---- eval=TRUE, warning=FALSE, cache=TRUE------------------------------------
# Inferred deletion summary table
del_binom_db <- deletionsByBinom(clip_db, chain = "IGH", 
                                 nonReliable_Vgenes = nonReliable_Vgenes)

## ---- eval=TRUE, warning=FALSE------------------------------------------------
# Using the deleted_genes and nonRelaible_Vgenes flags 
# to infer haplotype for a partial V coverage sequence dataset
haplo_db <- createFullHaplotype(clip_db, toHap_col=c("V_CALL","D_CALL"),
                        hapBy_col="J_CALL", hapBy="IGHJ6", 
                        toHap_GERM=c(HVGERM, HDGERM),
                        deleted_genes = del_binom_db, 
                        nonReliable_Vgenes = nonReliable_Vgenes)

## ---- eval=TRUE, warning=FALSE------------------------------------------------
# Generate interactive haplotype plot
p <- plotHaplotype(haplo_db, html_output = TRUE)

## ---- cache=TRUE, eval=FALSE, warning=FALSE-----------------------------------
#  # Saving the plot to html output
#  htmlwidgets::saveWidget(p, "haplotype.html", selfcontained = T)

## ---- cache=TRUE, eval=TRUE, warning=FALSE, message=F-------------------------
# Plotting the interactive haplotype inference
p

## ---- eval=TRUE, warning=FALSE,fig.height=12,fig.width=15---------------------
# Detecting non reliable genes
nonReliable_Vgenes <- nonReliableVGenes(samples_db)

# Infering double chromosome deletion
del_binom_db <- deletionsByBinom(samples_db, nonReliable_Vgenes = nonReliable_Vgenes)

# Inferred haplotype summary table for multiple subjects
haplo_db <- createFullHaplotype(samples_db, toHap_col=c("V_CALL","D_CALL"), 
                                hapBy_col="J_CALL", hapBy="IGHJ6", 
                                toHap_GERM=c(HVGERM, HDGERM), 
                                deleted_genes = del_binom_db, 
                                nonReliable_Vgenes = nonReliable_Vgenes)
# plot deletion heatmap
deletionHeatmap(haplo_db)

## ---- eval=TRUE, warning=FALSE, cache=T---------------------------------------
# Inferred deletion summary table
del_db <- deletionsByVpooled(samples_db, nonReliable_Vgenes = nonReliable_Vgenes)
head(del_db)

## ---- eval=TRUE, warning=FALSE,fig.height=4,fig.width=8-----------------------
# Plot the deletion heatmap
plotDeletionsByVpooled(del_db)

