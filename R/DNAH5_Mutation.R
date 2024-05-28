#      Created 04/11/2024             #
#                                     #
#      Author: Nicholas Hadas         #
#   DNAH5 JCI - Original Object       #


####-------Install Packages-------####
install.packages("ggraph")
install.packages("igraph")
install.packages('BiocManager')
install.packages("clustree")
BiocManager::install('multtest')
install.packages('SingleCellExperiment')
install.packages(c('dplyr', 'Seurat', 'patchwork', 'devtools', 'tidyverse', 
                   'gridExtra', 'harmony', 'metap', 'pals'))
devtools::install_github('satijalab/seurat-data')
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")
install.packages('hdf5r')
BiocManager::install("zellkonverter")

library(DoubletFinder)
library(Seurat)
library(SeuratData)
library(patchwork)
library(devtools)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(harmony)
library(sctransform)
library(BiocManager)
library(multtest)
library(metap)
library(ggraph)
library(clustree)
library(pals)
library(hdf5r)
library(SeuratDisk)
library(zellkonverter)
library(glue)


####-------Set a Working Directory-------####
setwd('/Users/labuser/Documents/DNAH5_JCI')


####-------Load Raw Data-------####
object <- readRDS("Objects/Rpca_13_SINGLETS_final.rds")


####-------Subset Out Mother-------####
print(unique(object@meta.data$stim2)) # View unique file names

# Subset #
object_subset <- subset(object, subset = stim2 %in% c("HNEC102", "HNEC103", "HNEC104", 
                                             "HNEC105", "HNEC106", "PCD111", 
                                             "PCD117", "PCD129", "PCD138"))

# Save #
saveRDS(object_subset, "Rpca_13_SINGLETS_final_subset.rds")


####-------Label Clusters-------####

# Set Idents to celltypefinal column #
Idents(object_subset) <- object_subset@meta.data$celltypefinal

# Rename Labels to Major Cell Groups #
object_subset <- RenameIdents(object_subset, 'Cil1' = "Ciliated")

object_subset <- RenameIdents(object_subset, 'Cil2' = "Ciliated")

object_subset <- RenameIdents(object_subset, 'Cil3' = "Ciliated")

object_subset <- RenameIdents(object_subset, 'Cil4' = "Ciliated")

object_subset <- RenameIdents(object_subset, 'Cil5' = "Ciliated")

object_subset <- RenameIdents(object_subset, 'Bas1' = "Basal")

object_subset <- RenameIdents(object_subset, 'Bas2' = "Basal")

object_subset <- RenameIdents(object_subset, 'Bas3' = "Basal")

object_subset <- RenameIdents(object_subset, 'Bas4' = "Basal")

object_subset <- RenameIdents(object_subset, 'Bas5' = "Basal")

object_subset <- RenameIdents(object_subset, 'Bas6' = "Basal")

object_subset <- RenameIdents(object_subset, 'Sec1' = "Secretory")

object_subset <- RenameIdents(object_subset, 'Sec2' = "Secretory")

object_subset <- RenameIdents(object_subset, 'Sec3' = "Secretory")

object_subset <- RenameIdents(object_subset, 'Sec4' = "Secretory")

object_subset <- RenameIdents(object_subset, 'Sec5' = "Secretory")

object_subset <- RenameIdents(object_subset, 'Iono' = "Ionocyte")

object_subset$final_cluster <- Idents(object_subset)

# Save Object #

saveRDS(object_subset, "Orig_Object_Clustered.rds")

####-------DE Analysis (Control vs PCD) - Set-Up-------####
DefaultAssay(object_subset)

# Add Column Discerning Control v PCD #
object_subset$stim4 <- paste(object_subset$final_cluster, object_subset$stim, sep='_')

# Set Idents #
Idents(object_subset) <- object_subset$stim4

# Create List of Cell Types #
ctrl_list <- levels(object_subset)[1:5]
ctrl_list <- ctrl_list[order(names(setNames(ctrl_list, ctrl_list)))]

pcd_list <- levels(object_subset)[6:10]
pcd_list <- pcd_list[order(names(setNames(pcd_list, pcd_list)))]


####-------DE Analysis (Control vs PCD) - RNA Assay-------####

# Run FindMarkers #
for (i in 1:5) {
  name = FindMarkers(object_subset, ident.1 = pcd_list[[i]], ident.2 = ctrl_list[[i]], verbose = TRUE)
  file = paste0(unlist((strsplit(ctrl_list[[i]], split = "_")))[1], "_RNA_", "orig", ".csv")
  write.csv(name, file)
}


####-------DE Analysis (Control vs PCD) - DESEQ2-------####

# Install #
BiocManager::install("DESeq2")

# Copy #
object_subset_DS2 <- object_subset

# Run FindMarkers #
for (i in 1:5) {
  name = FindMarkers(object_subset_DS2, ident.1 = pcd_list[[i]], ident.2 = ctrl_list[[i]], verbose = TRUE, test.use = "DESeq2")
  file = paste0(unlist((strsplit(ctrl_list[[i]], split = "_")))[1], "_DS2_", "orig", ".csv")
  write.csv(name, file)
}


####-------DE Analysis (Control vs PCD) - MAST-------####

# Install #
BiocManager::install("MAST")

# Run FindMarkers #
for (i in 1:5) {
  name = FindMarkers(object_subset, ident.1 = pcd_list[[i]], ident.2 = ctrl_list[[i]], verbose = TRUE, test.use = "MAST")
  file = paste0(unlist((strsplit(ctrl_list[[i]], split = "_")))[1], "_MAST_", "orig", ".csv")
  write.csv(name, file)
}


####-------DE Analysis (All) - Set-Up-------####

# Create List of Cell Types #
cell.types <- c("Basal", "Ciliated", "Div", "Ionocyte", "Secretory")

# Set Idents to Normal Cell Types #
Idents(object_subset) <- object_subset$final_cluster

####-------DE Analysis (All) - RNA-------####

# Run FindMarkers #
for  (i in 1:5) {
  name = FindMarkers(object_subset, ident.1 = cell.types[[i]], verbose = TRUE)
  file = paste0(unlist((strsplit(cell.types[[i]], split = "_")))[1], "_RNA_", "orig_all", ".csv")
  write.csv(name, file)
}

####-------DE Analysis (All) - DESEQ2-------####

# Copy #
object_subset_DS2 <- object_subset

# Run FindMarkers #
for  (i in 4:5) {
  name = FindMarkers(object_subset_DS2, ident.1 = cell.types[[i]], verbose = TRUE, test.use = "DESeq2")
  file = paste0(unlist((strsplit(cell.types[[i]], split = "_")))[1], "_DS2_", "orig_all", ".csv")
  write.csv(name, file)
} 


####-------DE Analysis (All) - MAST-------####

# Run FindMarkers #
for  (i in 1:5) {
  name = FindMarkers(object_subset, ident.1 = cell.types[i], verbose = TRUE, test.use = "MAST")
  file = paste0(unlist((strsplit(cell.types[[i]], split = "_")))[1], "_MAST_", "orig_all", ".csv")
  write.csv(name, file)
}

####-------PseudoBulk Analysis - DESeq2 - CTRL v PCD-------####

# Create PseudoBulk Object #
pseudo_object_subset <- AggregateExpression(object_subset, assays = "RNA", return.seurat = T, group.by = c("final_cluster", "stim2", "stim"))

# View PseudoBulk Profile #
tail(Cells(pseudo_object_subset))

# Create New PseudoBulk Ident column and set to Idents #
pseudo_object_subset$celltype.stim <- paste(pseudo_object_subset$final_cluster, pseudo_object_subset$stim, sep = "_")
Idents(pseudo_object_subset) <- "celltype.stim"

 # Run FindMarkers #
for (i in 1:5) {
  name = FindMarkers(pseudo_object_subset, ident.1 = pcd_list[[i]], ident.2 = ctrl_list[[i]], verbose = TRUE, test.use = "DESeq2")
  file = paste0(unlist((strsplit(ctrl_list[[i]], split = "_")))[1], "_DS2_", "orig_", "pseudo", ".csv")
  write.csv(name, file)
}

####-------PseudoBulk Analysis - DESeq2 - All-------####

# Create PseudoBulk Object #
pseudo_object_subset_all <- AggregateExpression(object_subset, assays = "RNA", return.seurat = T, group.by = c("final_cluster", "stim2"))

# View PseudoBulk Profile #
tail(Cells(pseudo_object_subset_all))

# Create New PseudoBulk Ident column and set to Idents #
pseudo_object_subset_all$celltype <- paste(pseudo_object_subset_all$final_cluster)
Idents(pseudo_object_subset_all) <- "celltype"

# Run FindMarkers #
for (i in 1:5) {
  name = FindMarkers(pseudo_object_subset_all, ident.1 =cell.types[i], verbose = TRUE, test.use = "DESeq2")
  file = paste0(unlist((strsplit(ctrl_list[[i]], split = "_")))[1], "_DS2_", "orig_", "pseudo", ".csv")
  write.csv(name, file)
}


####-------MAST with Random Effects-------####

# Install MAST and SCE #
BiocManager::install("MAST")
BiocManager::install("SingleCellExperiment")
library(MAST)
library(SingleCellExperiment)
library(lme4)

# Preview ColData #
head(colData(object_subset_sce))

# Convert Seurat Object to Single Cell Experiment #
object_subset_sce <- as.SingleCellExperiment(object_subset, assay = 'RNA')
object_subset_sca <- SceToSingleCellAssay(object_subset_sce)


##---Ciliated---##

# Make Ciliated Subset #
object_subset_ciliated <- subset(object_subset, subset = final_cluster %in% "Ciliated")

# Convert Seurat Object to Single Cell Experiment #
object_subset_ciliated_sce <- as.SingleCellExperiment(object_subset_ciliated, assay = 'RNA')
object_subset_ciliated_sca <- SceToSingleCellAssay(object_subset_ciliated_sce)

# Run Mast with random effects #
results1 <- zlm(~  stim + (1|stim2), object_subset_sca, method = 'glmer', ebayes = F, strictConvergence = FALSE)

# Filter Datatable #
ciliated_mast_re <- na.omit(summary(results1)$datatable[contrast=='stimPCD' & component=='logFC' ,. (primerid, contrast, component, z)])

# Order #
ciliated_mast_re <- ciliated_mast_re[order(-z)]

# Save as .csv # 
write.csv(ciliated_mast_re, "ciliated_mast_re.csv")


##---Secretory---##

# Make Secretory Subset #
object_subset_secretory <- subset(object_subset, subset = final_cluster %in% "Secretory")

# Convert Seurat Object to Single Cell Experiment #
object_subset_secretory_sce <- as.SingleCellExperiment(object_subset_secretory, assay = 'RNA')
object_subset_secretory_sca <- SceToSingleCellAssay(object_subset_secretory_sce)

# Run Mast with random effects #
results2 <- zlm(~  stim + (1|stim2), object_subset_secretory_sca, method = 'glmer', ebayes = F, strictConvergence = FALSE)

# Filter Datatable #
secretory_mast_re <- na.omit(summary(results2)$datatable[contrast=='stimPCD' & component=='logFC' ,. (primerid, contrast, component, z)])

# Order #
secretory_mast_re <- secretory_mast_re[order(-z)]

# Save as .csv # 
write.csv(secretory_mast_re, "secretory_mast_re.csv")

##---Basal---##

# Make Basal Subset #
object_subset_basal <- subset(object_subset, subset = final_cluster %in% "Basal")

# Convert Seurat Object to Single Cell Experiment #
object_subset_basal_sce <- as.SingleCellExperiment(object_subset_basal, assay = 'RNA')
object_subset_basal_sca <- SceToSingleCellAssay(object_subset_basal_sce)

# Run Mast with random effects #
results3 <- zlm(~  stim + (1|stim2), object_subset_basal_sca, method = 'glmer', ebayes = F, strictConvergence = FALSE)

# Filter Datatable #
basal_mast_re <- na.omit(summary(results3)$datatable[contrast=='stimPCD' & component=='logFC' ,. (primerid, contrast, component, z)])

# Order #
basal_mast_re <- basal_mast_re[order(-z)]

# Save as .csv # 
write.csv(basal_mast_re, "basal_mast_re.csv")


####-------MAST with Random Effects- Seurat Implementation-------####

# Load Object #
object <- readRDS("Objects/Orig_Object_Clustered.rds")

# Install MAST #
BiocManager::install("MAST")

# Add Column Discerning Control v PCD #
object$stim4 <- paste(object$final_cluster, object$stim, sep='_')

# Set Idents #
Idents(object) <- object$stim4

# Create List of Cell Types #
ctrl_list <- levels(object)[1:5]
ctrl_list <- ctrl_list[-c(4, 5)]
ctrl_list <- ctrl_list[order(names(setNames(ctrl_list, ctrl_list)))]

pcd_list <- levels(object)[6:10]
pcd_list <- pcd_list[-c(4, 5)]
pcd_list <- pcd_list[order(names(setNames(pcd_list, pcd_list)))]


# Change FindMarkers Code to Incorporate MAST with RE #

# We want a different behavior from MASTDETest, where the latent variable is used as a random effect in the model #
rem_MASTDETest <- function(
    data.use,
    cells.1,
    cells.2,
    latent.vars = NULL,
    verbose = TRUE,
    # New option - random effect variable (should be included in latent.vars)
    re.var = NULL,
    ...
) {
  # Check for MAST #
  if (!PackageCheck('MAST', error = FALSE)) {
    stop("Please install MAST - learn more at https://github.com/RGLab/MAST")
  }
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  latent.vars <- latent.vars %||% group.info
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])
  latent.vars.names <- c("condition", colnames(x = latent.vars))
  latent.vars <- cbind(latent.vars, group.info)
  latent.vars$wellKey <- rownames(x = latent.vars)
  fdat <- data.frame(rownames(x = data.use))
  colnames(x = fdat)[1] <- "primerid"
  rownames(x = fdat) <- fdat[, 1]
  sca <- MAST::FromMatrix(
    exprsArray = as.matrix(x = data.use),
    check_sanity = FALSE,
    cData = latent.vars,
    fData = fdat
  )
  cond <- factor(x = SummarizedExperiment::colData(sca)$group)
  cond <- relevel(x = cond, ref = "Group1")
  SummarizedExperiment::colData(sca)$condition <- cond
  
  # This is the main change in the code - we want ~ ( 1 | re.var) in the formula:~ condition + lat.vars + (1 | re.var) #
  if (!is.null(re.var)) {
    if (!re.var %in% latent.vars.names) {
      stop("Random effect variable should be included in latent variables!")
    }
    latent.vars.names <- latent.vars.names[!latent.vars.names %in% re.var]
    fmla <- as.formula(
      object = paste0(
        " ~ ", paste(latent.vars.names, collapse = "+"), glue("+ (1|{re.var})")
      )
    )
    # print(fmla) We need glmer to make this work#
    method <-  "glmer" # trying to troubleshoot this - it can clash with the already existing method var in the function call    
    zlmCond <- MAST::zlm(formula = fmla, sca = sca, method = "glmer", ...)
  } else {
    # Original code #
    fmla <- as.formula(
      object = paste0(" ~ ", paste(latent.vars.names, collapse = "+"))
    )
    zlmCond <- MAST::zlm(formula = fmla, sca = sca, ...)
  }
  summaryCond <- MAST::summary(object = zlmCond, doLRT = 'conditionGroup2')
  summaryDt <- summaryCond$datatable
  
  # The output format is slightly different, so we need adapt the code #
  if(!is.null(re.var)) {
    p_val <- summaryDt[summaryDt$"component" == "H", 4]$`Pr(>Chisq)`
    genes.return <- summaryDt[summaryDt$"component" == "H", 1]$primerid
  } else {
    p_val <- summaryDt[summaryDt[, "component"] == "H", 4]
    genes.return <- summaryDt[summaryDt[, "component"] == "H", 1]
  }
  
  to.return <- data.frame(p_val, row.names = genes.return)
  return(to.return)
}

# This is the zlm function with a minor change to fix a buggy variable being passed possibly due to namespace clash #
my_zlm <- function (formula, sca, method = "bayesglm", silent = TRUE, ebayes = TRUE, 
                    ebayesControl = NULL, force = FALSE, hook = NULL, parallel = TRUE, 
                    LMlike, onlyCoef = FALSE, exprs_values = MAST::assay_idx(sca)$aidx, 
                    ...) 
{
  dotsdata = list(...)$data
  if (!is.null(dotsdata)) {
    if (!missing(sca)) 
      stop("Cannot provide both `sca` and `data`")
    sca = dotsdata
  }
  if (!inherits(sca, "SingleCellAssay")) {
    if (inherits(sca, "data.frame")) {
      if (!is.null(dotsdata)) {
        return(.zlm(formula, method = method, silent = silent, 
                    ...))
      }
      else {
        return(.zlm(formula, data = sca, method = method, 
                    silent = silent, ...))
      }
    }
    else {
      stop("`sca` must inherit from `data.frame` or `SingleCellAssay`")
    }
  }
  if (missing(LMlike)) {
    method <- match.arg(method, MAST:::methodDict[, keyword])
    method <- MAST:::methodDict[keyword == method, lmMethod]
    if (!is(sca, "SingleCellAssay")) 
      stop("'sca' must be (or inherit) 'SingleCellAssay'")
    if (!is(formula, "formula")) 
      stop("'formula' must be class 'formula'")
    Formula <- MAST:::removeResponse(formula)
    priorVar <- 1
    priorDOF <- 0
    if (ebayes) {
      # if (!methodDict[lmMethod == method, implementsEbayes]) 
      if (!all(MAST:::methodDict[lmMethod == method, implementsEbayes]))
        stop("Method", method, " does not implement empirical bayes variance shrinkage.")
      ebparm <- MAST:::ebayes(t(SummarizedExperiment:::assay(sca, exprs_values)), ebayesControl, 
                              model.matrix(Formula, SummarizedExperiment::colData(sca)))
      priorVar <- ebparm[["v"]]
      priorDOF <- ebparm[["df"]]
      stopifnot(all(!is.na(ebparm)))
    }
    obj <- MAST:::new_with_repaired_slots(classname = method, design = SummarizedExperiment::colData(sca), 
                                          formula = Formula, priorVar = priorVar, priorDOF = priorDOF, 
                                          extra = list(...))
  }
  else {
    if (!missing(formula)) 
      warning("Ignoring formula and using model defined in 'objLMLike'")
    if (!inherits(LMlike, "LMlike")) 
      stop("'LMlike' must inherit from class 'LMlike'")
    obj <- LMlike
  }
  ee <- t(SummarizedExperiment:::assay(sca, exprs_values))
  genes <- colnames(ee)
  ng <- length(genes)
  MM <- MAST:::model.matrix(obj)
  coefNames <- colnames(MM)
  listEE <- setNames(seq_len(ng), genes)
  obj <- MAST:::fit(obj, ee[, 1], silent = silent)
  nerror <- totalerr <- 0
  pb = progress::progress_bar$new(total = ng, format = " Completed [:bar] :percent with :err failures")
  .fitGeneSet <- function(idx) {
    hookOut <- NULL
    tt <- try({
      obj <- MAST:::fit(obj, response = ee[, idx], silent = silent, 
                        quick = TRUE)
      if (!is.null(hook)) 
        hookOut <- hook(obj)
      nerror <- 0
    })
    if (is(tt, "try-error")) {
      obj@fitC <- obj@fitD <- NULL
      obj@fitted <- c(C = FALSE, D = FALSE)
      nerror <- nerror + 1
      totalerr = totalerr + 1
      if (nerror > 5 & !force) {
        stop("We seem to be having a lot of problems here...are your tests specified correctly?  \n If you're sure, set force=TRUE.", 
             tt)
      }
    }
    pb$tick(tokens = list(err = totalerr))
    if (onlyCoef) 
      return(cbind(C = coef(obj, "C"), D = coef(obj, "D")))
    summaries <- MAST:::summarize(obj)
    structure(summaries, hookOut = hookOut)
  }
  if (!parallel || getOption("mc.cores", 1L) == 1) {
    listOfSummaries <- lapply(listEE, .fitGeneSet)
  }
  else {
    listOfSummaries <- parallel::mclapply(listEE, .fitGeneSet, 
                                          mc.preschedule = TRUE, mc.silent = silent)
  }
  if (onlyCoef) {
    out <- do.call(abind, c(listOfSummaries, rev.along = 0))
    return(aperm(out, c(3, 1, 2)))
  }
  cls <- sapply(listOfSummaries, function(x) class(x))
  complain <- if (force) 
    warning
  else stop
  if (mean(cls == "try-error") > 0.5) 
    complain("Lots of errors here..something is amiss.")
  hookOut <- NULL
  if (!is.null(hook)) 
    hookOut <- lapply(listOfSummaries, attr, which = "hookOut")
  message("\nDone!")
  summaries <- MAST:::collectSummaries(listOfSummaries)
  summaries[["LMlike"]] <- obj
  summaries[["sca"]] <- sca
  summaries[["priorVar"]] <- obj@priorVar
  summaries[["priorDOF"]] <- obj@priorDOF
  summaries[["hookOut"]] <- hookOut
  summaries[["exprs_values"]] <- exprs_values
  summaries[["Class"]] <- "ZlmFit"
  zfit <- do.call(new, as.list(summaries))
  zfit
}

# We will replace the original function with ours inside the Seurat namespace #
assignInNamespace('MASTDETest', rem_MASTDETest, asNamespace("Seurat"))
# getFromNamespace("MASTDETest", "Seurat")

assignInNamespace('zlm', my_zlm, asNamespace("MAST"))
# getFromNamespace("zlm", "MAST")

# Run FindMarkers #
for (i in 1:3) {
  name = FindMarkers(object, ident.1 = pcd_list[[i]], ident.2 = ctrl_list[[i]], only.pos = TRUE, test.use = "MAST", latent.vars = "stim2", re.var = "stim2", ebayes = FALSE, min.pct = 0.1, verbose = TRUE)
  file = paste0(unlist((strsplit(ctrl_list[[i]], split = "_")))[1], "_MAST_", "RE", ".csv")
  write.csv(name, file)
}


####-------GSTA2 Expression Patterns-------####

# Load Raw Data #
object <- readRDS("Objects/Rpca_13_SINGLETS_final.rds")

# Subset Out Mother #
object_subset <- subset(object, subset = stim2 %in% c("HNEC102", "HNEC103", "HNEC104", 
                                                      "HNEC105", "HNEC106", "PCD111", 
                                                      "PCD117", "PCD129", "PCD138"))

# Set Idents to celltypefinal column #
Idents(object_subset) <- object_subset@meta.data$celltypefinal

# Rename Labels to Major Cell Groups #
object_subset <- RenameIdents(object_subset, 'Cil1' = "Ciliated 1")

object_subset <- RenameIdents(object_subset, 'Cil2' = "Ciliated 2")

object_subset <- RenameIdents(object_subset, 'Cil3' = "Ciliated 3")

object_subset <- RenameIdents(object_subset, 'Cil4' = "Ciliated 4")

object_subset <- RenameIdents(object_subset, 'Cil5' = "Ciliated 5")

object_subset <- RenameIdents(object_subset, 'Bas1' = "Basal 1")

object_subset <- RenameIdents(object_subset, 'Bas2' = "Basal 2")

object_subset <- RenameIdents(object_subset, 'Bas3' = "Basal 3")

object_subset <- RenameIdents(object_subset, 'Bas4' = "Basal 4")

object_subset <- RenameIdents(object_subset, 'Bas5' = "Basal 5")

object_subset <- RenameIdents(object_subset, 'Bas6' = "Basal 6")

object_subset <- RenameIdents(object_subset, 'Sec1' = "Secretory 1")

object_subset <- RenameIdents(object_subset, 'Sec2' = "Secretory 2")

object_subset <- RenameIdents(object_subset, 'Sec3' = "Secretory 3")

object_subset <- RenameIdents(object_subset, 'Sec4' = "Secretory 4")

object_subset <- RenameIdents(object_subset, 'Sec5' = "Secretory 5")

object_subset <- RenameIdents(object_subset, 'Iono' = "Ionocyte")

object_subset$final_cluster <- Idents(object_subset)


# GSTA2 Feature Plot #
FeaturePlot(object_subset, "GSTA2",label = TRUE, label.size = 5)

# FOXJ1 Feature Plot #
FeaturePlot(object_subset, "FOXJ1",label = TRUE, label.size = 5)


####-------MAST with Random Effects-------####

# Install MAST and SCE #
BiocManager::install("MAST")
BiocManager::install("SingleCellExperiment")
library(MAST)
library(SingleCellExperiment)
library(lme4)


# Convert Seurat Object to Single Cell Experiment #
object_subset_sce <- as.SingleCellExperiment(object_subset, assay = 'RNA')

# Preview ColData #
head(colData(object_subset_sce))

# Convert to Single Cell Assay #
object_subset_sca <- SceToSingleCellAssay(object_subset_sce)


##---Ciliated---##

# Make Ciliated Subset #
object_subset_ciliated <- subset(object_subset, subset = final_cluster %in% "Ciliated")

# Convert Seurat Object to Single Cell Assay #
object_subset_ciliated_sce <- as.SingleCellExperiment(object_subset_ciliated, assay = 'RNA')
object_subset_ciliated_sca <- SceToSingleCellAssay(object_subset_ciliated_sce)

# Filter Gene Expression #
object_subset_ciliated_sca <- filterLowExpressedGenes(object_subset_ciliated_sca, threshold = 0.1)

# Run Mast with random effects #
results1 <- zlm(~  stim + (1|stim2), object_subset_sca, method = 'glmer', ebayes = F, strictConvergence = FALSE)

# Filter Datatable #
ciliated_mast_re <- na.omit(summary(results1)$datatable[contrast=='stimPCD' & component=='logFC' ,. (primerid, contrast, component, z)])

# Order #
ciliated_mast_re <- ciliated_mast_re[order(-z)]

# Save as .csv # 
write.csv(ciliated_mast_re, "ciliated_mast_re.csv")


##---Secretory---##

# Make Secretory Subset #
object_subset_secretory <- subset(object_subset, subset = final_cluster %in% "Secretory")

# Convert Seurat Object to Single Cell Experiment #
object_subset_secretory_sce <- as.SingleCellExperiment(object_subset_secretory, assay = 'RNA')
object_subset_secretory_sca <- SceToSingleCellAssay(object_subset_secretory_sce)

# Filter Gene Expression #
object_subset_secretory_sca <- filterLowExpressedGenes(object_subset_secretory_sca, threshold = 0.1)

# Run Mast with random effects #
results2 <- zlm(~  stim + (1|stim2), object_subset_secretory_sca, method = 'glmer', ebayes = F, strictConvergence = FALSE)

# Filter Datatable #
secretory_mast_re <- na.omit(summary(results2)$datatable[contrast=='stimPCD' & component=='logFC' ,. (primerid, contrast, component, z)])

# Order #
secretory_mast_re <- secretory_mast_re[order(-z)]

# Save as .csv # 
write.csv(secretory_mast_re, "secretory_mast_re.csv")


##---Basal---##

# Make Basal Subset #
object_subset_basal <- subset(object_subset, subset = final_cluster %in% "Basal")

# Convert Seurat Object to Single Cell Experiment #
object_subset_basal_sce <- as.SingleCellExperiment(object_subset_basal, assay = 'RNA')
object_subset_basal_sca <- SceToSingleCellAssay(object_subset_basal_sce)

# Filter Gene Expression #
object_subset_basal_sca <- filterLowExpressedGenes(object_subset_basal_sca, threshold = 0.1)

# Run Mast with random effects #
results3 <- zlm(~  stim + (1|stim2), object_subset_basal_sca, method = 'glmer', ebayes = F, strictConvergence = FALSE)

# Filter Datatable #
basal_mast_re <- na.omit(summary(results3)$datatable[contrast=='stimPCD' & component=='logFC' ,. (primerid, contrast, component, z)])

# Order #
basal_mast_re <- basal_mast_re[order(-z)]

# Save as .csv # 
write.csv(basal_mast_re, "basal_mast_re.csv")


####-------Nebula DE Analysis-------####

# Install Nebula #
install_github("lhe17/nebula")
library(nebula)

# Subset Ciliated Cells #
object_subset_ciliated <- subset(object_subset, subset = final_cluster %in% "Ciliated")

# Convert to Nebula Object #
obj_neb <- scToNeb(object_subset_ciliated, assay = 'RNA', id = "stim2", pred = "stim")

# Create Design Matrix for Predictors #
df = model.matrix(~ stim, data = obj_neb$pred)

# Group Cells #
data_g <- group_cell(count = obj_neb$count,id = obj_neb$id, pred = df)

# Run Nebula #
nebula_ciiated = nebula(obj_neb$count, obj_neb$id, pred = df, verbose = TRUE)

nebula_ciliated_df <- nebula_ciiated$summary[c("logFC_stimPCD", "p_stimPCD", "gene")]

# Save as .csv # 
write.csv(nebula_ciliated_df, "ciliated_nebula.csv")

# Subset Secretory Cells #
object_subset_secretory <- subset(object_subset, subset = final_cluster %in% "Secretory")

# Convert to Nebula Object #
obj_neb <- scToNeb(object_subset_secretory, assay = 'RNA', id = "stim2", pred = "stim")

# Create Design Matrix for Predictors #
df = model.matrix(~ stim, data = obj_neb$pred)

# Group Cells #
data_g <- group_cell(count = obj_neb$count,id = obj_neb$id, pred = df)

# Run Nebula #
nebula_secretory = nebula(obj_neb$count, obj_neb$id, pred = df, verbose = TRUE)

nebula_secretory_df <- nebula_secretory$summary[c("logFC_stimPCD", "p_stimPCD", "gene")]

# Save as .csv # 
write.csv(nebula_secretory_df, "secretory_nebula.csv")

# Subset Basal Cells #
object_subset_basal <- subset(object_subset, subset = final_cluster %in% "Basal")

# Convert to Nebula Object #
obj_neb <- scToNeb(object_subset_basal, assay = 'RNA', id = "stim2", pred = "stim")

# Create Design Matrix for Predictors #
df = model.matrix(~ stim, data = obj_neb$pred)

# Group Cells #
data_g <- group_cell(count = obj_neb$count,id = obj_neb$id, pred = df)

# Run Nebula #
nebula_basal = nebula(obj_neb$count, obj_neb$id, pred = df, verbose = TRUE)

nebula_basal_df <- nebula_basal$summary[c("logFC_stimPCD", "p_stimPCD", "gene")]

# Save as .csv # 
write.csv(nebula_basal_df, "basal_nebula.csv")


####-------Muscat DE Analysis-------####

# Install Muscat #
BiocManager::install("muscat")
library(muscat)

# Make Copy of Object #
object_muscat <- object_subset

# Change Column Names #
colnames(object_muscat@meta.data)[6] <- "sample_id"
colnames(object_muscat@meta.data)[4] <- "group_id"
colnames(object_muscat@meta.data)[17] <- "cluster_id"

# Convert to Single Cell Experiment #
object_muscat_sce <- as.SingleCellExperiment(object_muscat, assay = 'RNA')

# Preview ColData #
head(colData(object_muscat_sce))

# Remove Lowly Expressed Genes #
dim(object_muscat_sce)
object_muscat_sce <- object_muscat_sce[rowSums(counts(object_muscat_sce) > 1) >= 10, ]
dim(object_muscat_sce)

# Make "group_id" Column Factorable #
object_muscat_sce$group_id <- as.factor(object_muscat_sce$group_id)

# Run Muscat #
mm <- mmDS(object_muscat_sce, method = "dream")

# Save Ciliated Dataframe #
write.csv(mm$Ciliated, "Ciliated_Muscat.csv")

# Save Basal Dataframe #
write.csv(mm$Basal, "Basal_Muscat.csv")

# Save Secretory Dataframe #
write.csv(mm$Secretory, "Secretory_Muscat.csv")


####-------Pseudobulk BoxPlot-------####

# Install SCpubr #
install.packages("SCpubr", dependencies = TRUE)
library(SCpubr)
library(ggsignif)

# Set Idents #
Idents(object_subset) <- object_subset$stim

# Subset Ciliated Cells #
object_subset_ciliated <- subset(object_subset, subset = final_cluster %in% "Ciliated")

# Run BoxPlot #
plot1 <- do_BoxPlot(object_subset_ciliated, "GSTA2", font.size = 11, outlier.alpha = 0, plot.title = "GSTA2", legend.position = "right", legend.title = "Condition", xlab = "Condition", ylab = "Expression Level", map_signif_level = TRUE, plot.grid = FALSE, use_test = TRUE, test = "wilcox.test", comparisons = list(c("CTRL", "PCD"))) +
  theme(plot.title = element_text(hjust = 0.39))

plot2 <- do_BoxPlot(object_subset_ciliated, "GSTA1", font.size = 11, outlier.alpha = 0, plot.title = "GSTA1", legend.position = "right", legend.title = "Condition", xlab = "Condition", ylab = "Expression Level", map_signif_level = TRUE, plot.grid = FALSE, use_test = TRUE, test = "wilcox.test", comparisons = list(c("CTRL", "PCD"))) +
  theme(plot.title = element_text(hjust = 0.39))

plot3 <- do_BoxPlot(object_subset_ciliated, "ALDH3A1", font.size = 11, outlier.alpha = 0, plot.title = "ALDH3A1", legend.position = "right", legend.title = "Condition", xlab = "Condition", ylab = "Expression Level", map_signif_level = TRUE, plot.grid = FALSE, use_test = TRUE, test = "wilcox.test", comparisons = list(c("CTRL", "PCD"))) +
  theme(plot.title = element_text(hjust = 0.39))

plot4 <- do_BoxPlot(object_subset_ciliated, "CBR1", font.size = 11, outlier.alpha = 0, plot.title = "CBR1", legend.position = "right", legend.title = "Condition", xlab = "Condition", ylab = "Expression Level", map_signif_level = TRUE, plot.grid = FALSE, use_test = TRUE, test = "wilcox.test", comparisons = list(c("CTRL", "PCD"))) +
  theme(plot.title = element_text(hjust = 0.39))

plot5 <- do_BoxPlot(object_subset_ciliated, "ALDH3B1", font.size = 11, outlier.alpha = 0, plot.title = "ALDH3B1", legend.position = "right", legend.title = "Condition", xlab = "Condition", ylab = "Expression Level", map_signif_level = TRUE, plot.grid = FALSE, use_test = TRUE, test = "wilcox.test", comparisons = list(c("CTRL", "PCD"))) +
  theme(plot.title = element_text(hjust = 0.39))

plot6 <- do_BoxPlot(object_subset_ciliated, "UCP2", font.size = 11, outlier.alpha = 0, plot.title = "UCP2", legend.position = "right", legend.title = "Condition", xlab = "Condition", ylab = "Expression Level", map_signif_level = TRUE, plot.grid = FALSE, use_test = TRUE, test = "wilcox.test", comparisons = list(c("CTRL", "PCD"))) +
  theme(plot.title = element_text(hjust = 0.39))

plot7 <- do_BoxPlot(object_subset_ciliated, "SLC7A2", font.size = 11, outlier.alpha = 0, plot.title = "SLC7A2", legend.position = "right", legend.title = "Condition", xlab = "Condition", ylab = "Expression Level", map_signif_level = TRUE, plot.grid = FALSE, use_test = TRUE, test = "wilcox.test", comparisons = list(c("CTRL", "PCD"))) +
  theme(plot.title = element_text(hjust = 0.39))

plot8 <- do_BoxPlot(object_subset_ciliated, "CDKN2A", font.size = 11, outlier.alpha = 0, plot.title = "CDKN2A", legend.position = "right", legend.title = "Condition", xlab = "Condition", ylab = "Expression Level", map_signif_level = TRUE, plot.grid = FALSE, use_test = TRUE, test = "wilcox.test", comparisons = list(c("CTRL", "PCD"))) +
  theme(plot.title = element_text(hjust = 0.39))

plot9 <- do_BoxPlot(object_subset_ciliated, "CBR3", font.size = 11, outlier.alpha = 0, plot.title = "CBR3", legend.position = "right", legend.title = "Condition", xlab = "Condition", ylab = "Expression Level", map_signif_level = TRUE, plot.grid = FALSE, use_test = TRUE, test = "wilcox.test", comparisons = list(c("CTRL", "PCD"))) +
  theme(plot.title = element_text(hjust = 0.39))

# Visualize #
plot1+plot2+plot3+plot4+plot5+plot6+plot7+plot8+plot9


####-------Bulk RNA Seq BoxPlot with P-Values-------####

# Install Packages #
library(ggplot2)

install.packages("ggpubr")
library(ggpubr)

##---GSTA2---##

# Add Expression Levels #
Expression <- c(4.928632654,4.28654947, 7.816881325, 5.518480975, 6.42110915, 5.815485347, 5.289505233, 6.033010953, 
                4.200619831,6.32007065,5.6973716,5.67821268,5.22265812,5.0514205,5.21177599,9.98914457,10.09819722,
                9.62400062,9.08375603,9.0028511,8.78810692,5.83847558,7.237522743,7.60411577, 9.4247855,9.27529608,
                9.320873663,8.131633277,8.22495267,8.152066832,5.911647249,7.041339043,7.307644461)

# Add Conditions #
Condition <- c("Control",  "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", 
               "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", 
               "Control", "Control", "PCD", "PCD", "PCD", "PCD", "PCD", "PCD", "PCD", "PCD", "PCD")

# Create DataFrame #
GSTA2_Expression <- data.frame(Expression, Condition)

# Create Comparisons #
my_comparisons <- list(c("Control", "PCD"))

# Set Condition as Factor #
GSTA2_Expression$Condition <- as.factor(GSTA2_Expression$Condition)

# Create Boxplot #
p1 <- ggboxplot(GSTA2_Expression, x = "Condition", y = "Expression", color = "black", palette = c("brown", "aquamarine4"), fill = "Condition") + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
  theme(
    legend.position = "right",
    axis.title.x = element_text(size=14, face="bold"),    
    axis.title.y = element_text(size=14, face="bold"),
    legend.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face = "bold")) + 
  ggtitle("GSTA2")

# Visualize #
p1


##---GSTA1---##

# Add Expression Levels #
Expression <- c(8.075297124,7.68539746,8.964583738,6.992218361,7.55351856,7.458259725,8.196369925,8.653491557,7.908795494,9.75772725,9.87409235,
                10.00932266,9.33069,9.1648972,9.40924988,10.38036962,10.46863151,10.1403429,10.28693478,10.57655819,10.2619898,9.20088788,9.377969359,
                9.67497817,11.1517231,11.11223193,11.00340362,9.306811028,9.33706693,9.384270968,9.256027057,9.57424925,9.758903928)

# Add Conditions #
Condition <- c("Control",  "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", 
               "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", 
               "Control", "Control", "PCD", "PCD", "PCD", "PCD", "PCD", "PCD", "PCD", "PCD", "PCD")

# Create DataFrame #
GSTA1_Expression <- data.frame(Expression, Condition)

# Create Comparisons #
my_comparisons <- list(c("Control", "PCD"))

# Set Condition as Factor #
GSTA1_Expression$Condition <- as.factor(GSTA1_Expression$Condition)

# Create Boxplot #
p2 <- ggboxplot(GSTA1_Expression, x = "Condition", y = "Expression", color = "black", palette = c("brown", "aquamarine4"), fill = "Condition") + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
  theme(
    legend.position = "right",
    axis.title.x = element_text(size=14, face="bold"),    
    axis.title.y = element_text(size=14, face="bold"),
    legend.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face = "bold")) + 
  ggtitle("GSTA1")

# Visualize #
p2


##---ALDH3A1---##

# Add Expression Levels #
Expression <- c(9.481587006,9.31849593,11.04168368,6.898535702,7.05049245,5.383743126,5.898995874,7.065630287,4.344460861,10.45233954,9.7256683,
                9.50535183,8.86183574,9.97469837,10.51338781,11.01397243,11.38797525,10.99568542,10.85288746,10.24885734,10.41356043,9.01344057,
                9.959393948,10.30712863,11.35929439,11.65265038,11.71884067,7.60209352,8.59636078,8.918207579,7.900998667,8.255399422,7.982577067)

# Add Conditions #
Condition <- c("Control",  "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", 
               "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", 
               "Control", "Control", "PCD", "PCD", "PCD", "PCD", "PCD", "PCD", "PCD", "PCD", "PCD")

# Create DataFrame #
ALDH3A1_Expression <- data.frame(Expression, Condition)

# Create Comparisons #
my_comparisons <- list(c("Control", "PCD"))

# Set Condition as Factor #
ALDH3A1_Expression$Condition <- as.factor(ALDH3A1_Expression$Condition)

# Create Boxplot #
p3 <- ggboxplot(ALDH3A1_Expression, x = "Condition", y = "Expression", color = "black", palette = c("brown", "aquamarine4"), fill = "Condition") + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
  theme(
    legend.position = "right",
    axis.title.x = element_text(size=14, face="bold"),    
    axis.title.y = element_text(size=14, face="bold"),
    legend.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face = "bold")) + 
  ggtitle("ALDH3A1")

# Visualize #
p3


##---CBR1---##

# Add Expression Levels #
Expression <- c(5.221459172,6.24258308,6.120424532,3.333773242,4.42715624,3.369508379,3.366300683,3.177368246,3.161931534,6.61739182,7.00894951,6.86815802,
                4.82053761,6.39732003,7.12448246,7.02074656,7.20951157,6.76007504,6.49929674,4.54310333,5.32450916,6.57997796,5.992400681,6.14603355,6.56688382
                ,6.99085274,6.779838917,4.723046764,5.680679709,4.448446041)

# Add Conditions #
Condition <- c("Control",  "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", 
               "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", 
               "Control", "Control", "PCD", "PCD", "PCD", "PCD", "PCD", "PCD")

# Create DataFrame #
CBR1_Expression <- data.frame(Expression, Condition)

# Create Comparisons #
my_comparisons <- list(c("Control", "PCD"))

# Set Condition as Factor #
CBR1_Expression$Condition <- as.factor(CBR1_Expression$Condition)

# Create Boxplot #
p4 <- ggboxplot(CBR1_Expression, x = "Condition", y = "Expression", color = "black", palette = c("brown", "aquamarine4"), fill = "Condition") + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
  theme(
    legend.position = "right",
    axis.title.x = element_text(size=14, face="bold"),    
    axis.title.y = element_text(size=14, face="bold"),
    legend.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face = "bold")) + 
  ggtitle("CBR1")

# Visualize #
p4


##---ALDH3B1---##

# Add Expression Levels #
Expression <- c(6.646237864,6.70820693,7.570301601,4.086386441,4.24955175,4.034586387,3.982843895,4.791774525,3.067144175,6.9394896,6.83341596,6.76968252,6.35001567,6.98858614,7.40408239,6.63997341,6.89425643
                ,6.4890521,6.42611922,5.34906636,5.95297574,5.05273095,5.760489539,5.78210636,5.70793736,6.26150751,6.063580707,5.701647556,5.7519623,5.805630136)

# Add Conditions #
Condition <- c("Control",  "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", 
               "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", 
               "Control", "Control", "PCD", "PCD", "PCD", "PCD", "PCD", "PCD")

# Create DataFrame #
ALDH3B1_Expression <- data.frame(Expression, Condition)

# Create Comparisons #
my_comparisons <- list(c("Control", "PCD"))

# Set Condition as Factor #
ALDH3B1_Expression$Condition <- as.factor(ALDH3B1_Expression$Condition)

# Create Boxplot #
p5 <- ggboxplot(ALDH3B1_Expression, x = "Condition", y = "Expression", color = "black", palette = c("brown", "aquamarine4"), fill = "Condition") + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
  theme(
    legend.position = "right",
    axis.title.x = element_text(size=14, face="bold"),    
    axis.title.y = element_text(size=14, face="bold"),
    legend.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face = "bold")) + 
  ggtitle("ALDH3B1")

# Visualize #
p5


##---UCP2---##

# Add Expression Levels #
Expression <- c(8.487304554,8.53565379,8.27684531,8.659390661,9.161922102,8.291806305,10.30370656,10.28546926,10.30964058,9.892525785,9.948161698,10.00354351)

# Add Conditions #
Condition <- c("Control",  "Control", "Control", "Control", "Control", "Control", "PCD", "PCD", "PCD", "PCD", "PCD", "PCD")

# Create DataFrame #
UCP2_Expression <- data.frame(Expression, Condition)

# Create Comparisons #
my_comparisons <- list(c("Control", "PCD"))

# Set Condition as Factor #
UCP2_Expression$Condition <- as.factor(UCP2_Expression$Condition)

# Create Boxplot #
p6 <- ggboxplot(UCP2_Expression, x = "Condition", y = "Expression", color = "black", palette = c("brown", "aquamarine4"), fill = "Condition", outlier.shape = NA) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
  theme(
    legend.position = "right",
    axis.title.x = element_text(size=14, face="bold"),    
    axis.title.y = element_text(size=14, face="bold"),
    legend.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face = "bold")) + 
  ggtitle("UCP2")

# Visualize #
p6


##---SLC7A2---##

# Add Expression Levels #
Expression <- c(1.544699838,1.3701287,4.6762498,1.585752459,2.42281731,2.013705419,0.41982032,1.029582222,-0.985967161,4.25647734,3.39366247,3.72163622,5.22512743,
                5.03953942,4.70039861,5.83717064,5.97248537,5.7170515,4.96024941,4.38754051,5.07810755,4.48122783,5.613827848,5.44462914,4.75079575,4.98346621,4.53101556
                ,3.237168483,3.63924389,2.614513692,1.558872832,2.157030515,2.502310068)

# Add Conditions #
Condition <- c("Control",  "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", 
               "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", 
               "Control", "Control", "PCD", "PCD", "PCD", "PCD", "PCD", "PCD", "PCD", "PCD", "PCD")

# Create DataFrame #
SLC7A2_Expression <- data.frame(Expression, Condition)

# Create Comparisons #
my_comparisons <- list(c("Control", "PCD"))

# Set Condition as Factor #
SLC7A2_Expression$Condition <- as.factor(SLC7A2_Expression$Condition)

# Create Boxplot #
p7 <- ggboxplot(SLC7A2_Expression, x = "Condition", y = "Expression", color = "black", palette = c("brown", "aquamarine4"), fill = "Condition", outlier.shape = NA) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
  theme(
    legend.position = "right",
    axis.title.x = element_text(size=14, face="bold"),    
    axis.title.y = element_text(size=14, face="bold"),
    legend.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face = "bold")) + 
  ggtitle("SLC7A2")

# Visualize #
p7


##---CDKN2A---##

# Add Expression Levels #
Expression <- c(2.799057008,2.66458924,2.514188132,2.515944021,2.38766998,2.972453735,1.665524458,1.778318225,1.569804922)

# Add Conditions #
Condition <- c("Control",  "Control", "Control", "Control", "Control", "Control", "PCD", "PCD", "PCD")

# Create DataFrame #
CDKN2A_Expression <- data.frame(Expression, Condition)

# Create Comparisons #
my_comparisons <- list(c("Control", "PCD"))

# Set Condition as Factor #
CDKN2A_Expression$Condition <- as.factor(CDKN2A_Expression$Condition)

# Create Boxplot #
p8 <- ggboxplot(CDKN2A_Expression, x = "Condition", y = "Expression", color = "black", palette = c("brown", "aquamarine4"), fill = "Condition", outlier.shape = NA) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
  theme(
    legend.position = "right",
    axis.title.x = element_text(size=14, face="bold"),    
    axis.title.y = element_text(size=14, face="bold"),
    legend.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face = "bold")) + 
  ggtitle("CDKN2A")

# Visualize #
p8


##---CBR3---##

# Add Expression Levels #
Expression <- c(-1.625225163,-1.7997963,-0.848041343,0.03215413,0.84775327,-0.241250122,0.396919918,0.418757658,0.287955561,0.18347049,1.37658896,0.37975853
                ,-1.21166121,-1.16226488,-2.29241101,1.71275766,2.09615609,1.59673818,1.57390095,-1.56820625,-2.1445776,-1.29459901,-1.509895644,-2.29394571,0.13897534
                ,2.30209923,1.870318937,1.483840881,1.10908041,0.881545277)

# Add Conditions #
Condition <- c("Control",  "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", 
               "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", 
               "Control", "Control", "PCD", "PCD", "PCD", "PCD", "PCD", "PCD")

# Create DataFrame #
CBR3_Expression <- data.frame(Expression, Condition)

# Create Comparisons #
my_comparisons <- list(c("Control", "PCD"))

# Set Condition as Factor #
CBR3_Expression$Condition <- as.factor(CBR3_Expression$Condition)

# Create Boxplot #
p9 <- ggboxplot(CBR3_Expression, x = "Condition", y = "Expression", color = "black", palette = c("brown", "aquamarine4"), fill = "Condition", outlier.shape = NA) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
  theme(
    legend.position = "right",
    axis.title.x = element_text(size=14, face="bold"),    
    axis.title.y = element_text(size=14, face="bold"),
    legend.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face = "bold")) + 
  ggtitle("CBR3")

# Visualize #
p9

# Visualize All
p1+p2+p3+p4+p5+p6+p7+p8+p9


####-------Bulk RNA Seq BoxPlot with P-Values - PCD 101-------####


##---GSTA2---##

# Add Expression Levels #
Expression <- c(4.928632654,4.28654947,7.816881325,6.32007065,5.6973716,5.67821268,5.22265812,5.0514205,5.21177599,9.98914457,10.09819722,9.62400062,9.08375603,9.0028511,8.78810692,5.83847558,7.237522743
                ,7.60411577,9.4247855,9.27529608,9.320873663)

# Add Conditions #
Condition <- c("Control",  "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", 
               "Control", "Control", "Control", "Control", "Control", "Control", "Control", "PCD", "PCD", "PCD")

# Create DataFrame #
GSTA2_Expression <- data.frame(Expression, Condition)

# Create Comparisons #
my_comparisons <- list(c("Control", "PCD"))

# Set Condition as Factor #
GSTA2_Expression$Condition <- as.factor(GSTA2_Expression$Condition)

# Create Boxplot #
p1 <- ggboxplot(GSTA2_Expression, x = "Condition", y = "Expression", color = "black", palette = c("brown", "aquamarine4"), fill = "Condition", outlier.shape = NA) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
  theme(
    legend.position = "right",
    axis.title.x = element_text(size=14, face="bold"),    
    axis.title.y = element_text(size=14, face="bold"),
    legend.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face = "bold")) + 
  ggtitle("GSTA2")

# Visualize #
p1

##---GSTA1---##

# Add Expression Levels #
Expression <- c(8.075297124,7.68539746,8.964583738,9.75772725,9.87409235,10.00932266,9.33069,9.1648972,9.40924988,10.38036962,10.46863151,10.1403429,10.28693478,10.57655819,10.2619898,9.20088788
                ,9.377969359,9.67497817,11.1517231,11.11223193,11.00340362)

# Add Conditions #
Condition <- c("Control",  "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", 
               "Control", "Control", "Control", "Control", "Control", "Control", "Control", "PCD", "PCD", "PCD")

# Create DataFrame #
GSTA1_Expression <- data.frame(Expression, Condition)

# Create Comparisons #
my_comparisons <- list(c("Control", "PCD"))

# Set Condition as Factor #
GSTA1_Expression$Condition <- as.factor(GSTA1_Expression$Condition)

# Create Boxplot #
p2 <- ggboxplot(GSTA1_Expression, x = "Condition", y = "Expression", color = "black", palette = c("brown", "aquamarine4"), fill = "Condition", outlier.shape = NA) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
  theme(
    legend.position = "right",
    axis.title.x = element_text(size=14, face="bold"),    
    axis.title.y = element_text(size=14, face="bold"),
    legend.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face = "bold")) + 
  ggtitle("GSTA1")

# Visualize #
p2


##---ALDH3A1---##

# Add Expression Levels #
Expression <- c(9.481587006,9.31849593,11.04168368,10.45233954,9.7256683,9.50535183,8.86183574,9.97469837,10.51338781,11.01397243,11.38797525,10.99568542,10.85288746,10.24885734,10.41356043,9.01344057
                ,9.959393948,10.30712863,11.35929439,11.65265038,11.71884067)

# Add Conditions #
Condition <- c("Control",  "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", 
               "Control", "Control", "Control", "Control", "Control", "Control", "Control", "PCD", "PCD", "PCD")

# Create DataFrame #
ALDH3A1_Expression <- data.frame(Expression, Condition)

# Create Comparisons #
my_comparisons <- list(c("Control", "PCD"))

# Set Condition as Factor #
ALDH3A1_Expression$Condition <- as.factor(ALDH3A1_Expression$Condition)

# Create Boxplot #
p3 <- ggboxplot(ALDH3A1_Expression, x = "Condition", y = "Expression", color = "black", palette = c("brown", "aquamarine4"), fill = "Condition", outlier.shape = NA) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
  theme(
    legend.position = "right",
    axis.title.x = element_text(size=14, face="bold"),    
    axis.title.y = element_text(size=14, face="bold"),
    legend.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face = "bold")) + 
  ggtitle("ALDH3A1")

# Visualize #
p3


##---CBR1---##

# Add Expression Levels #
Expression <- c(5.221459172,6.24258308,6.120424532,6.61739182,7.00894951,6.86815802,4.82053761,6.39732003,7.12448246,7.02074656,7.20951157,6.76007504,
                6.49929674,4.54310333,5.32450916,6.57997796,5.992400681,6.14603355,6.56688382,6.99085274,6.779838917)

# Add Conditions #
Condition <- c("Control",  "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", 
               "Control", "Control", "Control", "Control", "Control", "Control", "Control", "PCD", "PCD", "PCD")

# Create DataFrame #
CBR1_Expression <- data.frame(Expression, Condition)

# Create Comparisons #
my_comparisons <- list(c("Control", "PCD"))

# Set Condition as Factor #
CBR1_Expression$Condition <- as.factor(CBR1_Expression$Condition)

# Create Boxplot #
p4 <- ggboxplot(CBR1_Expression, x = "Condition", y = "Expression", color = "black", palette = c("brown", "aquamarine4"), fill = "Condition", outlier.shape = NA) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
  theme(
    legend.position = "right",
    axis.title.x = element_text(size=14, face="bold"),    
    axis.title.y = element_text(size=14, face="bold"),
    legend.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face = "bold")) + 
  ggtitle("CBR1")

# Visualize #
p4


##---ALDH3B1---##

# Add Expression Levels #
Expression <- c(6.646237864,6.70820693,7.570301601,6.9394896,6.83341596,6.76968252,6.35001567,6.98858614,7.40408239,6.63997341,6.89425643
                ,6.4890521,6.42611922,5.34906636,5.95297574,5.05273095,5.760489539,5.78210636,5.70793736,6.26150751,6.063580707)

# Add Conditions #
Condition <- c("Control",  "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", 
               "Control", "Control", "Control", "Control", "Control", "Control", "Control", "PCD", "PCD", "PCD")

# Create DataFrame #
ALDH3B1_Expression <- data.frame(Expression, Condition)

# Create Comparisons #
my_comparisons <- list(c("Control", "PCD"))

# Set Condition as Factor #
ALDH3B1_Expression$Condition <- as.factor(ALDH3B1_Expression$Condition)

# Create Boxplot #
p5 <- ggboxplot(ALDH3B1_Expression, x = "Condition", y = "Expression", color = "black", palette = c("brown", "aquamarine4"), fill = "Condition", outlier.shape = NA) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
  theme(
    legend.position = "right",
    axis.title.x = element_text(size=14, face="bold"),    
    axis.title.y = element_text(size=14, face="bold"),
    legend.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face = "bold")) + 
  ggtitle("ALDH3B1")

# Visualize #
p5


##---SLC7A2---##

# Add Expression Levels #
Expression <- c(1.544699838,1.3701287,4.6762498,4.25647734,3.39366247,3.72163622,5.22512743,5.03953942,4.70039861,5.83717064,5.97248537,5.7170515,4.96024941,4.38754051,5.07810755,4.48122783,5.613827848
                ,5.44462914,4.75079575,4.98346621,4.53101556)

# Add Conditions #
Condition <- c("Control",  "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", 
               "Control", "Control", "Control", "Control", "Control", "Control", "Control", "PCD", "PCD", "PCD")

# Create DataFrame #
SLC7A2_Expression <- data.frame(Expression, Condition)

# Create Comparisons #
my_comparisons <- list(c("Control", "PCD"))

# Set Condition as Factor #
SLC7A2_Expression$Condition <- as.factor(SLC7A2_Expression$Condition)

# Create Boxplot #
p6 <- ggboxplot(SLC7A2_Expression, x = "Condition", y = "Expression", color = "black", palette = c("brown", "aquamarine4"), fill = "Condition", outlier.shape = NA) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
  theme(
    legend.position = "right",
    axis.title.x = element_text(size=14, face="bold"),    
    axis.title.y = element_text(size=14, face="bold"),
    legend.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face = "bold")) + 
  ggtitle("SLC7A2")

# Visualize #
p6


##---CBR3---##

# Add Expression Levels #
Expression <- c(-1.625225163,-1.7997963,-0.848041343,0.18347049,1.37658896,0.37975853,-1.21166121,-1.16226488,-2.29241101,1.71275766,2.09615609,1.59673818
                ,1.57390095,-1.56820625,-2.1445776,-1.29459901,-1.509895644,-2.29394571,0.13897534,2.30209923,1.870318937)

# Add Conditions #
Condition <- c("Control",  "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", 
               "Control", "Control", "Control", "Control", "Control", "Control", "Control", "PCD", "PCD", "PCD")

# Create DataFrame #
CBR3_Expression <- data.frame(Expression, Condition)

# Create Comparisons #
my_comparisons <- list(c("Control", "PCD"))

# Set Condition as Factor #
CBR3_Expression$Condition <- as.factor(CBR3_Expression$Condition)

# Create Boxplot #
p7 <- ggboxplot(CBR3_Expression, x = "Condition", y = "Expression", color = "black", palette = c("brown", "aquamarine4"), fill = "Condition", outlier.shape = NA) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
  theme(
    legend.position = "right",
    axis.title.x = element_text(size=14, face="bold"),    
    axis.title.y = element_text(size=14, face="bold"),
    legend.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face = "bold")) + 
  ggtitle("CBR3")

# Visualize #
p7

# Visualize All #
p1+p2+p3+p4+p5+p6+p7

####-------Bulk RNA Seq BoxPlot with P-Values - PCD 111-------####

##---GSTA2---##

# Add Expression Levels #
Expression <- c(5.518480975,6.42110915,5.815485347,5.289505233,6.033010953,4.200619831,8.131633277,8.22495267,8.152066832)

# Add Conditions #
Condition <- c("Control",  "Control", "Control", "Control", "Control", "Control", "PCD", "PCD", "PCD")

# Create DataFrame #
GSTA2_Expression <- data.frame(Expression, Condition)

# Create Comparisons #
my_comparisons <- list(c("Control", "PCD"))

# Set Condition as Factor #
GSTA2_Expression$Condition <- as.factor(GSTA2_Expression$Condition)

# Create Boxplot #
p1 <- ggboxplot(GSTA2_Expression, x = "Condition", y = "Expression", color = "black", palette = c("brown", "aquamarine4"), fill = "Condition", outlier.shape = NA) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
  theme(
    legend.position = "right",
    axis.title.x = element_text(size=14, face="bold"),    
    axis.title.y = element_text(size=14, face="bold"),
    legend.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face = "bold")) + 
  ggtitle("GSTA2")

# Visualize #
p1


##---GSTA1---##

# Add Expression Levels #
Expression <- c(6.992218361,7.55351856,7.458259725,8.196369925,8.653491557,7.908795494,9.306811028,9.33706693,9.384270968)

# Add Conditions #
Condition <- c("Control",  "Control", "Control", "Control", "Control", "Control", "PCD", "PCD", "PCD")

# Create DataFrame #
GSTA1_Expression <- data.frame(Expression, Condition)

# Create Comparisons #
my_comparisons <- list(c("Control", "PCD"))

# Set Condition as Factor #
GSTA1_Expression$Condition <- as.factor(GSTA1_Expression$Condition)

# Create Boxplot #
p2 <- ggboxplot(GSTA1_Expression, x = "Condition", y = "Expression", color = "black", palette = c("brown", "aquamarine4"), fill = "Condition", outlier.shape = NA) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
  theme(
    legend.position = "right",
    axis.title.x = element_text(size=14, face="bold"),    
    axis.title.y = element_text(size=14, face="bold"),
    legend.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face = "bold")) + 
  ggtitle("GSTA1")

# Visualize #
p2


##---ALDH3A1---##

# Add Expression Levels #
Expression <- c(6.898535702,7.05049245,5.383743126,5.898995874,7.065630287,4.344460861,7.60209352,8.59636078,8.918207579)

# Add Conditions #
Condition <- c("Control",  "Control", "Control", "Control", "Control", "Control", "PCD", "PCD", "PCD")

# Create DataFrame #
ALDH3A1_Expression <- data.frame(Expression, Condition)

# Create Comparisons #
my_comparisons <- list(c("Control", "PCD"))

# Set Condition as Factor #
ALDH3A1_Expression$Condition <- as.factor(ALDH3A1_Expression$Condition)

# Create Boxplot #
p3 <- ggboxplot(ALDH3A1_Expression, x = "Condition", y = "Expression", color = "black", palette = c("brown", "aquamarine4"), fill = "Condition", outlier.shape = NA) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
  theme(
    legend.position = "right",
    axis.title.x = element_text(size=14, face="bold"),    
    axis.title.y = element_text(size=14, face="bold"),
    legend.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face = "bold")) + 
  ggtitle("ALDH3A1")

# Visualize #
p3


##---ALDH3B1---##

# Add Expression Levels #
Expression <- c(4.086386441,4.24955175,4.034586387,3.982843895,4.791774525,3.067144175,5.701647556,5.7519623,5.805630136)

# Add Conditions #
Condition <- c("Control",  "Control", "Control", "Control", "Control", "Control", "PCD", "PCD", "PCD")

# Create DataFrame #
ALDH3B1_Expression <- data.frame(Expression, Condition)

# Create Comparisons #
my_comparisons <- list(c("Control", "PCD"))

# Set Condition as Factor #
ALDH3B1_Expression$Condition <- as.factor(ALDH3B1_Expression$Condition)

# Create Boxplot #
p4 <- ggboxplot(ALDH3B1_Expression, x = "Condition", y = "Expression", color = "black", palette = c("brown", "aquamarine4"), fill = "Condition", outlier.shape = NA) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
  theme(
    legend.position = "right",
    axis.title.x = element_text(size=14, face="bold"),    
    axis.title.y = element_text(size=14, face="bold"),
    legend.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face = "bold")) + 
  ggtitle("ALDH3B1")

# Visualize #
p4


##---UCP2 ---##

# Add Expression Levels #
Expression <- c(8.487304554,8.53565379,8.27684531,8.659390661,9.161922102,8.291806305,10.30370656,10.28546926,10.30964058)

# Add Conditions #
Condition <- c("Control",  "Control", "Control", "Control", "Control", "Control", "PCD", "PCD", "PCD")

# Create DataFrame #
UCP2_Expression <- data.frame(Expression, Condition)

# Create Comparisons #
my_comparisons <- list(c("Control", "PCD"))

# Set Condition as Factor #
UCP2_Expression$Condition <- as.factor(UCP2_Expression$Condition)

# Create Boxplot #
p5 <- ggboxplot(UCP2_Expression, x = "Condition", y = "Expression", color = "black", palette = c("brown", "aquamarine4"), fill = "Condition", outlier.shape = NA) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
  theme(
    legend.position = "right",
    axis.title.x = element_text(size=14, face="bold"),    
    axis.title.y = element_text(size=14, face="bold"),
    legend.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face = "bold")) + 
  ggtitle("UCP2")

# Visualize #
p5


##---SLC7A2 ---##

# Add Expression Levels #
Expression <- c(1.585752459,2.42281731,2.013705419,0.41982032,1.029582222,-0.985967161,3.237168483,3.63924389,2.614513692)

# Add Conditions #
Condition <- c("Control",  "Control", "Control", "Control", "Control", "Control", "PCD", "PCD", "PCD")

# Create DataFrame #
SLC7A2_Expression <- data.frame(Expression, Condition)

# Create Comparisons #
my_comparisons <- list(c("Control", "PCD"))

# Set Condition as Factor #
SLC7A2_Expression$Condition <- as.factor(SLC7A2_Expression$Condition)

# Create Boxplot #
p6 <- ggboxplot(SLC7A2_Expression, x = "Condition", y = "Expression", color = "black", palette = c("brown", "aquamarine4"), fill = "Condition", outlier.shape = NA) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
  theme(
    legend.position = "right",
    axis.title.x = element_text(size=14, face="bold"),    
    axis.title.y = element_text(size=14, face="bold"),
    legend.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face = "bold")) + 
  ggtitle("SLC7A2")

# Visualize #
p6


##---CBR3 ---##

# Add Expression Levels #
Expression <- c(0.03215413,0.84775327,-0.241250122,0.396919918,0.418757658,0.287955561,1.483840881,1.10908041,0.881545277)

# Add Conditions #
Condition <- c("Control",  "Control", "Control", "Control", "Control", "Control", "PCD", "PCD", "PCD")

# Create DataFrame #
CBR3_Expression <- data.frame(Expression, Condition)

# Create Comparisons #
my_comparisons <- list(c("Control", "PCD"))

# Set Condition as Factor #
CBR3_Expression$Condition <- as.factor(CBR3_Expression$Condition)

# Create Boxplot #
p7 <- ggboxplot(CBR3_Expression, x = "Condition", y = "Expression", color = "black", palette = c("brown", "aquamarine4"), fill = "Condition", outlier.shape = NA) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
  theme(
    legend.position = "right",
    axis.title.x = element_text(size=14, face="bold"),    
    axis.title.y = element_text(size=14, face="bold"),
    legend.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face = "bold")) + 
  ggtitle("CBR3")

# Visualize #
p7

# Visualize All #
p1+p2+p3+p4+p5+p6+p7

####-------Bulk RNA Seq BoxPlot with P-Values - PCD 114-------####

##---GSTA2---##

# Add Expression Levels #
Expression <- c(5.518480975,6.42110915,5.815485347,5.289505233,6.033010953,4.200619831,5.911647249,7.041339043,7.307644461)

# Add Conditions #
Condition <- c("Control",  "Control", "Control", "Control", "Control", "Control", "PCD", "PCD", "PCD")

# Create DataFrame #
GSTA2_Expression <- data.frame(Expression, Condition)

# Create Comparisons #
my_comparisons <- list(c("Control", "PCD"))

# Set Condition as Factor #
GSTA2_Expression$Condition <- as.factor(GSTA2_Expression$Condition)

# Create Boxplot #
p1 <- ggboxplot(GSTA2_Expression, x = "Condition", y = "Expression", color = "black", palette = c("brown", "aquamarine4"), fill = "Condition", outlier.shape = NA) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
  theme(
    legend.position = "right",
    axis.title.x = element_text(size=14, face="bold"),    
    axis.title.y = element_text(size=14, face="bold"),
    legend.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face = "bold")) + 
  ggtitle("GSTA2")

# Visualize #
p1


##---GSTA1---##

# Add Expression Levels #
Expression <- c(6.992218361,7.55351856,7.458259725,8.196369925,8.653491557,7.908795494,9.256027057,9.57424925,9.758903928)

# Add Conditions #
Condition <- c("Control",  "Control", "Control", "Control", "Control", "Control", "PCD", "PCD", "PCD")

# Create DataFrame #
GSTA1_Expression <- data.frame(Expression, Condition)

# Create Comparisons #
my_comparisons <- list(c("Control", "PCD"))

# Set Condition as Factor #
GSTA1_Expression$Condition <- as.factor(GSTA1_Expression$Condition)

# Create Boxplot #
p2 <- ggboxplot(GSTA1_Expression, x = "Condition", y = "Expression", color = "black", palette = c("brown", "aquamarine4"), fill = "Condition", outlier.shape = NA) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
  theme(
    legend.position = "right",
    axis.title.x = element_text(size=14, face="bold"),    
    axis.title.y = element_text(size=14, face="bold"),
    legend.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face = "bold")) + 
  ggtitle("GSTA1")

# Visualize #
p2


##---ALDH3A1---##

# Add Expression Levels #
Expression <- c(6.898535702,7.05049245,5.383743126,5.898995874,7.065630287,4.344460861,7.900998667,8.255399422,7.982577067)

# Add Conditions #
Condition <- c("Control",  "Control", "Control", "Control", "Control", "Control", "PCD", "PCD", "PCD")

# Create DataFrame #
ALDH3A1_Expression <- data.frame(Expression, Condition)

# Create Comparisons #
my_comparisons <- list(c("Control", "PCD"))

# Set Condition as Factor #
ALDH3A1_Expression$Condition <- as.factor(ALDH3A1_Expression$Condition)

# Create Boxplot #
p3 <- ggboxplot(ALDH3A1_Expression, x = "Condition", y = "Expression", color = "black", palette = c("brown", "aquamarine4"), fill = "Condition", outlier.shape = NA) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
  theme(
    legend.position = "right",
    axis.title.x = element_text(size=14, face="bold"),    
    axis.title.y = element_text(size=14, face="bold"),
    legend.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face = "bold")) + 
  ggtitle("ALDH3A1")

# Visualize #
p3


##---CBR1---##

# Add Expression Levels #
Expression <- c(3.333773242,4.42715624,3.369508379,3.366300683,3.177368246,3.161931534,4.723046764,5.680679709,4.448446041)

# Add Conditions #
Condition <- c("Control",  "Control", "Control", "Control", "Control", "Control", "PCD", "PCD", "PCD")

# Create DataFrame #
CBR1_Expression <- data.frame(Expression, Condition)

# Create Comparisons #
my_comparisons <- list(c("Control", "PCD"))

# Set Condition as Factor #
CBR1_Expression$Condition <- as.factor(CBR1_Expression$Condition)

# Create Boxplot #
p4 <- ggboxplot(CBR1_Expression, x = "Condition", y = "Expression", color = "black", palette = c("brown", "aquamarine4"), fill = "Condition", outlier.shape = NA) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
  theme(
    legend.position = "right",
    axis.title.x = element_text(size=14, face="bold"),    
    axis.title.y = element_text(size=14, face="bold"),
    legend.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face = "bold")) + 
  ggtitle("CBR1")

# Visualize #
p4


##---UCP2---##

# Add Expression Levels #
Expression <- c(8.487304554,8.53565379,8.27684531,8.659390661,9.161922102,8.291806305,9.892525785,9.948161698,10.00354351)

# Add Conditions #
Condition <- c("Control",  "Control", "Control", "Control", "Control", "Control", "PCD", "PCD", "PCD")

# Create DataFrame #
UCP2_Expression <- data.frame(Expression, Condition)

# Create Comparisons #
my_comparisons <- list(c("Control", "PCD"))

# Set Condition as Factor #
UCP2_Expression$Condition <- as.factor(UCP2_Expression$Condition)

# Create Boxplot #
p5 <- ggboxplot(UCP2_Expression, x = "Condition", y = "Expression", color = "black", palette = c("brown", "aquamarine4"), fill = "Condition", outlier.shape = NA) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
  theme(
    legend.position = "right",
    axis.title.x = element_text(size=14, face="bold"),    
    axis.title.y = element_text(size=14, face="bold"),
    legend.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face = "bold")) + 
  ggtitle("UCP2")

# Visualize #
p5


##---SLC7A2---##

# Add Expression Levels #
Expression <- c(1.585752459,2.42281731,2.013705419,0.41982032,1.029582222,-0.985967161,1.558872832,2.157030515,2.502310068)

# Add Conditions #
Condition <- c("Control",  "Control", "Control", "Control", "Control", "Control", "PCD", "PCD", "PCD")

# Create DataFrame #
SLC7A2_Expression <- data.frame(Expression, Condition)

# Create Comparisons #
my_comparisons <- list(c("Control", "PCD"))

# Set Condition as Factor #
SLC7A2_Expression$Condition <- as.factor(SLC7A2_Expression$Condition)

# Create Boxplot #
p6 <- ggboxplot(SLC7A2_Expression, x = "Condition", y = "Expression", color = "black", palette = c("brown", "aquamarine4"), fill = "Condition", outlier.shape = NA) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
  theme(
    legend.position = "right",
    axis.title.x = element_text(size=14, face="bold"),    
    axis.title.y = element_text(size=14, face="bold"),
    legend.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face = "bold")) + 
  ggtitle("SLC7A2")

# Visualize #
p6


##---CDKN2A---##

# Add Expression Levels #
Expression <- c(2.799057008,2.66458924,2.514188132,2.515944021,2.38766998,2.972453735,1.665524458,1.778318225,1.569804922)

# Add Conditions #
Condition <- c("Control",  "Control", "Control", "Control", "Control", "Control", "PCD", "PCD", "PCD")

# Create DataFrame #
CDKN2A_Expression <- data.frame(Expression, Condition)

# Create Comparisons #
my_comparisons <- list(c("Control", "PCD"))

# Set Condition as Factor #
CDKN2A_Expression$Condition <- as.factor(CDKN2A_Expression$Condition)

# Create Boxplot #
p7 <- ggboxplot(CDKN2A_Expression, x = "Condition", y = "Expression", color = "black", palette = c("brown", "aquamarine4"), fill = "Condition", outlier.shape = NA) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
  theme(
    legend.position = "right",
    axis.title.x = element_text(size=14, face="bold"),    
    axis.title.y = element_text(size=14, face="bold"),
    legend.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face = "bold")) + 
  ggtitle("CDKN2A")

# Visualize #
p7

# Visualize All #
p1+p2+p3+p4+p5+p6+p7
