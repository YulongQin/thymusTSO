#' Find medulla clusters and corresponding edges.
#'
#' This function identifies medulla and cortex spots and their corresponding boundaries based on thymus spatial transcriptomics (ST) data.
#'
#' @param obj.st.lst A list of thymus spatial seurat objects.
#' @param medulla.genes A vector of genes associated with the medulla.
#' @param call.xgb Call trained XGBoost model or not. Default: FALSE.
#' @param p.cut Cutoff p-value to determine significant spots. Default is 1e-5.
#' @param module.size Minimum module size for modules. Default is 10.
#' @param remove.spots Optional list of spots to be removed from analysis.
#' @param out.figs.dir Directory to save output figures. If provided, a density plot of medulla scores will be saved.
#' @return Updated obj.st.lst with identified medulla clusters and edges.
#' @export tsoHis
#'
#' @examples
#' # Example usage:
#' sp.obj <- system.file("data/thymus_T2.RDS", package = "thymusTSO") %>% readRDS()
#' sp.obj <- tsoHis(sp.obj)
tsoHis <- function(obj.st.lst, medulla.genes = NULL, call.xgb = FALSE,cells_pct = 0.05, module.size = 10, remove.spots = NULL) {
    message(paste0("TSO-his Start at ", Sys.time()))
    if (!inherits(obj.st.lst, "list")) obj.st.lst <- list(TSOhis = obj.st.lst)
    if (is.null(medulla.genes)) medulla.genes <- c("EBI3", "CCL17", "CCR7", "CSF2RB", "CCL21", "CCL22", "TNFRSF18", "CCL27", "CXCL10", "CXCL9", "MS4A1", "LAMP3")
    obj.st.lst <- medullaScore(obj.st.lst, medulla.genes)
    if (length(cells_pct) == 1) cells_pct <- rep(cells_pct, length(obj.st.lst)) %>% `names<-`(names(obj.st.lst))
    if (length(module.size) == 1) module.size <- rep(module.size, length(obj.st.lst)) %>% `names<-`(names(obj.st.lst))
    module.size <- module.size[names(obj.st.lst)]

    obj.st.lst <- lapply(names(obj.st.lst), function(sn) {
        obj <- obj.st.lst[[sn]]
        message(sprintf("Running analysis sn: %s", sn))
        obj <- AddMetaData(obj, metadata = "Cortex", col.name = "HE.Labels")
        
        message(paste0("Step1: calcSpotsDist start at ", Sys.time()))
        dist.sig <- calcSpotsDist(obj, cells_pct = cells_pct[sn], call.xgb = call.xgb) # !!!: 修改了其中的代码
        sig.spots <- dist.sig$sig.spots
        elu.dist <- dist.sig$dist
        while (TRUE) {
            message(paste0("Step2: unpdateSigSpots start at ", Sys.time()))
            sig.spots.new <- unpdateSigSpots(elu.dist, sig.spots)
            bool.val <- length(setdiff(sig.spots.new, sig.spots)) > 0
            if (bool.val) {
                sig.spots <- sig.spots.new
            } else {
                break
            }
        }
        if (sn %in% names(remove.spots)) {
            cus.rm <- remove.spots[[sn]]
            idx <- which(sig.spots %in% cus.rm)
            if (length(idx) > 0) sig.spots <- sig.spots[-idx]
        }
        elu.dist.sig <- dist.sig$dist[sig.spots, ]
        message(paste0("Step3: removeLowConfModules start at ", Sys.time()))
        # debug(removeLowConfModules)
        sig.spots.classes <- removeLowConfModules(elu.dist.sig, module.size = module.size[sn]) # !!!: 修改了其中的代码
        cand.sig.dist <- dist.sig$dist[sig.spots.classes$remain %>% unlist(), ]

        min.dist <- apply(cand.sig.dist, 1, function(obj) {
            obj[order(obj)][2]
        })
        rm.spots <- min.dist[which(min.dist > mean(min.dist) + 3 * sd(min.dist))] %>% names()
        sig.spots.classes$remain <- lapply(sig.spots.classes$remain, function(obj) {
            idx <- which(obj %in% rm.spots)
            if (length(idx) > 0) {
                obj[-idx]
            } else {
                obj
            }
        })
        sig.spots.classes$remove <- c(sig.spots.classes$remove, rm.spots)
        message(paste0("Step4: candEdgeSpots start at ", Sys.time()))
        edge.spots <- candEdgeSpots(elu.dist.sig[sig.spots.classes$remain %>% unlist(., use.names = FALSE), ], dist.sig$dist)
        
        message(paste0("Step5: updateModuleCenter start at ", Sys.time()))
        sig.spots.classes$remain <- updateModuleCenter(obj, sig.spots.classes$remain, edge.spots)
        remain_list <- sig.spots.classes$remain
        tmp <- stack(remain_list)
        
        obj@meta.data[sig.spots.classes$remove, "HE.Labels"] <- "Medulla_lo"
        obj@meta.data[sig.spots.classes$remain %>% unlist(., use.names = FALSE), "HE.Labels"] <- "Medulla_hi"
        obj@meta.data["Cluster_adj_matrix"] <- NA
        obj@meta.data[tmp$values, "Cluster_adj_matrix"] <- tmp$ind %>% as.character()
        obj@meta.data[edge.spots, "HE.Labels"] <- "Medulla_edge"
        obj@meta.data[sig.spots.classes$remain %>% names(), "HE.Labels"] <- "Medulla_centric"
        obj@meta.data$Centric <- ifelse(obj@meta.data$HE.Labels == "Medulla_centric", "Y", "N")
        obj
    }) %>% `names<-`(names(obj.st.lst))
    # obj.st.lst[[1]] <- obj
    message(paste0("Step6: calcSpot2ModuleDist start at ", Sys.time()))
    obj.st.lst <- calcSpot2ModuleDist(obj.st.lst)
    obj.st.lst[[1]]@meta.data$Assign.Centric <- ifelse(is.na(obj.st.lst[[1]]@meta.data$Cluster_adj_matrix), 
                                                       obj.st.lst[[1]]@meta.data$Assign.Centric,obj.st.lst[[1]]@meta.data$Cluster_adj_matrix )
    # meta_data <- obj.st.lst[[1]]@meta.data
    print(paste0("All done at ", Sys.time()))
    return(obj.st.lst)
}
