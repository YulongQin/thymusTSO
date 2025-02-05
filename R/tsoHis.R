#' Find medulla clusters and corresponding edges.
#'
#' This function identifies medulla and cortex spots and their corresponding boundaries based on thymus spatial transcriptomics (ST) data.
#'
#' @param obj.st.lst A list of thymus spatial seurat objects.
#' @param medulla.genes A vector of genes associated with the medulla.
#' @param call.xgb Call trained XGBoost model or not. Default: FALSE.
#' @param cells_pct Percentage of cells to be used for clustering. Default is 0.05.
#' @param module.size Minimum module size for modules. Default is 10.
#' @param remove.spots Optional list of spots to be removed from analysis.
#' @param seq_method Sequencing method used for the data. Default is "stereo".
#' @return Updated obj.st.lst with identified medulla clusters and edges.
#' @export tsoHis
#'
#' @examples
#' # Example usage:
#' sp.obj <- system.file("data/thymus_T2.RDS", package = "thymusTSO") %>% readRDS()
#' sp.obj <- tsoHis(sp.obj)
tsoHis <- function(obj.st.lst, medulla.genes = NULL, call.xgb = FALSE,cells_pct = 0.05, 
                   module.size = 10, remove.spots = NULL,
                   seq_method = "stereo") {
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
        dist.sig <- calcSpotsDist(obj, cells_pct = cells_pct[sn], call.xgb = call.xgb, seq_method = seq_method) # !!!: 修改了其中的代码
        gc()
        sig.spots <- dist.sig$sig.spots
        elu.dist <- dist.sig$dist
        if(seq_method == "stereo") { # !!!: 添加代码
            spots.diff.thres <- cells_pct*ncol(obj)*0.0015
        }else{
            spots.diff.thres <- 0
        }
        message(paste0("spots.diff.thres: ", spots.diff.thres))
        while (TRUE) {
            message(paste0("Step2: updateSigSpots start at ", Sys.time()))
            sig.spots.new <- updateSigSpots(elu.dist, sig.spots, seq_method = seq_method) # !!!: 修改了其中的代码，耗时步骤
            gc()
            
            len_setdiff <- length(setdiff(sig.spots.new, sig.spots))
            bool.val <- len_setdiff > spots.diff.thres # !!!: 修改了其中的代码
            message(paste0("length of setdiff(sig.spots.new, sig.spots): ", len_setdiff))
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
        sig.spots.classes <- removeLowConfModules(elu.dist.sig, module.size = module.size[sn], seq_method = seq_method) # !!!: 修改了其中的代码
        gc()
        cand.sig.dist <- dist.sig$dist[sig.spots.classes$remain %>% unlist(), ]

        min.dist <- apply(cand.sig.dist, 1, function(obj) {
            obj[order(obj)][2]
        })
        if(seq_method == "stereo") { # !!!: 添加代码
            rm.spots <- min.dist[which(min.dist > 1)] %>% names() 
        }else{
            rm.spots <- min.dist[which(min.dist > mean(min.dist) + 3 * sd(min.dist))] %>% names() 
        }
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
        edge.spots <- candEdgeSpots(elu.dist.sig[sig.spots.classes$remain %>% unlist(., use.names = FALSE), ],
                                    dist.sig$dist, seq_method = seq_method) # !!!: 修改了其中的代码
        gc()
        message(paste0("Step5: updateModuleCenter start at ", Sys.time()))
        sig.spots.classes$remain <- updateModuleCenter(obj, sig.spots.classes$remain, edge.spots)
        gc()
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
    gc()
    obj.st.lst[[1]]@meta.data$Assign.Centric <- ifelse(is.na(obj.st.lst[[1]]@meta.data$Cluster_adj_matrix), 
                                                       obj.st.lst[[1]]@meta.data$Assign.Centric,obj.st.lst[[1]]@meta.data$Cluster_adj_matrix )
    # meta_data <- obj.st.lst[[1]]@meta.data
    print(paste0("All done at ", Sys.time()))
    return(obj.st.lst)
}
