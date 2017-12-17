options(warn=-1)


#' Main function
#'
#' @param file_name path to file
#' @param file_type list or file
#' @param save_files boolean
#' @param skip_frames number how many first frames should be skipped
#' @param plots boolean
#' @param protein_representation CA or all
#' @param chain chain ID
#' @param cluster_quantity number
#' @param receptor_peptide_mode boolean
#' @param reference_structure_pdb_code code from PDB web repository
#' @param reference_structure_file path to file
#' @param reference_chain_id chain ID
#' @param kmeans_algorithm Hartigan-Wong / MacQueen / Lloyd / Forgy
#' @param output_directory path to directory

#' @return return object with clusters data
#' @export
run <- function(file_name = "", file_type = "list", save_files = F,
                skip_frames = 0, plots = T, protein_representation = "CA",
                chain = "", cluster_quantity = 10, receptor_peptide_mode = F,
                reference_structure_pdb_code = '',
                reference_structure_file = '',
                reference_chain_id = '',
                kmeans_algorithm = "Hartigan-Wong",
                output_directory = '') {
  file_directory <<- dirname(file_name)
  file_name <- basename(file_name)
  file_type <<- file_type
  save_files <<- save_files
  skip_frames <<- skip_frames
  plots <<- plots
  protein_representation <<- protein_representation
  chain <<- chain
  cluster_quantity <<- cluster_quantity
  kmeans_algorithm <<- kmeans_algorithm
  r_p_mode <<-receptor_peptide_mode
  reference_structure_file <<- reference_structure_file
  reference_structure_pdb_code <<- reference_structure_pdb_code
  reference_chain_id <<- reference_chain_id

  if (output_directory == '') {
    output_directory = file_directory
  }

  if(!save_files){
    plots <<- FALSE
  }

  if (file.exists(paste(file_directory, file_name, sep = "/"))) {
    ptm <- proc.time()
    set_directory(file_directory)
    prepare_data(file_name)
    set_directory(output_directory)
    create_folder_structure()
    make_centroids_medoids()
    write_medoids_centroid()
    pca_reduction()
    k_means()
    k_means_cluster()
    print(proc.time() - ptm)
  } else {
    stop("No data src")
  }

  results <- list(medoids = klaMed, kmeans = kMeansResult, seq = SEQ)
}

set_directory <- function(directory) {
  setwd(directory)
}

create_folder_structure <<- function() {
  if(save_files){
    centroid_folder_name <<- "centroids/"
    medoid_folder_name <<- "medoids/"
    plots_folder_name <<- "plots/"
    reports_folder_name <<- "reports/"

    dir.create(centroid_folder_name, showWarnings = F)
    dir.create(medoid_folder_name, showWarnings = F)
    dir.create(plots_folder_name, showWarnings = F)
    dir.create(reports_folder_name, showWarnings = F)
  }
}

prepare_data <- function(file_name) {
  set.seed(1990)
  nN <<- 1
  nC <<- -1

  if(r_p_mode){
    if(reference_chain_id == ''){
      reference_chain_id <<- NULL
    }

    if(reference_structure_pdb_code != '' & reference_structure_file != ''){
      stop("Choose one reference source file. We can pass either pdb_code or file.")
    }

    if(reference_structure_pdb_code != ''){
      result = tryCatch({
          pdb <-read.pdb( get.pdb(reference_structure_pdb_code, URLonly = T) )
      }, warning = function(w) {
          stop("Wrong pdb code")
      }, error = function(e) {
          stop("Wrong pdb code")
      })
      a.inds <- atom.select.pdb(pdb,"calpha", chain = reference_chain_id)
      ligand_xyz <<- pdb$atom[c(a.inds$atom), c("x", "y", "z")]
      ligand_xyz <<- data.frame(ligand_xyz["x"], ligand_xyz["y"], ligand_xyz["z"])
    }
    if(reference_structure_file != ''){
      pdb <- read.pdb(reference_structure_file, ATOM.only = T)
      a.inds <- atom.select.pdb(pdb,"calpha", chain = reference_chain_id)
      ligand_xyz <<- pdb$atom[c(a.inds$atom), c("x", "y", "z")]
      ligand_xyz <<- data.frame(ligand_xyz["x"], ligand_xyz["y"], ligand_xyz["z"])
    }
  }

  if (file_type == "list") {
    file_list <<- readLines(file_name)
    for (file in file_list) {
      if (!(file.exists(file))) {
        stop("Missing input file: ", file)
      }
      if (!exists("dataset_from_file")) {
        dataset_from_file <- readPdbCoord(file)$selected
        dataset_from_file_all <- readPdbCoord(file)$all
        first_lenght <<- dim(dataset_from_file)[1]
        next
      }
      if (exists("dataset_from_file")) {
        temp_dataset <- readPdbCoord(file)$selected
        temp_dataset_all <- readPdbCoord(file)$all

        if (dim(temp_dataset)[1] != first_lenght) {
          stop("All the files should be of equal length: ", file)
        }
        dataset_from_file <- rbind(dataset_from_file, temp_dataset)
        dataset_from_file_all <- rbind(dataset_from_file_all, temp_dataset_all)
        rm(temp_dataset)
        rm(temp_dataset_all)
      }
    }
    temp_trajectory <- dataset_from_file
    temp_all_trajectory <- dataset_from_file_all
  } else {
    temp_trajectory <- readPdbCoord(file_name)$selected
    temp_all_trajectory <- readPdbCoord(file_name)$all
  }

  if (nrow(temp_trajectory) < 1) {
    stop("Empty data input.")
  }

  IncrementalTable <- c()
  control <- 1
  for (i in 1:nrow(temp_trajectory)) {
    if (control == 0) {
      next
    }
    if (i == nrow(temp_trajectory) ||
        temp_trajectory[i, 2] > temp_trajectory[i + 1, 2]) {
      IncrementalTable <- rbind(IncrementalTable, (temp_trajectory[i, ]))
      row_count <- i
      control <- 0
    } else {
      row_count <- i
      IncrementalTable <- rbind(IncrementalTable, (temp_trajectory[i, ]))
    }
  }

  if (skip_frames > 0) {
    temp_trajectory <- head(temp_trajectory, -(skip_frames * row_count))
  }

  trajectory <<- array(NA,
                       dim = c(row_count, 3, nrow(temp_trajectory)/row_count))
  j <- 1
  for (i in seq(from = 1, to = nrow(temp_trajectory), by = row_count)) {
    end <- i + row_count - 1
    trajectory[, , j] <<- as.matrix(temp_trajectory[i:end, 6:8])
    j <- j + 1
  }

  if (nrow(which(is.na(trajectory), arr.ind = TRUE)) > 0) {
    stop("Wrong PDB structure", which(is.na(trajectory), arr.ind = TRUE))
  }

  SEQ <<- data.frame(c(1:length(IncrementalTable[, 4])), IncrementalTable[, 4],
                     IncrementalTable[, 3], IncrementalTable[, 5],
                     IncrementalTable[, 9])


  IncrementalTable <- c()
  control <- 1
  for (i in 1:nrow(temp_all_trajectory)) {
    if (control == 0) {
      next
    }
    if (i == nrow(temp_all_trajectory) || temp_all_trajectory[i, 2] > temp_all_trajectory[i + 1, 2]) {
      IncrementalTable <- rbind(IncrementalTable, (temp_all_trajectory[i, ]))
      row_count <- i
      control <- 0
    } else {
      row_count <- i
      IncrementalTable <- rbind(IncrementalTable, (temp_all_trajectory[i, ]))
    }
  }

  if (skip_frames > 0) {
    temp_all_trajectory <- head(temp_all_trajectory, -(skip_frames * row_count))
  }

  trajectory_all <<- array(NA,
                       dim = c(row_count, 3, nrow(temp_all_trajectory)/row_count))
  j <- 1
  for (i in seq(from = 1, to = nrow(temp_all_trajectory), by = row_count)) {
    end <- i + row_count - 1
    trajectory_all[, , j] <<- as.matrix(temp_all_trajectory[i:end, 6:8])
    j <- j + 1
  }

  SEQ_all <<- data.frame(c(1:length(IncrementalTable[, 4])), IncrementalTable[, 4],
                     IncrementalTable[, 3], IncrementalTable[, 5],
                     IncrementalTable[, 9])
  if (nC == -1) {
    nC <<- dim(trajectory)[1]
  }
}

rmsd <- function(X, Y) {
  X <- t(t(X) - apply(X, 2, mean))
  Y <- t(t(Y) - apply(Y, 2, mean))
  A <- svd(t(X) %*% Y)
  Q <- A$u %*% t(A$v)
  if (det(Q) < 0) {
    A$u[, 3] <- -A$u[, 3]
    Q <- A$u %*% t(A$v)
  }
  sqrt(sum((Y - X %*% Q)^2)/nrow(Y))
}

rmsd_without_superimposition <- function(X, Y) {
  round(sqrt(sum((Y - X)^2)/nrow(Y)),2)
}

rmsdQ <- function(X, Y) {
  X <- t(t(X) - apply(X, 2, mean))
  Y <- t(t(Y) - apply(Y, 2, mean))
  A <- svd(t(X) %*% Y)
  Q <- A$u %*% t(A$v)
  if (det(Q) < 0) {
    A$u[, 3] <- -A$u[, 3]
    Q <- A$u %*% t(A$v)
  }
  Q
}

readPdbCoord <- function(plik) {
  Y <- scan(plik, what = "", sep = "\n", quiet = TRUE)
  if (protein_representation == "CA") {
    if (chain != "") {
      YY <- Y[substr(Y, 1, 4) == "ATOM" &
                substr(Y, 14, 15) == "CA" &
                substr(Y, 22, 22) == chain]
      rebuildPDB <- Y[substr(Y, 1, 4) == "ATOM" & (substr(Y, 14, 15) == "CA")]
    } else {
      YY <- Y[substr(Y, 1, 4) == "ATOM" & (substr(Y, 14, 15) == "CA")]
      rebuildPDB <- YY
    }
  } else {
    if (chain != "") {
      YY <- Y[substr(Y, 1, 4) == "ATOM" &
                substr(Y, 22, 22) == chain]
      rebuildPDB <- Y[substr(Y, 1, 4) == "ATOM"]
    } else {
      YY <- Y[substr(Y, 1, 4) == "ATOM"]
      rebuildPDB <- YY
    }
  }
  X <- matrix(NA, length(YY), 8)
  X1 <- as.character(substr(YY, 1, 4))
  X2 <- as.numeric(substr(YY, 7, 11))
  X3 <- trimws(as.character(substr(YY, 14, 16)))
  X4 <- as.character(substr(YY, 18, 20))
  X5 <- as.numeric(substr(YY, 23, 26))
  X6 <- as.numeric(substr(YY, 31, 38))
  X7 <- as.numeric(substr(YY, 39, 46))
  X8 <- as.numeric(substr(YY, 47, 54))
  X9 <- as.character(substr(YY, 22, 22))

  selected = data.frame(X1, X2, X3, X4, X5, X6, X7, X8, X9,
                        stringsAsFactors = FALSE)

  X <- matrix(NA, length(rebuildPDB), 8)
  X1 <- as.character(substr(rebuildPDB, 1, 4))
  X2 <- as.numeric(substr(rebuildPDB, 7, 11))
  X3 <- trimws(as.character(substr(rebuildPDB, 14, 16)))
  X4 <- as.character(substr(rebuildPDB, 18, 20))
  X5 <- as.numeric(substr(rebuildPDB, 23, 26))
  X6 <- as.numeric(substr(rebuildPDB, 31, 38))
  X7 <- as.numeric(substr(rebuildPDB, 39, 46))
  X8 <- as.numeric(substr(rebuildPDB, 47, 54))
  X9 <- as.character(substr(rebuildPDB, 22, 22))

  all = data.frame(X1, X2, X3, X4, X5, X6, X7, X8, X9,
                        stringsAsFactors = FALSE)

  list(selected = selected, all = all)
}

model2pdb <- function(SEQ, X, bFact, plik_out) {
  n <- nrow(SEQ)
  Y <- as.data.frame(matrix(NA, n, 11))
  Y[, 1] <- rep("ATOM", n)
  Y[, 2] <- SEQ[, 1]
  if (protein_representation == "all") {
    Y[, 3] <- SEQ[, 3]
  } else {
    Y[, 3] <- rep("CA", n)
  }
  Y[, 4] <- SEQ[, 2]
  Y[, 5] <- SEQ[, 5]
  Y[, 6] <- SEQ[, 4]
  Y[, 7:9] <- round(X[1:nrow(X) %in% SEQ[, 1], ], 3)
  Y[, 10] <- 1
  Y[, 11] <- 1#round(bFact[1:nrow(X) %in% SEQ[, 1]], 2)

  xx <- sprintf("%4s%7d%4s%5s%2s%4d%12.3f%8.3f%8.3f%6.2f%6.2f", Y[, 1], Y[, 2],
                Y[, 3], Y[, 4], Y[, 5], Y[, 6], Y[, 7], Y[, 8], Y[, 9], Y[, 10],
                Y[, 11])
  cat(noquote(xx), file = plik_out, sep = "\n")
}

klasterMed <- function(kol, n0, n1, k_count) {
  Yref <- 0 * k_count
  Ymed <- list()
  Ycen <- Ymed
  nKol <- sort(unique(kol))
  nKlast <- 0 * nKol
  mRMSD <- 0 * nKol
  files <- array(NA, dim = c(k_count))
  files_rebuild <- array(NA, dim = c(k_count))

  for (k in nKol) {
    iDecK <- which(kol == k)
    nDecK <- length(iDecK)
    nKlast[k] <- nDecK
    Ycen0 <- trajectory[, , iDecK[1]]
    Ymed0 <- Ycen0
    r1 <- 0
    if (nDecK > 1) {
      for (i in 2:nDecK) {
        Ycen0 <- ((i - 1) * Ycen0 + trajectory[, , iDecK[i]] %*% rmsdQ(
          trajectory[n0:n1, , iDecK[i]], Ycen0[n0:n1, ]))/i
      }

      R1 = apply(trajectory[, , iDecK], 3, function(Y)
        rmsd(Ycen0[n0:n1,], Y[n0:n1,]))

      r1 = mean(R1)

      Ymed0 <- trajectory[, , iDecK[which.min(R1)]]
      if (file_type == "list") {
        file_med <- file_list[iDecK[which.min(R1)] - 1]
        file_med_rebuild <-iDecK[which.min(R1)]
      } else {
        file_med <- iDecK[which.min(R1)]
      }
    }
    Ycen[[k]] <- Ycen0
    Ymed[[k]] <- Ymed0
    mRMSD[k] <- r1
    if(exists("file_med")){
      files[k] <- file_med
    }
    if(exists("file_med_rebuild")){
      files_rebuild[k] <- file_med_rebuild
    }
  }

  if(r_p_mode & exists("ligand_xyz")){
    for(i in 1:cluster_quantity){
      Yref[i] = rmsd_without_superimposition(Ymed[[i]], ligand_xyz)
    }
  } else {
    for(i in 1:cluster_quantity){
      Yref[i] = "n.a"
    }
  }

  list(cen = Ycen, med = Ymed, nClust = nKlast, mRMS = mRMSD,
       files_names = files, files_rebuild_indexes = files_rebuild, rmsdRef = Yref)
}

myLog <- function(x) ifelse(x > 0, log(x), 0)

mutInf <- function(kol, traj) {
  P <- table(kol, traj)
  P <- P/sum(P)
  Py <- apply(P, 1, sum)
  Px <- apply(P, 2, sum)
  PP <- outer(Py, Px, "*")
  Htraj <- -sum(myLog(Px) * Px)
  Hkol <- -sum(myLog(Py) * Py)
  infKL <- sum(myLog(P/PP) * P)
  list(redun = infKL/(Htraj + Hkol), infNorm = infKL/min(Htraj, Hkol),
       Cxy = infKL/Hkol)
}

gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

get_colors <- function(groups, group.col = palette()) {
  groups <- as.factor(groups)
  ngrps <- length(levels(groups))
  if (ngrps > length(group.col)) {
    group.col <- rep(group.col, ngrps)
  }

  color <- group.col[as.numeric(groups)]
  names(color) <- as.vector(groups)
  return(color)
}

plotPCA <- function(x, decoys, K) {
  groups <- x$cluster
  groups2 <- c(1:dim(Decoys)[1])
  selected_colors <- get_colors(groups, gg_color_hue(K))

  svg(file = paste(plots_folder_name, "pca_cluster_plot.svg", sep = ""),
      width = 8, height = 8, pointsize = 10)
  scatter3D(decoys[, 1], decoys[, 2], decoys[, 3], bty = "g", pch = 20,
            cex = 0.9, main = "PCA clusters", phi = 10, theta = 30,
            colkey = FALSE, xlab = "PC1", ylab = "PC2", zlab = "PC3",
            col.var = as.integer(groups), labels = c(1:K),
            col = selected_colors)
  legend("topright", inset = 0.001, bty = "n", cex = 1, title = "Clusters",
         legend = c(1:K), fill = gg_color_hue(K))
  dev.off()

  svg(file = paste(plots_folder_name, "pca_time_plot.svg", sep = ""), width = 8,
      height = 8, pointsize = 10)
  scatter3D(decoys[, 1], decoys[, 2], decoys[, 3], bty = "g", pch = 20,
            cex = 0.8, main = "PCA time", phi = 10, theta = 30, expand = 0.9,
            colvar = groups2, clab = c("Frame no"),
            colkey = list(side = 4, length = 0.4), xlab = "PC1", ylab = "PC2",
            zlab = "PC3", col = matlab.like2(dim(Decoys)[1]),
            col.var = as.integer(1:dim(Decoys)[1]),
            labels = c(1:dim(Decoys)[1]))
  dev.off()
}

klasterOut <- function(K) {
  kMeansResult <<- kmeans(Decoys, K, iter.max = 100, nstart = 50,
                          algorithm = kmeans_algorithm)
  kol <<- kMeansResult$cluster
  klaMed <<- klasterMed(kol, nN, nC, K)

  if (plots) {
    selected_colors <- get_colors(kMeansResult$cluster, gg_color_hue(K))
    plotPCA(kMeansResult, Decoys[, 1:3], K)
    mypalette <- matlab.like2(length(Decoys[, 1]))
    cl <- kMeansResult
    m <- as.data.frame(Decoys)
    m$cluster <- factor(cl$cluster)
    m$time <- seq(1, length(Decoys[, 1]))
    centers <- as.data.frame(cl$centers)
    pca_plot1 <- ggplot(data = m, aes(x = V1, y = V2, color = time)) +
      geom_point(size = 1, alpha = 1/2) +
      theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
            axis.ticks = element_blank()) +
      scale_color_gradientn(colours = mypalette, breaks = c(0, 1),
                            labels = NULL) +
      labs(title = "PCA 2d time", x = "PC1", y = "PC2")
    savePlot(pca_plot1, "PCA_2dplot_PC1_PC2_time.pdf")

    pca_plot2 <- ggplot(data = m, aes(x = V2, y = V3, color = time)) +
      geom_point(size = 1, alpha = 1/2) +
      theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
            axis.ticks = element_blank()) +
      scale_color_gradientn(colours = mypalette, breaks = c(0, 1),
                            labels = NULL) +
      labs(title = "PCA 2d time", x = "PC2", y = "PC3")
    savePlot(pca_plot2, "PCA_2dplot_PC2_PC3_time.pdf")

    pca_plot3 <- ggplot(data = m, aes(x = V1, y = V3, color = time)) +
      geom_point(size = 1, alpha = 1/2) +
      theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
            axis.ticks = element_blank()) +
      scale_color_gradientn(colours = mypalette, breaks = c(0, 1),
                            labels = NULL) +
      labs(title = "PCA 2d time", x = "PC1", y = "PC3")
    savePlot(pca_plot3, "PCA_2dplot_PC1_PC3_time.pdf")

    pca_kmeans1 <- ggplot(data = m, aes(x = V1, y = V2, color = cluster)) +
      geom_point(size = 2, alpha = 1/2, col = selected_colors) +
      labs(title = "Kmeans plot", x = "PC1", y = "PC2")
    savePlot(pca_kmeans1, "PCA_2dplot_PC1_PC2_kmeans_groups.pdf")

    pca_kmeans2 <- ggplot(data = m, aes(x = V1, y = V3, color = cluster)) +
      geom_point(size = 2, alpha = 1/2, col = selected_colors) +
      labs(title = "Kmeans plot", x = "PC1", y = "PC3")
    savePlot(pca_kmeans2, "PCA_2dplot_PC1_PC3_kmeans_groups.pdf")

    pca_kmeans3 <- ggplot(data = m, aes(x = V2, y = V3, color = cluster)) +
      geom_point(size = 2, alpha = 1/2, col = selected_colors) +
      labs(title = "Kmeans plot", x = "PC2", y = "PC3")
    savePlot(pca_kmeans3, "PCA_2dplot_PC2_PC3_kmeans_groups.pdf")
  }

  if(save_files){
    if (file_type == "list") {
      write.table(data.frame(1:K, klaMed$files_names),
                  paste(reports_folder_name, "medoid_file.csv", sep = ""),
                  sep = "\t", eol = "\r\n", row.names = FALSE,
                  col.names = c("medoid", "file_name"))
      write.table(data.frame(as.data.frame(file_list), kMeansResult$cluster),
                  paste(reports_folder_name, "file_cluster_assignment.csv",
                        sep = ""),
                  sep = "\t", eol = "\r\n", row.names = FALSE,
                  col.names = c("file_name", "cluster"))
    } else {
      write.table(data.frame(1:K, klaMed$files_names),
                  paste(reports_folder_name, "medoid_frame.csv", sep = ""),
                  sep = "\t", eol = "\r\n", row.names = FALSE,
                  col.names = c("medoid", "frame"))
      write.table(data.frame(1:dim(trajectory)[3], kMeansResult$cluster),
                  paste(reports_folder_name, "frame_cluster_assignment.csv",
                        sep = ""),
                  sep = "\t", eol = "\r\n",
                  row.names = FALSE, col.names = c("frame", "cluster"))
    }

    for (i in 1:K) {
      if (file_type == "list") {
        model2pdb(SEQ_all, trajectory_all[ , , klaMed$files_rebuild_indexes[[i]]], bFactM,
                  paste(medoid_folder_name, "medoid_", i, ".pdb", sep = ""))
      } else {
        model2pdb(SEQ_all, trajectory_all[ , , klaMed$files_names[[i]]], bFactM,
                  paste(medoid_folder_name, "medoid_", i, ".pdb", sep = ""))
      }
      model2pdb(SEQ, klaMed$cen[[i]], bFact,
                paste(centroid_folder_name, "centroid_", i, ".pdb", sep = ""))
    }
  }

  if (plots) {
    pdf(paste(plots_folder_name, "trajectory.pdf", sep = ""))
    par(las = 1)
    plot(kol, pch = ".", cex = 3, xlab = "decoy's number", ylab = "cluster",
         main = "Mixing trajectories in the clusters")
    for (i in 0:1) {
      abline(v = i * nDec0, col = 2)
    }
    graphics.off()
  }

  Cmed <- matrix(0, K + 1, K + 1)
  for (i in 1:(K - 1)) {
    for (j in (i + 1):K) {
      Cmed[i, j] <- rmsd(klaMed$med[[i]][nN:nC, ], klaMed$med[[j]][nN:nC, ])
    }
  }
  for (i in 1:K) {
    Cmed[i, K + 1] <- rmsd(klaMed$med[[i]][nN:nC, ], YcenRMS[nN:nC, ])
  }
  Cmed <- round(Cmed + t(Cmed), 2)
  Ccen <- matrix(0, K + 1, K + 1)

  for (i in 1:(K - 1)) {
    for (j in (i + 1):K) {
      Ccen[i, j] <- rmsd(klaMed$cen[[i]][nN:nC, ], klaMed$cen[[j]][nN:nC, ])
    }
  }
  for (i in 1:K) {
    Ccen[i, K + 1] <- rmsd(klaMed$cen[[i]][nN:nC, ], YcenRMS[nN:nC, ])
  }
  Ccen <- round(Ccen + t(Ccen), 2)
  Ccen_med <- matrix(0, K, K)

  for (i in 1:K) {
    for (j in 1:K) {
      Ccen_med[i, j] <- rmsd(klaMed$cen[[i]][nN:nC, ], klaMed$med[[j]][nN:nC, ])
    }
  }
  Ccen_med <- round(Ccen_med, 2)

  MI <- mutInf(kol, traj)
  MM <- cbind(c(1:K), round(klaMed$nClust/klaMed$mRMS, 1), klaMed$nClust,
              round(klaMed$mRMS, 2), Cmed[1:K, K + 1], Ccen[1:K, K + 1],
              diag(Ccen_med), klaMed$rmsdRef)
  MM <- as.data.frame(MM)
  colnames(MM) <- c("cluster", "density", "cardinality", "<RMSD>",
                    "(med_i,cen)", "(cen_i,cen)", "(cen_i,med_i)", "<refRMDS>")
  raport_name_csv <- "report.csv"
  raport_name_txt <- "report_with_legend.txt"
  if(save_files){
    write.table(MM[order(-MM$density), ],
                paste(reports_folder_name, raport_name_csv, sep = ""),
                sep = "\t", eol = "\r\n", row.names = FALSE, col.names = TRUE)

    sink(paste(reports_folder_name, raport_name_txt, sep = ""))
    cat("Global information")
    cat("\nrmsd(cen,med)=", round(rmsd(YcenRMS[nN:nC, ], YmedRMS[nN:nC, ]), 2),
        " <RMSD>_cen=", round(mean(rCen), 2))
    cat("\nseparability(", K, ")=", round(sepK[K], 2), " Cxy=", round(MI$Cxy, 2),
        "\n")
    cat("\nCluster info - (x,y) is a rmsd shortcut(x,y)\n")
    print(MM)
    sink()
    }
  }

savePlot <- function(myPlot, filename) {
  pdf(paste(plots_folder_name, filename, sep = ""))
  print(myPlot)
  dev.off()
}

make_centroids_medoids <- function() {
  YcenRMS <<- trajectory[, , 1]
  nDec <<- dim(trajectory)[3]
  nDec0 <<- dim(trajectory)[3]
  traj <<- rep(1:1, each = nDec0)

  for (iDec in 2:nDec) {
    YcenRMS <<- ((iDec - 1) * YcenRMS + trajectory[, , iDec] %*%
                   rmsdQ(trajectory[nN:nC, , iDec], YcenRMS[nN:nC, ]))/iDec
  }
  rCen <<- apply(trajectory, 3, function(Y) rmsd(YcenRMS[nN:nC, ], Y[nN:nC, ]))
  YmedRMS <<- trajectory[, , which.min(rCen)]
  rMed <- apply(trajectory, 3, function(Y) rmsd(YmedRMS[nN:nC, ], Y[nN:nC, ]))

  bFact <<- 0
  bFactM <<- 0
  for (iDec in 1:nDec) {
    Y <- trajectory[, , iDec]
    bFact <<- bFact + apply(
      (YcenRMS - Y %*% rmsdQ(Y[nN:nC, ], YcenRMS[nN:nC, ]))^2, 1, sum)
    bFactM <<- bFactM + apply(
      (YmedRMS - Y %*% rmsdQ(Y[nN:nC, ], YmedRMS[nN:nC, ]))^2, 1, sum)
  }
  YmedRMS <<- trajectory_all[, , which.min(rCen)]
  bFact <- sqrt(bFact/nDec)
  bFactM <- sqrt((3*bFact)/(8*pi*pi))

  bFactM <- sqrt(bFactM/nDec)
  bFactM <- sqrt((3*bFactM)/(8*pi*pi))
}

write_medoids_centroid <- function() {
  if(save_files){
    model2pdb(SEQ, YcenRMS, bFact,paste(centroid_folder_name, "centroid.pdb",
                                        sep = ""))
    model2pdb(SEQ_all, YmedRMS, bFactM, paste(medoid_folder_name, "medoid.pdb",
                                          sep = ""))
    col_set<-c("blue", "red")
    if (plots) {
      pdf(paste(plots_folder_name, "RMSF.pdf", sep = ""))
      par(mar = c(6, 6, 6, 6), mgp = c(4.5, 1, 0), las = 1)
      matplot(SEQ[, 1], cbind(bFact, bFactM), type = "l", lty = 1,
              xlab = "Residue number", ylab = "Prediction of RMSF",
              main = "RMSF plot", cex.lab = 1, cex.main = 1, cex.axis = 1, col = col_set)
      legend("top", legend = c("centroid", "medoid"), lty=c(1,1), lwd=c(2.5,2.5), col = col_set, cex = 0.7, box.lty=0)
      graphics.off()
    }
  }
}

pca_reduction <- function() {
  Decoys <- matrix(NA, (nC - nN + 1) * 4, nDec)
  for (iDec in 1:nDec) {
    Y <- trajectory[nN:nC, , iDec]
    Y <- prcomp(Y)$x
    r2 <- apply(Y^2, 1, sum)
    Decoys[, iDec] <- c(r2, as.vector(t(Y) * sign(apply(YcenRMS[nN:nC, ] * Y,
                                                        2, sum))))
  }
  Decoys <- t(Decoys - apply(Decoys, 1, mean))
  ss <<- svd(Decoys, nu = 100, nv = 0)
  rm(Decoys)
  sig2 <- ss$d^2
  cumSig2 <- cumsum(sig2/sum(sig2))
  nu <- min(c(which(cumSig2 > 0.99), 100))
  Decoys <<- t(t(ss$u[, 1:nu]) * ss$d[1:nu])
}

k_means <- function() {
  if (plots) {
    t0 <- sum(diag(var(Decoys) * (nDec - 1)/nDec))
    kN <<- 40
    for (k in 2:kN) {
      kM <- kmeans(Decoys, k, iter.max = 50, nstart = 25,
                   algorithm = kmeans_algorithm)
      t0 <- c(t0, sum(kM$with/nDec))
    }
    sepK <<- 1 - t0/t0[1]
    pdf(paste(plots_folder_name, "separability.pdf", sep = ""))
    plot(1:kN, sepK, type = "l", lty = 1, xlab = "Cluster number",
         ylab = "Separability", main = "Separability plot", cex.lab = 1,
         cex.main = 1, cex.axis = 1)
    ss <- smooth.spline(1:kN, sepK, df = 30)
    graphics.off()
  }
}

k_means_cluster <- function() {
  d_sepK <- (ss$y[-1] - ss$y[-kN]) * dim(trajectory)[3]
  KL <- min(c(which(d_sepK < 10), 20))
  KU <- min(c(which(d_sepK < 1), 20))
  K <- ceiling((KL + KU)/2)
  if (cluster_quantity) {
    K <- cluster_quantity
  }
  klasterOut(K)
}
