# Recursive Bi-partitioning
# Laplacian function
#' Title
#'
#' @param x is as square matrix from which we compute the Laplacian
#'
#' @return
#' @export
#'
#' @examples
lp_fn <- function(x) {
  Dinv <- Diagonal(nrow(x), 1 / Matrix::rowSums(x))

  lp_chr1 <- Diagonal(nrow(x), 1) - Dinv %*% x
  if (dim(lp_chr1)[1] > 10000) {
    return(eigs_sym(lp_chr1, k = 2, sigma = 0, which = "LM", maxitr = 10000))
  } else {
    temp <- eigen(lp_chr1)
    return(list(vectors = temp[["vectors"]][, c(length(temp$values) - 1, length(temp$values))], values = temp[["values"]][c(length(temp$values) - 1, length(temp$values))]))
  }
}
# Bipartition function
# delineate bi-partitions using kmeans and compute expansion
#' Title
#'
#' @param x is the eigen-decomposition step output
#' @param reff_g is the original graph on which we computed the eigen decomposition
#' @param tmp_res the resolution at which we are performing spectral clustering
#'
#' @return
#' @export
#'
#' @examples
part_cond_calc <- function(x, reff_g, tmp_res) {
  # perform kmeans with 2 clusters in first 2 smallest eigen vector space
  res <- kmeans(x$vectors, 2, nstart = 5)
  # calculate the expansion of the resulting k clusters
  l_temp_exp <- c()
  sub_g_list <- list()
  for (j in 1:2) {
    # create the subnetwork
    sub_g_temp <- induced_subgraph(reff_g, which(res$cluster == j))

    # save the members of considered cluster
    # create cluster label in considered subnetwork
    temp_name <- paste(tmp_res, length(V(sub_g_temp)$name), length(E(sub_g_temp)), min(as.numeric(unlist(lapply(strsplit(V(sub_g_temp)$name, ","), "[", 1)))), max(as.numeric(unlist(lapply(strsplit(V(sub_g_temp)$name, ","), "[", 1)))), sep = "_")
    sub_g_list[[temp_name]] <- V(sub_g_temp)$name
    rm(sub_g_temp)
  }
  return(sub_g_list)
}

# loop through all elements in nested list
ff <- function(x) {
  if (class(x) == "list" & length(x) > 0) {
    lapply(x, ff)
  } else {
    TRUE
  }
}

# recursive bi-partitioning
spec_bipart <- function(chr1_mat, g_chr1, res_set, res_num, chr_dat_l) {
  options(scipen = 999999999)
  # initialisation
  # container to save cluster hierarchy as list of lists
  chr1_tree_list <- list()
  # containers for cluster member list, cluster conductance/expansion
  chr1_tree_cl <- list()
  # compute the best resolution for considered chr given sparsity criteria
  tmp_res <- res_set[max(which(unlist(lapply(chr_dat_l, function(x) {
    nbin <- length(unique(c(x$X1, x$X2)))
    # 0.8 sparsity criteria
    return(nrow(x) - nbin > nbin * (nbin - 1) / 2 * 0.99)
  }))))]
  # if even the coarcest resolution doesn't fulffil the cirteria assign the starting resolution to it
  if (is.na(tmp_res)) {
    tmp_res <- res_set[1]
  }

  # whole chromosome laplacian
  lpe_chr1 <- lp_fn(chr1_mat)
  # spectral clusters
  res_chr1 <- part_cond_calc(lpe_chr1, g_chr1, tmp_res)

  # save cluster membership and expansion
  chr1_tree_cl <- list(chr1_tree_cl, res_chr1)
  chr1_tree_cl <- unlist(chr1_tree_cl, recursive = FALSE)

  # temporary list of candidate cluster to further partition
  ok_part <- names(chr1_tree_cl)

  # initiate the tree
  for (i in ok_part) {
    chr1_tree_list[[i]] <- list()
  }

  # save path to all considered leaves
  lnames <- names(unlist(lapply(chr1_tree_list, ff)))
  names(lnames) <- names(chr1_tree_cl)

  # given starting best resolution set the further resolutions to consider at later iterations
  tmp_res_set <- res_set[which(res_set == tmp_res):length(res_set)]
  # build a temporary list containing the subset of Hi-C interactions constitutive of the found clusters at all
  # higher resolutions for the consdiered cluster
  cl_dat_l <- vector("list", length(ok_part))
  for (cl in ok_part) {
    for (r in tmp_res_set) {
      # when considered resolution equals the original resolution
      if (r == tmp_res) {
        tmp_dat <- chr_dat_l[[r]] %>%
          filter(X1 %in% chr1_tree_cl[[cl]]) %>%
          filter(X2 %in% chr1_tree_cl[[cl]])
        cl_dat_l[[cl]][[r]] <- tmp_dat
      }
      # when resolution higher than the original resolution
      else {
        # create the higher resolution bins expecting a 0-start counting of bins
        r_bin <- unique(as.numeric(sapply(chr1_tree_cl[[cl]], function(x) {
          tmp <- seq(as.numeric(x), as.numeric(x) + res_num[tmp_res], by = res_num[r])
          return(tmp[-length(tmp)])
        })))
        # extract corresponding edgelist from Hi-C data
        tmp_dat <- chr_dat_l[[r]] %>%
          filter(X1 %in% as.character(r_bin)) %>%
          filter(X2 %in% as.character(r_bin))
        cl_dat_l[[cl]][[r]] <- tmp_dat
      }
    }
  }
  rm(lpe_chr1, res_chr1, r, r_bin, cl, tmp_res, tmp_res_set)
  # recursive looping
  while (length(ok_part) > 0) {
    temp_part <- c()
    tmp_cl_dat_l <- list()
    for (i in ok_part) {
      if (length(chr1_tree_cl[[i]]) < 3 & strsplit(i, split = "_")[[1]][[1]] == "5kb") {
        next
      }
      # set the original resolution to the resolution at the which the considered cluster was found
      tmp_res <- unlist(strsplit(i, split = "_"))[1]
      # extract the subset of Hi-C data at all required resolution for the considered cluster
      tmp_chr_dat_l <- cl_dat_l[[i]]
      # Only evaluate higher resolution if the considered cluster was found at coarser resolution than 5kb(highest resolution)
      if (tmp_res != "5kb") {
        print(paste(which(ok_part == i), "out of", length(ok_part)))
        # if cluster resolution > 5kb
        # check which resolution best for this cluster starting from the considered cluster resolution
        # vector containing the resolution selection criteria
        # only consider resolutions equal or higher than the considered cluster
        tmp_res_set <- res_set[which(res_set == tmp_res):length(res_set)]
        # create vector indicating which resolution adhere to the sparsity criteria
        res_select <- vector("logical", length(tmp_res_set))
        names(res_select) <- tmp_res_set
        # loop through selected resolutions
        for (r in tmp_res_set) {
          # when considered resolution equals the original resolution
          if (r == tmp_res) {
            tmp_dat <- tmp_chr_dat_l[[r]]
            nbin <- length(unique(c(tmp_dat$X1, tmp_dat$X2)))
            res_select[r] <- nrow(tmp_dat) - nbin > nbin * (nbin - 1) / 2 * 0.99
          }
          # when resolution higher than the original resolution
          else {
            # create the higher resolution bins expecting a 0-start counting of bins
            r_bin <- unique(as.numeric(sapply(chr1_tree_cl[[i]], function(x) {
              tmp <- seq(as.numeric(x), as.numeric(x) + res_num[tmp_res], by = res_num[r])
              return(tmp[-length(tmp)])
            })))
            # extract corresponding edgelist from Hi-C data
            tmp_dat <- tmp_chr_dat_l[[r]]
            nbin <- length(r_bin)
            res_select[r] <- nrow(tmp_dat) - nbin > nbin * (nbin - 1) / 2 * 0.99
          }
        }
        # if no higher or equal resolution is appropriate for further clustering keep original resolution
        if (sum(res_select) == 0) {
          tmp_res <- unlist(strsplit(i, split = "_"))[1]
        }
        # set best Hi-C resolution
        if (sum(res_select) > 0) {
          tmp_res <- names(res_select[max(which(res_select))])
        }
        # rm(res_select,tmp_dat,nbin,tmp_res_set,res_select,r)
      }
      # build the cluster edgelist at best resolution
      tmp_chr_dat <- tmp_chr_dat_l[[tmp_res]]

      if (length(unique(c(tmp_chr_dat$X1, tmp_chr_dat$X2))) < 3) {
        next
      }

      # create the subnetwork of considered cluster at best resolution
      sub_g1 <- graph_from_data_frame(tmp_chr_dat, directed = F)
      # eleminate self loop
      sub_g1 <- delete.edges(sub_g1, E(sub_g1)[which(which_loop(sub_g1))])
      # create the corresponding adjacency matrix
      sub_g1_adj <- get.adjacency(sub_g1, type = "both", attr = "weight")
      if (any(colSums(sub_g1_adj) == 0)) {
        out <- which(colSums(sub_g1_adj) == 0)
        sub_g1_adj <- sub_g1_adj[-out, ]
        sub_g1_adj <- sub_g1_adj[, -out]
      }
      if (nrow(sub_g1_adj) < 1) {
        next
      }

      print(paste("eigen decomposition:", i))
      lpe_sub_g1 <- lp_fn(sub_g1_adj)
      if (nrow(lpe_sub_g1$vectors) < 3) {
        next
      }
      print("kmeans")
      # find actual sub-structures using kmeans on fiedler vector
      res_subg1 <- part_cond_calc(lpe_sub_g1, sub_g1, tmp_res)

      print("cl append")
      for (k in names(res_subg1)) {
        chr1_tree_cl[[k]] <- res_subg1[[k]]
      }


      # Only consider future cluster partition if their expansion is majoritarily inside the cluster
      # ok_part_temp<-names(res_subg1[[1]])[which(res_subg1[[2]]<1)]
      ok_part_temp <- names(res_subg1)

      print("tree growth")
      for (j in ok_part_temp) {
        chr1_tree_list[[c(unlist(strsplit(lnames[i], split = "\\.")), j)]] <- list()
      }
      # update the temporary list of Hi-C tables containing the subset of interactions for the found clusters from
      # their respective parent cluster Hi-C table (much faster subsetting)
      tmp_res_set <- res_set[which(res_set == tmp_res):length(res_set)]
      for (cl in ok_part_temp) {
        for (r in tmp_res_set) {
          # when considered resolution equals the original resolution
          if (r == tmp_res) {
            tmp_dat <- tmp_chr_dat_l[[r]] %>%
              filter(X1 %in% chr1_tree_cl[[cl]]) %>%
              filter(X2 %in% chr1_tree_cl[[cl]])
            tmp_cl_dat_l[[cl]][[r]] <- tmp_dat
          }
          # when resolution higher than the original resolution
          else {
            # create the higher resolution bins expecting a 0-start counting of bins
            r_bin <- unique(as.numeric(sapply(chr1_tree_cl[[cl]], function(x) {
              tmp <- seq(as.numeric(x), as.numeric(x) + res_num[tmp_res], by = res_num[r])
              return(tmp[-length(tmp)])
            })))
            # extract corresponding edgelist from Hi-C data
            tmp_dat <- tmp_chr_dat_l[[r]] %>%
              filter(X1 %in% as.character(r_bin)) %>%
              filter(X2 %in% as.character(r_bin))
            tmp_cl_dat_l[[cl]][[r]] <- tmp_dat
          }
        }
      }
      temp_part <- c(temp_part, ok_part_temp)
      rm(sub_g1, sub_g1_adj, lpe_sub_g1, res_subg1, tmp_chr_dat)
    }
    print(length(temp_part))
    ok_part <- temp_part
    cl_dat_l <- tmp_cl_dat_l
    # gather updated paths
    lnames <- names(unlist(lapply(chr1_tree_list, ff)))
    # name each path according to leaf of considered path
    names(lnames) <- lapply(lnames, function(x) unlist(strsplit(x, split = "\\."))[length(unlist(strsplit(x, split = "\\.")))])

    rm(tmp_cl_dat_l, temp_part, k, nbin, tmp_chr_dat_l, tmp_dat, r_bin, tmp_res_set, tmp_res, r, res_select, ok_part_temp, i, j, cl)
  }

  return(list(cl_member = chr1_tree_cl, part_tree = chr1_tree_list))
}
