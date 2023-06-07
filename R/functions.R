get_coexpr = function(mat, metric) {
  # set seed for reproducibilty of some metrics, e.g., GENIE3
  set.seed(0)
  
  if (metric %in% c("pearson", "spearman", "kendall")) {
    cor = stats::cor(mat, method = metric, use = 'pairwise.complete.obs')
  } else if (metric == 'bicor') {
    cor = WGCNA::bicor(mat, use = 'pairwise.complete.obs')
  } else if (metric == 'binomial') {
    cor = dismay::binomial(mat)
  } else if (metric == 'cosine') {
    cor = lsa::cosine(mat)
  } else if (metric == 'jaccard') {
    cor = dismay::jaccard(mat)
  } else if (metric == 'canberra') {
    cor = -1.0 * as.matrix(dist(t(mat), method = 'canberra'))
  } else if (metric == 'euclidean') {
    cor = -1.0 * as.matrix(dist(t(mat), method = 'euclidean'))
  } else if (metric == 'manhattan') {
    cor = -1.0 * as.matrix(dist(t(mat), method = 'manhattan'))
  } else if (metric == 'weighted_rank') {
    cor = dismay::wtd_rank(mat)
  } else if (metric == 'hamming') {
    cor = -1.0 * dismay::hamming(mat)
  } else if (metric == 'dice') {
    cor = -1.0 * dismay::dice(mat)
  } else if (metric == 'phi_s') {
    cor = -1.0 * propr::phis(mat, select = colnames(mat))@matrix
  } else if (metric == 'rho_p') {
    cor = propr::perb(mat, select = colnames(mat))@matrix
  } else if (metric == 'zi_kendall') {
    cor = dismay::kendall_zi(mat)
  } else if (metric == 'MI') {
    cor = WGCNA::mutualInfoAdjacency(mat)$AdjacencyUniversalVersion1
    rownames(cor) = colnames(cor) = colnames(mat)
  } else if (metric == 'count') {
    cor = crossprod(!is.na(mat)) / nrow(mat)
  } else if (metric == 'bayes_cor') {
    source("R/util/BayesianCorrelation.R")
    cor = Bayes_Corr_Prior3(t(mat))
  } else if (metric == 'bray_curtis') {
    cor = -1.0 * as.matrix(vegan::vegdist(t(mat), method = 'bray',
                                          na.rm = T))
  } else if (metric == 'wccor') {
    permutations = gtools::permutations(n = ncol(mat), r = 2, 
                                        v = seq_len(ncol(mat)))
    cor = matrix(NA, nrow = ncol(mat), ncol = ncol(mat))
    cor[permutations] = map2_dbl(permutations[, 1], permutations[, 2],
                                 ~ wccsom::wcc(mat[, .x], mat[, .y], 
                                               trwdth = 1))
    rownames(cor) = colnames(cor) = colnames(mat)
  } else if (metric == 'distance_cor') {
    cor = Pigengene::dcor.matrix(mat)
  } else if (metric == 'profile_cor') {
    cor1 = stats::cor(mat, method = 'p', use = 'pairwise.complete.obs')
    cor = stats::cor(cor1, use = 'pairwise.complete.obs')
  } else if (metric == 'partial_cor') {
    cor = ppcor::pcor(mat)$estimate
    rownames(cor) = colnames(cor) = colnames(mat)
  } else if (metric == 'GENIE3') {
    cor = GENIE3::GENIE3(t(mat))
  } else if (metric == 'treeClust') {
    input = mat %>%
      t() %>%
      as.data.frame() %>%
      set_colnames(paste0('fraction', seq_len(ncol(.))))
    dist = treeClust::treeClust(dfx = input, d.num = 2, 
                                control = treeClust::treeClust.control(
                                  return.dists = T))$dists
    cor = -1.0 * as.matrix(dist)
  } else {
    stop("don't know what to do for metric: ", metric)
  }
  
  return(cor)
}

# define column names for a GOA file in GAF format
gaf_colnames = c('DB',
                 'DB Object ID',
                 'DB Object Symbol',
                 'Qualified',
                 'GO ID',
                 'DB:Reference',
                 'Evidence Code',
                 'With (or) From',
                 'Aspect',
                 'DB Object Name',
                 'DB Object Synonym',
                 'DB Object Type', 
                 'Taxon',
                 'Date',
                 'Assigned By',
                 'Annotation Extension',
                 'Gene Product Form ID')
