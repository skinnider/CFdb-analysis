#' Calculate a random set of features for a given CF-MS matrix, and return them
#' as a pairwise data frame.
get_features = function(mat0, 
                        n_features,
                        matrix_dir,
                        quant_mode = 'iBAQ',
                        feature_grid = NULL
                        ) {
  if (is.null(feature_grid)) {
    # define grid of possible feature combinations
    feature_grid = get_feature_grid()
    
    print(str(feature_grid))
    
    # sample a random set of them
    feature_grid %<>% 
      sample_frac(1)
    
    print(head(feature_grid, n_features))
    print(n_features)
  } else {
    # get features for a predefined grid, passed as argument
    print(str(feature_grid))
  }
  
  # get the matrices
  matrices = list()
  for (grid_idx in seq_len(nrow(feature_grid))) {
    # check if we can stop
    if (length(matrices) == n_features) {
      message("-> got ", n_features, " features: breaking loop ...")
      break
    }
    
    message("  creating matrix: ", paste0(unlist(feature_grid[grid_idx, ]),
                                          collapse = ' / '))
    metric = feature_grid$metric[grid_idx]
    transform = feature_grid$transform[grid_idx]
    normalize = feature_grid$normalize[grid_idx]
    missing = feature_grid$missing[grid_idx]
    
    # make a copy of the matrix
    mat = mat0
    
    key = paste(metric, transform, normalize, missing, sep = '_')
    tryCatch({
      if (metric %in% c('MI', 'wccor', 'zi_kendall', 'GENIE3', 'distance_cor') &
          !is.null(matrix_dir)) {
        # read from file 
        matrix_file = file.path(matrix_dir, 
                                paste0(quant_mode,
                                       '-metric=', metric,
                                       "-scale=", FALSE,
                                       "-transform=", transform,
                                       "-normalize=", normalize,
                                       "-missing=", missing, 
                                       ".rds"))
        cor = readRDS(matrix_file)
        
        # optionally, we may need to map protein groups to genes
        id_map = attr(mat, 'id_map')
        if ('gene_map' %in% names(attributes(mat))) 
          id_map = attr(mat, 'gene_map')
        if (!is.null(id_map)) {
          cor %<>% extract(id_map$protein_group, id_map$protein_group)
          dimnames(cor) = list(id_map[[1]], id_map[[1]])
        }
      } else {
        ## calculate matrix on the fly
        # step 1: impute missing values (optional)
        if (missing == 'zero') {
          mat[is.na(mat)] = 0
        } else if (missing == 'noise') {
          mat[mat == 0] = NA
          mat %<>% clean_profiles(impute_NA = TRUE, smooth = FALSE)
        } else if (missing == 'NA') {
          mat[mat == 0] = NA
        }
        
        # step 2: log-transform (optional)
        if (transform == "log") {
          mat %<>% log()
          if (missing == 'zero') {
            mat[mat == -Inf] = 0
          }
        }
        # skip step 3 (scaling)
        
        # step 4: quantile normalize (optional)
        if (normalize == 'quantile') {
          dims = dimnames(mat)
          mat %<>% preprocessCore::normalize.quantiles()
          dimnames(mat) = dims
        }
        
        # step 5: calculate pairwise interactions
        cor = try(get_coexpr(t(mat), metric))
      }
      
      # throw error for treeClust
      if ("try-error" %in% class(cor))
        stop("cor has class try-error")
      
      # append to list
      matrices[[key]] = cor
    }, error = function(e) {
      message("error generating matrix ", key, ": ", e)
    })
  }
  
  # melt and join the matrices
  features = map2(matrices, names(matrices), ~
                    reshape2::melt(.x, varnames = c('protein_A', 'protein_B'),
                                   value.name = .y, as.is = TRUE)) %>%
    Reduce(function(x, y) full_join(x, y, by = c('protein_A', 'protein_B')), .) 
  
  # filter to unique pairs
  features %<>% filter(protein_A < protein_B)
  
  # return data frame
  return(features)
}

#' Define the grid of possible feature combinations to calculate
get_feature_grid = function() {
  # establish grid of analyses
  source("R/functions.R") ## contains metrics used to predict interactions
  opts = list(
    metric = metrics %>% setdiff(c('partial_cor', 'count')),
    scale = FALSE, # c(T, F),
    transform = c('none', 'log'),
    normalize = c('quantile', 'none'),
    missing = c('NA', 'zero', 'noise')
  )
  grid = do.call(expand.grid, c(opts, stringsAsFactors = F)) %>%
    # pick only one transformation
    filter(
      # log-transform
      (scale == FALSE & transform == 'log' & normalize == 'none') |
        # quantile normalize
        (scale == FALSE & transform == 'none' & normalize == 'quantile') |
        # scale
        (scale == TRUE & transform == 'none' & normalize == 'none') |
        # ... or none of them
        (scale == FALSE & transform == 'none' & normalize == 'none')) %>%
    # filter some combinations that cannot deal with NAs
    filter(!(missing == 'NA' & 
               metric %in% c('phi_s', 'rho_p', 'partial_cor', 'GENIE3',
                             'distance_cor', 'cosine', 'wccor'))) %>%
    # filter some combinations for which zeros and NAs are identical
    filter(!(missing == 'NA' &
               metric %in% c('dice', 'hamming', 'jaccard', 'zi_kendall',
                             'binomial', 'bayes_cor'))) %>%
    # don't impute missing values with noise for co-occurrence metrics
    filter(!(missing == 'noise' & 
               metric %in% c('dice', 'hamming', 'jaccard', 'zi_kendall',
                             'binomial'))) %>%
    # proportionality does not work with log-transform and near-zero noise
    filter(!(metric %in% c("phi_s", 'rho_p') & transform == 'log' &
               missing == 'noise')) %>%
    # ignore scale
    dplyr::select(-scale)
  return(grid)
}
