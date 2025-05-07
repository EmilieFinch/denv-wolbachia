#### Script with code to output chain binomials for odin model 

## For S -> I

n_strains <- 60

for (j in 1:n_strains) {
  if (j == 1) {
    cat(sprintf(
      "n_SI[, %d] <- if (sum(lambda[%d:%d]) > 0) Binomial(S_out[i], lambda[%d] / sum(lambda[%d:%d])) else 0\n",
      j, j, n_strains, j, j, n_strains))
  } else if (j < n_strains) {
    cat(sprintf(
      "n_SI[, %d] <- if (sum(lambda[%d:%d]) > 0) Binomial(S_out[i] - sum(n_SI[i, 1:%d]), lambda[%d] / sum(lambda[%d:%d])) else 0\n",
      j, j, n_strains, j - 1, j, j, n_strains))
  } else {
    cat(sprintf(
      "n_SI[, %d] <- S_out[i] - sum(n_SI[i, 1:%d])\n",
      j, j - 1))
  }
}

## For R -> I

## Generate R_out
serotypes <- 1:4

# Immune histories (no naive state)
immune_histories <- list(
  c(1),         # 1
  c(2),         # 2
  c(3),         # 3
  c(4),         # 4
  c(1,2),       # 5
  c(1,3),       # 6
  c(1,4),       # 7
  c(2,3),       # 8
  c(2,4),       # 9
  c(3,4),       #10
  c(1,2,3),     #11
  c(1,2,4),     #12
  c(1,3,4),     #13
  c(2,3,4)      #14
)

code_lines <- c()

for (j in seq_along(immune_histories)) {
  immune_to <- immune_histories[[j]]
  susceptible <- setdiff(serotypes, immune_to)
  
  if (length(susceptible) == 0) {
    prob <- "0"
  } else {
    lambda_sum <- paste0("lambda_", susceptible, collapse = " + ")
    prob <- sprintf("1 - exp(- (%s)* dt)", lambda_sum)
  }
  
  line <- sprintf("R_out[, %d] <- Binomial(R[i, %d], %s)", j, j, prob)
  code_lines <- c(code_lines, line)
}

cat(paste(code_lines, collapse = "\n"))

## n_RI Generate chain binomial ## 

# Map serotypes to strain indices
serotype_strains <- list(
  `1` = 1:20,
  `2` = 21:40,
  `3` = 41:50,
  `4` = 51:60
)

# Add one more level to immune histories (15 total)
immune_histories[[15]] <- c(1, 2, 3, 4)  # 1234

for (h_idx in seq_along(immune_histories)) {
  immune <- immune_histories[[h_idx]]
  remaining_serotypes <- setdiff(1:4, immune)
  if (length(remaining_serotypes) == 0) next  # Fully immune
  
  eligible_strains <- unlist(serotype_strains[as.character(remaining_serotypes)])
  k <- length(eligible_strains)
  
  for (j_idx in seq_along(eligible_strains)) {
    j <- eligible_strains[j_idx]
    lambda_range <- sprintf("lambda[%d:%d]", eligible_strains[1], eligible_strains[k])
    
    if (j_idx == 1) {
      cat(sprintf(
        "n_RI[, %d, %d] <- if (sum(%s) > 0) Binomial(R_out[i, %d], lambda[%d] / sum(%s)) else 0\n",
        h_idx, j, lambda_range, h_idx, j, lambda_range
      ))
    } else if (j_idx < k) {
      cat(sprintf(
        "n_RI[, %d, %d] <- if (sum(%s) > 0) Binomial(R_out[i, %d] - sum(n_RI[i, %d, %d:%d]), lambda[%d] / sum(%s)) else 0\n",
        h_idx, j, lambda_range, h_idx, h_idx, eligible_strains[1], eligible_strains[j_idx - 1], j, lambda_range
      ))
    } else {
      cat(sprintf(
        "n_RI[, %d, %d] <- R_out[i, %d] - sum(n_RI[i, %d, %d:%d])\n",
        h_idx, j, h_idx, h_idx, eligible_strains[1], eligible_strains[j_idx - 1]
      ))
    }
  }
}

## For n_IC 

# I histories (with naive)
I_histories <- list(
  c(), c(1), c(2), c(3), c(4),
  c(1,2), c(1,3), c(1,4),
  c(2,3), c(2,4), c(3,4),
  c(1,2,3), c(1,2,4), c(1,3,4), c(2,3,4),
  c(1,2,3,4)
)

# C histories (no naive)
C_histories <- list(
  c(1), c(2), c(3), c(4),
  c(1,2), c(1,3), c(1,4),
  c(2,3), c(2,4), c(3,4),
  c(1,2,3), c(1,2,4), c(1,3,4), c(2,3,4),
  c(1,2,3,4)
)

key <- function(set) paste(sort(set), collapse = "")
C_map <- setNames(seq_along(C_histories), sapply(C_histories, key))

sero_ranges <- list(
  `1` = "1:20",
  `2` = "21:40",
  `3` = "41:50",
  `4` = "51:60"
)

C_contributions <- vector("list", length(C_histories))

for (i in seq_along(I_histories)) {
  existing <- I_histories[[i]]
  possible_serotypes <- setdiff(1:4, existing)
  
  for (sero in possible_serotypes) {
    new_set <- union(existing, sero)
    new_key <- key(new_set)
    if (length(new_set) == 0) next
    
    c_index <- C_map[[new_key]]
    strain_slice <- sero_ranges[[as.character(sero)]]
    expr <- sprintf("sum(I_out[i, %d, %s])", i, strain_slice)
    
    C_contributions[[c_index]] <- c(C_contributions[[c_index]], expr)
  }
}

code_lines <- mapply(function(exprs, c_index) {
  if (length(exprs) == 0) return(NULL)
  joined <- paste(exprs, collapse = " + ")
  sprintf("n_IC[, %d] <- Binomial(%s, p_IC)", c_index, joined, c_index)
}, C_contributions, seq_along(C_contributions))

cat(paste(code_lines[!sapply(code_lines, is.null)], collapse = "\n"))
