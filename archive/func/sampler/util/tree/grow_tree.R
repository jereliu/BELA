# grow tree
init_tree_rrt <- function(U_old){
  root_tree <- as.data.frame(t(U_old))
  root_tree$id <- 1
  root_tree$parent <- 0
  
  root_tree
}

prop_dir_rrt <- 
  function(U_old){
    # sample randomly in space to propose a new direction
    dim1 <- 1
    dim2 <- length(U_old)
    U_new <- 
      rmvnorm(1, as.matrix(U_old), #rep(0, dim2), 
              diag(dim2)*1e10)
    U_new
  }

prop_node_line <- 
  function(new_node_count, tree, tree_dir, step_dir){
    root_idx <- 
      ifelse(tree_dir == 1, nrow(tree), 1)
    root_node <- tree[root_idx, ]
      
    # propose node 
    step_list <- 
      sapply(1:new_node_count, 
             function(count) count*step_dir)
    
    node_prop <- (root_node + tree_dir*step_list) %>% 
      t %>% matrix(nrow = new_node_count)
    
    node_prop
  }

prop_node_rrt <- 
  function(U_new, tree, edge_max = 0.1, random = FALSE){
    dim2 <- dim(U_new)[2]
    
    if (random){
      parent_id <- sample(nrow(tree), 1)
      parent_coord <- tree[parent_id, 1:dim2]
      grow_dir <- 
        matrix(rnorm(k), nrow = 1) %>% 
        (function(x) x/sqrt(sum(x^2)))
    } else {
      # compute cost
      cost <- 
        apply(as.matrix(tree[, 1:dim2]), 1, 
              function(x) sqrt(sum((x - U_new)^2))
        )
      
      # choose node to grow
      parent_id <- which.min(cost)
      parent_coord <- tree[parent_id, 1:dim2]
      grow_dir <- 
        (parent_coord - U_new) %>% 
        (function(x) x/sqrt(sum(x^2)))
    }
    grow_dist <- runif(1) * edge_max
    
    # propose node 
    node_prop <- parent_coord - grow_dir * grow_dist
    
    out <- data.frame(node_prop, 
                      id = nrow(tree) + 1, 
                      parent = parent_id)
    names(out) <- names(tree)
    out
  }

flex_node <- 
  function(new_node, parent_node, prop_grad){
    coord_id <- 1:(length(parent_node)-2)
    
    grow_grad <- 
      (parent_node - new_node)[coord_id]
    grow_dist <- sqrt(sum(grow_grad^2))
    grow_dir <- as.matrix(grow_grad/grow_dist)
    prop_dir <- prop_grad/sqrt(sum(prop_grad^2))
    
    flex_dir <- 
      grow_dir - 2 * prop_dir * sum(grow_dir * prop_dir) 
    #
    flex_node_prop <- new_node
    flex_node_prop[coord_id] <- 
      new_node[coord_id] - flex_dir * grow_dist
    flex_node_prop
  }
