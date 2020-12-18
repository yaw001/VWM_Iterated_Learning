#Silhouette coefficient 
#a(o): the average distance between object o and all other objects in the cluster to which o belongs.
#b(0): the minimum average distance from o to all cluster to which o does not belong.
#Silhouette coefficient = (b(o) - a(o))/ max(a(o),b(o))

#TODO S_O on color iterated learning exp
dist_same <- function(object_same){
  dist_mat <- as.matrix(dist(object_same))
  a_o = apply(dist_mat, MARGIN = 1, sum)/(nrow(object_same)-1)
  return(a_o)
}

dist_diff <- function(object_diff) {
  dist_mat <- as.matrix(dist(object_diff))
  b_o = mean(dist_mat[,1])
  return(b_o)
}

sil_coef <- function(data) {
  num_cluster = data$assignment %>% unique() %>% length()
  dist_mat = dist(data %>% select(x,y)) %>% as.matrix()
  index_list = list()
  for(i in 1:num_cluster){
    index_list[[i]] = which(data$assignment == i)
  }
  a_o_list = list()
  b_o_list = list()
  dist_diff = list()
  for(j in 1:num_cluster){
    a_o_list[[j]] = apply(dist_mat[index_list[[j]],index_list[[j]]], MARGIN = 2, sum)/(length(index_list[[j]])-1)
    other = c(1:num_cluster)[-j]
    for(k in other){
      dist_diff[[which(other==k)]] = apply(dist_mat[which(data$assignment == k),index_list[[j]]], MARGIN = 2, mean)
    }
    b_o_list[[j]] = apply(matrix(unlist(dist_diff),ncol = length(other)),MARGIN = 1, min)
  }
  ab = matrix(c(unlist(a_o_list),unlist(b_o_list)),ncol = 2)
  ab_diff = ab[,2] - ab[,1]
  min_ab = apply(ab,MARGIN = 1,max)
  s_o = mean(ab_diff/min_ab)
  return(s_o)
}

#gibbs sampling assignment consistency rate 
