Overlap_Rec = function(rec1, rec2) 
{
  upmat <- rec1
  upmat[, 1] <- rec2[, 2]
  lowmat <- rec1
  lowmat[, 2] <- rec2[, 1]
  rec <- rec1
  rec[, 1] <- apply(lowmat, 1, max)
  rec[, 2] <- apply(upmat, 1, min)
  return(rec)
}



surface_funs <- function(rec, label, reclst0, labellst0, epsilon=10^(-12))
{
  # Processing all overlapping cells
  d <- dim(rec)[1]
  V <- prod(rec[,2] - rec[,1])
  S_faces <- numeric(length = d)
  
  # the overlapping surface between rec and other rectangles 
  # that are labeled 1, at two faces of feature j
  overlap <- matrix(0, d, 2)
  
  # sub_overlap is a list, with each element as [start, end, sub overlapping surface]
  sub_overlap <- vector("list", d) ## not same 
  
  for (i in 1:d)
  {
    sub_overlap[[i]] <- list()
    S_faces[i] <- V / (rec[i,2] - rec[i,1])
  }
  S <- sum(S_faces) * 2
  ans <- vector("list", d+1) 
  
  # If reclst is empty
  if (length(labellst0) == 0)
  {
    for (i in 1:d)
    {
      intercepts10 <- 2*S_faces[i]
      slopes10 <- (S - S_faces[i]*2) / (rec[i,2] - rec[i,1])
      ans[[i]] <- c(rec[i,1], slopes10, intercepts10, S_faces[i])
    }
    ans[[d+1]] <- c(0, S)
    return(ans)
  }
  else
  {
    for(i in 1:length(labellst0)) #67
    {
      if(labellst0[i] == 0)
      {next}
      recnow <- reclst0[i] 
      recnow <- recnow[[1]]
      contact_feat <- -1
      for(j in 1:d) #72
      {
        if(rec[j,1] == recnow[j,2])
        {
          contact_feat <- j
          contact_direct <- 1
          break
        }
        else if (rec[j,2] == recnow[j,1])
        {
          contact_feat <- j
          contact_direct <- 2
          break
        }
      }
      if(contact_feat == -1)
      {next}
      overlap_rec <- Overlap_Rec(rec, recnow)
      overlap_rec_del <- overlap_rec[-contact_feat,]
      if(min(overlap_rec_del[,2]-overlap_rec_del[,1]) <= 0)
      {next}
      overlap_V <- prod(overlap_rec_del[,2] - overlap_rec_del[,1])
      overlap[contact_feat, contact_direct] <-  overlap_V + overlap[contact_feat, contact_direct] # 88
      feats <- 1:d
      feats <- feats[-contact_feat]
      
      for(j in feats)
      {
        sub_overlap[[j]][1] <- overlap_rec[j,1]
        sub_overlap[[j]][2] <- overlap_rec[j,2]
        sub_overlap[[j]][3] <- overlap_V/(overlap_rec[j,2]-overlap_rec[j,1])
      }
    }
  }
  
  # Compute piecewise linear functions with overlapping information
  s_0 <- sum(overlap)
  s_1 = S - s_0
  ans[[d+1]] <- c(s_0, s_1)
  
  for(j in 1:d)
  {
    sub_overlap_j = sub_overlap[[j]]
    if(length(sub_overlap_j) == 0)
    {
      intercepts10 <- (s_0 - overlap[j,1] + S_faces[j] - overlap[j,1] + S_faces[j])
      slopes10 = ((S - S_faces[j]*2) / (rec[j,2] - rec[j,1]))
      ans[[j]] = c(rec[j,1], slopes10, intercepts10, S_faces[j]) 
      next
    }
    
    # both slopes_changes and slopes depicts slopes 
    # overlapping with elements of reclst with label 1
    slopes_changes <- matrix(0, 2*length(sub_overlap_j), 2)
    sidelen_j = rec[j,2] - rec[j,1]
    
    for(i in 1:length(sub_overlap_j))
    {
      slopes_changes[i+i-2,] <- c(sub_overlap_j[i][1], sub_overlap_j[i][3])
      slopes_changes[i+i+1-2,] <- c(sub_overlap_j[i][2], -sub_overlap_j[i][3])
    }
    slopes_changes <- slopes_changes[order(slopes_changes[,1]),]
    checkpoints <- list()
    slopes <- list()
    value[j] <- rec[j,1]
    slope_all <- (S - S_faces[j]*2) / sidelen_j
    sl <- 0
    
    for(k in 1:length(slopes_changes)) #115
    {
      if((abs(slopes_changes[k,1]-value[k])) < epsilon)
      {
        sl <- sl + slopes_changes[k,2]
      }
      else
      {
        checkpoints <- c(checkpoints,value[k])
        value[k] <- slopes_changes[k,1]
        slopes <- c(slopes, sl)
        sl <- sl + slopes_changes[k,2]
        
        if(abs(slopes_changes[k,1]-rec[j,2]) < epsilon)
        {break}
      }
    }
    
    tryCatch(
      {
        if(length(checkpoints) == 0)
        {
          intercepts10 = (s_0 - overlap[j,1] + S_faces[j] - overlap[j,1] + S_faces[j])
          slopes10 = ((S - S_faces[j]*2) / (rec[j,2] - rec[j,1]))
          ans[[j]] = c(rec[j,1], slopes10, intercepts10, S_faces[j])   
          next
        }
        if(abs(checkpoints[length(checkpoints)]-value) >= epsilon)
        {
          checkpoints <- c(checkpoints,value)
          slopes <- c(slopes, sl)
        }
      }, 
      error= function(e) {
        browser()
        debug <- checkpoints}
    )
    slopes10 = slope_all - 2*as.numeric(slopes)
    intercepts10 = numeric(length = length(checkpoints))
    intercepts10[1] = s_0 - overlap[j,1] + S_faces[j] - overlap[j,1] + S_faces[j]
    
    for(k in 2:length(checkpoints))
    {
      if(checkpoints[k] < checkpoints[k-1])
      {
        print(paste("Error: invalid checkpoints: ",(checkpoints[k])))
      }
      intercepts10[k] = intercepts10[k-1] + slopes10[k-1]*(checkpoints[k]-checkpoints[k-1])
    }
    ans[[j]] = c(checkpoints, slopes10, intercepts10, S_faces[j]) 
  }
  return(ans)
}

data_standardize <- function(X)
{
  di <- dim(X)
  border <- matrix(rep(0, di[2]*2), di[2], 2)
  feat_min <- numeric(di[2])
  feat_max <- numeric(di[2])
  border_dist <- numeric(di[2])
  for(j in 1:di[2])
  {
    feat_min[j] = min(X[,j])
    feat_max[j] = max(X[,j])
    
    if(feat_max[j] == feat_min[j])
    {
      print(paste("feature",j,"has only one value"))
    }
    else
    {
      border_dist[j] <- (feat_max[j]-feat_min[j])/(di[1]-1)
      border[j,1] <- feat_min[j]-border_dist[j]
      border[j,2] <- feat_max[j]+border_dist[j]
    }
  } 
  shifts <- -border[,1]
  multipliers <- diag((1/(border[,2]-border[,1])), di[2], di[2])
  new <- matrix(rep(shifts,di[1]), di[1], di[2], byrow=TRUE)
  standardize_para <- list(shifts, multipliers)
  rtn <- (as.matrix(X) + new)%*% multipliers
  return(rtn)
}

sv_regular <- function(surface, volume, d)
{
  return (surface/volume)
}

Compute_Impu <- function(wy, w)
{
  rtn <- 1 - (wy/w)^2 - ((w-wy)/w)^2
  return(rtn)
}

Compute_SignImpu <- function(wy, w, label=1)
{
  if(as.integer(wy/w>=0.5) == label)
  {
    return (1 - (wy/w)^2 - ((w-wy)/w)^2)
  }
  else{return (wy/w)^2 + ((w-wy)/w)^2} 
}

Compute_NodeImpu <- function(wyleft, wleft, wy, w)
{
  rtn <- (1- ((wyleft/wleft)^2 + ((wleft-wyleft)/wleft)^2)*wleft/w 
          - (((wy-wyleft)/(w-wleft))^2 + ((w-wleft-wy+wyleft)/(w-wleft))^2)*(w-wleft)/w) 
  return(rtn)
}

Compute_SignNodeImpu <- function(wyleft, wleft, wy, w, child_labels)
{
  impu_left = 1 - (wyleft/wleft)^2 - ((wleft-wyleft)/wleft)^2
  impu_right = 1 - ((wy-wyleft)/(w-wleft))^2 - ((w-wleft-wy+wyleft)/(w-wleft))^2
  
  if(as.integer(wyleft/wleft >= 0.5) == child_labels[1])
  {impu_left_sign = impu_left}
  else
  {impu_left_sign = 1 - impu_left}
  
  if(as.integer((wy-wyleft)/(w-wleft)>=0.5) == child_labels[2])
  {impu_right_sign = impu_right}
  else
  {impu_right_sign = 1 - impu_right}
  
  return(impu_left_sign*wleft/w + impu_right_sign*(w-wleft)/w)
}



fit_sv <- function(X, Y, pen, c0=1, weight=1, border=NULL, standardize=FALSE, 
                   criterion='gini', min_split_weight=NULL, min_leaf_weight=NULL, tol=10^(-10), maximal_leaves=NULL)
{
  n <- dim(X)[1]
  d <- dim(X)[2]
  if(is.null(border) == TRUE)
  {
    border <- matrix(0,d,2)
    border[,2] <- rep(1,d)
  }
  if(standardize==FALSE)
  {
    X <- data_standardize(X)
  }
  if(is.null(min_split_weight) == TRUE)
  {
    min_split_weight <- weight+1
  }
  if(is.null(min_leaf_weight) == TRUE)
  {
    min_leaf_weight <- weight+1
  }
  if(is.null(maximal_leaves) == TRUE)
  {
    maximal_leaves <- floor(sqrt(n))
  }
  wn_all <- length(Y) + (weight-1)*sum(Y)
  wn <- wn_all
  wy <- weight*sum(Y)
  impu <- Compute_Impu(wy, wn)
  class_label <- 1
  sign_impu <- Compute_SignImpu(wy, wn, class_label)
  tree_impu <- impu
  tree_sign_impu <- sign_impu
  
  volume <- prod(border[,2] - border[,1]) 
  surface <- 0
  for(j in 1:d)
  {
    surface <- surface + (2 * volume / (border[j,2]-border[j,1]))
  }
  sv_reg_min <- sv_regular(surface, volume, d)
  risk <- tree_impu + (pen * sv_reg_min)
  class_label <- as.numeric((wy/wn) >=0.5)
  rec <- border
  
  node_que <- deque()
  #list()
  rec_que <- deque()
  rec_que <- rec_que$push(border) #list(border) #check plz
  
  label_que <- deque() 
  label_que <- label_que$push(1)
    
  reclst_leg <- list()
  labellst_leg <- list()
  
  feats_usage <- vector(length=d) # false false
  n_operate_nodes = 1
  
  while(length(node_que) > 0 && n_operate_nodes < maximal_leaves)
  {
    n_operate_nodes <- n_operate_nodes + 1
    
    node <- node_que$popleft() # look
    rec <- rec_que$popleft()
    reclst <- rec_que$as_list()
    reclst <- c(reclst, unlist(reclst_leg))
    label <- label_que$popleft()
    
    labellst <- label_que$as_list()
    labellst <- c(labellst, unlist(labellst_leg))
    
    ans <- surface_funs(rec, label, reclst, labellst)
    s_0 <- ans[[d+1]][1]
    s_1 <- ans[[d+1]][2]
    
    if(label == 1){s_origin <- s_1}
    else{s_origin <- s_0}
    volume0 <- volume - label * prod(rec[,2] - rec[,1])
    # print c((volume0, rec))#class_label
    
    if(volume0 < -tol)
    {
      browser()
      print(paste('Negative volume0: ',volume0))
      if (volume0 < 0) 
      {
        stop(paste("Negative volume0: ", volume0))
      }
      
    }
    surface0 <- surface
    tree_impu0 <- tree_impu - impu * wn
    tree_sign_impu0 <- tree_sign_impu - sign_impu * wn/wn_all
    
    featureid <- -1 # featureid = -1 means no better partition is found
    feats_reorder <- c(which(feats_usage != 0), which(1-feats_usage != 0))
    node_impu_selected <- impu
    S_faces <- numeric(d)
    
    for(j in feats_reorder)
    {
      checkpoints <- ans[[j]][1]
      slope10 <- ans[[j]][2]
      intercept10 <- ans[[j]][3]
      S_faces[j] <- ans[[j]][4]
      
      loc = 0          ## loc is the largest index of checkpoints that are no greater than thre
      wleft = 0
      wyleft = 0
      dat <- data.frame(feature = X[,j], label=Y)
      dat <- dat[order(dat$feature),]
      
      for (sa in 1:(length(Y)-1))
      {
        wyleft <- wyleft + weight*dat[sa,]$label
        wleft <- wleft + 1 + (weight-1)*dat[sa,]$label
        
        tryCatch({
          if (wleft < min_leaf_weight || dat[sa,]$feature - rec[j,1] < tol) { } # pass
          else if (wn - wleft < min_leaf_weight || rec[j,2]-dat[sa+1,]$feature < tol) {}
        }, 
        error = function(e) 
        {
          browser()
          cat(dat[sa,]$feature, dat[sa+1,]$feature, rec[j,1], rec[j,2])
        })
        
        if(wleft < min_leaf_weight || dat[sa][1]-rec[j,1] < tol)
        {next}
        else if(wn - wleft < min_leaf_weight || rec[j,2]-dat[sa+1][1]< tol)
        {break}
        
        if(dat[sa+1][1] != dat[sa][1])
        {
          
          thre_new <- (dat[sa+1][1]+dat[sa][1]) / 2
          node_impu_new <- Compute_NodeImpu(wyleft, wleft, wy, wn)
          while (loc < len(checkpoints)-1 && checkpoints[loc+1] <= thre_new)
          {loc <- loc + 1}
          
          tree_impu_new <- node_impu_new * wn / wn_all + tree_impu0
          
          
            
          tree_sign_impu_new_lst <- rep(tree_sign_impu0, times = 4)
          
          surface_new_lst <- vector(length =  4)
          volume_new_lst <- vector(length = 4)
          risk_new_lst <- vector(length = 4)
          
          child_labels_lst <- vector("list", 4)
          child_labels_lst[[1]] <- c(1,1)
          child_labels_lst[[2]] <- c(0,0)
          child_labels_lst[[3]] <- c(0,1)
          child_labels_lst[[1]] <- c(1,0)
          
          # If both child nodes are labeled 1
          surface_new_lst[1] <- surface0 + s_1 - s_origin
          volume_new_lst[1] <- prod(rec[,2] - rec[,1]) + volume0
          tree_sign_impu_new_lst[1] <- tree_sign_impu_new_lst[1] + wn / wn_all * Compute_SignNodeImpu(wyleft, wleft, wy, wn, c(1,1))
          if (volume_new_lst[1] <= 0 || surface_new_lst[1] <= 0)
          {svr <- sv_reg_min}
          else
          {
            svr <- sv_regular(surface_new_lst[1], volume_new_lst[1], d)
          }
          risk_new_lst[1] <- tree_sign_impu_new_lst[1] + pen*svr
          
          # If both child nodes are labeled 0
          surface_new_lst[2] = surface0 + s_0 - s_origin
          volume_new_lst[2] = volume0
          tree_sign_impu_new_lst[2] = tree_sign_impu_new_lst[2] + wn / wn_all * Compute_SignNodeImpu(wyleft, wleft, wy, wn, c(0,0))
          
          if (volume_new_lst[2] <= 0 || surface_new_lst[2] <= 0)
          {svr <- sv_reg_min}
          else
          {svr <- sv_regular(surface_new_lst[2], volume_new_lst[2], d)}
          risk_new_lst[2] = tree_sign_impu_new_lst[2] + pen*svr
          
          # If left child is labeled 0 and right child is labeled 1'''
          surface_new_lst[3] <- surface0 + s_0 + s_1 + 2*S_faces[j] - (intercept10[loc] + slope10[loc]*(thre_new-checkpoints[loc])) - s_origin
          volume_new_lst[3] <- volume0 + prod((rec[-j,2]) - (rec[-j,1])) * (rec[j,2]-thre_new)
          tree_sign_impu_new_lst[[3]] <- tree_sign_impu_new_lst[[3]] + wn / wn_all * Compute_SignNodeImpu(wyleft, wleft, wy, wn, c(0,1))
          if (volume_new_lst[2] <= 0 || surface_new_lst[2] <= 0)
          {svr = sv_reg_min}
          else
          {svr <- sv_regular(surface_new_lst[3], volume_new_lst[3], d)}
          risk_new_lst[3] <- tree_sign_impu_new_lst[3] + pen*svr 
          
          # If left child is labeled 1 and right child is labeled 0
          surface_new_lst[4] = surface0 + intercept10[loc] + slope10[loc]*(thre_new-checkpoints[loc]) - s_origin
          volume_new_lst[4] = volume0 + prod((rec[-j,2])-np.delete(rec[-j,1])) * (thre_new-rec[j,1])
          tree_sign_impu_new_lst[4] = tree_sign_impu_new_lst[4] + wn / wn_all * Compute_SignNodeImpu(wyleft, wleft, wy, wn, c(1,0))
          if (volume_new_lst[4] <= 0 || surface_new_lst[4] <= 0)
          {svr = sv_reg_min}
          else
          {svr = sv_regular(surface_new_lst[4], volume_new_lst[4], d)}
          risk_new_lst[4] = tree_sign_impu_new_lst[4] + pen*svr
          
          argmin = which.min(risk_new_lst)
          
          if (min(surface_new_lst) < -tol)
          {
            print(paste('reclst:', reclst))
            print(paste('rec:', rec, 'length(reslst):', length(reclst)))                            
            print(paste('slope10:', slope10, 'intercept10:', intercept10))
            print(paste('Negative surface: ',(min(surface_new_lst)),'  pen: ',(pen),'  type: ',(which.min(surface_new_lst))))
            print(paste('volume0:', volume0, 'surface0:', surface0))
            print(paste('featureid_now:', j, 'thre_now:', thre_new))
            
            browser()
            print(paste('Negative surface: ',(min(surface_new_lst)),'  pen:',(pen)))
          }
          
          if (min(tree_sign_impu_new_lst) < - tol)
          {print(paste('Negative tree signed impurity: ',(min(tree_sign_impu_new_lst))))}
          
          if (risk_new_lst[argmin] < risk)
          {
            thre <- thre_new
            featureid <- j
            child_labels <- child_labels_lst[[argmin]]
            surface <- surface_new_lst[argmin]
            volume <- volume_new_lst[argmin]
            tree_impu <- tree_impu_new   
            tree_sign_impu <- tree_sign_impu_new_lst[argmin]
            risk <- risk_new_lst[argmin]
            
            if (risk < -tol)
            {
              print(paste('Negative risk: ',(risk),'  pen:',(pen)))
              print(paste('signed impu: ',(tree_sign_impu)))
              print(paste('volume: ',(volume)))
              print(paste('surface: ',(surface)))
              browser()
              print(paste('Negative risk: ',(risk),'  pen:'(pen)))
            }
          }
          
        }
        
      }
    }
    if (featureid >= 0)
    {
      leaf = False
      feats_usage[featureid] = True
      split = c(featureid, thre)
      
      left$standardize_para = standardize_para
      leftind = leftind <- which(X[,featureid] <= thre)
      left$X = X[leftind,]
      left$Y = Y[leftind]
      left$wn = length(left$Y) + (weight-1) * sum(left$Y)
      left$wy = weight * sum(left$Y)
      left$impu = Compute_Impu(left$wy, left$wn)
      left$class_label = child_labels[[1]]
      left$sign_impu = Compute_SignImpu(left$wy, left$wn, left$class_label)
      left$rec = rec
      left$rec[featureid,1] = thre
      
      if (left$wy == 0 || left$wy == left$wn || left$wn < min_split_weight)
      {
        left$leaf = TRUE
        if(left$class_label == 1)
        {
          reclst_leg <- c(reclst_leg, left$rec)
          labellst_leg <- c(labellst_leg, 1)
        }
      }
      else
      {
        node_que <- c(node_que, left) # in the form of dqeue
        rec_que <- c(rec_que, left$rec)
        label_que <- c(label_que, left$class_label)
      }
      
      
      right$standardize_para = standardize_para
      rightind = which(X[, featureid] > thre)
      right$X = X[rightind,]
      right$Y = Y[rightind]
      right$wn = length(right$Y) + (weight-1) * sum(right$Y)
      right$wy = weight * sum(right$Y)
      right$impu = Compute_Impu(right$wy, right$wn)
      right$class_label = child_labels[[2]]
      right$sign_impu = Compute_SignImpu(right$wy, right$wn, right$class_label)
      right$rec = rec
      right$rec[featureid,1] = thre
      
      if (right$wy == 0 || right$wy == right$wn || right$wn < min_split_weight)
      {
        right$leaf = TRUE
        if(right$class_label == 1)
        {
          reclst_leg <- c(reclst_leg, right$rec)
          labellst_leg <- c(labellst_leg, 1)
        }
      }
      else
      {
        node_que <- c(node_que, right)
        rec_que <- c(rec_que, right$rec)
        label_que <- c(label_que, right$class_label)
      }
    }
    else
    {
      if (class_label == 1)
      {
        reclst_leg <- c(reclst_leg, node$rec)
        labellst_leg <- c(labellst_leg, 1)
      }
    }
    
    
  }
  feats_usage <- feats_usage
  return()
}

tree_predict <- function( X) 
{
  d <- dim(X)[2]
  
  if (!is.null(standardize_para)) 
  {
    X <- data_standardize(X)
  }
  
  return(tree_localpredict(X))
}

localpredict <- function(X) 
{
  if (leaf) 
  {
    return(class_label*rep(class_label, nrow(X)))
  } 
  else 
  {
    Y <- rep(0, nrow(X))
    featureid <- as.integer(split[1])
    thre <- split[2]
    leftind <- which(X[, featureid] <= thre)
    Y[leftind] <- left$localpredict(X[leftind, ]) #look
    rightind <- which(X[, featureid] > thre)
    Y[rightind] <- right$localpredict(X[rightind, ])
    return(Y)
  }
}

print_tree <- function( leaf = TRUE, init = TRUE, print_weight = FALSE, print_impu = FALSE, class_label = NULL) 
{
  if(init) 
  {
    codename <- "root"
  }
  if (leaf) 
  {
    cat(codename, ":", class_label, "\n")
    if (print_weight) 
    {
      cat("class 1 weight, total weight:", wy, wn, "\n")
    }
    if (print_impu) 
    {
      cat("impurity:", impu, "\n")
    }
  } 
  else 
  {
    cat(codename, ": feature", split[1], "<=", split[2], "\n")
    if (print_weight) 
    {
      cat("class 1 weight and total weight:", wy, wn, "\n")
    }
    if (print_impu) 
    {
      cat("impurity, impurity_decr:", impu, impu_decr, tot_impudecr, alpha, "\n")
    }
    left$codename <- paste(codename, "left", sep = ".") #look
    print_tree(left, FALSE, print_weight, print_impu)
    right$codename <- paste(codename, "right", sep = ".")
    print_tree(right, FALSE, print_weight, print_impu)
  }
}


da0 <- read.csv("Indian Liver Patient Dataset (ILPD).csv",header = FALSE)
da0[2] <- as.numeric(da0[2]=="Male")
da0 <- as.data.frame(da0)

di <- dim(da0)
n <- di[1]

Xall <- da0
n1_ind <- which(da0[11]==2)

Yall <- rep(0,di[1])
Yall[n1_ind] <- 1
times <- floor((di[1]-length(n1_ind))/length(n1_ind)) - 1

train_test_split <- function(Xall, Yall, test_size=0.2)
{
  set.seed(221258)
  si <- floor(test_size*di[1]) + 1
  rs_test <- sample(1:di[1], si)
  rs_train <- (1:di[1])[-rs_test]
  
  X_train <- Xall[rs_train,-di[2]] 
  X_test <- Xall[rs_test,-di[2]] 
  y_train <- Yall[rs_train]
  y_test <-  Yall[rs_test]
  
  return(list(X_train, X_test, y_train, y_test))
}

X_train <- na.omit(train_test_split(Xall,Yall,test_size = 0.2)[[1]])
X_test <- na.omit(train_test_split(Xall,Yall,test_size = 0.2)[[2]])
y_train <- na.omit(train_test_split(Xall,Yall,test_size = 0.2)[[3]])
y_test <- na.omit(train_test_split(Xall,Yall,test_size = 0.2)[[4]])


n1 <- sum(Yall==1)
n0 <- di[1]-n1
train_ratio <- 0.8

n1train <- floor(n1*train_ratio)
n0train <- floor(n0*train_ratio)
pen_lst <- c(0, 1, 1.4, 2, 2.8, 4, 5.7, 8, 11, 16, 22, 32, 44, 64, 89, 128, 179, 256, 358, 512, 716, 1024) * 10^(-3) * (n0train+n1train)^(-1/3)


bootstrap_sample <- function(X, y)
{
  n_samples <- dim(X)[1]
  idxs <- sample(1:n_samples, n_samples, replace = TRUE)
  return(list(X[idxs,], y[idxs]))
}

most_common_label <- function(y)
{
  sum0 <- sum(y==0)
  sum1 <- sum(y==1)
  
  if(sum0 > sum1){return(0)}
  else{return(1)}
}

fit <- function(n_trees=40, min_samples_split=2,
                max_depth=100, n_feats=NULL, weight=1)
{
  trees <- list()
  for(i in 1:n_trees)
  {
    X_samp = bootstrap_sample(X, y)[[1]]
    y_samp = bootstrap_sample(X, y)[[2]]
    F_lst <- numeric(length = pen_lst)
    
    for(j in 1:length(pen_lst))
    {
      TP=0
      FP=0
      fit_sv(X_samp, y_samp, pen_lst[j], weight=weight, maximal_leaves=2*sqrt(n*2/3))
      Y_pred_temp = tree_predict(X_samp)
      TP = TP + sum(Y_pred_temp[which(y_samp != 0)])
      FP = FP + sum(Y_pred_temp[which(y_samp == 0)])
    }
    if (TP > 0) 
    {
      tpr <- TP / n1train
      precision <- TP / (TP + FP)
      F_lst[j] <- 2 * tpr * precision / (tpr + precision)
    }
    para_id <- which.max(F_lst)
    trees[[i]] = fit_sv(X_samp, y_samp, pen_lst[para_id], weight=self.weight, maximal_leaves=2*np.sqrt(n*2/3))
    return(trees)
  }
}

rv_predict <- function(X)
{
  tree_preds <- matrix(0, nrow = nrow(X), ncol = n_trees)
  for (i in 1:n_trees) 
  {
    tree_preds[, i] <- trees[[i]]$predict(X)
  }
  y_pred <- numeric(tree_pred)
  for( tree_pred in tree_preds)
  {
    y_pred[tree_pred] <-  most_common_label(tree_pred)
  }
  return(y_pred)
}



X_samp = bootstrap_sample(X_train, y_train)[[1]]
y_samp = bootstrap_sample(X_train, y_train)[[2]]
fit_sv(X_samp, y_samp, pen_lst[1], weight=times +1, maximal_leaves=2*sqrt(n*2/3))


model = RandomForest(weight=times+1)
fit(X_train, y_train)
preds = .predict(X_test)

print(accuracy_score(y_test, preds))
print(roc_auc_score(y_test, preds))
print(f1_score(y_test, preds))


