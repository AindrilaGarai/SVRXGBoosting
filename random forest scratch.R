library(datasets)
data(iris)
data(ToothGrowth)

entropy <- function(y){
  # assumes y is a factor with all levels
  if(length(y)==0) return(0)
  p <- table(y)/length(y)
  sum(-p*log2(p+1e-9))
}

gini_impurity <- function(y){
  # assumes y if a factor with all levels
  if(length(y) == 0) return(0)
  p <- table(y)/length(y)
  1-sum(p^2)
}

variance <- function(y){
  if(length(y) <= 1) return(0)
  var(y)
}


entropy(iris[,'Species'])
gini_impurity(iris[,'Species'])
variance(ToothGrowth[,'len'])

information_gain <- function(y,mask,func=entropy){
  s1 = sum(mask)
  s2 = length(mask)-s1
  if ( s1 == 0 | s2 == 0) return(0)
  func(y)-s1/(s1+s2)*func(y[mask])-s2/(s1+s2)*func(y[!mask])
}

information_gain(iris[,'Species'],iris[,'Petal.Width'] < 1.5)
information_gain(iris[,'Species'],iris[,'Petal.Width'] < 1.5,gini_impurity)
information_gain(ToothGrowth[,'len'],ToothGrowth[,'dose']<= 0.5,variance)

max_information_gain_split <- function(y,x,func=gini_impurity){
  best_change = NA
  split_value = NA
  is_numeric = !(is.factor(x)|is.logical(x)|is.character(x))
  for( val in sort(unique(x))){
    mask <- x == val
    if (is_numeric) mask <- x < val
    change <- information_gain(y,mask,func) 
    if(is.na(best_change) | change > best_change){
      best_change = change
      split_value = val
    }
  }
  return(list("best_change"=best_change,
              "split_value"=split_value,
              "is_numeric"=is_numeric))
}

print(unlist(max_information_gain_split(iris$Species,iris$Sepal.Length)))

print(unlist(max_information_gain_split(iris$Species, iris$Petal.Width, entropy)))

print(unlist(max_information_gain_split(ToothGrowth$len,ToothGrowth$dose,variance)))

sapply(iris[,1:4],function(x) max_information_gain_split(iris[,5],x))

sapply(ToothGrowth[,2:3],function(x) max_information_gain_split(ToothGrowth[,1],x))

best_feature_split <- function(X,y){
  results <- sapply(X,function(x) max_information_gain_split(y,x))
  best_name <- names(which.max(results['best_change',]))
  best_result <- results[,best_name]
  best_result[["name"]] <- best_name
  best_result
} 

as.data.frame(best_feature_split(iris[,1:4],iris[,5]))
as.data.frame(best_feature_split(ToothGrowth[,2:3],ToothGrowth[,1]))


get_best_mask <- function(X,best_feature_list){
  best_mask <- X[,best_feature_list$name] == best_feature_list$split_value
  
  if(best_feature_list$is_numeric){
    best_mask <- X[,best_feature_list$name] < best_feature_list$split_value
  }
  return(best_mask)
}

best_first_split <- best_feature_split(ToothGrowth[,2:3],ToothGrowth[,1])
best_tooth_mask <- get_best_mask(ToothGrowth[,2:3],best_first_split)
leftDf = ToothGrowth[best_tooth_mask,]
rightDf = ToothGrowth[!best_tooth_mask,]

print.DecisionTreeNode <- function(node,level=0){
  response <- paste("|->",node$split_description)
  if(level < node$max_depth){
    if(!is.null(node$branches$left)){
      
      response <- paste0(response,"\n",paste(rep(" ",2*(level+1)),collapse=" "),print(node$branches$left,level+1))
      
    }
    if(!is.null(node$branches$right)){
      
      response <- paste0(response,"\n",paste(rep(" ",2*(level+1)),collapse=" "),print(node$branches$right,level+1))
      
    }
  }
  
  if(level==0) {
    cat(response)
  } else {
    return(response)
  }
}

DecisionTreeNode <- setRefClass("DecisionTreeNode",
                                fields = list(x = "data.frame",
                                              y = "ANY",
                                              is_leaf="logical",
                                              split_description="character",
                                              best_split="list",
                                              branches="list",
                                              depth="numeric",
                                              minimize_func="function",
                                              min_information_gain="numeric",
                                              min_leaf_size="numeric",
                                              max_depth="numeric"),
                                methods = list(
                                  initialize = function(...){
                                    defaults <- list(x = data.frame(),
                                                     y=c(),
                                                     depth=0,
                                                     minimize_func=gini_impurity,
                                                     min_information_gain=1e-3,
                                                     min_leaf_size=20,
                                                     max_depth=3,
                                                     is_leaf=T,
                                                     split_description="root",
                                                     best_split=NULL,
                                                     branches=list("left"=NULL,"right"=NULL))
                                    params <- list(...)
                                    fields <- names(getRefClass()$fields())
                                    for( field in fields){
                                      if (!(field %in% names(params))) {
                                        params[[field]] <- defaults[[field]]
                                      }
                                    }
                                    for( param in names(params)){
                                      do.call("<<-",list(param, params[[param]]))
                                    }
                                    
                                  },
                                  information_gain = function(mask){
                                    
                                    s1 = sum(mask)
                                    s2 = length(mask)-s1
                                    if ( s1 == 0 | s2 == 0) return(0)
                                    minimize_func(y)-s1/(s1+s2)*minimize_func(y[mask])-s2/(s1+s2)*minimize_func(y[!mask])
                                  },
                                  max_information_gain_split = function(feature){
                                    
                                    best_change = NA
                                    split_value = NA
                                    is_numeric = !(is.factor(feature)|is.logical(feature)|is.character(feature))
                                    
                                    previous_val <- NA
                                    for( val in sort(unique(feature))){
                                      
                                      mask <- feature == val
                                      
                                      if (is_numeric) mask <- feature < val
                                      change <- information_gain(mask) 
                                      
                                      s1 = sum(mask)
                                      s2 = length(mask)-s1
                                      if(is.na(best_change) & s1 >= min_leaf_size & s2 >= min_leaf_size){
                                        best_change = change
                                        split_value = ifelse(is.na(previous_val),
                                                             val,
                                                             mean(c(val,previous_val)))
                                      } else if( change > best_change & s1 >= min_leaf_size & s2 >= min_leaf_size ){
                                        best_change = change
                                        split_value = ifelse(is_numeric,
                                                             mean(c(val,previous_val)),
                                                             val)
                                      }
                                      previous_val <- val
                                      
                                    }
                                    return(list("best_change"=best_change,
                                                "split_value"=split_value,
                                                "is_numeric"=is_numeric))
                                  },
                                  best_feature_split = function(){
                                    results <- sapply(x,function(feature) max_information_gain_split(feature))
                                    if (!all(is.na(unlist(results['best_change',])))) {
                                      best_name <- names(which.max(results['best_change',]))
                                      best_result <- results[,best_name]
                                      best_result[["name"]] <- best_name
                                      best_split <<- best_result
                                    }
                                    
                                  },
                                  best_mask = function(){
                                    best_mask <- x[,best_split$name] == best_split$split_value
                                    if(best_split$is_numeric){
                                      best_mask <- x[,best_split$name] < best_split$split_value
                                    }
                                    return(best_mask)
                                    
                                  },
                                  split_node = function() {
                                    if(depth < max_depth){ 
                                      best_feature_split() 
                                      if(!is.null(best_split) & best_split$best_change > min_information_gain ){
                                        
                                        mask = best_mask()
                                        if(sum(mask) >= min_leaf_size && length(mask)-sum(mask) >= min_leaf_size){
                                          is_leaf <<- F
                                          
                                          branches$left <<- .self$copy()
                                          branches$left$is_leaf <<- T
                                          branches$left$x <<-  branches$left$x[mask,]
                                          branches$left$y <<-  branches$left$y[mask]
                                          
                                          branches$left$split_description <<- ifelse(best_split$is_numeric,
                                                                                     paste(c(best_split$name,
                                                                                             "<",
                                                                                             as.numeric(as.character(best_split$split_value))),
                                                                                           collapse = " "),
                                                                                     paste(c(best_split$name,
                                                                                             "=",
                                                                                             best_split$split_value),
                                                                                           collapse = " "))
                                          
                                          branches$left$depth <<-  branches$left$depth+1
                                          branches$left$branches <<- list("left"=NULL,"right"=NULL)
                                          branches$left$split_node()
                                          
                                          branches$right <<- .self$copy()
                                          branches$right$is_leaf <<- T
                                          branches$right$x <<-  branches$right$x[!mask,]
                                          branches$right$y <<-  branches$right$y[!mask]
                                          
                                          branches$right$split_description <<- ifelse(best_split$is_numeric, 
                                                                                      paste(c(best_split$name, ">=",
                                                                                              best_split$split_value),
                                                                                            collapse = " "),
                                                                                      paste(c(best_split$name,
                                                                                              "!=",
                                                                                              best_split$split_value),
                                                                                            collapse = " "))
                                          
                                          branches$right$depth <<-  branches$right$depth+1
                                          branches$right$branches <<- list("left"=NULL,"right"=NULL)
                                          branches$right$split_node()
                                        }
                                      }
                                    }
                                    if(is_leaf){
                                      split_description <<- ifelse(identical(minimize_func,variance),
                                                                   paste(c(split_description,
                                                                           ":",
                                                                           "predict - ",
                                                                           mean(y)),
                                                                         collapse=" "),
                                                                   paste(c(split_description,
                                                                           ":",
                                                                           "predict - ",
                                                                           names(which.max(table(y)))),
                                                                         collapse=" "))
                                    }
                                    
                                  },
                                  predict_row = function(row){
                                    if(is_leaf){
                                      predict_value <- ifelse(identical(minimize_func,variance),
                                                              mean(y),
                                                              names(which.max(table(y))))
                                    } else {
                                      if(best_split$is_numeric){
                                        left = row[best_split$name] < best_split$split_value
                                      } else{
                                        left = row[best_split$name] == best_split$split_value
                                      }
                                      if(left){
                                        predict_value = branches$left$predict_row(row)
                                      } else {
                                        predict_value = branches$right$predict_row(row)
                                      }
                                    }
                                    return(predict_value)
                                  },
                                  
                                  predict = function(features){
                                    pred <- character(length=dim(features)[1])
                                    if(identical(minimize_func,variance)) pred <- numeric(length=dim(features)[1])
                                    for(i in 1:dim(features)[1]){
                                      pred[i] = predict_row(features[i,])
                                    }
                                    pred
                                  }
                                ))

dtn <- DecisionTreeNode$new(x=iris[,1:4],y=iris[,5],max_depth=3)
dtn$split_node()
print(dtn)

dtn <- DecisionTreeNode$new(x=ToothGrowth[,2:3],
                            y=ToothGrowth[,1],
                            minimize_func=variance,
                            min_leaf_size=5,
                            max_depth=5,
                            min_information_gain=1e-3)
dtn$split_node()
print(dtn)

DecisionTree <- setRefClass("DecisionTree",
                            fields = list(minimize_func="function",
                                          min_information_gain="numeric",
                                          min_leaf_size="numeric",
                                          max_depth="numeric",
                                          root = "DecisionTreeNode"),
                            methods = list(
                              initialize = function(...){
                                defaults <- list(minimize_func=gini_impurity,
                                                 min_information_gain=1e-3,
                                                 min_leaf_size=20,
                                                 max_depth=3,
                                                 root=NULL)
                                
                                params <- list(...)
                                fields <- names(getRefClass()$fields())
                                for( field in fields){
                                  if (!(field %in% names(params))) {
                                    params[[field]] <- defaults[[field]]
                                  }
                                }
                                for( param in names(params)){
                                  do.call("<<-",list(param, params[[param]]))
                                }
                                
                              },
                              fit = function(features,target){
                                root <<- DecisionTreeNode$new(x=features,
                                                              y=target,
                                                              minimize_func=minimize_func,
                                                              min_information_gain=min_information_gain,
                                                              min_leaf_size=min_leaf_size,
                                                              max_depth=max_depth
                                )
                                root$split_node()
                                
                              },
                              predict = function(features){
                                root$predict(features)
                              }
                            ))

print.DecisionTree <- function(tree){
  print(tree$root)
}

dt_bts <- DecisionTree(max_depth=5,min_leaf_size=20,min_information_gain=1e-7)
dt_bts$fit(iris[,1:4],iris[,5])
print(dt_bts)

dt_bts <- DecisionTree(max_depth=5,min_leaf_size=5,minimize_func=variance)
dt_bts$fit(ToothGrowth[,2:3],ToothGrowth[,1])
print(dt_bts)

rf <- function(no_of_trees,x_train,y_train,x_test)
{
  pred <- matrix(0,dim(x_train)[1],no_of_trees)
  for(i in 1:no_of_trees)
  {
    boot_data <- sample(1:dim(x_train)[1], dim(x_train)[1], replace = TRUE)
    x <- x_train[boot_data,]
    y <- y_train[boot_data]
    
    dt_bts <- DecisionTree(max_depth=3,min_leaf_size=20,min_information_gain=1e-7)
    dt_bts$fit(x,y)
    dt_pred <- dt_bts$predict(x_test) 
    for(j in 1:dim(x_train)[1])
    {
      pred[j,i] <- (dt_pred[j]) 
    }
  }
  rvf <- NULL
  new <- numeric(no_of_trees)
  for (k in 1:dim(x)[1])
  {
    new <- pred[k,]
    rvf <- c(rvf,which.max(table(new)))
  }
  rtn <- rvf 
  return(rtn)
  
}

sam <- sample(1:150, 100)
x_train <- iris[sam,1:4]
x_test <- iris[-sam,1:4]
y_train <- iris[sam,5]
y_test <- iris[-sam,5]

fp <- as.integer(as.factor(names(rf(100,x_train,y_train,x_test))))
oy <- as.integer(as.factor(y_test))
table(fp==oy)


tsam <- sample(1:60,40)
x_train <- ToothGrowth[tsam,2:3]
x_test <- ToothGrowth[-tsam,2:3]
y_train <- ToothGrowth[tsam,1]
y_test <- ToothGrowth[-tsam,1]


rf <- function(no_of_trees,x_train,y_train,x_test)
{
  pred <- matrix(0,dim(x_train)[1],no_of_trees)
  for(i in 1:no_of_trees)
  {
    boot_data <- sample(1:dim(x_train)[1], dim(x_train)[1], replace = TRUE)
    x <- x_train[boot_data,]
    y <- y_train[boot_data]
    
    dt_bts <- DecisionTree(max_depth=5,min_leaf_size=5,minimize_func=variance)
    dt_bts$fit(x,y)
    dt_pred <- dt_bts$predict(x_test) 
    for(j in 1:dim(x_train)[1])
    {
      pred[j,i] <- (dt_pred[j]) 
    }
  }
  rvf <- NULL
  new <- numeric(no_of_trees)
  for (k in 1:dim(x)[1])
  {
    new <- pred[k,]
    rvf <- c(rvf,which.max(table(new)))
  }
  rtn <- rvf 
  return(rtn)
  
}
fpt <- as.numeric(names(rf(100,x_train,y_train,x_test)))
oyt <- (y_test)
cov(oyt,fpt)/sqrt(var(oyt)*var(fpt))

