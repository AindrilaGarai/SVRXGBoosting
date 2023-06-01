library(datasets)
data(iris)

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

sv_ratio <- function(x, feature, val) # x= whole data, feature= which feature , mask= same feature region
{
  # new_data <- x[,!names(iris) %in% feature]
  new_data <- x[,-feature]
  d <- dim(new_data)[2]
  
  ran <- numeric(d+1)
  ran[d+1] <- val-min(x[,feature]) # length of feature
  
  for(i in 1:d)
  { ran[i] <- max(x[,i]) - min(x[,i]) }
  
  if(ran[d+1] == 0)
  { 
    volume <- prod(ran[-(d+1)])
    
    sur <- numeric(length=d-1)
    for(j in 1:d-1)
    {
      term <- 0
      for(k in (j+1):d) 
      {
        term <- term + ran[j]*ran[k]
        sur[j] <- term
      }
    }
    surface <- 2*sum(sur)
  }
  else { volume <- prod(ran)
  
  sur <- numeric(length=d)
  for(j in 1:d)
  {
    term <- 0
    for(k in (j+1):(d+1)) 
    {
      term <- term + ran[j]*ran[k]
      sur[j] <- term
    }
  }
  surface <- 2*sum(sur)
  
  } 
  rtn <- surface/volume
  return(rtn)
}
sv_ratio(iris[1:4],1,5.5)

 a <- unique(sort(iris$Sepal.Width))
l <- length(a)
vec <- numeric(l)
for (i in 1:l) {
  vec[i] <- sv_ratio(iris[1:4],2,a[i])
}
vec # decreasing


a <- which(iris$Sepal.Length<6)
table(iris$Species[a])

impurity <- function(y,alpha,mask) # y= species, alpha=vector of weights of size how many different class of y for balanced alpha=1,1
{
  len <- length(unique(y))
  wght_prob <- (alpha * table(y)) / sum(alpha * table(y))
  
  imp <- numeric(len)
  for (i in 1:len)
  {
    index <- which(y==names(table(y)[i]))
    impu_prob <- sum(mask[index]) / table(y)[i]
    
    if(impu_prob>0.5)
    {
      imp[i] <- 1- (impu_prob^2) - ((1-impu_prob)^2)
    }
    else
    {
      imp[i] <- (impu_prob^2) + ((1-impu_prob)^2)
    }
  }
  
  rtn <- sum(wght_prob*imp) 
  return(rtn)
}

y <- iris$Species
alpha <- c(1,1,1)
mask <- iris$Sepal.Length<5.8
len <- length(unique(y))
wght_prob <- (alpha * table(y)) / sum(alpha * table(y))
imp <- numeric(len)
i=1
index <- which(y==names(table(y)[i]))
impu_prob <- sum(mask[index]) / table(y)[i]
imp[i] <- 1- (impu_prob^2) - ((1-impu_prob)^2)
imp[i]

a <- unique(sort(iris$Sepal.Width))
l <- length(a)
vec <- numeric(l)
for (i in 1:l) {
  vec[i] <- impurity(iris[,5],1,iris$Sepal.Length<a[i])
}
vec

information_gain <- function(x, y, feature, mask, val, pen,alpha) # x= whole data, y= 1 2 3 ei gulo, mask= region of a feature
{
  rtn <- impurity(y,alpha,mask) + (pen*sv_ratio(x, feature, val))
  return(rtn)
}
information_gain(x=iris[,1:4],y=iris[,5],feature=1,mask=iris$Sepal.Length<5.5,val=5.5,pen=1,alpha=c(1,1,1))



de <- sum(table(iris$Species)/sum(table(iris$Species)))
de
A <- iris$Sepal.Length<5.5
y_A <- table(iris$Species[A])/sum(table(iris$Species[A]))
y_A
prob <- y_A/de
prob
#pm <- which.max(prob)
#pm
rest <- prob[-1]
rt <- sum(rest)
rt
gi_prob <- c(prob[which.max(prob)],sum(prob[-which.max(prob)]))
gi_prob
max(gi_prob)==1
gi_imp <- 1-(sum(gi_prob^2))
gi_imp
1-gi_imp
p <- sum(A)/length(A)
p*gi_imp
p*(1-gi_imp)

impurity <- function(mask,y)
{
  p <- sum(mask)/length(mask)
  de <- sum(table(y)/sum(table(y)))
  
  if(p==0){y_mask <- rep(0,length(unique(y)))}
  else{
  y_mask <- table(y[mask])/sum(table(y[mask]))}
  prob <- y_mask/de
  
  mx <- max(prob)
  
  gini <- 1- (mx^2) - ((1-mx)^2)
  if(mx > .5){gini_imp <- gini}
  else{gini_imp <- 1-gini}
  return(gini_imp)
}
impurity(iris$Sepal.Length<5.5,iris$Species)
impurity(iris$Sepal.Length<9,iris$Species)

a <- unique(sort(iris$Sepal.Width))
vec1<- numeric(length(a))
for (i in 1:length(a)) {
  vec1[i] <- impurity(iris$Sepal.Width<a[i],iris$Species)
  
}

pen

vec1
min(vec1+vec*pen[22])
vec1+vec*pen[11]



a[1]
length(a)


> iris <- iris[1:100,]
> iris
Sepal.Length Sepal.Width Petal.Length
1            5.1         3.5          1.4
2            4.9         3.0          1.4
3            4.7         3.2          1.3
4            4.6         3.1          1.5
5            5.0         3.6          1.4
6            5.4         3.9          1.7
7            4.6         3.4          1.4
8            5.0         3.4          1.5
9            4.4         2.9          1.4
10           4.9         3.1          1.5
11           5.4         3.7          1.5
12           4.8         3.4          1.6
13           4.8         3.0          1.4
14           4.3         3.0          1.1
15           5.8         4.0          1.2
16           5.7         4.4          1.5
17           5.4         3.9          1.3
18           5.1         3.5          1.4
19           5.7         3.8          1.7
20           5.1         3.8          1.5
21           5.4         3.4          1.7
22           5.1         3.7          1.5
23           4.6         3.6          1.0
24           5.1         3.3          1.7
25           4.8         3.4          1.9
26           5.0         3.0          1.6
27           5.0         3.4          1.6
28           5.2         3.5          1.5
29           5.2         3.4          1.4
30           4.7         3.2          1.6
31           4.8         3.1          1.6
32           5.4         3.4          1.5
33           5.2         4.1          1.5
34           5.5         4.2          1.4
35           4.9         3.1          1.5
36           5.0         3.2          1.2
37           5.5         3.5          1.3
38           4.9         3.6          1.4
39           4.4         3.0          1.3
40           5.1         3.4          1.5
41           5.0         3.5          1.3
42           4.5         2.3          1.3
43           4.4         3.2          1.3
44           5.0         3.5          1.6
45           5.1         3.8          1.9
46           4.8         3.0          1.4
47           5.1         3.8          1.6
48           4.6         3.2          1.4
49           5.3         3.7          1.5
50           5.0         3.3          1.4
51           7.0         3.2          4.7
52           6.4         3.2          4.5
53           6.9         3.1          4.9
54           5.5         2.3          4.0
55           6.5         2.8          4.6
56           5.7         2.8          4.5
57           6.3         3.3          4.7
58           4.9         2.4          3.3
59           6.6         2.9          4.6
60           5.2         2.7          3.9
61           5.0         2.0          3.5
62           5.9         3.0          4.2
63           6.0         2.2          4.0
64           6.1         2.9          4.7
65           5.6         2.9          3.6
66           6.7         3.1          4.4
67           5.6         3.0          4.5
68           5.8         2.7          4.1
69           6.2         2.2          4.5
70           5.6         2.5          3.9
71           5.9         3.2          4.8
72           6.1         2.8          4.0
73           6.3         2.5          4.9
74           6.1         2.8          4.7
75           6.4         2.9          4.3
76           6.6         3.0          4.4
77           6.8         2.8          4.8
78           6.7         3.0          5.0
79           6.0         2.9          4.5
80           5.7         2.6          3.5
81           5.5         2.4          3.8
82           5.5         2.4          3.7
83           5.8         2.7          3.9
84           6.0         2.7          5.1
85           5.4         3.0          4.5
86           6.0         3.4          4.5
87           6.7         3.1          4.7
88           6.3         2.3          4.4
89           5.6         3.0          4.1
90           5.5         2.5          4.0
91           5.5         2.6          4.4
92           6.1         3.0          4.6
93           5.8         2.6          4.0
94           5.0         2.3          3.3
95           5.6         2.7          4.2
96           5.7         3.0          4.2
97           5.7         2.9          4.2
98           6.2         2.9          4.3
99           5.1         2.5          3.0
100          5.7         2.8          4.1
Petal.Width    Species
1           0.2     setosa
2           0.2     setosa
3           0.2     setosa
4           0.2     setosa
5           0.2     setosa
6           0.4     setosa
7           0.3     setosa
8           0.2     setosa
9           0.2     setosa
10          0.1     setosa
11          0.2     setosa
12          0.2     setosa
13          0.1     setosa
14          0.1     setosa
15          0.2     setosa
16          0.4     setosa
17          0.4     setosa
18          0.3     setosa
19          0.3     setosa
20          0.3     setosa
21          0.2     setosa
22          0.4     setosa
23          0.2     setosa
24          0.5     setosa
25          0.2     setosa
26          0.2     setosa
27          0.4     setosa
28          0.2     setosa
29          0.2     setosa
30          0.2     setosa
31          0.2     setosa
32          0.4     setosa
33          0.1     setosa
34          0.2     setosa
35          0.2     setosa
36          0.2     setosa
37          0.2     setosa
38          0.1     setosa
39          0.2     setosa
40          0.2     setosa
41          0.3     setosa
42          0.3     setosa
43          0.2     setosa
44          0.6     setosa
45          0.4     setosa
46          0.3     setosa
47          0.2     setosa
48          0.2     setosa
49          0.2     setosa
50          0.2     setosa
51          1.4 versicolor
52          1.5 versicolor
53          1.5 versicolor
54          1.3 versicolor
55          1.5 versicolor
56          1.3 versicolor
57          1.6 versicolor
58          1.0 versicolor
59          1.3 versicolor
60          1.4 versicolor
61          1.0 versicolor
62          1.5 versicolor
63          1.0 versicolor
64          1.4 versicolor
65          1.3 versicolor
66          1.4 versicolor
67          1.5 versicolor
68          1.0 versicolor
69          1.5 versicolor
70          1.1 versicolor
71          1.8 versicolor
72          1.3 versicolor
73          1.5 versicolor
74          1.2 versicolor
75          1.3 versicolor
76          1.4 versicolor
77          1.4 versicolor
78          1.7 versicolor
79          1.5 versicolor
80          1.0 versicolor
81          1.1 versicolor
82          1.0 versicolor
83          1.2 versicolor
84          1.6 versicolor
85          1.5 versicolor
86          1.6 versicolor
87          1.5 versicolor
88          1.3 versicolor
89          1.3 versicolor
90          1.3 versicolor
91          1.2 versicolor
92          1.4 versicolor
93          1.2 versicolor
94          1.0 versicolor
95          1.3 versicolor
96          1.2 versicolor
97          1.3 versicolor
98          1.3 versicolor
99          1.1 versicolor
100         1.3 versicolor
> de <- sum(table(iris$Species)/sum(table(iris$Species)))
> de
[1] 1
> A <- iris$Sepal.Length<4.7
> y_A <- table(iris$Species[A])/sum(table(iris$Species[A]))
> y_A

setosa versicolor  virginica 
1          0          0 
> table(iris$Species)

setosa versicolor  virginica 
50         50          0 
> A <- iris$Sepal.Length<5.5
> y_A <- table(iris$Species[A])/sum(table(iris$Species[A]))
> y_A

setosa versicolor  virginica 
0.8823529  0.1176471  0.0000000 
> prob <- y_A/de
> prob

setosa versicolor  virginica 
0.8823529  0.1176471  0.0000000 
> A <- iris$Sepal.Length<5.5
> A
[1]  TRUE  TRUE  TRUE  TRUE
[5]  TRUE  TRUE  TRUE  TRUE
[9]  TRUE  TRUE  TRUE  TRUE
[13]  TRUE  TRUE FALSE FALSE
[17]  TRUE  TRUE FALSE  TRUE
[21]  TRUE  TRUE  TRUE  TRUE
[25]  TRUE  TRUE  TRUE  TRUE
[29]  TRUE  TRUE  TRUE  TRUE
[33]  TRUE FALSE  TRUE  TRUE
[37] FALSE  TRUE  TRUE  TRUE
[41]  TRUE  TRUE  TRUE  TRUE
[45]  TRUE  TRUE  TRUE  TRUE
[49]  TRUE  TRUE FALSE FALSE
[53] FALSE FALSE FALSE FALSE
[57] FALSE  TRUE FALSE  TRUE
[61]  TRUE FALSE FALSE FALSE
[65] FALSE FALSE FALSE FALSE
[69] FALSE FALSE FALSE FALSE
[73] FALSE FALSE FALSE FALSE
[77] FALSE FALSE FALSE FALSE
[81] FALSE FALSE FALSE FALSE
[85]  TRUE FALSE FALSE FALSE
[89] FALSE FALSE FALSE FALSE
[93] FALSE  TRUE FALSE FALSE
[97] FALSE FALSE  TRUE FALSE
> sum(A)
[1] 51
> table(A)
A
FALSE  TRUE 
49    51 
> length(A)
[1] 100
> A <- iris$Sepal.Length<5.5
> y_A <- table(iris$Species[A])/sum(table(iris$Species[A]))
> y_A

setosa versicolor  virginica 
0.8823529  0.1176471  0.0000000 
> max(y_A)=="setosa"
[1] FALSE
> max(y_A)==0.8823529
[1] FALSE
> max(y_A)
[1] 0.8823529
> y_A <- table(iris$Species[A])/sum(table(iris$Species[A]))
> y_A

setosa versicolor  virginica 
0.8823529  0.1176471  0.0000000 
> max(y_A)
[1] 0.8823529
> y_A[1]
setosa 
0.8823529 
> y_A[1]>.5
setosa 
TRUE 
> impurity <- function(mask,y)
  + {
    +   p <- sum(mask)/length(mask)
    +   z <- 1
    +   de <- sum(table(y)/sum(table(y)))
    +   
      +   y_mask <- table(y[mask])/sum(table(y[mask])) 
      +   prob <- y_mask/de
      +   gini <- 1- (sum(prob^2))
      +   if(prob[1]>.5){gini_imp <- gini}
      +   else{gini_imp <- 1-gini}
      +   return(gini_imp)
      + }
> iris
Sepal.Length Sepal.Width
1            5.1         3.5
2            4.9         3.0
3            4.7         3.2
4            4.6         3.1
5            5.0         3.6
6            5.4         3.9
7            4.6         3.4
8            5.0         3.4
9            4.4         2.9
10           4.9         3.1
11           5.4         3.7
12           4.8         3.4
13           4.8         3.0
14           4.3         3.0
15           5.8         4.0
16           5.7         4.4
17           5.4         3.9
18           5.1         3.5
19           5.7         3.8
20           5.1         3.8
21           5.4         3.4
22           5.1         3.7
23           4.6         3.6
24           5.1         3.3
25           4.8         3.4
26           5.0         3.0
27           5.0         3.4
28           5.2         3.5
29           5.2         3.4
30           4.7         3.2
31           4.8         3.1
32           5.4         3.4
33           5.2         4.1
34           5.5         4.2
35           4.9         3.1
36           5.0         3.2
37           5.5         3.5
38           4.9         3.6
39           4.4         3.0
40           5.1         3.4
41           5.0         3.5
42           4.5         2.3
43           4.4         3.2
44           5.0         3.5
45           5.1         3.8
46           4.8         3.0
47           5.1         3.8
48           4.6         3.2
49           5.3         3.7
50           5.0         3.3
51           7.0         3.2
52           6.4         3.2
53           6.9         3.1
54           5.5         2.3
55           6.5         2.8
56           5.7         2.8
57           6.3         3.3
58           4.9         2.4
59           6.6         2.9
60           5.2         2.7
61           5.0         2.0
62           5.9         3.0
63           6.0         2.2
64           6.1         2.9
65           5.6         2.9
66           6.7         3.1
67           5.6         3.0
68           5.8         2.7
69           6.2         2.2
70           5.6         2.5
71           5.9         3.2
72           6.1         2.8
73           6.3         2.5
74           6.1         2.8
75           6.4         2.9
76           6.6         3.0
77           6.8         2.8
78           6.7         3.0
79           6.0         2.9
80           5.7         2.6
81           5.5         2.4
82           5.5         2.4
83           5.8         2.7
84           6.0         2.7
85           5.4         3.0
86           6.0         3.4
87           6.7         3.1
88           6.3         2.3
89           5.6         3.0
90           5.5         2.5
91           5.5         2.6
92           6.1         3.0
93           5.8         2.6
94           5.0         2.3
95           5.6         2.7
96           5.7         3.0
97           5.7         2.9
98           6.2         2.9
99           5.1         2.5
100          5.7         2.8
Petal.Length Petal.Width
1            1.4         0.2
2            1.4         0.2
3            1.3         0.2
4            1.5         0.2
5            1.4         0.2
6            1.7         0.4
7            1.4         0.3
8            1.5         0.2
9            1.4         0.2
10           1.5         0.1
11           1.5         0.2
12           1.6         0.2
13           1.4         0.1
14           1.1         0.1
15           1.2         0.2
16           1.5         0.4
17           1.3         0.4
18           1.4         0.3
19           1.7         0.3
20           1.5         0.3
21           1.7         0.2
22           1.5         0.4
23           1.0         0.2
24           1.7         0.5
25           1.9         0.2
26           1.6         0.2
27           1.6         0.4
28           1.5         0.2
29           1.4         0.2
30           1.6         0.2
31           1.6         0.2
32           1.5         0.4
33           1.5         0.1
34           1.4         0.2
35           1.5         0.2
36           1.2         0.2
37           1.3         0.2
38           1.4         0.1
39           1.3         0.2
40           1.5         0.2
41           1.3         0.3
42           1.3         0.3
43           1.3         0.2
44           1.6         0.6
45           1.9         0.4
46           1.4         0.3
47           1.6         0.2
48           1.4         0.2
49           1.5         0.2
50           1.4         0.2
51           4.7         1.4
52           4.5         1.5
53           4.9         1.5
54           4.0         1.3
55           4.6         1.5
56           4.5         1.3
57           4.7         1.6
58           3.3         1.0
59           4.6         1.3
60           3.9         1.4
61           3.5         1.0
62           4.2         1.5
63           4.0         1.0
64           4.7         1.4
65           3.6         1.3
66           4.4         1.4
67           4.5         1.5
68           4.1         1.0
69           4.5         1.5
70           3.9         1.1
71           4.8         1.8
72           4.0         1.3
73           4.9         1.5
74           4.7         1.2
75           4.3         1.3
76           4.4         1.4
77           4.8         1.4
78           5.0         1.7
79           4.5         1.5
80           3.5         1.0
81           3.8         1.1
82           3.7         1.0
83           3.9         1.2
84           5.1         1.6
85           4.5         1.5
86           4.5         1.6
87           4.7         1.5
88           4.4         1.3
89           4.1         1.3
90           4.0         1.3
91           4.4         1.2
92           4.6         1.4
93           4.0         1.2
94           3.3         1.0
95           4.2         1.3
96           4.2         1.2
97           4.2         1.3
98           4.3         1.3
99           3.0         1.1
100          4.1         1.3
Species
1       setosa
2       setosa
3       setosa
4       setosa
5       setosa
6       setosa
7       setosa
8       setosa
9       setosa
10      setosa
11      setosa
12      setosa
13      setosa
14      setosa
15      setosa
16      setosa
17      setosa
18      setosa
19      setosa
20      setosa
21      setosa
22      setosa
23      setosa
24      setosa
25      setosa
26      setosa
27      setosa
28      setosa
29      setosa
30      setosa
31      setosa
32      setosa
33      setosa
34      setosa
35      setosa
36      setosa
37      setosa
38      setosa
39      setosa
40      setosa
41      setosa
42      setosa
43      setosa
44      setosa
45      setosa
46      setosa
47      setosa
48      setosa
49      setosa
50      setosa
51  versicolor
52  versicolor
53  versicolor
54  versicolor
55  versicolor
56  versicolor
57  versicolor
58  versicolor
59  versicolor
60  versicolor
61  versicolor
62  versicolor
63  versicolor
64  versicolor
65  versicolor
66  versicolor
67  versicolor
68  versicolor
69  versicolor
70  versicolor
71  versicolor
72  versicolor
73  versicolor
74  versicolor
75  versicolor
76  versicolor
77  versicolor
78  versicolor
79  versicolor
80  versicolor
81  versicolor
82  versicolor
83  versicolor
84  versicolor
85  versicolor
86  versicolor
87  versicolor
88  versicolor
89  versicolor
90  versicolor
91  versicolor
92  versicolor
93  versicolor
94  versicolor
95  versicolor
96  versicolor
97  versicolor
98  versicolor
99  versicolor
100 versicolor
> impurity(iris$Sepal.Length<5.5,iris$Species)
[1] 0.2076125
> impurity(iris$Sepal.Length<5.8,iris$Species)
[1] 0.42
> impurity(iris$Sepal.Length<4.7,iris$Species)
[1] 0
> impurity(iris$Sepal.Length<7.8,iris$Species)
[1] 0.5
> a <- unique(sort(iris$Sepal.Length))
> a <- unique(sort(iris$Sepal.Length))
> vec<- numeric(length(a))
> for (i in 1:length(a)) {
  +   vec[i] <- impurity(iris$Sepal.Length<a[i],iris$Species)
  +   
    + }
Error in if (prob[1] > 0.5) { : missing value where TRUE/FALSE needed
  > vec
  [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0
  [15] 0 0 0 0 0 0 0 0 0 0 0 0 0 0
  > y
  Error: object 'y' not found
  > y <- iris$Species
  > y
  [1] setosa     setosa    
  [3] setosa     setosa    
  [5] setosa     setosa    
  [7] setosa     setosa    
  [9] setosa     setosa    
  [11] setosa     setosa    
  [13] setosa     setosa    
  [15] setosa     setosa    
  [17] setosa     setosa    
  [19] setosa     setosa    
  [21] setosa     setosa    
  [23] setosa     setosa    
  [25] setosa     setosa    
  [27] setosa     setosa    
  [29] setosa     setosa    
  [31] setosa     setosa    
  [33] setosa     setosa    
  [35] setosa     setosa    
  [37] setosa     setosa    
  [39] setosa     setosa    
  [41] setosa     setosa    
  [43] setosa     setosa    
  [45] setosa     setosa    
  [47] setosa     setosa    
  [49] setosa     setosa    
  [51] versicolor versicolor
  [53] versicolor versicolor
  [55] versicolor versicolor
  [57] versicolor versicolor
  [59] versicolor versicolor
  [61] versicolor versicolor
  [63] versicolor versicolor
  [65] versicolor versicolor
  [67] versicolor versicolor
  [69] versicolor versicolor
  [71] versicolor versicolor
  [73] versicolor versicolor
  [75] versicolor versicolor
  [77] versicolor versicolor
  [79] versicolor versicolor
  [81] versicolor versicolor
  [83] versicolor versicolor
  [85] versicolor versicolor
  [87] versicolor versicolor
  [89] versicolor versicolor
  [91] versicolor versicolor
  [93] versicolor versicolor
  [95] versicolor versicolor
  [97] versicolor versicolor
  [99] versicolor versicolor
  3 Levels: setosa ... virginica
  > mask <- iris$Sepal.Length <5.5
  >   p <- sum(mask)/length(mask)
  >   z <- 1
  >   de <- sum(table(y)/sum(table(y)))
  >   y_mask <- table(y[mask])/sum(table(y[mask])) 
  >   prob <- y_mask/de
  >   gini <- 1- (sum(prob^2))
  > prob
  
  setosa versicolor  virginica 
  0.8823529  0.1176471  0.0000000 
  > prob[1]
  setosa 
  0.8823529 
  > prob[1]>0.5
  setosa 
  TRUE 
  >   if(prob[1]>.5){gini_imp <- gini}
  >   else{gini_imp <- 1-gini}
  Error: unexpected 'else' in "  else"
  > gini_imp
  [1] 0.2076125
  > impurity <- function(mask,y)
    + {
      +   p <- sum(mask)/length(mask)
      +   z <- 1
      +   de <- sum(table(y)/sum(table(y)))
      +   
        +   y_mask <- table(y[mask])/sum(table(y[mask])) 
        +   prob <- y_mask/de
        +   gini <- 1- (sum(prob^2))
        +   if(prob[1]>.5){gini_imp <- gini}
        +   else{gini_imp <- 1-gini}
        +   return(gini_imp)
        + }
  > impurity(iris$Sepal.Length<7.8,iris$Species)
  [1] 0.5
  > a
  [1] 4.3 4.4 4.5 4.6 4.7 4.8 4.9
  [8] 5.0 5.1 5.2 5.3 5.4 5.5 5.6
  [15] 5.7 5.8 5.9 6.0 6.1 6.2 6.3
  [22] 6.4 6.5 6.6 6.7 6.8 6.9 7.0
  > impurity(iris$Sepal.Length<6.9,iris$Species)
  [1] 0.4997918
  > a <- unique(sort(iris$Sepal.Length))
  > vec<- numeric(length(a))
  > for (i in 1:length(a)) {
    +   vec[i] <- impurity(iris$Sepal.Length<a[i],iris$Species)
    +   
      + }
  Error in if (prob[1] > 0.5) { : missing value where TRUE/FALSE needed
    > impurity(iris$Sepal.Length<a,iris$Species)
    [1] 0.4108642
    Warning message:
      In iris$Sepal.Length < a :
      longer object length is not a multiple of shorter object length
    > a[1]
    [1] 4.3
    > length(a)
    [1] 28
    > impurity(iris$Sepal.Length<4.3,iris$Species)
    Error in if (prob[1] > 0.5) { : missing value where TRUE/FALSE needed
      > impurity(iris$Sepal.Length<4.3,iris$Species)
      Error in if (prob[1] > 0.5) { : missing value where TRUE/FALSE needed
        > iris$Sepal.Length<4.3
        [1] FALSE FALSE FALSE FALSE
        [5] FALSE FALSE FALSE FALSE
        [9] FALSE FALSE FALSE FALSE
        [13] FALSE FALSE FALSE FALSE
        [17] FALSE FALSE FALSE FALSE
        [21] FALSE FALSE FALSE FALSE
        [25] FALSE FALSE FALSE FALSE
        [29] FALSE FALSE FALSE FALSE
        [33] FALSE FALSE FALSE FALSE
        [37] FALSE FALSE FALSE FALSE
        [41] FALSE FALSE FALSE FALSE
        [45] FALSE FALSE FALSE FALSE
        [49] FALSE FALSE FALSE FALSE
        [53] FALSE FALSE FALSE FALSE
        [57] FALSE FALSE FALSE FALSE
        [61] FALSE FALSE FALSE FALSE
        [65] FALSE FALSE FALSE FALSE
        [69] FALSE FALSE FALSE FALSE
        [73] FALSE FALSE FALSE FALSE
        [77] FALSE FALSE FALSE FALSE
        [81] FALSE FALSE FALSE FALSE
        [85] FALSE FALSE FALSE FALSE
        [89] FALSE FALSE FALSE FALSE
        [93] FALSE FALSE FALSE FALSE
        [97] FALSE FALSE FALSE FALSE
        > mask <- iris$Sepal.Length<4.3
        >   p <- sum(mask)/length(mask)
        > p
        [1] 0
        >   de <- sum(table(y)/sum(table(y)))
        >   y_mask <- table(y[mask])/sum(table(y[mask])) 
        > y_mask
        
        setosa versicolor  virginica 
        
        >   prob <- y_mask/de
        > prob
        
        setosa versicolor  virginica 
        
        >   gini <- 1- (sum(prob^2))
        > gini
        [1] NaN
        > a <- unique(sort(iris$Sepal.Length))
        > vec<- numeric(length(a))
        > for (i in 2:length(a)) {
          +   vec[i] <- impurity(iris$Sepal.Length<a[i],iris$Species)
          +   
            + }
        > vec
        [1] 0.00000000 0.00000000
        [3] 0.00000000 0.00000000
        [5] 0.00000000 0.00000000
        [7] 0.00000000 0.09070295
        [9] 0.17481790 0.18000000
        [11] 0.20144628 0.19753086
        [13] 0.20761246 0.30737218
        [15] 0.37893676 0.42000000
        [17] 0.43827611 0.45013850
        [19] 0.46875000 0.48185941
        [21] 0.48674959 0.49236208
        [23] 0.49510929 0.49621928
        [25] 0.49796288 0.49952173
        [27] 0.49979175 0.49994898
        > sv_ratio <- function(x, feature, val) # x= whole data, feature= which feature , mask= same feature region
          + {
            +   # new_data <- x[,!names(iris) %in% feature]
              +   new_data <- x[,-feature]
              +   d <- dim(new_data)[2]
              +   
                +   ran <- numeric(d+1)
                +   ran[d+1] <- val-min(x[,feature]) # length of feature
                +   
                  +   for(i in 1:d)
                    +   { ran[i] <- max(x[,i]) - min(x[,i]) }
                +   
                  +   if(ran[d+1] == 0)
                    +   { 
                      +     volume <- prod(ran[-(d+1)])
                      +     
                        +     sur <- numeric(length=d-1)
                        +     for(j in 1:d-1)
                          +     {
                            +       term <- 0
                            +       for(k in (j+1):d) 
                              +       {
                                +         term <- term + ran[j]*ran[k]
                                +         sur[j] <- term
                                +       }
                            +     }
                        +     surface <- 2*sum(sur)
                        +   }
                +   else { volume <- prod(ran)
                +   
                  +   sur <- numeric(length=d)
                  +   for(j in 1:d)
                    +   {
                      +     term <- 0
                      +     for(k in (j+1):(d+1)) 
                        +     {
                          +       term <- term + ran[j]*ran[k]
                          +       sur[j] <- term
                          +     }
                      +   }
                  +   surface <- 2*sum(sur)
                  +   
                    +   } 
                +   rtn <- surface/volume
                +   return(rtn)
                + }
        > a <- unique(sort(iris$Sepal.Length))
        > l <- length(a)
        > vec1 <- numeric(l)
        > for (i in 2:l) {
          +     vec1[i] <- sv_ratio(iris[1:4],1,a[i])
          + }
        > vec
        [1] 0.00000000 0.00000000 0.00000000
        [4] 0.00000000 0.00000000 0.00000000
        [7] 0.00000000 0.09070295 0.17481790
        [10] 0.18000000 0.20144628 0.19753086
        [13] 0.20761246 0.30737218 0.37893676
        [16] 0.42000000 0.43827611 0.45013850
        [19] 0.46875000 0.48185941 0.48674959
        [22] 0.49236208 0.49510929 0.49621928
        [25] 0.49796288 0.49952173 0.49979175
        [28] 0.49994898
        > vec1
        [1]  0.000000 21.311352 11.001957
        [4]  7.565492  5.847260  4.816320
        [7]  4.129027  3.638104  3.269911
        [10]  2.983539  2.754441  2.566998
        [13]  2.410795  2.278623  2.165333
        [16]  2.067148  1.981237  1.905432
        [19]  1.838051  1.777762  1.723502
        [22]  1.674410  1.629780  1.589032
        [25]  1.551679  1.517314  1.485593
        [28]  1.456221
        > vec+vec1
        [1]  0.000000 21.311352 11.001957
        [4]  7.565492  5.847260  4.816320
        [7]  4.129027  3.728807  3.444729
        [10]  3.163539  2.955888  2.764529
        [13]  2.618407  2.585995  2.544270
        [16]  2.487148  2.419513  2.355571
        [19]  2.306801  2.259621  2.210252
        [22]  2.166772  2.124889  2.085251
        [25]  2.049642  2.016836  1.985385
        [28]  1.956170
        > vec+(vec1*10^(-3))
        [1] 0.000000000 0.021311352 0.011001957
        [4] 0.007565492 0.005847260 0.004816320
        [7] 0.004129027 0.094341052 0.178087809
        [10] 0.182983539 0.204200722 0.200097862
        [13] 0.210023252 0.309650799 0.381102093
        [16] 0.422067148 0.440257351 0.452043937
        [19] 0.470588051 0.483637172 0.488473096
        [22] 0.494036485 0.496739067 0.497808313
        [25] 0.499514558 0.501039049 0.501277346
        [28] 0.501405206
        > vec+(vec1*10^(-3)*256)
        [1] 0.0000000 5.4557061 2.8165011
        [4] 1.9367660 1.4968985 1.2329780
        [7] 1.0570310 1.0220575 1.0119152
        [10] 0.9437860 0.9065833 0.8546823
        [13] 0.8247760 0.8906997 0.9332620
        [16] 0.9491900 0.9454727 0.9379292
        [19] 0.9392910 0.9369665 0.9279661
        [22] 0.9210109 0.9123330 0.9030114
        [25] 0.8951926 0.8879541 0.8801035
        [28] 0.8727417
        > vec+(vec1*10^(-3)*128)
        [1] 0.0000000 2.7278531 1.4082505
        [4] 0.9683830 0.7484493 0.6164890
        [7] 0.5285155 0.5563802 0.5933665
        [10] 0.5618930 0.5540148 0.5261066
        [13] 0.5161942 0.5990359 0.6560994
        [16] 0.6845950 0.6918744 0.6940339
        [19] 0.7040205 0.7094129 0.7073578
        [22] 0.7066865 0.7037212 0.6996153
        [25] 0.6965778 0.6937379 0.6899476
        [28] 0.6863453
        > min(vec+(vec1*10^(-3)*128))
        [1] 0
        > sort(vec+(vec1*10^(-3)*128))
        [1] 0.0000000 0.5161942 0.5261066
        [4] 0.5285155 0.5540148 0.5563802
        [7] 0.5618930 0.5933665 0.5990359
        [10] 0.6164890 0.6560994 0.6845950
        [13] 0.6863453 0.6899476 0.6918744
        [16] 0.6937379 0.6940339 0.6965778
        [19] 0.6996153 0.7037212 0.7040205
        [22] 0.7066865 0.7073578 0.7094129
        [25] 0.7484493 0.9683830 1.4082505
        [28] 2.7278531
        > a[14]
        [1] 5.6
