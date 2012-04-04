##################################################" nouvelle classe


validPairedObject<-function(object){
test1<-(dim(object)[2]==2)
if(test1) {test<-TRUE} else 
{
return("paired data must have two columns")
}
test2<-class(object[,1])==class(object[,2])
if(test2) {test<-TRUE} else 
{
return("paired data must have the same classes")
}
if(is.factor(object[,1])) {test3<-all(levels(object[,1])==levels(object[,2]))
if(test3) {test<-TRUE} else 
{
return("paired factor data must have the same levels")
}
}
return(test)
}


### Les constructeurs

setClass(Class="paired",contains="data.frame",validity=validPairedObject)

paired<-function(x,y){
object<-new(Class="paired",data.frame(x,y))
name.x<-deparse(substitute(x))
name.y<-deparse(substitute(y))
colnames(object)<-c(name.x,name.y)
object
}



setMethod("summary",
          signature(object = "paired"),
          function(object){
if(is.numeric(object[,1])){paired.summary(object[,1],object[,2])} else
{summary.data.frame(object)}
          }
)


setMethod(f="effect.size", 
signature(object = "paired"),
  definition=function(object){
if(!is.numeric(object[,1])){return("only suitable to numeric paired data")}
paired.effect.size(object[,1],object[,2])
}
)

setMethod("plot", signature="paired",
  function(x,groups=NULL,facet=TRUE,...){
df<-x
if(is.null(colnames(df))){colnames(df)<-c("C1","C2")}
conditions<-colnames(df)
if(!is.null(groups)){
df<-data.frame(df,groups)	
paired.plotCor(df,conditions[1],conditions[2],groups=colnames(df)[3],facet=facet,...)
}
else{
paired.plotCor(df,conditions[1],conditions[2],groups=NULL,facet=facet,...)
}
}
)

setMethod(f="plotBA", 
signature(object = "paired"),
  definition=function(object,groups=NULL,facet=TRUE,...){
if(!is.numeric(object[,1])){return("only suitable to numeric paired data")}
if(is.null(groups)){
df<-object
condition1<-colnames(df)[1]
condition2<-colnames(df)[2]
return(paired.plotBA(df, condition1, condition2, groups = NULL, facet = facet,...))
}
else{
df<-data.frame(object,groups)
condition1<-colnames(df)[1]
condition2<-colnames(df)[2]
colnames(df)[3]<-deparse(substitute(groups))
groups<-colnames(df)[3]
return(paired.plotBA(df, condition1, condition2, groups, facet = facet,...))
}
}
)

setMethod(f="plotMcNeil", 
signature(object = "paired"),
  definition=function(object,groups=NULL,subjects=NULL,facet=TRUE,...){
if(!is.numeric(object[,1])){return("only suitable to numeric paired data")}

df<-object
condition1<-colnames(df)[1]
condition2<-colnames(df)[2]

if(is.null(subjects)){
n<-dim(df)[1]
subjects<-factor(paste("S",1:n,sep=""))
df<-data.frame(df,subjects)
}
else{
df<-data.frame(df,subjects)
}

if(is.null(groups)){
return(paired.plotMcNeil(df, condition1, condition2, groups = NULL,subjects=colnames(df)[3], facet = facet,...))
}
else{
df<-data.frame(df,groups)
colnames(df)[4]<-deparse(substitute(groups))
groups<-colnames(df)[4]
return(paired.plotMcNeil(df, condition1, condition2, groups,subjects=colnames(df)[3], facet = facet,...))
}
}
)

setMethod(f="plotProfiles", 
signature(object = "paired"),
  definition=function(object,groups=NULL,subjects=NULL,facet=TRUE,...){
if(!is.numeric(object[,1])){return("only suitable to numeric paired data")}

df<-object
condition1<-colnames(df)[1]
condition2<-colnames(df)[2]

if(is.null(subjects)){
n<-dim(df)[1]
subjects<-factor(paste("S",1:n,sep=""))
df<-data.frame(df,subjects)
}
else{
df<-data.frame(df,subjects)
}

if(is.null(groups)){
return(paired.plotProfiles(df, condition1, condition2, groups = NULL,subjects=colnames(df)[3], facet = facet,...))
}
else{
df<-data.frame(df,groups)
colnames(df)[4]<-deparse(substitute(groups))
groups<-colnames(df)[4]
return(paired.plotProfiles(df, condition1, condition2, groups,subjects=colnames(df)[3], facet = facet,...))
}
}
)



setMethod(f="plotSliding", 
signature="paired",
  definition=function(object,...){
if(!is.numeric(object[,1])){return("only suitable to numeric paired data")}

df<-object
condition1<-colnames(df)[1]
condition2<-colnames(df)[2]
return(paired.plotSliding(df, condition1, condition2,...))

}
)


#### Enfin les fonctions de simulation :

rpaired.contaminated<-
function (n, d1 = c(0.1, 10, 1), d2 = c(0.1, 10, 1), r = 0.5) 
{
    require(mvtnorm)
    eps1 <- d1[1]
    k1 <- d1[2]
    Sigma1 <- d1[3]
    eps2 <- d2[1]
    k2 <- d2[2]
    Sigma2 <- d2[3]
    X <- rmvnorm(n, mean = c(0, 0), sigma = matrix(c(1, r, r, 
        1), ncol = 2))
    u1 <- pnorm(X[, 1])
    b1 <- rbinom(n, size = 1, prob = eps1)
    SD1 <- Sigma1 * (b1 * k1 + (1 - b1) * 1)
    x <- qnorm(u1, mean = 0, sd = SD1)
    u2 <- pnorm(X[, 2])
    b2 <- rbinom(n, size = 1, prob = eps2)
    SD2 <- Sigma2 * (b2 * k2 + (1 - b2) * 1)
    y <- qnorm(u2, mean = 0, sd = SD2)
    return(paired(x,y))
}

rpaired.gld<-
function (n, d1=c(0.000,0.1974,0.1349,0.1349), d2=c(0.000,0.1974,0.1349,0.1349), r) 
{
    require(mvtnorm)
    require(gld)
    X <- rmvnorm(n, mean = c(0, 0), sigma = matrix(c(1, r, r, 
        1), ncol = 2))
    u1 <- pnorm(X[, 1])
    u2 <- pnorm(X[, 2])
    x <- qgl(u1, lambda1 = d1[1], lambda2 = d1[2], lambda3 = d1[3], 
        lambda4 = d1[4], param = "rs")
    y <- qgl(u2, lambda1 = d2[1], lambda2 = d2[2], lambda3 = d2[3], 
        lambda4 = d2[4], param = "rs")
    return(paired(x,y))
}

