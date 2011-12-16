

### Fonctions de simulation
rpaired.gld<-
function(n,d1,d2,r){
require(mvtnorm)
require(gld)
X<-rmvnorm(n,mean=c(0,0),sigma=matrix(c(1,r,r,1),ncol=2))
u1<-pnorm(X[,1])
u2<-pnorm(X[,2])
x<-qgl(u1,lambda1=d1[1],lambda2=d1[2],lambda3=d1[3],lambda4=d1[4],param="rs")
y<-qgl(u2,lambda1=d2[1],lambda2=d2[2],lambda3=d2[3],lambda4=d2[4],param="rs")
return(list(x=x,y=y))
}



### Simulations appariees par des lois normales contaminees
rpaired.contaminated<-
function(n,d1=c(0.1,10,1),d2=c(0.1,10,1),r=0.5){
require(mvtnorm)
eps1<-d1[1]
k1<-d1[2]
Sigma1<-d1[3]
eps2<-d2[1]
k2<-d2[2]
Sigma2<-d2[3]
X<-rmvnorm(n,mean=c(0,0),sigma=matrix(c(1,r,r,1),ncol=2))
u1<-pnorm(X[,1])
b1<-rbinom(n,size=1,prob=eps1)
SD1<-Sigma1*(b1*k1+(1-b1)*1)
x<-qnorm(u1,mean=0,sd=SD1)
u2<-pnorm(X[,2])
b2<-rbinom(n,size=1,prob=eps2)
SD2<-Sigma2*(b2*k2+(1-b2)*1)
y<-qnorm(u2,mean=0,sd=SD2)
return(list(x=x,y=y))
}

### Les fonctions de test

yuen1.test<-
function(x,tr=.2,conf.level = 0.95,mu=0,alternative=c("two.sided", "less", "greater")){
#
# inspired from yuend (WILCOX)
#
#  Test if the trimmed mean of x can be mu
#  Compute ac CI for the trimmed mean
#  The default amount of trimming is 0.2
#
#  Missing values (values stored as NA) are not allowed.
#
#
#  This function uses winvar.

    alternative <- match.arg(alternative)

    if(!missing(mu) && (length(mu) != 1 || is.na(mu)))
        stop("'mu' must be a single number")
    if(!missing(conf.level) &&
       (length(conf.level) != 1 || !is.finite(conf.level) ||
        conf.level < 0 || conf.level > 1))
        stop("'conf.level' must be a single number between 0 and 1")
   
alpha<-1-conf.level
n<-length(x)
g<-floor(tr*n)
df<-n-2*g-1
sw<-sqrt(winvar(x,tr))
se<-sw/((1-2*tr)*sqrt(n))
dif<-mean(x,tr)
test<-(dif-mu)/se

if (alternative == "less") {
crit<-qt(1-alpha,df)
low<- -Inf
up<-dif+crit*se
pval <- pt(test, df)
    }

if (alternative == "greater") {
crit<-qt(1-alpha,df)
low<- dif-crit*se
up<- Inf
pval <- 1-pt(test, df)
    }

if (alternative == "two.sided") {
crit<-qt(1-alpha/2,df)
low<-dif-crit*se
up<-dif+crit*se
pval <- 2*(1-pt(abs(test),df))
    }

 
estimate <- dif
tstat <- test
cint <- c(low,up)
method <- paste("One sample Yuen test, trim=",tr,sep="")
names(estimate) <- c("trimmed mean of x")
dname <- deparse(substitute(x))
names(tstat) <- "t"
names(df) <- "df"
names(mu) <- c("trimmed means")   


   
    attr(cint,"conf.level") <- conf.level
    rval <- list(statistic = tstat, parameter = df, p.value = pval,
       conf.int=cint, estimate=estimate, null.value = mu,
       alternative=alternative,
       method=method, data.name=dname)
    class(rval) <- "htest"
    return(rval)
}


yuenp.test <-
function(x, y = NULL, alternative = c("two.sided", "less", "greater"),
         mu = 0, conf.level = 0.95,tr=0.2)
{
# inspired from yuend (WILCOX)
#
#  Compare the trimmed means of two dependent random variables
#  using the data in x and y.
#  The default amount of trimming is 0.2

    alternative <- match.arg(alternative)

    if(!missing(mu) && (length(mu) != 1 || is.na(mu)))
        stop("'mu' must be a single number")
    if(!missing(conf.level) &&
       (length(conf.level) != 1 || !is.finite(conf.level) ||
        conf.level < 0 || conf.level > 1))
        stop("'conf.level' must be a single number between 0 and 1")
   
 if( !is.null(y) ) {
dname <- paste(deparse(substitute(x)),"and",deparse(substitute(y)))
}
   
    else {
stop("'y' is missing for paired test")
    }
    xok <- yok <- complete.cases(x,y)
    x <- x[xok]
    y <- y[yok]


h1<-length(x)-2*floor(tr*length(x))
q1<-(length(x)-1)*winvar(x,tr)
q2<-(length(y)-1)*winvar(y,tr)
q3<-(length(x)-1)*wincor(x,y,tr)$cov
df<-h1-1
se<-sqrt((q1+q2-2*q3)/(h1*(h1-1)))
dif<-mean(x,tr)-mean(y,tr)
test<-(dif-mu)/se
alpha<-1-conf.level

if (alternative == "less") {
crit<-qt(1-alpha,df)
low<- -Inf
up<-dif+crit*se
pval <- pt(test, df)
    }

if (alternative == "greater") {
crit<-qt(1-alpha,df)
low<- dif-crit*se
up<- Inf
pval <- 1-pt(test, df)
    }

if (alternative == "two.sided") {
crit<-qt(1-alpha/2,df)
low<-dif-crit*se
up<-dif+crit*se
pval <- 2*(1-pt(abs(test),df))
    }

estimate<-c(mean(x,tr), mean(y,tr))
cint<-c(low,up)
dif<-dif
tstat<-test



    
method <- paste("Paired Yuen test, trim=",tr,sep="")
names(estimate) <- c("trimmed mean of x","trimmed mean of y")
    

       
    names(tstat) <- "t"
    names(df) <- "df"
    names(mu) <- c("difference in trimmed means")
    attr(cint,"conf.level") <- conf.level
    rval <- list(statistic = tstat, parameter = df, p.value = pval,
       conf.int=cint, estimate=estimate, null.value = mu,
       alternative=alternative,
       method=method, data.name=dname)
    class(rval) <- "htest"
    return(rval)
}

### Autres fonctions necessaires a la precedente…



winvar<-
function(x,tr=.2,na.rm=FALSE){
#
#  Compute the gamma Winsorized variance for the data in the vector x.
#  tr is the amount of Winsorization which defaults to .2.
#
if(na.rm)x<-x[!is.na(x)]
y<-sort(x)
n<-length(x)
ibot<-floor(tr*n)+1
itop<-length(x)-ibot+1
xbot<-y[ibot]
xtop<-y[itop]
y<-ifelse(y<=xbot,xbot,y)
y<-ifelse(y>=xtop,xtop,y)
winvar<-var(y)
winvar
}

wincor<-
function(x,y,tr=.2){
#   Compute the Winsorized correlation between x and y.
#
#   tr is the amount of Winsorization
#   This function also returns the Winsorized covariance
#

g<-floor(tr*length(x))
xvec<-winval(x,tr)
yvec<-winval(y,tr)
wcor<-cor(xvec,yvec)
wcov<-var(xvec,yvec)
if(sum(x==y)!=length(x)){
test<-wcor*sqrt((length(x)-2)/(1.-wcor^2))
sig<-2*(1-pt(abs(test),length(x)-2*g-2))
}
list(cor=wcor,cov=wcov,siglevel=sig)
}

winval<-
function(x,tr=.2){
#
#  Winsorize the data in the vector x.
#  tr is the amount of Winsorization which defaults to .2.
#
#  This function is used by several other functions that come with this book.
#
y<-sort(x)
n<-length(x)
ibot<-floor(tr*n)+1
itop<-length(x)-ibot+1
xbot<-y[ibot]
xtop<-y[itop]
winval<-ifelse(x<=xbot,xbot,x)
winval<-ifelse(winval>=xtop,xtop,winval)
winval
}

pitman.morgan.test <-
function(x, y = NULL, alternative = c("two.sided", "less", "greater"),
         omega = 1, conf.level = 0.95)
{
    alternative <- match.arg(alternative)

    if(!missing(omega) && (length(omega) != 1 || is.na(omega)))
        stop("'omega' must be a single number")
    if(!missing(conf.level) &&
       (length(conf.level) != 1 || !is.finite(conf.level) ||
        conf.level < 0 || conf.level > 1))
        stop("'conf.level' must be a single number between 0 and 1")
   
 if( !is.null(y) ) {
dname <- paste(deparse(substitute(x)),"and",deparse(substitute(y)))
}
   
    else {
stop("'y' is missing for paired test")
    }
    xok <- yok <- complete.cases(x,y)
    x <- x[xok]
    y <- y[yok]

 
alpha=1-conf.level
n<-length(x)
df=n-2
r<-cor(x,y)
Var1<-var(x)
Var2<-var(y)
w<-var(x)/var(y)
stat.t<-((w-omega)*sqrt(n-2))/sqrt(4*(1-r^2)*w*omega)

if (alternative == "two.sided") {
k<-qt(1-alpha/2,df=n-2)
K<-1+(2*(1-r^2)*k^2)/(n-2)
low<-w*(K-sqrt(K^2-1))
up<- w*(K+sqrt(K^2-1))
pval<-2*(1-pt(abs(stat.t),df=df))
    }

if (alternative == "less") {
k<-qt(1-alpha,df=n-2)
K<-1+(2*(1-r^2)*k^2)/(n-2)
low<-0
up<- w*(K+sqrt(K^2-1))
pval<-pt(stat.t,df=df)
    }

if (alternative == "greater") {
k<-qt(1-alpha,df=n-2)
K<-1+(2*(1-r^2)*k^2)/(n-2)
low<-w*(K-sqrt(K^2-1))
up<- Inf
pval<-1-pt(stat.t,df=df)
    }

estimate<-c(Var1, Var2)
cint<-c(low,up)
tstat<-stat.t


method <- c("Paired Pitman-Morgan test")
names(estimate) <- c("variance of x","variance of y")
    


    names(tstat) <- "t"
    names(df) <- "df"
    names(omega) <- c("ratio of variances")
    attr(cint,"conf.level") <- conf.level
    rval <- list(statistic = tstat, parameter = df, p.value = pval,
       conf.int=cint, estimate=estimate, null.value = omega,
       alternative=alternative,
       method=method, data.name=dname)
    class(rval) <- "htest"
    return(rval)
}


grambsch.test <-
function(x, y = NULL, alternative = c("two.sided", "less", "greater"))
{
    alternative <- match.arg(alternative)
   
 if( !is.null(y) ) {
dname <- paste(deparse(substitute(x)),"and",deparse(substitute(y)))
}
   
    else {
stop("'y' is missing for paired test")
    }
    xok <- yok <- complete.cases(x,y)
    x <- x[xok]
    y <- y[yok]

S<-x+y
D<-x-y
Z<-(D-mean(D))*(S-mean(S))
statistique<-sqrt(length(Z))*mean(Z)/sd(Z)
if(alternative=="two.sided"){probabilite<-2*pnorm(abs(statistique),lower.tail=FALSE)}
if(alternative=="less"){probabilite<-pnorm(statistique,lower.tail=TRUE)}
if(alternative=="greater"){probabilite<-pnorm(statistique,lower.tail=FALSE)}

tstat <- statistique
method <- c("Paired Grambsch test")
pval <- probabilite
    names(tstat) <- "z"
omega<-1
    names(omega) <- c("ratio of variances")

    rval <- list(statistic = tstat, p.value = pval,
       null.value = omega,
       alternative=alternative,
       method=method, data.name=dname)
    class(rval) <- "htest"
    return(rval)
}



bonett.seier.test <-
function(x, y = NULL, alternative = c("two.sided", "less", "greater"),
         omega = 1, conf.level = 0.95)
{
    alternative <- match.arg(alternative)

    if(!missing(omega) && (length(omega) != 1 || is.na(omega)))
        stop("'omega' must be a single number")
    if(!missing(conf.level) &&
       (length(conf.level) != 1 || !is.finite(conf.level) ||
        conf.level < 0 || conf.level > 1))
        stop("'conf.level' must be a single number between 0 and 1")
   
 if( !is.null(y) ) {
dname <- paste(deparse(substitute(x)),"and",deparse(substitute(y)))
}
   
    else {
stop("'y' is missing for paired test")
    }
    xok <- yok <- complete.cases(x,y)
    x <- x[xok]
    y <- y[yok]

alpha <- 1 - conf.level


n<-length(x)
tauX<-mean(abs(x-median(x)))
tauY<-mean(abs(y-median(y)))
ratio<-tauX/tauY

muX<-mean(x)
sdX<-sd(x)
gammaX<-(sdX^2)/(tauX^2)
deltaX<-(muX-median(x))/tauX
var.ln.tauX<-(gammaX+deltaX^2-1)/n

muY<-mean(y)
sdY<-sd(y)
gammaY<-(sdY^2)/(tauY^2)
deltaY<-(muY-median(y))/tauY
var.ln.tauY<-(gammaY+deltaY^2-1)/n

dX<-abs(x-median(x))
dY<-abs(y-median(y))
corD<-cor(dX,dY)
var.ln.ratio<- var.ln.tauX+ var.ln.tauY-2*corD*sqrt(var.ln.tauX* var.ln.tauY)



stat.t<-(log(ratio)-log(omega))/sqrt(var.ln.ratio)

if (alternative == "two.sided") {
z<-qnorm(1-alpha/2)
low<-exp(log(ratio)-z*sqrt(var.ln.ratio))
up<- exp(log(ratio)+z*sqrt(var.ln.ratio))
pval<-2*(1-pnorm(abs(stat.t)))
    }

if (alternative == "less") {
z<-qnorm(1-alpha)
low<-0
up<- exp(log(ratio)+z*sqrt(var.ln.ratio))
pval<-pnorm(stat.t,lower.tail=TRUE)
    }

if (alternative == "greater") {
z<-qnorm(1-alpha)
low<-exp(log(ratio)-z*sqrt(var.ln.ratio))
up<- Inf
pval<-pnorm(stat.t,lower.tail=FALSE)
    }


estimate<-c(tauX, tauY)
tstat<-stat.t
cint<-c(low,up)



    

method <- c("Paired Bonett-Seier test")
names(estimate) <- c("mean abs. dev. of x"," mean abs. dev. of y")
    


     

    names(tstat) <- "z"
    names(omega) <- c("ratio of means absolute deviations")
    attr(cint,"conf.level") <- conf.level
    rval <- list(statistic = tstat, p.value = pval,
       conf.int=cint, estimate=estimate, null.value = omega,
       alternative=alternative,
       method=method, data.name=dname)
    class(rval) <- "htest"
    return(rval)
}



### statistiques classiques pour donnees appariees

paired.summary<-
function(x,y){
tt<-matrix(numeric(44),nrow=4)
X<-x
Y<-y
# calculs pour x
x<-X[!is.na(X)]
tt[1,1]<-length(x)
tt[1,2]<-mean(x)
tt[1,3]<-median(x)
tt[1,4]<-mean(x,0.2)
tt[1,5]<-sd(x)
tt[1,6]<-IQR(x)/1.349
tt[1,7]<-mad(x)
tt[1,8]<-mean(abs(x-median(x)))*sqrt(pi/2)
tt[1,9]<-sqrt(winvar(x))/0.6421
tt[1,10]<-min(x)
tt[1,11]<-max(x)

# calculs pour y
y<-Y[!is.na(Y)]
tt[2,1]<-length(y)
tt[2,2]<-mean(y)
tt[2,3]<-median(y)
tt[2,4]<-mean(y,0.2)
tt[2,5]<-sd(y)
tt[2,6]<-IQR(y)/1.349
tt[2,7]<-mad(y)
tt[2,8]<-mean(abs(y-median(y)))*sqrt(pi/2)
tt[2,9]<-sqrt(winvar(y))/0.6421
tt[2,10]<-min(y)
tt[2,11]<-max(y)

xok <- yok <- complete.cases(X,Y)
    x <- X[xok]
    y <- Y[yok]

# calculs pour x-y
tt[3,1]<-length(x-y)
tt[3,2]<-mean(x-y)
tt[3,3]<-median(x-y)
tt[3,4]<-mean(x-y,0.2)
tt[3,5]<-sd(x-y)
tt[3,6]<-IQR(x-y)/1.349
tt[3,7]<-mad(x-y)
tt[3,8]<-mean(abs(x-y-median(x-y)))*sqrt(pi/2)
tt[3,9]<-sqrt(winvar(x-y))/0.6421
tt[3,10]<-min(x-y)
tt[3,11]<-max(x-y)

# calculs pour (x+y)/2
tt[4,1]<-length((x+y)/2)
tt[4,2]<-mean((x+y)/2)
tt[4,3]<-median((x+y)/2)
tt[4,4]<-mean((x+y)/2,0.2)
tt[4,5]<-sd((x+y)/2)
tt[4,6]<-IQR((x+y)/2)/1.349
tt[4,7]<-mad((x+y)/2)
tt[4,8]<-mean(abs((x+y)/2-median((x+y)/2)))*sqrt(pi/2)
tt[4,9]<-sqrt(winvar((x+y)/2))/0.6421
tt[4,10]<-min((x+y)/2)
tt[4,11]<-max((x+y)/2)

colnames(tt)<-c("n","mean","median","trim","sd","IQR (*)","median ad (*)","mean ad (*)","sd(w)","min","max")
rownames(tt)<-c("x","y","x-y","(x+y)/2")
tt
}


paired.effect.size<-
function(x,y){
tt<-matrix(numeric(8),nrow=2)

### Classical effect sizes
tt[1,1]<-(mean(x)-mean(y))/sqrt((sd(x)^2+sd(y)^2)/2)

tt[1,2]<-(mean(x)-mean(y))/sd(x)

tt[1,3]<-(mean(x)-mean(y))/sd(y)

d<-x-y

tt[1,4]<-mean(d)/sd(d)


### Robust effect sizes (trim=0.2)

Vx<-winvar(x,tr=0.2)
Vy<-winvar(y,tr=0.2)
Vd<-winvar(d,tr=0.2)

tt[2,1]<-0.642*(mean(x,tr=0.2)-mean(y,tr=0.2))/sqrt((Vx+Vy)/2)

tt[2,2]<-0.642*(mean(x,tr=0.2)-mean(y,tr=0.2))/sqrt(Vx)

tt[2,3]<-0.642*(mean(x,tr=0.2)-mean(y,tr=0.2))/sqrt(Vy)

tt[2,4]<-0.642*mean(d,tr=0.2)/sqrt(Vd)

colnames(tt)<-c("Average","Single (x)","Single (y)","Difference")
rownames(tt)<-c("OLS","Robust")
tt
}












### Les fonctions graphiques

plotCor<-function(df,condition1,condition2,groups=NULL,facet=TRUE,...){
plotP<-ggplot(data=df)+aes_string(x=condition1,y=condition2)
if(is.null(groups)){
plotP+geom_point()+coord_equal()+ geom_abline(intercept =0, slope =1)
}
else{
if(facet){
formula<-paste(groups,"~.",sep="")
plotP+geom_point()+coord_equal()+ geom_abline(intercept =0, slope =1)+facet_grid(formula)
}
else{
plotP+geom_point()+coord_equal()+ geom_abline(intercept =0, slope =1)+aes_string(colour=groups)
}
}
}


plotBA<-function(df,condition1,condition2,groups=NULL,facet=TRUE,...){
df2<-df
df2$M<-(df[,condition1]+df[,condition2])/2
df2$D<-(df[,condition2]-df[,condition1])
plotP<-ggplot(data=df2)+aes(x=M,y=D)
if(is.null(groups)){
plotP+geom_point()+geom_abline(intercept=0,slope=0)+geom_smooth(method="lm",formula="y~1")+xlab("Mean")+ylab("Difference")
}
else{
if(facet){
formula<-paste(groups,"~.",sep="")
plotP+geom_point()+geom_abline(intercept=0,slope=0)+geom_smooth(method="lm",formula="y~1",fullrange=TRUE)+xlab("Mean")+ylab("Difference")+facet_grid(formula)
}
else{
plotP+geom_point()+geom_abline(intercept=0,slope=0)+geom_smooth(method="lm",formula="y~1",fullrange=TRUE)+aes_string(colour=groups,fill=groups)+xlab("Mean")+ylab("Difference")
}
}
}

unroll<-function(df,c1,c2,subjects){
df2<-data.frame(c(df[,c1],df[,c2]),rep(c(c1,c2),rep(length(df[,c1]),2)),rep(df[,subjects],2))
names(df2)<-c("Measurements","Conditions","Subjects")
df2
}
unrollG<-function(df,c1,c2,subjects,groups){
df2<-data.frame(c(df[,c1],df[,c2]),rep(c(c1,c2),rep(length(df[,c1]),2)),rep(df[,subjects],2),rep(df[,groups],2))
names(df2)<-c("Measurements","Conditions","Subjects","Groups")
df2
}

plotMcNeil<-function(df,condition1,condition2,subjects,groups=NULL,facet=TRUE,...){
if(is.null(groups)){
df2<-unroll(df,condition1,condition2,subjects)
}
else{
df2<-unrollG(df,condition1,condition2,subjects,groups)
}
df2[,"Subjects"]<-reorder(df2[,"Subjects"],df2[,"Measurements"])
plotP<-ggplot(data=df2)+aes_string(x="Measurements",y="Subjects",group="Subjects")
if(is.null(groups)){
plotP+geom_point()+aes_string(colour="Conditions")
}
else{
if(facet){
formula<-"Groups~."
plotP+geom_point()+aes_string(colour="Conditions")+facet_grid(formula,scales="free_y")
}
else{
plotP+geom_point()+aes_string(colour="Conditions")+aes_string(shape="Groups")
}
}
}


plotProfiles<-function(df,condition1,condition2,subjects,groups=NULL,facet=TRUE,...){
if(is.null(groups)){
df2<-unroll(df,condition1,condition2,subjects)
}
else{
df2<-unrollG(df,condition1,condition2,subjects,groups)
}
plotP<-ggplot(data=df2)+aes_string(x="Conditions",y="Measurements",group="Subjects")
if(is.null(groups)){
plotP+geom_boxplot(aes(group=NULL) ,width=0.1)+geom_line()
}
else{
if(facet){
formula<-"Groups~."
plotP+geom_boxplot(aes(group=NULL) ,width=0.1)+geom_line()+facet_grid(formula)
}
else{
plotP+geom_boxplot(aes(group=NULL) ,width=0.1)+geom_line(aes_string(colour="Groups"))+aes_string(fill="Groups")
}
}
}

plotSliding<-
function (df,condition1,condition2,xlab = "", ylab = "", ...) 
{
x<-df[,condition1]
y<-df[,condition2]
    require(MASS)
    par(mar = c(2, 2, 5, 4))
    mX <- min(x)
    MX <- max(x)
    mY <- min(y)
    MY <- max(y)
    addX <- (max(x) - min(x))/25
    addY <- (max(y) - min(y))/25
    addXY <- min(addX, addY)
    minX <- mX - (MY - mY)/2
    maxX <- MX + (MY - mY)/2
    minY <- mY - (MX - mX)/2
    maxY <- MY
    eqscplot(x, y, xlim = c(minX - addX - addXY, maxX + addX + 
        addXY), ylim = c(minY - addY - addXY, maxY + addY + addXY), 
        axes = FALSE, xlab = "", ylab = "", ...)
    polygon(c(mX - addX, mX - addX, MX + addX, MX + addX), c(mY - 
        addY, MY + addY, MY + addY, mY - addY), lty = 3)
    rug(x, side = 3, ticksize = 0.02)
    rug(y, side = 4, ticksize = 0.02)
    lines(c(minX, (mX + MX)/2 + addXY) - addX, c((mY + MY)/2, 
        minY - addXY) - addY, col = "red", lwd = 0.25)
    lines(c((mX + MX)/2 - addXY, maxX) + addX, c(minY - addXY, 
        (mY + MY)/2) - addY, col = "blue", lwd = 0.25)
    axis(3, at = pretty(x))
    if (xlab == "") {
        xlab <- condition1
    }
    mtext(xlab, side = 3, line = 2)
    axis(4, at = pretty(y))
    if (ylab == "") {
        ylab <- condition2
    }
    mtext(ylab, side = 4, line = 2)
    xy <- (x + y)/2
    lines(c(max(mX - addX, mY - addY), min(MX + addX, MY + addY)), 
        c(max(mX - addX, mY - addY), min(MX + addX, MY + addY)), 
        lty = 2)
    xx <- mX - addX
    yy <- mY - addY
    c <- xx + yy
    points(c/2 - addXY + (x - y)/2, c/2 - addXY + (y - x)/2, 
        col = "red")
    x1 <- MX + addX
    y1 <- mY - addY
    b <- y1 - x1
    e <- b/2
    points(xy - e + addXY, xy + e - addXY, col = "blue")
}
























