# Fonksiyon 7.1: İşaret testi fonksiyonu
sign.test <- function(x, y=NULL, mu=0, alternative="two.sided",
  conf.level=0.95){
  if(is.null(y)) y <- mu
  n <- sum((y - x)!=0)
  sneg <- sum(x < y)
  spos <- sum(x > y)
  s <- spos
  if(alternative=="greater")
    s <- sneg
  if(alternative=="two.sided")
    s <- min(sneg, spos)
  sonuc <- signtest.prob(n=n, s=s, alternative=alternative)
  return(sonuc)
}

# Fonksiyon 7.2: İşaret testi olasılık hesaplama fonksiyonu
# Binomiyal dağılış olasılık hesaplama fonksiyonu
signtest.prob <-function(n, s, alternative="two.sided"){
  if (alternative=="less") 
    p.value <- pbinom(q=s, size=n, prob=0.5)
  else if (alternative=="greater")
    p.value <- 1-pbinom(q=s, size=n, prob=0.5, 
      lower.tail=FALSE)
  else if (alternative=="two.sided"){
    p1 <- pbinom(q=s, size=n, prob=0.5)
     p2 <- 1-pbinom(q=s, size=n, prob=0.5)
     p.value <- 2*min(p1, p2)
  }	
  list(hipotez=alternative, n=n, s=s, p.value=p.value)
}

# Fonksiyon 7.3: Pratt yöntemi ile sıra değeri hesaplama
pratt.ranks <- function(d){
  d <- abs(d)
  zidx <- which(d==0)
  xm <- replace(x,zidx, NA)
  rx <- rank(xm, ties.method="average", na.last=FALSE)  
  xr <- replace(rx,zidx, 0)
  return(xr)
}

# Fonksiyon 7.4: WSRT W dağılışı için kesin test
wsrt.exact <- function(x, y=NULL, mu=0, 
  alternative="two.sided", conf.level=0.95){
  if(!is.null (y)){
    d <- x-y
    ad <- abs(d - mu)
  }else{
    d <- x-mu
    ad <- abs(d)
  }
  n <- length(d)
  rnk <- ad
  zidx <- which(ad!=0)
  rx <- rank(ad[zidx])
  rnk[zidx] <- rx
  rnk[-zidx] <- NA
  if(!is.null(y))
     atable <- cbind(x, y, d, sgn=sign(d), rnk)
  else
     atable <- cbind(x, d, sgn=sign(d), rnk)
  Wneg <- sum(rnk[d < 0])
  Wpos <- sum(rnk[d > 0])
  n0 <- length(rnk[d == 0])
  n <- length(rnk[d != 0])
  if(alternative=="less"){
    W <- Wpos
    p.value <- psignrank(W, n)
  }
  if(alternative=="greater"){
    W <- Wneg
    p.value <- 1-psignrank(W, n, lower.tail=FALSE) 
  }
  if(alternative=="two.sided"){
    W <- min(Wneg, Wpos)
    p.value <- 2*min(psignrank(W, n), 1-psignrank(W, n,
    lower.tail=FALSE))
  }
  walmat <- outer(d, d, "+") / 2
  walvec <- walmat[lower.tri(walmat, diag = TRUE)]
  n.wal <- length(walvec)
  est.med <- median(walvec)
  sig.level <- 1 - conf.level
  k <- qsignrank(sig.level/2, n)
  if (k == 0) k <- k + 1
  achieved.CI <-  1 - 2 * psignrank(k-1, n)
  pop.med.CI <- c(walvec[k], walvec[n.wal + 1 - k])
  sonuc <- list(atable=atable, Wneg=Wneg, Wpos=Wpos, 
    n=n, n0=n0, W=W, est.med=est.med, 
    achieved.CI= achieved.CI,
    pop.med.CI=pop.med.CI, p.value=p.value)
  return(sonuc)
}

# Fonksiyon 7.5: WSRT W dağılışı için yaklaşma testi
wsrt.approx <- function(x, y=NULL, mu=0, 
  method=c("wilcoxon", "pratt"), approx=c("normal","t"),
  alternative="two.sided", conf.level=0.95,
  c.cont=c(0,1/2,1/4,3/8), c.ties=TRUE){
  if(missing(method)) method <- "wilcoxon"
  if(missing(approx)) approx <- "normal"
  if(missing(c.cont)) c.cont <- 1/2
  if(!is.null (y)){
    d <- x-y
    ad <- abs(d - mu)
  }else{
    d <- x-mu
    ad <- abs(d)
  }
  n <- length(d)
  rnk1 <- ad
  zidx <- which(ad!=0)
  rx <- rank(ad[zidx])
  rnk1[zidx] <- rx
  rnk1[-zidx] <- NA
  xm <- replace(ad,-zidx, NA) #rank with zero
  rx <- rank(xm, ties.method="average", na.last=FALSE)
  rnk2 <- replace(rx,-zidx, 0)
  rnk <- rnk1
  if(method=="pratt") rnk <- rnk2
  trnk <- tabulate(rnk[which(rnk>0)])     	#ties
  sum.ties <- sum(trnk*(trnk-1)*(trnk+1))
  if(!is.null(y))
    atable <- cbind(x,y,d,sgn.d=sign(d),abs.d=ad,rnk1, rnk2)
  else
    atable <- cbind(x,d,sgn.d=sign(d),abs.d=ad,rnk1,rnk2)
  Wneg <- sum(rnk[which(d < 0)])
  Wpos <- sum(rnk[which(d > 0)])
  n0 <- length(which(d == 0))
  n <- length(d)-n0
  expw1 <- n*(n+1)/4
  expw2 <- (n*(n+1)-n0*(n0+1))/4
  varw1 <- n*(n+1)*(2*n+1)/24
  varw2 <- (n*(n+1)*(2*n+1)-n0*(n0+1)*(2*n0+1))/24
  expw <- expw1
  varw <- varw1
  if(method=="pratt") {
     expw <- expw2
     varw <- varw2
  }
  if(c.ties) 
     varw <- varw - 1/2*sum.ties/24
  if(alternative=="less"){
     W <- Wpos
     if(approx=="normal")
        p <- pnorm(W+c.cont, expw, sqrt(varw))
     else 
        p <- pt(q=W, df=n-1) 
  }
  if(alternative=="greater"){
     W <- Wneg
     if(approx=="normal")
        p <- 1-pnorm(W-c.cont, expw, sqrt(varw), 
             lower.tail=FALSE)
     else 
       p <- 1-pt(q=W, df=n-1, lower.tail=FALSE) 
  }
  if(alternative=="two.sided"){
     W <- min(Wneg, Wpos)
     if(approx=="normal"){
       p1 <- pnorm(W+c.cont, expw, sqrt(varw)) 
     p2 <- 1-pnorm(W-c.cont, expw, sqrt(varw),lower.tail=FALSE) 
       p <- 2*min(p1, p2) 
     }else{
        p1 <- pt(q=W, df=n-1) 
        p2 <- 1-pt(q=W, df=n-1,lower.tail=FALSE) 
        p <- 2*min(p1, p2) 
    }
  }
  walmat <- outer(d, d, "+") / 2
  walvec <- walmat[lower.tri(walmat, diag=TRUE)]
  n.wal <- length(walvec)
  est.med <- median(walvec)
  sig.level <- 1 - conf.level
  k <- qsignrank(sig.level/2, n)
  if (k == 0) k <- k + 1
  achieved.CI <-  1 - 2 * psignrank(k - 1, n)
  pop.med.CI <- c(walvec[n.wal + 1 - k], walvec[k])
  sonuc <- list(atable=atable, Wneg=Wneg, Wpos=Wpos, W=W,
    n=n, n0=n0, sum.ties=sum.ties, p.value=p, est.med=est.med,
    achieved.CI=achieved.CI, pop.med.CI=pop.med.CI)
  return(sonuc)
}

# Fonksiyon 8.1: Permütasyon testi
# Bağımlılık – Fonksiyon: 8.2
#
perm.test <- function(gozlem, grup, stat=c("mean", "median"),
   N=10000){
   if(missing(stat)) stat="mean"
   n <- length(gozlem)/2
   xidx <- grup == unique(grup)[1]
   yidx <- !xidx
   nx <- sum(xidx)
   ny <- sum(yidx)
   w <- replicate(N, stat.diff(gozlem, nx, ny, n, stat=stat))
   x <- gozlem[xidx]
   y <- gozlem[yidx]
   if(stat=="mean"){
      smean <- mean(x)-mean(y)
      (mean(abs(w) > abs(smean)))
   }else{
      smedian <- median(x)-median(y)
     (mean(abs(w) > abs(smedian)))
   }
}

# Fonksiyon 8.2: Örneklem istatistikleri farkını hesaplama
stat.diff <- function(xdata, nx, ny, n, stat="mean"){
   xidx <- sample.int(nx+ny, n)
   x <- xdata[xidx]
   y <- xdata[-xidx]
   if(stat=="mean"){
      w <- mean(x)-mean(y)
   }else{
      w <- median(x)-median(y)
   }
}

# Fonksiyon 9.1: İzotonik regresyon test istatistiği 
# Bağımlılık – Fonksiyon: 9.2
test.stat <- function(x, w) {
   mu <- pava(x, w) # w tartı faktörüdür
   mu0 <- sum(w * x) / sum(w)
   ss.alt <- sum(w * (x-mu)^2)
   ss.null <- sum(w * (x-mu0)^2)
   return(ss.null-ss.alt)
}  

# Fonksiyon 9.2: pava fonksiyonu
pava <- function(x, w) {
   if (any(w <= 0))
       stop("Ağırlıklar pozitif olmalıdır")
   if (length(x) != length(w))
       stop("Argümanlar aynı uzunlukta olmalıdır")
   n <- length(x)
   design <- diag(n)
   repeat {
       out <- lm(x ~ design + 0, weights = w)
       mu <- coefficients(out)
       dmu <- diff(mu)
       if (all(dmu >= 0)) break
       j <- min(seq(along = dmu)[dmu < 0])
       design[ , j] <- design[ , j] + design[ , j + 1]
       design <- design[, -(j + 1), drop=FALSE]
   }
   return(as.numeric(design %*% mu))
}

# Fonksiyon 10.1: Durbin testi fonksiyonu 
durbinD <- function(x){
 if(!is.matrix(x)) stop("x matris olmalı")
 y <- x
 for(i in 1:nrow(y)) 
   y[i,] <- rank(x[i,], na.last="keep")
 t <- ncol(y) ; b <- nrow(y)
 r <- length(which(y[1,]!="NA"))
 k <- length(which(y[,1]!="NA"))
 Rj <- colSums(y, na.rm=TRUE)
 Rjo <- Rj/r
 C <- 12*(t-1)*r/(t*(k^2-1))
 V <- (1/4)*(sum(y^2,na.rm=TRUE)/((r*t)-(k+1)^2)) 
 C1 <- (t-1)*r/(t*V)  
 C <- C-C1
 D <- C * sum((Rjo-(k+1)/2)^2)
 v1 <- t-1 						
 v2 <- b*k-t-b+1 
 p.value <- pchisq(q=D, df=v1,lower.tail=FALSE)
 sonuc <- list(y=y, D=D, C=C, C1=C1, 
   df.1=v1, df.2=v2, p.value=p.value)
 return(sonuc)
}

# Fonksiyon 11.1: Steel-Dwass testi fonksiyonu
sd.test <- function(x, grup){        
 xtam <- complete.cases(x, grup)      
 x <- x[xtam]    
 g <- grup[xtam]    
 nj <- table(g)           
 ng <- length(nj)  
 t <- combn(ng, 2, function(ij){
   i <- ij[1]      
   j <- ij[2]          
   r <- rank(c(x[g == i], x[g == j]))    
   R <- sum(r[1:nj[i]])                
   n <- nj[i]+nj[j]      
   xort <- nj[i]*(n+1)/2                
   xsd <- sqrt(nj[i]*nj[j]/(n*(n-1))*(sum(r^2)-n*(n+1)^2/4))  
   return(abs(R-xort)/xsd)           
  })    
  p <- ptukey(t*sqrt(2), ng, Inf, lower.tail=FALSE)  # P-değeri    
  sonuc <- cbind(t, p)                  
  rownames(sonuc) <- combn(ng, 2, paste, collapse=":")    
  return(sonuc)
}

# Fonksiyon 11.2: Düzeltilmiş p değerleri grafiği fonksiyonu
plot.pvalues <- function(pv){
 renkpalet <- c(1:6,"#9920AA")
  pv <- pv[order(pv[,1]),]
  plot(pv[,1], pv[,1], type="b", col=1, 
    lwd=2, ylim=c(0,1), pch=1,
    xlab="p değeri", ylab="Düzeltilmiş p-değeri")
  n <- nrow(pv)
  k <- ncol(pv)
  for(j in 2:k){
    lines(pv[,1], pv[,j], 
    lty=j, col=renkpalet[j], type="b", lwd=2, pch=j)
  }
  legend(0,0.95,legend=colnames(pv), 
    col=renkpalet[1:n], lty=1:n, cex=0.8)
}

# Fonksiyon 13.1: LOOCV çapraz doğrulama fonksiyonu
loocv <- function(x, y, h=1.0){
  n <- length(x)
  tahmin.hata <- rep(NA, n)
  for(i in 1:n){
    # Test verisi
    x.test <- x[i]
    y.test <- y[i]
    # Eğitim seti
    x.egitim <- x[-i]
    y.egitim <- y[-i]
    # Tahmin değeri
    y.tahmin <- ksmooth(x=x.egitim, y=y.egitim, 
       kernel = "normal", bandwidth=h, x.points = x.test)
    # Tahmin hatası
    tahmin.hata[i] <- (y.test - y.tahmin$y)^2
  }
  #Tahmin ortalaması
  return(mean(tahmin.hata, na.rm=TRUE))
}

# Fonksiyon 13.2: k-CV çapraz doğrulama fonksiyonu
kcv <- function(x, y, h=1, k=2){
   xy <- data.frame(x,y)
   xy <- xy[sample(nrow(xy)),]
   katlar <- cut(seq(1, nrow(xy)), breaks=k, labels=FALSE)
   hata <- 0
   for(i in 1:k){
    idx <- which(katlar==i, arr.ind=TRUE)
    xy.test <- xy[idx, ]
    xy.egitim <- xy[-idx, ]
    kreg <- ksmooth(x=xy.egitim$x, y=xy.egitim$y, 
        kernel="normal", bandwidth=h, x.points=xy.test$x)
     hata <- hata + mean((xy.test$y[order(xy.test$x)]-
        kreg$y)^2, na.rm=T)
   }
   hata <- hata / k
   return(hata)
}

# Fonksiyon 14.1: Kesin koşu testi
# Bağımlılık – Fonksiyon: 14.2
runs.exact.test <- function(rmin, rmax, m1, m2){
  plt <- 0 ; pgt <- 0
  for(i in rmin:m1){
   plt <- plt + prob.calc(n1,n2,i)
  }
  for(i in m2:rmax){
   pgt <- pgt + prob.calc(n1,n2,i)
  }
  p <- plt + pgt
  list(p=p, plt=plt, pgt=pgt)
}

# Fonksiyon 14.2: Olasılık hesaplama
# Bağımlılık – Fonksiyon: 14.3a,b
prob.calc <- function(n1, n2, i){
   k <- i%/%2
   if(i%%2==0) 
     p<- 2*comb(n1-1,k-1)* comb(n2-1,k-1)/comb(n1+n2, n1)
   else 
     p <- (comb(n1-1,k)*comb(n2-1,k-1)+
           comb(n1-1,k-1)*comb(n2-1,k))/comb(n1+n2, n1)
  return(p)
}

# Fonksiyon 14.3a: Kombinasyon hesaplama
comb <- function(n,r){
   (factorial(n)/(factorial(r)*factorial(n-r)))
}

# Fonksiyon 14.3b: Permütasyon hesaplama
perm <- function(n,r){
   (factorial(n)/factorial(n-r))
}

# Fonksiyon 14.4: Kategorik koşu testi fonksiyonu
kat.kosu.test <- function(x) {
  xf <- factor(x)
  n <- length(xf)
  ni <- table(xf)
  ni <- structure(.Data=as.vector(ni), .Names=names(ni))
  expr <- (n*(n+1)-sum(ni^2))/n
  vnom <- sum(ni^2*(sum(ni^2)+n*(n+1)))-2*n*sum(ni^3)-n^3
  vdenom <- n^2*(n-1)
  varr <- vnom/vdenom
  sdr <- sqrt(varr)
  r <- sum(diff(as.numeric(xf)) != 0)+1 #NCSS’de +1 yok
  cats <- names(ni)
  xr <- rle(x)
  runs <- c()
  for(cat in cats){
    runs[cat] <- length(which(xr$values==cat))
  }
  names(runs)<- cats
  z.stat <- (r-expr) / sdr
  p.value <- 2*(1-pnorm(abs(z.stat)))
  sonuc <- list(cat.size=ni,n=n, cat.runs=runs, 
    r=r, statistics=z.stat,p.value=p.value)
  return(sonuc)
}

# Fonksiyon 14.5: Kategorik dizilerde koşu sayısını saptama
kosu.say <- function(x){
  kosu <- 0
  for(i in 2:length(x)){
    if(x[i]!=x[i-1]) kosu <- kosu+1
  }
  kosu <- kosu+1
  return(kosu)
}

# Fonksiyon 14.6: Mann-Kendall Testi için işaretler tablosu
mk.isaret.tablosu <- function(x){
 n <- length(x)
 sayi <- 0
 for(i in 2:n){
   s <- x[i]-x[1:i-1]
   cat(sign(s), ":", sum(sign(s)), "\n")
   sayi <- sayi + sum(sign(s))
 }
 return(sayi)
}

# Fonksiyon 15.1: Fi-kare ölçüsü fonksiyonu
phisq <- function(x){
  n <- sum(x)
  phisq <- chisq.test(x, correct=FALSE)$statistic/n 
  print.noquote("Phi Sq.:")
  return(as.numeric(phisq))
}

# Fonksiyon 15.2: Cramér’in V katsayısı fonksiyonu
cramer.v <- function(x) {
  chi.sq <- chisq.test(x, correct=FALSE)$statistic
  m <- min(nrow(x),ncol(x))-1
  n <- sum(x)
  crmV <- sqrt(chi.sq /(n*m))
  print.noquote("Cramér V: ")
  return(as.numeric(crmV))
}

# Fonksiyon 15.3: Yule’un Q katsayısı
yule.q <- function(x){
  q <- (x[1,1]*x[2,2] - x[1,2]*x[2,1])/ (x[1,1]*x[2,2] +
    x[1,2]*x[2,1])
  print.noquote("Yule’s Q:")
  return(q)
}

