
#########################################################################################################
#######                                                                                      ############
###### copy file name's[skewnormal_var_ridgeline_May12_working.R] dated 05/06/2017-13:15   ##############
#########################################################################################################
#********************************************************************************************************
#********************************************************************************************************


##########################################################################################################
# This code to find the ridgeline of skew normal with different varince
# this code has been fixed by Dr.Ray on 16- May 2017
#
##########################################################################################################
########################### ridgeline of skwe normal #####################################################
##########################################################################################################
#-----------------------------------------------------------

# Now starting for bivariate skew mixtures with equal variance to start with
# For derivative
#  This function uses inverse of square root matrices to ease computational burden


#' Title
#'
#' @param x  test
#' @param alpha1 test
#' @param alpha2 test
#' @param xi1 test
#' @param xi2 test
#' @param inv.omega1 test
#' @param inv.omega2 test
#' @param inv.w1 test
#' @param inv.w2 test
#' @param a test
#'
#' @return test
#' @export test
#'
#' @examples test
ridge.skew.var<-function(x,alpha1,alpha2,xi1,xi2,inv.omega1,inv.omega2,inv.w1,inv.w2,a){
  # lambda skewnees paramater and xi's the mean

  z1=inv.omega1%*%(x-xi1)
  z2=inv.omega2%*%(x-xi2)

  u1= as.vector(alpha1%*%inv.w1%*%(x-xi1))
  u2=  as.vector(alpha2%*%inv.w2%*%(x-xi2))

  ##    cat(dnorm(u2),"\t",pnorm(u2),"\n")

  ## If else to avoid division by 0
  if(pnorm(u1)>0)          comp1= -z1 + inv.w1%*%alpha1 *dnorm(u1)/pnorm(u1)
  else  comp1= -z1
  if(pnorm(u2)>0)
    comp2= -z2 +  inv.w2%*%alpha2 * dnorm(u2)/pnorm(u2)
  else  comp2= -z2

  return( ((1-a)*comp1+a*comp2))

}



ridgeline.skew.var<-function(alpha1,alpha2,xi1,xi2,omega1,omega2,by=.01){


  ## Here we changed the x1i to d , to added extrea paremeter to get 4 mode
  d=length(xi1)
  avec <- seq(0,1,by=by)#alpha

  ridge=matrix(0,length(avec),d);
  dens=rep(0,length(avec));

  ind=1;

  mymode=xi1
  ind=1
  inv.omega1=solve(omega1)
  inv.omega2=solve(omega2)
  inv.w1=diag(1/sqrt(diag(omega1)))
  inv.w2=diag(1/sqrt(diag(omega2)))


  for ( a in avec){

    ss <- multiroot(f=function(x) ridge.skew.var(x,alpha1,alpha2,xi1,xi2,inv.omega1,inv.omega2,inv.w1,inv.w2,a), mymode)
    mymode<-ss$root
    ridge[ind,]=mymode;
    dens[ind]= (dmsn(mymode, xi1, omega1, alpha1)+dmsn(mymode, xi2, omega2, alpha2))/2
    ind=ind+1

    #points(mymode[1],mymode[2],pch="*",cex=2,col=2)


  }
  return(data.frame(alpha=avec,ridge=ridge,dens=dens))
}






contour.skew.var<-function(alpha1,alpha2,xi1,xi2,omega1,omega2,by=.01){
  x <- seq(-2,4, length.out = 500)
  y <- seq(-2, 4, length.out = 500)
  d1 <- expand.grid(x = x, y = y)



  #    xi1 <- c(0, 0)
  #     Omega[2,1] <- Omega[1,2] <- 0
  #     alpha1 <- c(0,a)
  pdf1 <- dmsn(d1, xi1, omega1, alpha1)
  #   contour(x,y,matrix(pdf1,500,500))


  pdf2 <- dmsn(d1, xi2, omega2, alpha2)
  pdf=0.5*(pdf1+pdf2)

  contour(x,y,matrix(pdf,500,500))
  title(main="Two dimensional equal variance example with three modes",sub="mean at (0,0) and (1, -1) unit variance and skewness=10")
}


contour.skew.dat<-function(data=data,alpha1,alpha2,xi1,xi2,omega1,omega2,by=.01,...){
  # get contour dimensions from data
  x <- seq(min(data[,1]),max(data[,1]), length.out = 500)
  y <- seq(min(data[,2]),max(data[,2]), length.out = 500)
  d1 <- expand.grid(x = x, y = y)



  #    xi1 <- c(0, 0)
  #     Omega[2,1] <- Omega[1,2] <- 0
  #     alpha1 <- c(0,a)
  pdf1 <- dmsn(d1, xi1, omega1, alpha1)
  #   contour(x,y,matrix(pdf1,500,500))


  pdf2 <- dmsn(d1, xi2, omega2, alpha2)
  pdf=0.5*(pdf1+pdf2)

  contour(x,y,matrix(pdf,500,500),...)
}





##########################################################################################################
###################### Mraging the ridgeline by the ration of maxmum , and low point #### ################
##########################################################################################################
########################################################################################################
##########################################################################################################

localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}

#######################################################################################################
#######################################################################################################

ridgeline.diagnosis.skew<-function (mixsmsn.out,data=NULL,
                                    ipairs = "all", compute.ratio = FALSE, by = 0.001, ratiocutoff = 0.8,
                                    ridgelineplot = "matrix")
{

#  ridgelineplot	one of "none", "matrix", "all"  "megred",
  # propvector<-mixsmsn.out
  if(!is(mixsmsn.out,"Skew.normal")) {
    cat("Object not of class Skew.Normal \n First perform a skew normal fit to the data\n")
      return(1)
    }


  mulist= mixsmsn.out$mu
  Sigmalist= mixsmsn.out$Sigma
  alphalist=mixsmsn.out$shape
  original.group=mixsmsn.out$group
  k=length(mixsmsn.out$mu)


  comat <- diag(k)
  ratiomatrix <- comat

 # if (compute.ratio)
#    ratiomatrix <- comat
#  else ratiomatrix <- NULL
  pairlist <- ipairs
  if (identical(ipairs, "all")) {
    pairlist <- list()
    m <- 1
    for (i in 1:(k - 1)) for (j in (i + 1):k) {
      pairlist[[m]] <- c(i, j)
      m <- m + 1
    }
  }
  m <- length(pairlist)
 # Create matrix to store the ridgeline dataframe for each [pair]
  bb.list<-rep(list(diag(2)), m)


  if (ridgelineplot == "matrix" | ridgelineplot == "all" )
    par(mfrow = c(3,3))

  for (q in 1:m) {
    if (q == 1)
      ia <- TRUE
    else ia <- !(pairlist[[q]][1] == pairlist[[q - 1]][1])
    i <- pairlist[[q]][1]
    j <- pairlist[[q]][2]


          bb=ridgeline.skew.var(as.vector(alphalist[[i]]),as.vector(alphalist[[j]]),as.vector(mulist[[i]]),as.vector(mulist[[j]]),round(Sigmalist[[i]],10),round(Sigmalist[[j]],10))  # CHANGE XI1 AND XI2 AND OMEGA'S
          bb.list[[q]]=bb
    #    contour.skew.dat(data=faithful,alpha1,alpha2,xi1,xi2,omega1,omega2)
    #       if(max(bb$den)!=max(bb$dens[1],bb$dens[length(bb$dens)]) )

    # Plots the ridgeline





    #Calculates the number of peaks
    peaks<-localMaxima(bb$dens)
    zerosij <-length(peaks)


      comat[i, j] <- comat[j, i] <- (zerosij <   2)
      ratiomatrix[i,j]=ratiomatrix[j,i]=min(bb$dens)/min(bb$dens[peaks])

      if(zerosij <2) ratiomatrix[i,j]=ratiomatrix[j,i]<-comat[i,j]



    if(ridgelineplot=="matrix" |  ridgelineplot == "all" ){
        plot(bb$alpha,bb$dens,type="l",xlab='alpha',ylab='Density',main=paste("Component",i," ",j))
        if(compute.ratio)
          legend("topleft",paste("ratio=",round(ratiomatrix[i,j],2)),bty="n",text.col="blue")
  }
    }




 if(ridgelineplot=="merged" |  ridgelineplot == "all" ){
     if(is.null(data)) cat ("No data to create contour plot")

      else{ par(mfrow=c(1,1))
        smoothScatter(data,nrpoints=Inf,pch=16,cex=.3)
        for (q in 1:m) {
          if (q == 1)
            ia <- TRUE
          else ia <- !(pairlist[[q]][1] == pairlist[[q - 1]][1])
          i <- pairlist[[q]][1]
          j <- pairlist[[q]][2]

          bb=bb.list[[q]]
          points(bb[c(1,nrow(bb)),2:3],pch="*",col="red",cex=2)

           if(compute.ratio){
             if(ratiomatrix[i,j]>ratiocutoff){ # Merging based on mathematical modes -> if 1 mode merge
         #      cat("Merging based ratio >", ratiocutoff,"\n")
               lines(bb$ridge.1,bb$ridge.2,col=2,lwd=2)
             }

           }

          else{
          if(comat[i,j]==1){ # Merging based on mathematical modes -> if 1 mode merge
        #    cat("Merging based on mathematical modes -> if 1 mode merge\n")
            bb=bb.list[[q]]

            points(bb[c(1,nrow(bb)),2:3],pch="*",col="red",cex=2)
            lines(bb$ridge.1,bb$ridge.2,col=2,lwd=2)
          }
          }



      }
  }
 }

  ratio.comat=0+(ratiomatrix>ratiocutoff)
  merged.clusters.ratio=con.comp(ratio.comat)


    # WE DO THE FINAL GROUPING HERE


      merged.clusters <- con.comp(comat)
   merged.group.ratio=original.group



  original.group=mixsmsn.out$group
  merged.group=original.group
  if (max(original.group)!=max(merged.clusters)){
  for(i in 1:max(original.group))
    merged.group[original.group == i] = merged.clusters[i]
}



  if (max(original.group)!=max(merged.clusters.ratio)){

  for(i in 1:max(original.group))
    merged.group.ratio[original.group == i] = merged.clusters.ratio[i]
  }




  out <- list(merged.clusters = merged.clusters, connection.matrix = comat,ratiomatrix=round(ratiomatrix,3),ratio.comat=ratio.comat,merged.clusters.ratio=merged.clusters.ratio,merged.group=merged.group,merged.group.ratio=merged.group.ratio)
  out
}
