# 'creates test statistic value and critical value
#' @export
#' @param x data matrix
#' @param alpha numeric variable
NBTVMF3=function(x,alpha)
{
  library(blockmatrix)
  p=3
  k=3
  sequence_of_n=c(nrow(x[[1]]),nrow(x[[2]]),nrow(x[[3]]))
  y=vector(mode = "list", length = k)
  k_hat_0=vector(mode = "list", length = k)
  k_hat_hat_0=vector(mode = "list", length = k)
  g=vector(mode = "list", length = k)
  z=vector(mode = "list", length = k)
  r_bar=array(0,k)
  g_bar=array(0,k)
  k_hat_0=vector(mode = "list", length = k)
  D1=diag(1,nrow = (p),ncol = (p))
  O=diag(0,nrow = (p),ncol = (p))
  l1=list(D1=D1,D2=-D1,D3=O)
  m=array(c("D1","D3","D2","D1","D3","D2"),c((k-1),k))
  C=blockmatrix(value=m,list=l1)
  C1=as.matrix(C)

  H=1000 # no of inner bootstrap simulation
  inner_loop_statistic_value_array=array(0,H)
  for(i in 1:k)
  {

    y[[i]]=c(sum(x[[i]][,1]),sum(x[[i]][,2]),sum(x[[i]][,3]))

    r_bar[i]=sqrt(sum(y[[i]]**2))/sequence_of_n[i]

    k_hat_0[[i]]=y[[i]]/sqrt(sum(y[[i]]**2))

  }
  X=matrix(c(unlist(k_hat_0)),nrow = k*(p),ncol = 1,byrow = FALSE)
  X1=C1%*%(X)
  N=sum(sequence_of_n)    #calculations of Q^
  Q=sqrt(N)*X1
  statistic_value<-sum(abs(Q)**2)


  for(h in 1:H) #inner loop
  {
    for(i in 1:k)
    {

      g[[i]] <- x[[i]][sample(nrow(x[[i]]),replace = TRUE),]

      z[[i]]=c(sum(g[[i]][,1]),sum(g[[i]][,2]),sum(g[[i]][,3]))

      g_bar[i]=(sqrt(sum(z[[i]]**2)))

      k_hat_hat_0[[i]]=z[[i]]/g_bar[i]
    }
    X2=matrix(c(unlist(k_hat_hat_0)),nrow = k*(p),ncol = 1,byrow = FALSE)
    X3=C1%*%(X2-X)
    N=sum(sequence_of_n)
    Q1=sqrt(N)*X3

    inner_loop_statistic_value=sum(abs(Q1)**2)
    inner_loop_statistic_value_array[h]<-inner_loop_statistic_value
  }

  critical_value=quantile(inner_loop_statistic_value_array,1-alpha)
  return(list(statistic_value,critical_value))

}
