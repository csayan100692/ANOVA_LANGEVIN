# 'creates test statistic value and critical value
#' @export
#' @param x data matrix
#' @param alpha numeric variable
####################################################
PBTVMF3=function(x,alpha)
{
p=3
k=3

sequence_of_n=c(nrow(x[[1]]),nrow(x[[2]]),nrow(x[[3]]))
H=1000
y=vector(mode = "list", length = k)
k_hat_0=vector(mode = "list", length = k)
k_hat_hat_0=vector(mode = "list", length = k)
g=vector(mode = "list", length = k)
z=vector(mode = "list", length = k)
r_bar=array(0,k)
r_bar1=array(0,k)
g_bar=array(0,k)
g_bar1=array(0,k)
a=array(0,k)
b=array(0,k)
inner_loop_statistic_value_array=array(0,H)
####################################################################
  for(i in 1:k)
  {
    y[[i]]=c(sum(x[[i]][,1]),sum(x[[i]][,2]),sum(x[[i]][,3]))

    r_bar[i]=sqrt(sum(y[[i]]**2))/sequence_of_n[i]
    r_bar1[i]=sqrt(sum(y[[i]]**2))
    k_hat_0[[i]]=y[[i]]/sqrt(sum(y[[i]]**2))
    a[i]=r_bar[i]*(p-r_bar[i]^2)/(1-r_bar[i]^2)

  }
  y1=(a[1]*r_bar1[1]*k_hat_0[[1]]+a[2]*r_bar1[2]*k_hat_0[[2]]+a[3]*r_bar1[3]*k_hat_0[[3]])
  mu_hat=y1/sqrt(sum(y1**2))
  N=sum(sequence_of_n)
  A1=a[1]*sqrt(sum(y[[1]]**2))+a[2]*sqrt(sum(y[[2]]**2))+a[3]*sqrt(sum(y[[3]]**2))
  A2=a[1]*y[[1]]+a[2]*y[[2]]+a[3]*y[[3]]
  A3=sqrt(sum(A2**2))
  statistic_value=2*(A1-A3)
for(h in 1:H)
  {
    for(i in 1:k)
    {

    g[[i]] <- rmovMF(sequence_of_n[i],a[i]*mu_hat)

    z[[i]]=c(sum(g[[i]][,1]),sum(g[[i]][,2]),sum(g[[i]][,3]))
    g_bar[i]=sqrt(sum(z[[i]]**2))/sequence_of_n[i]
    g_bar1[i]=sqrt(sum(z[[i]]**2))
    k_hat_hat_0[[i]]=z[[i]]/sqrt(sum(z[[i]]**2))
    b[i]=g_bar[i]*(p-g_bar[i]^2)/(1-g_bar[i]^2)
    }

    B1=b[1]*g_bar1[1]+b[2]*g_bar1[2]+b[3]*g_bar1[3]
    B2=b[1]*z[[1]]+b[2]*z[[2]]+b[3]*z[[3]]
    B3=sqrt(sum(B2**2))
    inner_loop_statistic_value_array[h]=2*(B1-B3)
  }

  critical_value<-quantile(inner_loop_statistic_value_array,1-alpha)
  return(list(statistic_value,critical_value))
}

