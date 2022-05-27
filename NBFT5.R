# 'creates test statistic value and critical value
#' @export
#' @param x data matrix
#' @param alpha numeric variable
##############initialization#####################
NBFT5=function(x,alpha)
{

  ##############initialization#####################

  rm(list = ls())

  sequence_of_n=c(nrow(x[[1]]),nrow(x[[2]]),nrow(x[[3]]),nrow(x[[4]]),nrow(x[[5]]))

  k=5 # 3 population anova
  p=3 # 3 dimensional case(sphere)

  mu_hat_0=vector(mode = "list", length = k)
  mu_hat_hat_0=vector(mode = "list", length = k)
  k_hat_0=vector(mode = "list", length = k) #initializing X^
  k_hat_hat_0=vector(mode = "list", length = k) #initializing X*
  H=1000 # no of inner bootstrap simulation

  alpha=0.05
  N=sum(sequence_of_n)

  inner_loop_statistic_value_array=array(0,H) # initializing inner loop array

  g=vector(mode = "list", length = k)
  y=vector(mode = "list", length = k)
  z=vector(mode = "list", length = k)

  r_bar=array(0,k)
  r_bar1=array(0,k)
  q_bar=array(0,k)
  q_bar1=array(0,k)

  l=array(0,k)
  l1=array(0,k)

  n11=vector(mode = "list", length = k)
  n1=vector(mode = "list", length = k)
  n1_hat=vector(mode = "list", length = k)
  gamma1=array(0,k)

  Q1=vector(mode = "list", length = k)

  A=vector(mode = "list", length = k)
  ###############################################################

  for(i in 1:k)
  {


    y[[i]]=c(sum(x[[i]][,1]),sum(x[[i]][,2]),sum(x[[i]][,3]))

    r_bar1[i]=sqrt(sum(y[[i]]**2))




    r_bar[i]=sqrt(sum(y[[i]]**2))/sequence_of_n[i]

    mu_hat_0[[i]]=y[[i]]/r_bar1[i]

    l[i]=(r_bar[i]*(p-r_bar[i]^2))/(1-r_bar[i]^2)

    k_hat_0[[i]]=matrix(c(mu_hat_0[[i]]),nrow = p,ncol = 1)




  }



  r=sqrt(sum((l[1]*y[[1]]+l[2]*y[[2]]+l[3]*y[[3]]+l[4]*y[[4]]+l[5]*y[[5]])^2))




  num=(sum(l*r_bar1)-r)/((k-1)*(p-1))
  den=(sum(sequence_of_n*l)-sum(l*r_bar1))/((N-k)*(p-1))


  statistic_value<-num/den


  mu1=l[1]*r_bar1[1]*mu_hat_0[[1]]+l[2]*r_bar1[2]*mu_hat_0[[2]]+l[3]*r_bar1[3]*mu_hat_0[[3]]+l[4]*r_bar1[4]*mu_hat_0[[4]]+l[5]*r_bar1[5]*mu_hat_0[[5]]
  mu_hat=mu1/sqrt(sum(mu1**2))
  mu_hat1=matrix(c(mu_hat),nrow = p,ncol = 1)


  M1= matrix(0,nrow = sequence_of_n[1],ncol = p)
  M2= matrix(0,nrow = sequence_of_n[2],ncol = p)
  M3= matrix(0,nrow = sequence_of_n[3],ncol = p)
  M4= matrix(0,nrow = sequence_of_n[4],ncol = p)
  M5= matrix(0,nrow = sequence_of_n[5],ncol = p)

  M=list(M1,M2,M3,M4,M5)

  for(i in 1:k)
  {
    n11[[i]]=(t(mu_hat1)%*%k_hat_0[[i]])
    gamma1[i]=acos(c(n11[[i]]))

    n1[[i]]=k_hat_0[[i]]-mu_hat1*c(n11[[i]])
    n1_hat[[i]]=n1[[i]]/sqrt(sum(n1[[i]]**2))

    A[[i]]=mu_hat1%*%t(n1_hat[[i]])-n1_hat[[i]]%*%t(mu_hat1)
    Q1[[i]]=diag(1,p)+sin(gamma1[i])*A[[i]]+(cos(gamma1[i])-1)*(mu_hat1%*%t(mu_hat1)+n1_hat[[i]]%*%t(n1_hat[[i]]))
    for(j in 1:sequence_of_n[i])
    {
      M[[i]][j,]=t(Q1[[i]]%*%matrix(c(x[[i]][j,]),nrow = p,ncol = 1))
    }
  }

  ############################################




  for(h in 1:H) #inner loop
  {
    for(i in 1:k)
    {

      g[[i]]=M[[i]][sample(nrow(M[[i]]),replace = TRUE),]
      z[[i]]=c(sum(g[[i]][,1]),sum(g[[i]][,2]),sum(g[[i]][,3]))

      q_bar[i]=sqrt(sum(z[[i]]**2))/sequence_of_n[i]
      q_bar1[i]=sqrt(sum(z[[i]]**2))

      mu_hat_hat_0[[i]]=z[[i]]/q_bar1[i]


      l1[i]=(q_bar[i]*(p-q_bar[i]^2))/(1-q_bar[i]^2)
    }

    r1=sqrt(sum((l1[1]*z[[1]]+l1[2]*z[[2]]+l1[3]*z[[3]]+l1[4]*z[[4]]+l1[5]*z[[5]])^2))



    num1=(sum(l1*q_bar1)-r1)/((k-1)*(p-1))
    den1=(sum(sequence_of_n*l1)-sum(l1*q_bar1))/((N-k)*(p-1))


    test_statistic=num1/den1

    inner_loop_statistic_value=test_statistic
    inner_loop_statistic_value_array[h]<-inner_loop_statistic_value
  }
  critical_value=quantile(inner_loop_statistic_value_array,1-alpha)
  return(list(statistic_value,critical_value))
}



