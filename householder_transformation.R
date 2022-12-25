# Householder transformation convert a symmertic to a tridiagonal matrix
A<-matrix(c(1,-1,2,2,-1,2,1,-1,2,1,3,2,2,-1,2,1),ncol = 4,byrow = TRUE)
A<-matrix(c(2,-1,1,-1,3,2,1,2,1),ncol = 3,byrow = TRUE)

Householder_trasformation<-function(A){
n<-ncol(A)

for (i in 1:(n-2)){

s=0
for (j in (i+1):n){
  s=s+(A[j,i])^2
  
}
s<-sqrt(s)

if (A[i+1,i]<0){
  a=-1
}else {
  a=1
  print(a)
}

r=0.5*(1+a*A[i+1,i]/s)

w<-matrix(0,n,1)
for (k in 1:i){
  w[k]=0
}
w[i+1]=sqrt(r)

for (z in (i+2):n){
  w[z]<-a*A[i,z]/(2*s*w[i+1])
}


P=diag(1,n,n)-2*w%*%t(w)

A<-P%*%A%*%P
print(A)
}


}