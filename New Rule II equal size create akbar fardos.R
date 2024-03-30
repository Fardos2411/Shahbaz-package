
################################################################################
# R-coded Algorithm Based on Cyclic Shifts (Rule II) to Generate Efficient 
# Classes of Circular Balanced Repeated Measurements Designs  for period of 
# equal size(P)
################################################################################

# Algorithm from paper:

# Muhammad Shahbaz, Muhammad Riaz, H. M. kashif Rasheed, Abid Khan and Rashid Ahmed (2024). 
# Classes of Minimal Circular balance RMDs for period of equal 
# size(P)  CBRMDs2 -designs Generated with R 


# Coded by shahbaz et al., 01-08-2023 to 01-01-2024
# Version 2.0  (2024-01-02)
################################################################################




################################################################
# Division of adjusted A in i-1 groups of size p and one
# of size p-2 to get the set(s) of shifts
################################################################
grouping1<-function(A,p,v,i){
  bs<-c()
  z=0;f=1
  A1=A
  while(f<=(i-1)){
    
    for(y in 1:5000){
      comp<-sample(1:length(A1),p)
      com<-A1[comp]
      cs<-sum(com)
      if(cs%%v==0){
        bs<-rbind(bs,com)
        A1<-A1[-comp]
        z<-z+1
        f=f+1
      }
      if(z==(i-1)) break
    }
    if(z<(i-1)) {bs<-c();z=0;f=1;A1=A}  
    
  }
  
  
  bs1<-t(apply(bs,1,sort))
  bs1<-cbind(bs1,rowSums(bs),rowSums(bs)/(v))
  rownames(bs1)<-paste("G",1:(i-1), sep="")
  colnames(bs1)<-c(paste(1:p, sep=""),"sum" ,"sum/(v-1)")
  
  bs2<-t(apply(as.matrix(A1),1,sort))
  bs2<-cbind(bs2,rowSums(bs2),rowSums(bs2)/(v))
  rownames(bs2)<-paste("G",i, sep="")
  colnames(bs2)<-c(paste(1:(p-2), sep=""),"sum" ,"sum/(v-1)")
  
  
  fs1<-t(apply(bs,1,sort))
  fs1<-delmin(fs1)
  rownames(fs1)<-paste("S",1:(i-1), sep="")
  colnames(fs1)<-rep("",(p-1))
  
  fs2<-t(apply(as.matrix(A1),1,sort))
  rownames(fs2)<-paste("S",i, sep="")
  colnames(fs2)<-rep("",(p-2))
  
  
  list(B1=list(fs1,fs2),B2=list(bs1,bs2))
}


#######################################################################
# Obtain set(s) of shifts by deleting smallest value of group
#######################################################################

delmin<-function(z){
  fs<-c()
  n<-nrow(z)
  c<-ncol(z)-1
  for(i in 1:n){
    z1<-z[i,]
    z2<-z1[z1!=min(z1)]
    fs<-rbind(fs,z2)
  }
  return(fs)
}


################################################################################
# Selection of adjusted A and the set(s) of shifts to obtain Circular 
# BRMDs for period of equal size.
################################################################################

#   p: period sizes
#   i: Number of set of shifts for K


RMD_equalsize<-function(v,p,i,C=1){
  
  if(p<=2) stop("p= Block size: Block size must be greater than 2")
  if(i<=1) stop("i= Must be a greater than 1")
  
  setClass( "stat_test", representation("list"))
  
  setMethod("show", "stat_test", function(object) {
    
    row <- paste(rep("=", 50), collapse = "")
    cat(row, "\n")
    print(object$S[[1]])
    print(object$S[[2]])
  })
  
  if(v%%2==0){
  if(C==1){
    v=i*p; m=(v-2)
    A=c(1,2:m)
    A1<-grouping1(A,p,v=(v-1),i)
    A2<-c(v,p);names(A2)<-c("V","p")
    x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
  }
    if(C==2){
      v=i*p+1; m=(v-2)
      A<-c(1:((v-2)/2),((v+2)/2),((v+4)/2):m)
      A1<-grouping1(A,p,v=(v-1),i)
      A2<-c(v,p);names(A2)<-c("V","p")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
  if(C==3){
      v=i*p+2; m=(v-2)
      A<-c(1:((v-2)/2),((v+4)/2),((v+6)/2):m)
      A1<-grouping1(A,p,v=(v-1),i)
      A2<-c(v,p);names(A2)<-c("V","p")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
  }
    if(C==4){
      v=i*p-1; m=(v-2)
      A<-c(1:m,(v/2))
      A1<-grouping1(A,p,v=(v-1),i)
      A2<-c(v,p);names(A2)<-c("V","p")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
 }    
    if(C==5){
      v=i*p-2; m=(v-2)
      A<-c(1:m,(v/2),((v+2)/2))
      A1<-grouping1(A,p,v=(v-1),i)
      A2<-c(v,p);names(A2)<-c("V","p")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    
    if(C==6){
      v=i*p-1; m=(v-2)
      A=c(0,1,2:m)
      A1<-grouping1(A,p,v=(v-1),i)
      A2<-c(v,p);names(A2)<-c("V","p")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(C==7){
      v=i*p; m=(v-2)
      A<-c(0,1:((v-2)/2),((v+2)/2),((v+4)/2):m)
      A1<-grouping1(A,p,v=(v-1),i)
      A2<-c(v,p);names(A2)<-c("V","p")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(C==8){
      v=i*p+1; m=(v-2)
      A<-c(0,1:((v-2)/2),((v+4)/2),((v+6)/2):m)
      A1<-grouping1(A,p,v=(v-1),i)
      A2<-c(v,p);names(A2)<-c("V","p")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
    if(C==9){
      v=i*p-2; m=(v-2)
      A<-c(0,1:m,(v/2))
      A1<-grouping1(A,p,v=(v-1),i)
      A2<-c(v,p);names(A2)<-c("V","p")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }    
    if(C==10){
      v=i*p-3; m=(v-2)
      A<-c(0,1:m,(v/2),((v+2)/2))
      A1<-grouping1(A,p,v=(v-1),i)
      A2<-c(v,p);names(A2)<-c("V","p")
      x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
    }
  }
   if(v%%2!=0){
      if(C==1){
        v=i*p; m=(v-2)
        A=c(1,2:m)
        A1<-grouping1(A,p,v=(v-1),i)
        A2<-c(v,p);names(A2)<-c("V","p")
        x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
      }      
  if(C==2){
    v=i*p+1; m=(v-2)
    A<-c(1:((v-1)/2),((v+3)/2),((v+5)/2):m)
    A1<-grouping1(A,p,v=(v-1),i)
    A2<-c(v,p);names(A2)<-c("V","p")
    x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
  }
  
  if(C==3){
    v=i*p+2; m=(v-2)
    A<-c(1:((v-3)/2),((v+3)/2),((v+5)/2):m)
    A1<-grouping1(A,p,v=(v-1),i)
    A2<-c(v,p);names(A2)<-c("V","p")
    x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
  }
  
  if(C==4){
    v=i*p-1; m=(v-2)
    A<-c(1:m,((v+1)/2))
    A1<-grouping1(A,p,v=(v-1),i)
    A2<-c(v,p);names(A2)<-c("V","p")
    x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
  }
  
  if(C==5){
    v=i*p-2; m=(v-2)
    A<-c(1:m,((v-1)/2),((v+1)/2))
    A1<-grouping1(A,p,v=(v-1),i)
    A2<-c(v,p);names(A2)<-c("V","p")
    x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
  }
     
     if(C==6){
       v=i*p-1; m=(v-2)
       A=c(0,1,2:m)
       A1<-grouping1(A,p,v=(v-1),i)
       A2<-c(v,p);names(A2)<-c("V","p")
       x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
     }      
     if(C==7){
       v=i*p; m=(v-2)
       A<-c(0,1:((v-1)/2),((v+3)/2),((v+5)/2):m)
       A1<-grouping1(A,p,v=(v-1),i)
       A2<-c(v,p);names(A2)<-c("V","p")
       x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
     }
     
     if(C==8){
       v=i*p+1; m=(v-2)
       A<-c(0,1:((v-3)/2),((v+3)/2),((v+5)/2):m)
       A1<-grouping1(A,p,v=(v-1),i)
       A2<-c(v,p);names(A2)<-c("V","p")
       x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
     }
     
     if(C==9){
       v=i*p-2; m=(v-2)
       A<-c(0,1:m,((v+1)/2))
       A1<-grouping1(A,p,v=(v-1),i)
       A2<-c(v,p);names(A2)<-c("V","p")
       x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
     }
     
     if(C==10){
       v=i*p-3; m=(v-2)
       A<-c(0,1:m,((v-1)/2),((v+1)/2))
       A1<-grouping1(A,p,v=(v-1),i)
       A2<-c(v,p);names(A2)<-c("V","p")
       x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
     }
}
  
  new("stat_test", x)
  
}

##################################################################
# Generation of design using sets of cyclical shifts
###################################################################
# H is an output object from CBRMDs_equalsize
# The output is called using the design_CWBND to generate design

design<-function(H){
  
  setClass( "CRMD_design", representation("list"))
  setMethod("show", "CRMD_design", function(object) {
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    cat("Following is minimal CRMD for", "v=" ,object$R[1], "and","k=",object$R[2], "\n")
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    for(i in 1:length(ss)){
      W<-ss[[i]]
      nr<-dim(W)[1]
      for(j in 1:nr){
        print(object$Design[[i]][[j]])
        cat("\n\n")
      }}
  })  
  
  v<-H$R[1]
  p<-H$R[2]
  ss<-H$S  
  treat<-(1:v)-1
  fn<-(1:v)
  G<-list()
  
  
  for(j in 1:length(ss)){ 
    W<-ss[[j]]
    nr<-dim(W)[1]
    nc<-dim(W)[2]
    D<-list()
    
    for(i in 1:nr){
      dd<-c()
      d1<-matrix(treat,(nc+1),v,byrow = T)
      ss1<-cumsum(c(0,W[i,]))
      dd2<-d1+ss1
      dd<-rbind(dd,dd2)
      rr<-dd[which(dd>=v)]%%v
      dd[which(dd>=v)]<-rr
      colnames(dd)<-paste("B",fn, sep="")
      rownames(dd)<-rep("",(nc+1))
      fn<-fn+v
      D[[i]]<-dd
    }
    G[[j]]<-D
    
  }
  
  x<-list(Design=G,R=H$R)
  new("CRMD_design", x)
}


Create_Designs<-function(v,p,C,i){
  if((v)%%(p)==0 & (i=(v)/(p)) & (i%%1==0) &(i>0) ){
    
    cat("MCBRMDs is possible for","v=",v,"P=",p,"i=",i,rep("", 30))
    (H<-RMD_equalsize(v=v,p=p,i=i,C=1))
    
  }else
    if((v-1)%%(p)==0 & (i=(v-1)/(p))  & (i%%1==0) &(i>0)) {
      
      cat("MCPBRMDs-I is possible for","v=",v,"p=",p,"i=",i,rep("", 30))
      (H<-RMD_equalsize(v=v,p=p,i=i,C=2))        
    }else
      if((v-2)%%(p)==0 & (i=(v-2)/(p))  & (i%%1==0) &(i>0)){
        
        cat("MCPBRMDs-II is possible for","v=",v,"p=",p,"i=",i,rep("", 30))
        (H<-RMD_equalsize(v=v,p=p,i=i,C=3))        
        
      }else
        if((v+1)%%(p)==0 & (i=(v+1)/(p)) & (i%%1==0) &(i>0)) {
          
          cat("MCWBRMDs-I is possible for","v=",v,"p=",p,"i=",i,rep("", 30))
          (H<-RMD_equalsize(v=v,p=p,i=i,C=4))        
        }else
          if((v+2)%%(p)==0 & (i=(v+2)/(p))  & (i%%1==0) &(i>0)){
            
            cat("MCWBRMDs-II is possible for","v=",v,"P=",p,"i=",i,rep("", 30))
            (H<-RMD_equalsize(v=v,p=p,i=i,C=5))
            
          }else
            if((v+1)%%(p)==0 & (i=(v+1)/(p))  & (i%%1==0) &(i>0)) {
              
              cat("MCNSBRMDs is possible for","v=",v,"p=",p,"i=",i,rep("", 30))
              (H<-RMD_equalsize(v=v,p=p,i=i,C=6))        
            }else
              if((v)%%(p)==0 & (i=(v)/(p)) & (i%%1==0) &(i>0) ){
                
                cat("MCSPBRMDs-I is possible for","v=",v,"p=",p,"i=",i,rep("", 30))
                (H<-RMD_equalsize(v=v,p=p,i=i,C=7))        
                
              }else
                if((v-1)%%(p)==0 & (i=(v-1)/(p)) & (i%%1==0) &(i>0)) {
                  
                  cat("MCSPBRMDs-II is possible for","v=",v,"p=",p,"i=",i,rep("", 30))
                  (H<-RMD_equalsize(v=v,p=p,i=i,C=8))        
                }else
                  if((v+2)%%(p)==0 & (i=(v+2)/(p)) & (i%%1==0) &(i>0) ){
                    
                    cat("MCSBGRMDs-I is possible for","v=",v,"P=",p,"i=",i,rep("", 30))
                    (H<-RMD_equalsize(v=v,p=p,i=i,C=9))
                    
                  }else
                    if((v+3)%%(p)==0 & (i=(v+3)/(p)) & (i%%1==0) &(i>0) ) {
                      
                      cat("MCSBGRMDs-II is possible for","v=",v,"p=",p,"i=",i,rep("", 30))
                      (H<-RMD_equalsize(v=v,p=p,i=i,C=10))        
                    }else
                      cat("MCBRMDs, MCPBRMDs-I, MCPBRMDs-II, MCWBRMDs-I, MCWBRMDs-II,
MCNSBRMDs, MCSPBRMDs-I, MCSPBRMDs-II, MCSBGRMDs-I and MCSBGRMDs-II
are not possible for ","v=",v,"and", "p=",p,"i=",i,rep("", 30))
  
  
}















################################################################################
# Examples: Using CBRMDs2_equal size function to obtain the set(s) of shifts
# for construction of Circular BRMDs for equal period  
# sizes (p)
################################################################################

(H<-Create_Designs(v=24, p=8))

