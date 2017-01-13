####### F U N C T I O N S :   Data --> Linear Model

WspKierK = function(L,C,P) {
  dxCL=L[1]-C[1]; dxCP=P[1]-C[1]
  dyCL=L[2]-C[2]; dyCP=P[2]-C[2]
  dCL2=dxCL^2+dyCL^2; dCP2=dxCP^2+dyCP^2
  XL=dxCL/dCL2; YL=dyCL/dCL2
  XP=dxCP/dCP2; YP=dyCP/dCP2
  sinPminL=dxCL*dyCP-dxCP*dyCL
  cosPminL=dxCL*dxCP+dyCL*dyCP
  betaPrz=atan(sinPminL/cosPminL)           # beta w 1. cwiartce
  if (cosPminL<0) betaPrz=betaPrz+pi        # beta w 2. lub 3. cwiartce
  else if (sinPminL<0) betaPrz=betaPrz+2*pi # beta w 4. cwiartce
  c(YL,-XL,YP-YL,XL-XP,-YP,XP,betaPrz*200/pi)

}

WspKierA = function(C,P) {
  Kat=WspKierK(C+c(100,0),C,P);
  XP=Kat[6]; YP=-Kat[5]
  c(YP,-XP,-YP,XP,Kat[7])
}

WspKierZ = function(C,P,Zero) {
  Kat1=WspKierK(C+c(100,0),C,P);
  beta1=Kat1[7]; XP=Kat1[6]; YP=-Kat1[5]
  beta0=WspKierK(C+c(100,0),C,Zero)[7]
  if(beta1<beta0) beta1=beta1+400
  c(YP,-XP,-YP,XP,beta1)

}

WspKierD = function(J,K) {
  dxJK=K[1]-J[1]; dyJK=K[2]-J[2]
  dJK=sqrt(dxJK^2+dyJK^2)
  cosJK=dxJK/dJK; sinJK=dyJK/dJK
  c(-cosJK,-sinJK,cosJK,sinJK,dJK)
}

helmert = function(P, method="plain"){
  k=nrow(P); P_rot=cbind(-P[,2],P[,1])
  Helm = cbind( rep(1:0,k), rep(0:1,k), as.vector(t(P_rot)), as.vector(t(P)) )
  #Helm0=rbind(matrix(0,5,3),Helm); Helm0[1:5,3]=-1
  HelmNorm=Helm;
  HelmNorm[,3]=as.vector(t(P_rot)-apply(P_rot,2,mean))
  HelmNorm[,4]=as.vector(t(P)-apply(P,2,mean))
  HelmNorm=apply(HelmNorm,2,function(x) x/sqrt(sum(x^2)))
  if(method == "normalize") Helm=HelmNorm
  if(method == "nuisance") Helm=Helm0
  Helm
}

rowsMatCompl = function(A) {
  QR=qr(t(A)); r=QR$rank
  qr.Q(QR)[,(r+1):min(ncol(A),nrow(A))]
}


setModel = function( D, P0 ) {
  n=nrow(D);  pOrient=length(unique(D[D[,3]==-1,1]))
  pPts=length( union( unique(as.vector(D[D[,3]<1,1:2])),
                      unique(as.vector(D[D[,3]>0,1:3])) ) )
  p=2*pPts+pOrient;
  A=matrix(0,n,p); name0=rep(NA,n); p1=0
  for(i in 1:n)
    if(D[i,3]==-1)
      if(D[i,4]==0) {name0[i]=D[i,2]; p1=p1+1; A[i,p1]=-1
      }else {name0[i]=name0[i-1]; A[i,p1]=-1}

  dObs=rep(NA,n)
  for(i in 1:n){
    if(D[i,3]==-1){## directions
      j=which(P0[,1]==D[i,1]); k=which(P0[,1]==D[i,2]); z=which(P0[,1]==name0[i])
      rob = WspKierZ(P0[j,2:3], P0[k,2:3], P0[z,2:3])
      A[i, pOrient + c(2*j-1,2*j,2*k-1,2*k) ] = 200/pi*rob[-5]
      d0=D[i,4]-rob[5]; dObs[i]=ifelse(d0<0,d0,d0-400)
    }
    if(D[i,3] > 0){## angles
      jC=which(P0[,1]==D[i,1]); jL=which(P0[,1]==D[i,2]); jP=which(P0[,1]==D[i,3])
      rob = WspKierK(P0[jL,2:3], P0[jC,2:3], P0[jP,2:3])
      A[i, pOrient + c(2*jL-1,2*jL,2*jC-1,2*jC,2*jP-1,2*jP) ] = 200/pi*rob[-7]
      dObs[i]=D[i,4]-rob[7]
    }
    if(D[i,3]==0){## distances
      j=which(P0[,1]==D[i,1]); k=which(P0[,1]==D[i,2])
      rob = WspKierD(P0[j,2:3], P0[k,2:3])
      A[i, pOrient + c(2*j-1,2*j,2*k-1,2*k) ] = rob[-5]
      dObs[i]=D[i,4]-rob[5]
    }
  }
  list( modMat=A, delObs=dObs )
}

setCommonData = function(D1,D2){
  n1=nrow(D1); n2=nrow(D2); ind1=rep(0,n1)
  for(i1 in 1:n1)
    for(i2 in 1:n2)
      if( D1[i1,1]==D2[i2,1] & D1[i1,2]==D2[i2,2] & D1[i1,3]==D2[i2,3] ) ind1[i1]=i2
      ind1a=which(ind1>0); ind2a=ind1[ind1>0]
      cbind(D1[ind1a,1:3], D2[ind2a,4]-D1[ind1a,4], D1[ind1a,4], D2[ind2a,4])
}

setCommonModel = function( D, P0 ) {
  n=nrow(D);  pOrient=length(unique(D[D[,3]==-1,1]))
  pPts=length( union( unique(as.vector(D[D[,3]<=0,1:2])),
                      unique(as.vector(D[D[,3]> 0,1:3])) ) )
  p=2*pPts+pOrient;
  A=matrix(0,n,p); p1=0; station=0
  for(i in 1:n)
    if(D[i,3]==-1)
      if(station!=D[i,1]){station=D[i,1]; p1=p1+1; A[i,p1]=-1} else A[i,p1]=-1

  for(i in 1:n){
    if(D[i,3]==-1){## directions
      j=which(P0[,1]==D[i,1]); k=which(P0[,1]==D[i,2]);
      rob = WspKierA(P0[j,2:3], P0[k,2:3])
      A[i, pOrient + c(2*j-1,2*j,2*k-1,2*k) ] = 200/pi*rob[-5]
    }
    if(D[i,3] > 0){## angles
      jC=which(P0[,1]==D[i,1]); jL=which(P0[,1]==D[i,2]); jP=which(P0[,1]==D[i,3])
      rob = WspKierK(P0[jL,2:3], P0[jC,2:3], P0[jP,2:3])
      A[i, pOrient + c(2*jL-1,2*jL,2*jC-1,2*jC,2*jP-1,2*jP) ] = 200/pi*rob[-7]
    }
    if(D[i,3]==0){## distances
      j=which(P0[,1]==D[i,1]); k=which(P0[,1]==D[i,2])
      rob = WspKierD(P0[j,2:3], P0[k,2:3])
      A[i, pOrient + c(2*j-1,2*j,2*k-1,2*k) ] = rob[-5]
    }
  }
  list( modMat=A, delObs=D[,4])
}

numOrient = function(D) length(unique(D[D[,3]==-1,1]))

####### F U N C T I O N S :   Fitting Linear Model

restGenInv.2 = function(A,ref){
  A0=A[,-ref]; qr0=qr(A0); Q0=qr.Q(qr0)
  A1=A[, ref]; qr1=qr(A1); Q1=qr.Q(qr1)
  A_ref = matrix(NA,ncol(A),nrow(A))
  G = t(Q0) %*% Q1; Q12 = t(Q1) - t(G)%*%t(Q0)
  A_ref[-ref,] = backsolve( qr.R(qr0), t(Q0) - G %*% Q12 )
  A_ref[ ref,] = backsolve( qr.R(qr1), Q12 )
  A_ref
}

restGenInv.3 = function(A,ref){
  A0=A[,-ref]; qr0=qr(A0); Q0=qr.Q(qr0)
  A1=A[, ref]
  A_ref = matrix(NA,ncol(A),nrow(A))
  A11_ = ginv( A1 - Q0 %*% (t(Q0)%*%A1) )
  A_ref[-ref,] = backsolve( qr.R(qr0), t(Q0) - (t(Q0)%*%A1)%*%A11_ )
  A_ref[ ref,] = A11_
  A_ref
}

restGenInvTest = function(A,ref){
  A0=A[,-ref]; qr0=qr(A0); Q0=qr.Q(qr0)
  A1=A[, ref]; qr1=qr(A1); Q1=qr.Q(qr1)
  B1 = ginv( A1 - Q0 %*% (t(Q0)%*%A1) )
  B2 = backsolve( qr.R(qr1), t(Q1) - (t(Q1)%*%Q0)%*%t(Q0) )
  all.equal(B1,B2)
}

free.lm = function(y,X){
  S=svd(X);  tol = sqrt(.Machine$double.eps); nz = S$d > tol * S$d[1]
  U=S$u[,nz]; d=S$d[nz]; V=S$v[,nz]
  hi = apply(U^2,1,sum)
  yy = (U %*% (t(U) %*% y))
  X_ = t(t(V)/d) %*% t(U); beta= as.vector(X_ %*% y)
  list(fitted = yy, coefficients = beta, hatvalues = hi)
}


rest.free.lm = function(y,A,ref){
  Ao=A[,-ref]; Ar=A[,ref]; Q=qr.Q(qr(Ao));
  AAr = Ar - Q %*% (t(Q) %*% Ar);
  yr  = y - Q %*% (t(Q) %*% y)
  xr=free.lm(yr,AAr)$coeff; xo=lm(y-Ar%*%xr ~ Ao +0)$coef
  x=1:ncol(A); x[ref]=xr; x[-ref]=xo
  yy=A%*%x
  S=svd(A);  tol = sqrt(.Machine$double.eps); rA = S$d > tol * S$d[1]
  hat = apply(S$u[,rA]^2,1,sum); hat=ifelse(hat>=1,hat-tol,hat)

  list( coeff = x, fitted=yy, hatvalues = hat, rank=sum(rA) )
}

geolm = function(Dane){## Adjustment of Control Network ~~ Fitting Linear Model
  D=Dane$D; sObs=Dane$sObs; P0=Dane$P0; ref=Dane$ref
  Mod=setModel(D,P0); pOrient=numOrient(D)

  A0=Mod$modMat; n=nrow(A0); p=ncol(A0); delObs=Mod$delObs
  A=A0; A[,(pOrient+1):p]=A[,(pOrient+1):p]/sObs; l=delObs/sObs

  RFLM = rest.free.lm(l,A,ref);
  ll=RFLM$fitted; v = as.vector(l-ll);
  lAdj = D[,4] - delObs + ll*sObs; lAdj = ifelse(lAdj>0,lAdj,lAdj+400)
  defect=p-RFLM$rank; s0 = sqrt(sum(v^2)/(n-p+defect))
  hat=RFLM$hatval; reliab=sqrt(1-hat); vStand=v/(s0*reliab)
  OutObs=cbind(D[,1:4],lAdj,hat,vStand); #rownames(OutObs)=1:n
  colnames(OutObs)=c("C","L","P","lObs","lAdj","hat","vStand")

  del=RFLM$coeff; delp = del[(pOrient+1):p]
  Delp0 = matrix(delp,ncol=2,byrow=TRUE);
  Delp=10^3*Delp0; Delp = cbind( Delp, sqrt(apply(Delp^2,1,sum)) )
  colnames(Delp)=c("dX","dY","dP"); #rownames(Delp)=P0[,1]
  P1 = cbind(P0[,1], P0[,2:3] + Delp0); colnames(P1)=c("Label","x","y")

  orient=1:pOrient
  Apts=A[,-orient]; Q=qr.Q(qr(A[,orient])); Apts = Apts - Q %*% (t(Q) %*% Apts)
  C=matrix(0,p-pOrient,defect); C[ref-pOrient,]=helmert(P0[,2:3])[ref-pOrient,-4]
  Ac=rbind(Apts,t(C)); QR.Ac=qr(Ac); Qa=qr.Q(QR.Ac)[1:n,]
  A_c = backsolve( qr.R(QR.Ac), t(Qa));

  s_delp=s0*sqrt(apply(A_c^2,1,sum))
  Sdelp=10^3*cbind( matrix(s_delp,ncol=2,byrow=TRUE))
  Sdelp = cbind( Sdelp, sqrt(apply(Sdelp^2,1,sum)) )
  colnames(Sdelp)=c("sX","sY","sP")
  OutP=cbind(P1, Delp, Sdelp); #rownames(OutP)=1:nrow(P0)

  list(delp=delp, Ainv=A_c, sigma=s0, obs=OutObs, points=OutP)
}

####### F U N C T I O N S :   Robust Deformation Analysis

robSimilarity = function(dX0p,SimilP,ref,iter){
  dX0=dX0p[ref]; Simil=SimilP[ref,]
  c=1E-7; c2=c^2; w05=rep(1,nrow(Simil))
  for(i in 1:iter){
    QR=qr(Simil*w05); Q=qr.Q(QR)
    z =  t(Q) %*% as.vector(dX0*w05)
    dX = as.vector(dX0 - (Q %*% z)/w05)
    w05 = 1 / sqrt( c + sqrt( rep( apply(matrix(dX^2,ncol=2,byrow=TRUE),1,sum), each=2 ) ) )
    #w05 = 1 / sqrt( c + rep( apply(matrix(dX^2,ncol=2,byrow=TRUE),1,sum), each=2 ) )
    #w05 = 1 / sqrt(sqrt(dX^2+c2))
    #w05 = 1 /  sqrt(abs(dX)+c)
  }
  QR=qr(Simil*w05); Qw=qr.Q(QR)*w05; R=qr.R(QR); M = backsolve(R,t(Qw))
  p=nrow(SimilP); B=matrix(0,p,p); B[,ref] = SimilP %*% M
  ## C=matrix(0,p,ncol(Simil)); C[ref,]=Simil*(w05^2)
  ## B2 = SimilP %*% solve(t(C)%*%SimilP,t(C))
  ##print(all.equal(B,B2))

  list( displ = dX0p - B%*%as.vector(dX0p), trans = diag(p)-B )
}

redod = function(dl0,N1,iter){
  c=1E-15; c2=c^2; w05=rep(1,nrow(N1));
  for(i in 1:iter){
    QR=qr(N1/w05); Q=qr.Q(QR); R=qr.R(QR)
    dl = Q %*% forwardsolve(t(R),dl0)
    dl = as.vector(dl/w05)
    w05 = 1 / sqrt( sqrt( c + rep( apply(matrix(dl^2,ncol=2,byrow=TRUE),1,sum), each=2 ) ) )
    #w05 = 1 / sqrt(sqrt(dl^2+c2))
    #w05 = 1 / sqrt(abs(dl)+c)
    #w05 = 1 / sqrt( (1+sqrt(dl^2+c2))*sqrt(dl^2+c2) ) ### kara log(1+sqrt(dl^2+c2))
  }
  dl
}

rda = function(delp,Ainv, P,refP,method="par4"){
  ref = as.vector(t(cbind(2*refP-1,2*refP)))
  SimP = helmert(P[,2:3]); if(method == "par3") SimP=SimP[,-4]

  RobSim= robSimilarity(delp,SimP,ref,iter=20)

  Displ0= matrix(1000*delp,ncol=2,byrow=TRUE)
  displ = RobSim$displ
  Displ = matrix(1000*displ,ncol=2,byrow=TRUE)
  Displ = cbind( Displ, sqrt(apply(Displ^2,1,sum)) )
  Binv = RobSim$trans %*% Ainv
  nP=nrow(P); qua=rep(0,nP); qua0=qua; qua2=qua
  for(i in 1:nP){
    pts=(2*i-1):(2*i);
    qua0[i]=t(delp[pts])  %*% solve(Ainv[pts,] %*% t(Ainv[pts,]), delp[pts])
    qua[i] =t(displ[pts]) %*% solve(Ainv[pts,] %*% t(Ainv[pts,]), displ[pts])
    qua2[i]=t(displ[pts]) %*% solve(Binv[pts,] %*% t(Binv[pts,]), displ[pts])
  }
  Out=cbind(P[,1], P[,2:3]+Displ[,1:2], Displ0, Displ, qua0, qua, qua2);
  colnames(Out)=c(" label","  XRob","  YRob","    dX0","    dY0","  dXRob","  dYRob",
                  "  dPRob","   test0","    test"," testRob")
  Out
}

