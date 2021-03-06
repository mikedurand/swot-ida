function [Delta,DeltaA,B,C,thetauA0,thetauna,thetaux1,thetauq,R] = ...
               InitializeMetropolis (D,C,P,R)

Delta = CalcDelta(D.nR,D.nt,D.L);
DeltaA = CalcADelta(D.nR,D.nt);
B = CalcB(D.nR,D.nt);
              
%3.1) allocations, and initial state set
C.thetaA0=nan(D.nR,C.N); %this is just A0
C.thetaA0(:,1)=P.meanA0;
thetauA0=C.thetaA0(:,1);

C.thetana=nan(D.nR,C.N);
C.thetana(:,1)=P.meanna;
thetauna=C.thetana(:,1);

C.thetax1=nan(D.nR,C.N);
C.thetax1(:,1)=P.meanx1;
thetaux1=C.thetax1(:,1);

% C.thetaq=nan(D.nR*(D.nt-1),C.N);
% C.thetaq(:,1)=P.meanq;
thetauq=nan; 

% C.thetaQbar=nan(1,C.N);
% C.thetaQbar(:,1)=P.meanQbar*ones(1,1);
% thetauQbar=C.thetaQbar(:,1);

%Get Random numbers
rng(R.Seed)
% R.z1=randn(1,C.N);  %spatially invariant increments... not good
R.z1=randn(D.nR,C.N);  
R.z2=randn(D.nR,C.N);
R.z3=randn(D.nR,C.N);
% R.z3=randn(D.nR*(D.nt-1),C.N);

R.u1=rand(C.N,1); %used for acceptance of A0
R.u2=rand(C.N,1); %used for acceptance of na
R.u3=rand(C.N,1); %used for acceptance of x1

return