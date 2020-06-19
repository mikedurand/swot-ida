function [Estimate,C] = CalculateEstimates (C,D,Obs,Prior,DAll,AllObs,BjerklienOpt)

% 1) estimates on the chain A0, n, q: means & covariances
Estimate.A0hat=mean(C.thetaA0(:,C.Nburn+1:end),2);
Estimate.CA0=cov(C.thetaA0(:,C.Nburn+1:end)');
Estimate.stdA0Post=sqrt(diag(Estimate.CA0));

Estimate.nahat=mean(C.thetana(:,C.Nburn+1:end),2);
Estimate.Cna=cov(C.thetana(:,C.Nburn+1:end)');
Estimate.stdnaPost=sqrt(diag(Estimate.Cna));

Estimate.x1hat=mean(C.thetax1(:,C.Nburn+1:end),2);
Estimate.Cx1=cov(C.thetax1(:,C.Nburn+1:end)');
Estimate.stdx1Post=sqrt(diag(Estimate.Cx1));

c1=0.85;
for r=1:D.nR,
    Estimate.nhat(r,:)=calcnhat(Obs.w(r,:),Obs.h(r,:),AllObs.hmin(r),Estimate.A0hat(r)+Obs.dA(r,:),...
        Prior.Wa(r),Prior.Ha(r),c1,Estimate.x1hat(r),Estimate.nahat(r),BjerklienOpt);
    
    Estimate.nhatAll(r,:)=calcnhat(AllObs.w(r,:),AllObs.h(r,:),AllObs.hmin(r),Estimate.A0hat(r)+AllObs.dA(r,:)-AllObs.A0Shift(r),...
    Prior.Wa(r),Prior.Ha(r),c1,Estimate.x1hat(r),Estimate.nahat(r),BjerklienOpt);    
end


% Estimate.qhat=mean(C.thetaq(:,C.Nburn+1:end),2);
% Estimate.Cq=cov(C.thetaq(:,C.Nburn+1:end)');
% Estimate.stdqpost=sqrt(diag(Estimate.Cq));

%2) calculate the Q chain, and estimate mean and std
for i=1:C.N,
    for r=1:D.nR,
        nhat(r,:)=calcnhat(Obs.w(r,:),Obs.h(r,:),AllObs.hmin(r),Estimate.A0hat(r)+Obs.dA(r,:), ...
            Prior.Wa(r),Prior.Ha(r),c1,...
            C.thetax1(r,i),C.thetana(r,i),BjerklienOpt);    
        nhatAll(r,:)=calcnhat(AllObs.w(r,:),AllObs.h(r,:),AllObs.hmin(r),Estimate.A0hat(r)+AllObs.dA(r,:)-AllObs.A0Shift(r), ...
            Prior.Wa(r),Prior.Ha(r),c1,...
            C.thetax1(r,i),C.thetana(r,i),BjerklienOpt);    
    end
    
    C.thetaQ(:,:,i) = 1./nhat .* ...
        (C.thetaA0(:,i)*ones(1,D.nt)+Obs.dA).^(5/3).*Obs.w.^(-2/3).*sqrt(Obs.S) ;
    C.thetaAllQ(:,:,i) = 1./nhatAll .* ...
        ( (C.thetaA0(:,i)-AllObs.A0Shift)*ones(1,DAll.nt)+AllObs.dA).^(5/3).*AllObs.w.^(-2/3).*sqrt(AllObs.S);
end

Estimate.QhatPost=mean(C.thetaQ(:,:,C.Nburn+1:end),3);
Estimate.QstdPost=std(C.thetaQ(:,:,C.Nburn+1:end),[],3);

%3) Calculate Q prior estimate
for r=1:D.nR,
    nhat(r,:)=calcnhat(Obs.w(r,:),Obs.h(r,:),AllObs.hmin(r),Prior.meanA0(r)*ones(1,D.nt)+Obs.dA(r,:),Prior.Wa(r),Prior.Ha(r),c1,...
        Prior.meanx1(r),Prior.meanna(r),BjerklienOpt);            
end

Estimate.QhatPrior=1./nhat .* ...
    (Prior.meanA0*ones(1,D.nt)+Obs.dA).^(5/3).*Obs.w.^(-2/3).*sqrt(Obs.S);

clear nhat
for r=1:D.nR,
    nhat(r,:)=calcnhat(AllObs.w(r,:),AllObs.h(r,:),AllObs.hmin(r),Prior.meanA0(r)-AllObs.A0Shift(r)+AllObs.dA(r,:),Prior.Wa(r),Prior.Ha(r),c1,...
        Prior.meanx1(r),Prior.meanna(r),BjerklienOpt);                
end

Estimate.QhatAllPrior=1./nhat .* ...
    ( (Prior.meanA0-AllObs.A0Shift) *ones(1,DAll.nt)+AllObs.dA).^(5/3).*AllObs.w.^(-2/3).*sqrt(AllObs.S);


%4)   Discharge error budget: all done for Q(nr x nt)
%4.1) Uncertainty estimate of the dA term
Obs.sigdAv=sqrt(diag(Obs.CdA));
Obs.sigdA=reshape(Obs.sigdAv,D.nt,D.nR)';
%4.2) Uncertainty (variance) of the Manning terms
% Estimate.QhatUnc.n=(Estimate.stdnPost./Estimate.nhat).^2;
Estimate.QhatUnc.na=(Estimate.stdnaPost./Prior.meanna).^2;
Estimate.QhatUnc.x1=(log(Obs.w./(Estimate.A0hat*ones(1,D.nt)+Obs.dA)) .* ...
                    (Estimate.stdx1Post*ones(1,D.nt))).^2;

Estimate.QhatUnc.w=(2/3*Obs.sigw./Obs.w).^2;
% Estimate.QhatUnc.A0=(5/3.*Estimate.stdA0Post*ones(1,D.nt)./(Estimate.A0hat*ones(1,D.nt)+Obs.dA)).^2; % change
% Estimate.QhatUnc.dA=(5/3.*Obs.sigdA./(Estimate.A0hat*ones(1,D.nt)+Obs.dA)).^2;
Estimate.QhatUnc.S=(1/2.*Obs.sigS./Obs.S).^2;

Estimate.QhatUnc.A0=((5/3-Prior.meanx1*ones(1,D.nt)) .* ...
                    (Estimate.stdA0Post*ones(1,D.nt)) ./ ...
                    (Estimate.A0hat*ones(1,D.nt)+Obs.dA)).^2;
Estimate.QhatUnc.dA=((5/3-Prior.meanx1*ones(1,D.nt)) .* ...
                    Obs.sigdA ./ ...
                    (Estimate.A0hat*ones(1,D.nt)+Obs.dA)).^2;
     
% Need to revise this for variable n... 

%4.3) Estimate correlation coefficient between A0 & na, A0 & x1, na & x1
 for i=1:D.nR,
     R_A0na=corrcoef([C.thetaA0(i,:); C.thetana(i,:);]');
     Estimate.rho_A0na(i,1)=R_A0na(1,2);
     R_A0x1=corrcoef([C.thetaA0(i,:); C.thetax1(i,:);]');
     Estimate.rho_A0x1(i,1)=R_A0x1(1,2);
     R_nax1=corrcoef([C.thetana(i,:); C.thetax1(i,:);]');
     Estimate.rho_nax1(i,1)=R_nax1(1,2);
 end

%4.4) Estimate uncertainty of Manning's Q
%4.4.1) Estimate uncertainty in Q due to cross-correlation of A0 & na, na & x1, x1 & A0 
Estimate.QhatUnc.A0na=-2*(Estimate.rho_A0na*ones(1,D.nt)) .* ...
                      (5/3-Prior.meanx1*ones(1,D.nt)) .* ...
                      (Estimate.stdnaPost*ones(1,D.nt)) .* ...
                      (Estimate.stdA0Post*ones(1,D.nt)) ./ ...
                      (Prior.meanna*ones(1,D.nt)) ./ ...
                      (Estimate.A0hat*ones(1,D.nt)+Obs.dA);
              
Estimate.QhatUnc.nax1=-2*(Estimate.rho_nax1*ones(1,D.nt)) .* ...
                      (Estimate.stdx1Post*ones(1,D.nt)) .* ...
                      (Estimate.stdnaPost*ones(1,D.nt)) .* ...
                      log(Obs.w./(Estimate.A0hat*ones(1,D.nt)+Obs.dA)) ./ ...
                      (Prior.meanna*ones(1,D.nt));
                  
Estimate.QhatUnc.A0x1=2*(Estimate.rho_A0x1*ones(1,D.nt)) .* ...
                      (Estimate.stdA0Post*ones(1,D.nt)) .* ...
                      (Estimate.stdx1Post*ones(1,D.nt)) .* ...
                      (5/3-Prior.meanx1*ones(1,D.nt)) .* ...
                      log(Obs.w./(Estimate.A0hat*ones(1,D.nt)+Obs.dA)) ./ ...
                      (Estimate.A0hat*ones(1,D.nt)+Obs.dA);

% %4.4.2) Estimate total Q uncertinty
% Estimate.QhatUnc.Hat=sqrt( Estimate.QhatUnc.n+mean(Estimate.QhatUnc.w,2)+...
%     mean(Estimate.QhatUnc.A0,2)+mean(Estimate.QhatUnc.dA,2)+mean(Estimate.QhatUnc.S,2)+...
%     mean(Estimate.QhatUnc.A0n,2)  );
% 
% %4.4.2) Estimate total Q uncertinty wrt time
% Estimate.QhatUnc.HatAll=sqrt( Estimate.QhatUnc.n*ones(1,D.nt) + Estimate.QhatUnc.w+...
%     Estimate.QhatUnc.A0+Estimate.QhatUnc.dA+Estimate.QhatUnc.S+...
%     Estimate.QhatUnc.A0n  );
% 
% %4.4.3) Discharge error budget
% Estimate.QerrVarSum(1:D.nR,1)=Estimate.QhatUnc.n;
% Estimate.QerrVarSum(1:D.nR,2)=mean(Estimate.QhatUnc.A0,2);
% Estimate.QerrVarSum(1:D.nR,3)=mean(Estimate.QhatUnc.dA,2);
% Estimate.QerrVarSum(1:D.nR,4)=mean(Estimate.QhatUnc.w,2);
% Estimate.QerrVarSum(1:D.nR,5)=mean(Estimate.QhatUnc.S,2);
% 
% %5 all estimates
% 

%4.4.2) Estimate total Q uncertinty
Estimate.QhatUnc.Hat=sqrt( Estimate.QhatUnc.na+mean(Estimate.QhatUnc.x1,2)+mean(Estimate.QhatUnc.w,2)+...
    mean(Estimate.QhatUnc.A0,2)+mean(Estimate.QhatUnc.dA,2)+mean(Estimate.QhatUnc.S,2)+...
    mean(Estimate.QhatUnc.A0na,2)+mean(Estimate.QhatUnc.nax1,2)+mean(Estimate.QhatUnc.A0x1,2));

%4.4.2) Estimate total Q uncertinty wrt time
Estimate.QhatUnc.HatAll=sqrt( Estimate.QhatUnc.na*ones(1,D.nt)+Estimate.QhatUnc.x1+Estimate.QhatUnc.w+...
    Estimate.QhatUnc.A0+Estimate.QhatUnc.dA+Estimate.QhatUnc.S+...
    Estimate.QhatUnc.A0na+Estimate.QhatUnc.nax1+Estimate.QhatUnc.A0x1);

%4.4.3) Discharge error budget
Estimate.QerrVarSum(1:D.nR,1)=Estimate.QhatUnc.na;
Estimate.QerrVarSum(1:D.nR,2)=mean(Estimate.QhatUnc.x1,2);
Estimate.QerrVarSum(1:D.nR,3)=mean(Estimate.QhatUnc.A0,2);
Estimate.QerrVarSum(1:D.nR,4)=mean(Estimate.QhatUnc.dA,2);
Estimate.QerrVarSum(1:D.nR,5)=mean(Estimate.QhatUnc.w,2);
Estimate.QerrVarSum(1:D.nR,6)=mean(Estimate.QhatUnc.S,2);

%5 all estimates

Estimate.AllQ = 1./Estimate.nhatAll .* ...
        ( (Estimate.A0hat-AllObs.A0Shift)*ones(1,DAll.nt)+AllObs.dA).^(5/3).*AllObs.w.^(-2/3).*sqrt(AllObs.S) ;


return
