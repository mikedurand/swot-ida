function [D,Obs,AllObs,DAll,Truth,Lats]=SelObs(DAll,Obs,Exp,AllTruth,AllLats)

AllObs=Obs;

% for i=1:Exp.Est_nt,
%     iEst(i)=find(abs(DAll.t-Exp.tUse(i))<(1*Exp.tStep/2));
% end

iEst=find(DAll.t>=Exp.tUse(1)&DAll.t<=Exp.tUse(end));

Obs.h=AllObs.h(:,iEst);
Obs.h0=Obs.h(:,1);
Obs.S=AllObs.S(:,iEst);
Obs.w=AllObs.w(:,iEst);

D=DAll;
D.nt=length(iEst);
%D.nt=Exp.Est_nt;
%D.t=Exp.tUse;
D.t=DAll.t(iEst);
%reshape new data -- this is copied/modified from ReadObs.m
Obs.hv=reshape(Obs.h',D.nR*D.nt,1);
Obs.wv=reshape(Obs.w',D.nR*D.nt,1);
Obs.Sv=reshape(Obs.S',D.nR*D.nt,1);

%calculate new dt -- this is copied/modified from ReadObs.m
D.dt=reshape( (diff(D.t)'.*86400*ones(1,D.nR)),D.nR*(D.nt-1),1);

Truth=AllTruth;
Truth.Q=AllTruth.Q(:,iEst);
Truth.dA=AllTruth.dA(:,iEst);
Truth.h=AllTruth.h(:,iEst);
Truth.W=AllTruth.W(:,iEst);

%these lines snagged from ReadTruth
Truth.dAv=reshape(Truth.dA',D.nR*D.nt,1);
Truth.hv=reshape(Truth.h',D.nR*D.nt,1);
Truth.Wv=reshape(Truth.W',D.nR*D.nt,1);

Lats.q=AllLats.q(:,iEst(1:end-1));
Lats.qv=reshape(Lats.q',D.nR*(D.nt-1),1);

return
