% this script runs a single experiment, taking Run Directory as input

function RunExp(RunDir,ShowFigs,Laterals)

% RunDir='.';

disp(['Running ' RunDir])

ObsFile=[RunDir '/SWOTobs.txt'];
[DAll,Obs] = ReadObs(ObsFile);

ParamFile=[RunDir '/params.txt'];
[Chain,Prior,R,Exp] = ReadParams(ParamFile);

TruthFile=[RunDir '/truth.txt'];
AllTruth=ReadTruth (TruthFile,DAll);

if Laterals.UseMean
    LateralMeanFile=[RunDir '/LateralsMean.txt'];
    Prior.AllLats=ReadLats(LateralMeanFile,DAll);    
else
    Prior.AllLats=nan(size(Obs.h));
end

[D,Obs,AllObs,DAll,Truth,Prior.Lats]=SelObs(DAll,Obs,Exp,AllTruth,Prior.AllLats);

[Obs] = CalcdA(D,Obs);
[AllObs] = CalcdA(DAll,AllObs);

[Prior,jmp,AllObs]=ProcessPrior(Prior,AllObs,DAll,Obs,D,false,Exp.nOpt); 

[Obs,Prior] = GetCovMats(D,Obs,Prior);

%limit slopes to zero...
Obs.S(Obs.S<0)=0;
AllObs.S(AllObs.S<0)=0;

Chain=MetropolisCalculations(Prior,D,Obs,jmp,Chain,R,DAll,AllObs,Exp.nOpt);

[Estimate,Chain]=CalculateEstimates (Chain,D,Obs,Prior,DAll,AllObs,Exp.nOpt);

[Estimate] = FilterEstimate(Estimate,Chain,D,Obs);

Err=CalcErrorStats(AllTruth,Estimate,DAll);

Err=DispRMSEStats(Err,Truth,Prior,Estimate);

if ShowFigs,
    MakeFigs(D,Truth,Prior,Chain,Estimate,Err,AllTruth,DAll,AllObs);
end

% WriteSummary (R,Err,Estimate,RunDir);

% CheckBals(Truth,Obs,D,Prior,Chain,R,Estimate)

save([RunDir '/RunData.mat'])

return