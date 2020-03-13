function [Chain,Prior,R,Exp] = ReadParams(fname)

fid=fopen(fname,'r');
fgetl(fid); Chain.N=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); Chain.Nburn=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); R.Seed=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); tUse1=fscanf(fid,'%f',1); tUseEnd=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); Prior.meanQbar=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); Prior.covQbar=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); Exp.nOpt=fscanf(fid,'%f',1); fscanf(fid,'\n');
fgetl(fid); Exp.tStep=fscanf(fid,'%f',1); fscanf(fid,'\n'); %time step as fraction of day

if fgetl(fid)~=-1 %if extra line in param file for geomorph priors
    Prior.Geomorph.Use=true;
    Prior.Geomorph.loga_hat=fscanf(fid,'%f',1);   %variable names from Craig's spreadsheet
    Prior.Geomorph.loga_sigma=fscanf(fid,'%f',1); 
    Prior.Geomorph.b_hat=fscanf(fid,'%f',1); 
    Prior.Geomorph.b_sigma=fscanf(fid,'%f',1); 
    Prior.Geomorph.logA0_hat=fscanf(fid,'%f',1); 
    Prior.Geomorph.logA0_sigma=fscanf(fid,'%f',1); 
    fscanf(fid,'\n');
else
    Prior.Geomorph.Use=false;
end


Exp.tUse(1)=tUse1;
Exp.tUse(2)=tUseEnd;

fclose(fid);

return
