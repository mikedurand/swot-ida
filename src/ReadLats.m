function [Lats] = ReadLats(fname,D)

fid=fopen(fname,'r');
fgetl(fid);

for i=1:D.nR %read laterals
    Lats.q(i,:)=fscanf(fid,'%f',D.nt-1); fscanf(fid,'\n');
end

fclose(fid);

Lats.qv=reshape(Lats.q',D.nR*(D.nt-1),1);

return