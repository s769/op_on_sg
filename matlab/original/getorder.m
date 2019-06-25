function indices = getorder()
% helper function to reorder function value points from polydata_20.txt
%
% Calls functions:
% junctionindex

fid=fopen('indexorder','r');

for j=1:3282
   %v =  textscan(fid, '%c');
   v = fscanf(fid,'%c',1);
   r=[];
   while v~= ']'
       v=fscanf(fid,'%c',1);
       r=[r str2num(v)+1];
   end
   indices(j) = junctionindex(r,8);
   v=fscanf(fid,'%c',2);
    
end

fclose(fid);