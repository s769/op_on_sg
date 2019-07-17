function T = readPolys(level, deg)
%
% reads the monomials in polydata_20.txt
%
% Calls the functions:
% getorder
% alternateindex

fid=fopen('polydata_trial64new.txt');
m=level;
T=zeros(3^(m+1),1,3);
indices=getorder();
num = ((3^(level + 1) + 3)/2);
for j=1:deg
    for i=1:3
        for k=1:num
            index = indices(k);
            tline=fgetl(fid);
            T(index,j,i) = str2double(tline);
        end
        
    end
end

for j=1:deg
    for i=1:3
        for k=1:3^(m+1)
            if T(k,j,i)==0
                v=alternateindex(k,m);
                T(k,j,i)=T(v,j,i);
            end
        end
    end
end

T(1,1,1)=1;
T(3^(m+1),1,1)=1;
%gaskplot(T(:,5,2),m)
fclose(fid);