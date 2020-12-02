function [largest,smallest,maxmols,minmols] = max_mols(clus1,clus2,ch_info)
% Input:
%    names of two clusters and info on charges
% Output:
%    the name of the larger and smaller cluster and
%    the corresponding total molecule numbers (not including
%    an extra or a removed proton)

clusters={clus1 clus2};
totmols=zeros(1,2);

for i=1:2
    totmols(i)=calc_mols(clusters{i},ch_info);
end

[maxmols,nmax]=max(totmols);
if maxmols==0
    error(['Maximum number of total molecules is zero? ',clusters])
end
largest=clusters{nmax};

if nmax==1
    nmin=2;
else
    nmin=1;
end

smallest=clusters{nmin};
minmols=totmols(nmin);

end