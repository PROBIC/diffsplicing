function[matrixP matrixR vectorAP]=getPRcurves(experiment, files)

nFiles=length(files);
exp_file=[experiment,'_true'];
f0=fopen(exp_file)
t0=textscan(f0,'%s%f');
fclose(f0)
I=t0{2};

valuesMatrix=zeros(length(I),nFiles);

for i=1:nFiles
    filename=files{i};
    f0=fopen(filename)
    t0=textscan(f0,'%s%f')
    fclose(f0)
    valuesMatrix(:,i)=t0{2};
end

nans=isnan(valuesMatrix);
common_indices=find(((all(nans'==0))==1))';

valuesMatrix=valuesMatrix(common_indices,:);
I=I(common_indices);

matrixP=zeros(length(common_indices),nFiles);
matrixR=zeros(length(common_indices),nFiles);
vectorAP=zeros(1,nFiles);

for j=1:nFiles
    [matrixP(:,j) matrixR(:,j) vectorAP(j)] = getPrecRec(valuesMatrix(:,j),I);
end


end