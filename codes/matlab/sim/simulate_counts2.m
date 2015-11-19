%load('Homo_sapiens.GRCh37.73.cdna.chr1.ref.tr_workspace.mat')

function[type,gene_levels,T,trans_type,TT]=simulate_counts2(genes,c)

N_gene=length(genes); % total number of genes
N_tr=0;
for i=1:N_gene
    N_tr=N_tr+length(c{i,2}); % total number of transcripts
end

tr_len=zeros(N_tr,1);
for i=1:N_gene
    tr_len(c{i,2})=len{i};  % vector of the transcript lengths
end


%DE1=floor(10*N1/100);  % number of DE genes
n_DEgene=384;   % number of DE genes

n_nonDEgene=N_gene-n_DEgene;  % number of non DE genes
n_nonDEgene_DEtr=floor(10*N_gene/100); % number of nonDE genes with different levels of transcripts

n_tr=zeros(N_gene,1);
for i=1:N_gene
    n_tr(i)=length(c{i,2});  % number of transcripts of each gene
end

n_2tr=find(n_tr>=2); % indices of the genes with at least 2 transcripts


%ff1=1:N;
ind_DEgene=randsample(n_2tr,n_DEgene);  % indices for DE genes, which all have more than 2 transcripts
ind_nonDEgene=setdiff(1:N_gene,ind_DEgene);  % indices for nonDE genes

remaining=setdiff(n_2tr,ind_DEgene);  % remaining indices of genes with at least 2 transcripts
ind_nonDEgene_DEtr=randsample(remaining,n_nonDEgene_DEtr); % indices of the nonDE genes with DE transcripts
ind_nonDEgene_nonDEtr=setdiff(ind_nonDEgene, ind_nonDEgene_DEtr); % indices of nonDE genes with all nonDE transcripts 
%ind_DE1_1=randsample(ind_DE1,DE1/2);


gene_type=zeros(N_gene,1); % gene types: DE(1) or nonDE(0)
gene_type(ind_nonDEgene)=repmat(0,n_nonDEgene,1); % type=0 for nonDE genes
gene_levels=zeros(N_gene,10);

% Assign expressions in 4 different ranges (3-50, 51-100, 101-150, 151-200) to the all nonDE genes:
m=unifrnd(3,50,761,1);
gene_levels(ind_nonDEgene_nonDEtr(1:761),:)=repmat(m,1,10); % expression levels of the all constant genes
m=unifrnd(51,100,762,1);
gene_levels(ind_nonDEgene_nonDEtr(762:1523),:)=repmat(m,1,10);
m=unifrnd(101,150,761,1);
gene_levels(ind_nonDEgene_nonDEtr(1524:2284),:)=repmat(m,1,10);
m=unifrnd(151,200,762,1);
gene_levels(ind_nonDEgene_nonDEtr(2285:3046),:)=repmat(m,1,10);

% Assign expressions in 4 different ranges to the nonDE genes with DE
% transcripts:
m=unifrnd(3,50,95,1);
gene_levels(ind_nonDEgene_DEtr(1:95),:)=repmat(m,1,10);
m=unifrnd(51,100,95,1);
gene_levels(ind_nonDEgene_DEtr(96:190),:)=repmat(m,1,10);
m=unifrnd(101,150,95,1);
gene_levels(ind_nonDEgene_DEtr(191:285),:)=repmat(m,1,10);
m=unifrnd(151,200,96,1);
gene_levels(ind_nonDEgene_DEtr(286:381),:)=repmat(m,1,10);

t=(1:10)'; % time points

T=10000*ones(nT,10);
trans_type=zeros(nT,1);
[T, trans_type]=expr_nonDEgenes_nonDEtr(c,ind_nonDEgene_nonDEtr,gene_levels,trans_type,T);

[T, trans_type]=expr_nonDEgenes_DEtr(c,t,n_nonDEgene_DEtr,ind_nonDEgene_DEtr,gene_levels,trans_type,T);



ind_DE1_1=[randsample(ind_DEgene(1:48),n_DEgene/8);randsample(ind_DEgene(49:96),n_DEgene/8);randsample(ind_DEgene(97:144),n_DEgene/8);randsample(ind_DEgene(145:192),n_DEgene/8)];
ind_DE1_2=setdiff(ind_DEgene,ind_DE1_1);
[gene_levels,gene_type]=expr_DEgenes(gene_levels,gene_type,ind_DEgene,n_DEgene,ind_DE1_1,ind_DE1_2);

[T, trans_type]=expr_DEgenes_DEtr(c,n_DEgene,gene_levels,t,trans_type,T,ind_DE1_1,ind_DE1_2);


%TT=round(T); % get integer absolute counts


LEN=repmat(tr_len,1,10);
COUNTS=LEN.*T;
COUNTS=COUNTS./(10^2);

save('sim_base2.mat')

% V=T+((0.01)^2)*(T.^2);
% sim1=zeros(nT,10);
% sim2=zeros(nT,10);
% sim3=zeros(nT,10);
% for i=1:nT
%     i
%      for j=1:10
%         v=sqrt(V(i,j));
%         test=0;
%         while test==0
%             a=repmat(T(i,j),3,1)+normrnd(0,v,3,1);
%             if (length(find(a<0))==0)
%                 test=1;
%             end
%         end
%         sim1(i,j)=a(1);
%         sim2(i,j)=a(2);
%         sim3(i,j)=a(3);
%      end
% end
% sim1=ceil(sim1);
% sim2=ceil(sim2);
% sim3=ceil(sim3);

%

% normgene_levels=gene_levels./(repmat(sum(gene_levels,1),N,1));
% normgene_levels=normgene_levels*1500000;
% normT=T./(repmat(sum(T,1),nT,1));
% normT=normT*1500000;

close all
fclose all
clear all

load('sim_base2.mat')
normT=COUNTS;
ovd=0.05; %overdispersion parameter {0.05, 0.1, 0.2}
normV=normT+((ovd)^2)*(normT.^2); 
normsim1=zeros(nT,10);
normsim2=zeros(nT,10);
normsim3=zeros(nT,10);
for i=1:nT
    i
     for j=1:10
        v=normV(i,j);
        mn=normT(i,j);
        p=1-(mn/v);
        r=(mn*(1-p))/p;
        test=0;
        while test==0
            %a=repmat(normT(i,j),3,1)+normrnd(0,v,3,1);
            a=nbinrnd(r,(1-p),3,1);
            if (length(find(a<0))==0)
                test=1;
            end
        end
        normsim1(i,j)=a(1);
        normsim2(i,j)=a(2);
        normsim3(i,j)=a(3);
     end
end
normsim1=ceil(normsim1);
normsim2=ceil(normsim2);
normsim3=ceil(normsim3);

ffff=find(normsim1==0);
normsim1(ffff)=ones(length(ffff),1);
ffff=find(normsim2==0);
normsim2(ffff)=ones(length(ffff),1);
ffff=find(normsim3==0);
normsim3(ffff)=ones(length(ffff),1);

ind_not_expressed=find(trans_type==4);
ind_expressed=setdiff(1:nT,ind_not_expressed)';

str_ovd=num2str(ovd);
str_ovd(str_ovd=='.') = [];
fname=['simulated_counts_',str_ovd,'.mat'];
%save('simulated_counts_005.mat')

end

close all
fclose all
clear all

cd('02/')
load('simulated_counts_02.mat')
f0=fopen('Homo_sapiens.GRCh37.73.cdna.chr1.ref.tr')
t0=textscan(f0,'%s%s%f%f')
fclose(f0)
gn1=t0{1}(ind_expressed);
trn1=t0{2}(ind_expressed);
len1=len_(ind_expressed);
dlmwrite('genes',char(gn1),'delimiter','')
dlmwrite('transcripts',char(trn1),'delimiter','')
dlmwrite('lengths',len1)

normsim1_new=normsim1(ind_expressed,:);
normsim2_new=normsim2(ind_expressed,:);
normsim3_new=normsim3(ind_expressed,:);
trans_type_new=trans_type(ind_expressed);

for i=1:10
     result=[normsim1_new(:,i) trans_type_new];
     filename=['t',num2str(i),'_r1'];
     dlmwrite(filename,result,'delimiter', ' ')
 end
 
 for i=1:10
     result=[normsim2_new(:,i) trans_type_new];
     filename=['t',num2str(i),'_r2'];
     dlmwrite(filename,result,'delimiter', ' ')
 end
 
 for i=1:10
     result=[normsim3_new(:,i) trans_type_new];
     filename=['t',num2str(i),'_r3'];
     dlmwrite(filename,result,'delimiter', ' ')
 end
 
% paste genes transcripts lengths > base
% paste t1_r1 base > t1_r1.exp
% sed 's/\t/ /g' t1_r1.exp > t1_r1.exps

% sh generate_files.sh

% # M 14495

% python select.py Homo_sapiens.GRCh37.73.cdna.chr1.ref.fa transcripts > expressed.fa

% python getReads.py --ref expressed.fa t1_r1.exps

iii=find(type~=0);
hh=zeros(200,1);
for i=1:200
hh(i)=max(normgene_levels(iii(i),:))/min(normgene_levels(iii(i),:));
end
