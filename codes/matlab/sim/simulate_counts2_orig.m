%load('Homo_sapiens.GRCh37.73.cdna.chr1.ref.tr_workspace.mat')

function[type,gene_levels,T,trans_type,TT]=simulate_counts2_orig(genes,c)

N=length(genes); % total number of genes
N1=N;

len_=zeros(15530,1);
for i=1:N
len_(c{i,2})=len{i};
end


%DE1=floor(10*N1/100);  % number of DE genes
DE1=384;

CO1=N1-DE1;  % number of constants
CO1_tr=floor(10*N1/100); % number of constant genes with different levels of transcripts

s=zeros(N,1);
for i=1:N
s(i)=length(c{i,2});
end
d=find(s>=2); % indices of the genes with at least 2 transcripts
nT=sum(s); % total number of transcripts

ff1=1:N;
ind_DE1=randsample(d,DE1);  % indices for DE genes


ind_CO1=setdiff(ff1,ind_DE1);  % indices for constant genes

remaining=setdiff(d,ind_DE1);  % remaining indices of genes with at least 2 transcripts
ind_CO1_tr=randsample(remaining,CO1_tr);
ind_CO1_gene=setdiff(ind_CO1, ind_CO1_tr); % indices for all constant genes 
%ind_DE1_1=randsample(ind_DE1,DE1/2);
ind_DE1_1=[randsample(ind_DE1(1:48),DE1/8);randsample(ind_DE1(49:96),DE1/8);randsample(ind_DE1(97:144),DE1/8);randsample(ind_DE1(145:192),DE1/8)];
ind_DE1_2=setdiff(ind_DE1,ind_DE1_1);


type=zeros(N,1);
gene_levels=zeros(N,10);

ind_CO=[ind_CO1]; % indices of constant genes
CO=CO1; % total number of constant genes

m=unifrnd(3,50,761,1);
%gene_levels(ind_CO1_gene(1:25),:)=repmat(0,25,10);
gene_levels(ind_CO1_gene(1:761),:)=repmat(m,1,10); % expression levels of the all constant genes
m=unifrnd(3,50,95,1);
gene_levels(ind_CO1_tr(1:95),:)=repmat(m,1,10);
m=unifrnd(51,100,762,1);
gene_levels(ind_CO1_gene(762:1523),:)=repmat(m,1,10);
m=unifrnd(51,100,95,1);
gene_levels(ind_CO1_tr(96:190),:)=repmat(m,1,10);
m=unifrnd(101,150,761,1);
gene_levels(ind_CO1_gene(1524:2284),:)=repmat(m,1,10);
m=unifrnd(101,150,95,1);
gene_levels(ind_CO1_tr(191:285),:)=repmat(m,1,10);
m=unifrnd(151,200,762,1);
gene_levels(ind_CO1_gene(2285:3046),:)=repmat(m,1,10);
m=unifrnd(151,200,96,1);
gene_levels(ind_CO1_tr(286:381),:)=repmat(m,1,10);
type(ind_CO)=repmat(0,CO,1); % type=0 for constant
t=[1;2;3;4;5;6;7;8;9;10];

%% genes with changing expression levels

AA=[];
for i=1:DE1/2
    
    i
    test=0;
    while test==0
        %
        sigma1=unifrnd(1,3);
        %sigma2=unifrnd(1,5);
        sigma2=unifrnd(2,5);
        %c1=unifrnd(0,10);
        c1=unifrnd(3,8);
        c2=unifrnd(c1-sigma1,c1+sigma1); 
        %alpha1=unifrnd(1,5);
        %alpha2=unifrnd(0,1);
        alpha1=unifrnd(0.6,0.9);
        alpha2=1-alpha1;
        %fx_gene=3*(1+((alpha1*exp(-((t-c1)./sigma1).^2))+(alpha2*exp(-((t-c2)./sigma2).^2)))); 
        %fx_gene=10*(1+((alpha1*exp(-((t-c1)./sigma1).^2))+(alpha2*exp(-((t-c2)./sigma2).^2)))); 
        %fx_gene=30*(1+((alpha1*exp(-((t-c1)./sigma1).^2))+(alpha2*exp(-((t-c2)./sigma2).^2)))); 
        
        if i<=48
            fx_gene1=3*(1+((alpha1*exp(-((t-c1)./sigma1).^2))+(alpha2*exp(-((t-c2)./sigma2).^2)))); 
        elseif (i>48 & i<=96)
            fx_gene1=10*(1+((alpha1*exp(-((t-c1)./sigma1).^2))+(alpha2*exp(-((t-c2)./sigma2).^2)))); 
        elseif (i>96 & i<=144)
            fx_gene1=30*(1+((alpha1*exp(-((t-c1)./sigma1).^2))+(alpha2*exp(-((t-c2)./sigma2).^2)))); 
        elseif (i>144 & i<=192)
            fx_gene1=100*(1+((alpha1*exp(-((t-c1)./sigma1).^2))+(alpha2*exp(-((t-c2)./sigma2).^2)))); 
        end
        
        %fx_gene1=500+1000*((alpha1*exp(-((t-c1)./sigma1).^2))+(alpha2*exp(-((t-c2)./sigma2).^2)));
        %r=unifrnd(0.1,0.9);
        %fx_gene1=100./(0.1+exp(-(r*t)));
        mm1=max(fx_gene1)/min(fx_gene1);
        mm2=min(fx_gene1);
        %if (mm1>=1.2) & (mm1<=2.2) & (mm2>=10)
        if (mm1>=1.2) & (mm1<=2.2) 
            test=1;
            AA=[AA;mm1];
        end
    end      
    fx_gene2=flipud(fx_gene1);
    type_gene1=1;  
    type_gene2=2;
    
    gene_levels(ind_DE1_1(i),:)=fx_gene1;
    gene_levels(ind_DE1_2(i),:)=fx_gene2;
    type(ind_DE1_1(i))=type_gene1;
    type(ind_DE1_2(i))=type_gene2;
    
end

T=10000*ones(nT,10);
trans_type=zeros(nT,1);

%% All constant

for i=1:length(ind_CO1_gene)
    
    lg=length(c{ind_CO1_gene(i),2});
    %ratios=drchrnd(lg,1);
    dene=min(floor(gene_levels(ind_CO1_gene(i),1)),lg);
    
    ratios=repmat((1/dene),1,dene);
  
    for h=1:dene       
        T(c{ind_CO1_gene(i),2}(h),:)=gene_levels(ind_CO1_gene(i),:)*ratios(h);
        trans_type(c{ind_CO1_gene(i),2}(h))=0;       
    end
    
    if dene<lg
        for h=dene+1:lg
            T(c{ind_CO1_gene(i),2}(h),:)=repmat(0,1,10);
            trans_type(c{ind_CO1_gene(i),2}(h))=4; % not expressed
        end
    end
        
end

%% Constant genes but changing transcripts
AA1=[];
AA2=[];
for i=1:CO1_tr
    lg=length(c{ind_CO1_tr(i),2});
    dene=min(floor(gene_levels(ind_CO1_tr(i),1)),lg);
    cons=1/dene;
    cons_ch=cons*2;
    test=0;
  %  while test==0
        %
        sigma1=unifrnd(1,3);
        %sigma2=unifrnd(1,5);
        sigma2=unifrnd(2,5);
        %c1=unifrnd(0,10);
        c1=unifrnd(3,8);
        c2=unifrnd(c1-sigma1,c1+sigma1); 
        %alpha1=unifrnd(1,5);
        %alpha2=unifrnd(0,1);
        alpha1=unifrnd(0.6,0.9);
        alpha2=1-alpha1;
        fx_gene=3*(1+((alpha1*exp(-((t-c1)./sigma1).^2))+(alpha2*exp(-((t-c2)./sigma2).^2)))); 
        %fx_gene=10*(1+((alpha1*exp(-((t-c1)./sigma1).^2))+(alpha2*exp(-((t-c2)./sigma2).^2)))); 
        %fx_gene=30*(1+((alpha1*exp(-((t-c1)./sigma1).^2))+(alpha2*exp(-((t-c2)./sigma2).^2)))); 
        %fx_gene=100*(1+((alpha1*exp(-((t-c1)./sigma1).^2))+(alpha2*exp(-((t-c2)./sigma2).^2))));
       
        %
        %r=unifrnd(0.1,0.9);
        %fx_gene=0.1./(0.1+exp(-(r*t)));
        minmax=min(fx_gene)+max(fx_gene);
        fx_gene1=(cons_ch/minmax)*fx_gene;
        fx_gene2=repmat(cons_ch,10,1)-fx_gene1;
        %for h=1:2
        T(c{ind_CO1_tr(i),2}(1),:)=gene_levels(ind_CO1_tr(i),:).*fx_gene1';
        T(c{ind_CO1_tr(i),2}(2),:)=gene_levels(ind_CO1_tr(i),:).*fx_gene2';
        trans_type(c{ind_CO1_tr(i),2}(1))=1;
        trans_type(c{ind_CO1_tr(i),2}(2))=2;
        %end
        mm1=max(T(c{ind_CO1_tr(i),2}(1),:))/min(T(c{ind_CO1_tr(i),2}(1),:));
        mm2=max(T(c{ind_CO1_tr(i),2}(2),:))/min(T(c{ind_CO1_tr(i),2}(2),:));
        AA1=[AA1;max(fx_gene1)/min(fx_gene1)];
        AA2=[AA2;max(fx_gene2)/min(fx_gene2)];
      %  if (mm1>=1.2)  & (mm1<=2.2) & (mm2>=1.2) & (mm2<=2.2)
      %      test=1;
      %  end
    %end           
    if lg>2 & dene>2 
        for h=3:dene
            T(c{ind_CO1_tr(i),2}(h),:)=gene_levels(ind_CO1_tr(i),:)*cons;
            trans_type(c{ind_CO1_tr(i),2}(h))=0;
        end
        if  dene<lg
            for h=dene+1:lg
                T(c{ind_CO1_tr(i),2}(h),:)=repmat(0,1,10);
                trans_type(c{ind_CO1_tr(i),2}(h))=4; % not expressed
            end
        end
    end
end

AA_1=[];
AA_2=[];
%% Gene level changing
for ha=1:4
    for i=1:(DE1/4)/4
    lg=length(c{ind_DE1_1((ha-1)*48+i),2});
    %cons_ch=0.5;
    dene=min(floor(min(gene_levels(ind_DE1_1((ha-1)*48+i),:))),lg);
    cons=1/dene;
   % cons=1/lg;
    %cons=(1-cons_ch)/(lg-2);
    cons_ch=cons*2;
    test=0;
    %while test==0
        %
        sigma1=unifrnd(1,3);
        %sigma2=unifrnd(1,5);
        sigma2=unifrnd(2,5);
        %c1=unifrnd(0,10);
        c1=unifrnd(3,8);
        c2=unifrnd(c1-sigma1,c1+sigma1); 
        %alpha1=unifrnd(1,5);
        %alpha2=unifrnd(0,1);
        alpha1=unifrnd(0.6,0.9);
        alpha2=1-alpha1;
        fx_gene=3*(1+((alpha1*exp(-((t-c1)./sigma1).^2))+(alpha2*exp(-((t-c2)./sigma2).^2)))); 
        %
        %r=unifrnd(0.1,0.9);
        %fx_gene=0.1./(0.1+exp(-(r*t)));
        minmax=min(fx_gene)+max(fx_gene);
        fx_gene1=(cons_ch/minmax)*fx_gene;
        fx_gene2=repmat(cons_ch,10,1)-fx_gene1;
        
        mm1=max(fx_gene1)/min(fx_gene1);
        mm2=max(fx_gene2)/min(fx_gene2);
        AA_1=[AA_1;mm1];
        AA_2=[AA_2;mm2];
     %   if (mm1>=1.2)  & (mm1<=2.2) & (mm2>=1.2) & (mm2<=2.2)
     %       test=1;
     %   end
    %end     
        %for h=1:2
        T(c{ind_DE1_1((ha-1)*48+i),2}(1),:)=gene_levels(ind_DE1_1((ha-1)*48+i),:).*fx_gene1';
        T(c{ind_DE1_1((ha-1)*48+i),2}(2),:)=gene_levels(ind_DE1_1((ha-1)*48+i),:).*fx_gene2';
        trans_type(c{ind_DE1_1((ha-1)*48+i),2}(1))=1;
        trans_type(c{ind_DE1_1((ha-1)*48+i),2}(2))=2;
        %end
    %    mm1=max(T(c{ind_DE1_1(i),2}(1),:))/min(T(c{ind_DE1_1(i),2}(1),:));
    %    mm2=max(T(c{ind_DE1_1(i),2}(2),:))/min(T(c{ind_DE1_1(i),2}(2),:));
    %    if (mm1>=2) & (mm2>=2)
    %        test=1;
    %    end
    %end           
    if lg>2 & dene>2
        for h=3:dene
            T(c{ind_DE1_1((ha-1)*48+i),2}(h),:)=gene_levels(ind_DE1_1((ha-1)*48+i),:)*cons;
            trans_type(c{ind_DE1_1((ha-1)*48+i),2}(h))=0;
        end
        if  dene<lg
            for h=dene+1:lg
                T(c{ind_DE1_1((ha-1)*48+i),2}(h),:)=repmat(0,1,10);
                trans_type(c{ind_DE1_1((ha-1)*48+i),2}(h))=4; % not expressed
            end
        end
    end
    end
end

for ha=1:4
    for i=1:(DE1/4)/4
    lg=length(c{ind_DE1_2((ha-1)*48+i),2});
    %cons_ch=0.5;
    %cons_ch=0.5;
    dene=min(floor(min(gene_levels(ind_DE1_2((ha-1)*48+i),:))),lg);
    cons=1/dene;
    %cons=(1-cons_ch)/(lg-2);
    cons_ch=cons*2;
    test=0;
    %while test==0
        %
        sigma1=unifrnd(1,3);
        %sigma2=unifrnd(1,5);
        sigma2=unifrnd(2,5);
        %c1=unifrnd(0,10);
        c1=unifrnd(3,8);
        c2=unifrnd(c1-sigma1,c1+sigma1); 
        %alpha1=unifrnd(1,5);
        %alpha2=unifrnd(0,1);
        alpha1=unifrnd(0.6,0.9);
        alpha2=1-alpha1;
        fx_gene=3*(1+((alpha1*exp(-((t-c1)./sigma1).^2))+(alpha2*exp(-((t-c2)./sigma2).^2)))); 
        %
        %r=unifrnd(0.1,0.9);
        %fx_gene=0.1./(0.1+exp(-(r*t)));
        minmax=min(fx_gene)+max(fx_gene);
        fx_gene1=(cons_ch/minmax)*fx_gene;
        fx_gene2=repmat(cons_ch,10,1)-fx_gene1;
        
        mm1=max(fx_gene1)/min(fx_gene1);
        mm2=max(fx_gene2)/min(fx_gene2);
        AA_1=[AA_1;mm1];
        AA_2=[AA_2;mm2];
      %  if (mm1>=1.2) & (mm2>=1.2) & (mm1<=2.2) & (mm2<=2.2)
      %      test=1;
      %  end
    %end    
        %for h=1:2
        T(c{ind_DE1_2((ha-1)*48+i),2}(1),:)=gene_levels(ind_DE1_2((ha-1)*48+i),:).*fx_gene1';
        T(c{ind_DE1_2((ha-1)*48+i),2}(2),:)=gene_levels(ind_DE1_2((ha-1)*48+i),:).*fx_gene2';
        trans_type(c{ind_DE1_2((ha-1)*48+i),2}(1))=1;
        trans_type(c{ind_DE1_2((ha-1)*48+i),2}(2))=2;
        %end
      %  mm1=max(T(c{ind_DE1_2(i),2}(1),:))/min(T(c{ind_DE1_2(i),2}(1),:));
      %  mm2=max(T(c{ind_DE1_2(i),2}(2),:))/min(T(c{ind_DE1_2(i),2}(2),:));
      %  if (mm1>=2) & (mm2>=2)
      %      test=1;
      %  end
    %end           
    if lg>2 & dene>2
        for h=3:dene
            T(c{ind_DE1_2((ha-1)*48+i),2}(h),:)=gene_levels(ind_DE1_2((ha-1)*48+i),:)*cons;
            trans_type(c{ind_DE1_2((ha-1)*48+i),2}(h))=0;
        end
        if  dene<lg
            for h=dene+1:lg
                T(c{ind_DE1_2((ha-1)*48+i),2}(h),:)=repmat(0,1,10);
                trans_type(c{ind_DE1_2((ha-1)*48+i),2}(h))=4; % not expressed
            end
        end
    end
end
end

for ha=1:4
    for i=(DE1/4)/4+1:(DE1/2)/4
        lg=length(c{ind_DE1_1((ha-1)*48+i),2});
    %cons_ch=0.5;
    dene=min(floor(min(gene_levels(ind_DE1_1((ha-1)*48+i),:))),lg);
    cons=1/dene;
   % cons=1/lg;
    %cons=(1-cons_ch)/(lg-2);
    cons_ch=cons*2;
    test=0;
    %while test==0
        %
        sigma1=unifrnd(1,3);
        %sigma2=unifrnd(1,5);
        sigma2=unifrnd(2,5);
        %c1=unifrnd(0,10);
        c1=unifrnd(3,8);
        c2=unifrnd(c1-sigma1,c1+sigma1); 
        %alpha1=unifrnd(1,5);
        %alpha2=unifrnd(0,1);
        alpha1=unifrnd(0.6,0.9);
        alpha2=1-alpha1;
        fx_gene=3*(1+((alpha1*exp(-((t-c1)./sigma1).^2))+(alpha2*exp(-((t-c2)./sigma2).^2)))); 
        %
        %r=unifrnd(0.1,0.9);
        %fx_gene=0.1./(0.1+exp(-(r*t)));
        minmax=min(fx_gene)+max(fx_gene);
        fx_gene1=(cons_ch/minmax)*fx_gene;
        fx_gene2=repmat(cons_ch,10,1)-fx_gene1;
        
        mm1=max(fx_gene1)/min(fx_gene1);
        mm2=max(fx_gene2)/min(fx_gene2);
        AA_1=[AA_1;mm1];
        AA_2=[AA_2;mm2];
     %   if (mm1>=1.2)  & (mm1<=2.2) & (mm2>=1.2) & (mm2<=2.2)
     %       test=1;
     %   end
    %end     
        %for h=1:2
        T(c{ind_DE1_1((ha-1)*48+i),2}(1),:)=gene_levels(ind_DE1_1((ha-1)*48+i),:).*fx_gene1';
        T(c{ind_DE1_1((ha-1)*48+i),2}(2),:)=gene_levels(ind_DE1_1((ha-1)*48+i),:).*fx_gene2';
        trans_type(c{ind_DE1_1((ha-1)*48+i),2}(1))=1;
        trans_type(c{ind_DE1_1((ha-1)*48+i),2}(2))=2;
        %end
    %    mm1=max(T(c{ind_DE1_1(i),2}(1),:))/min(T(c{ind_DE1_1(i),2}(1),:));
    %    mm2=max(T(c{ind_DE1_1(i),2}(2),:))/min(T(c{ind_DE1_1(i),2}(2),:));
    %    if (mm1>=2) & (mm2>=2)
    %        test=1;
    %    end
    %end           
    if lg>2 & dene>2
        for h=3:dene
            T(c{ind_DE1_1((ha-1)*48+i),2}(h),:)=gene_levels(ind_DE1_1((ha-1)*48+i),:)*cons;
            trans_type(c{ind_DE1_1((ha-1)*48+i),2}(h))=0;
        end
        if  dene<lg
            for h=dene+1:lg
                T(c{ind_DE1_1((ha-1)*48+i),2}(h),:)=repmat(0,1,10);
                trans_type(c{ind_DE1_1((ha-1)*48+i),2}(h))=4; % not expressed
            end
        end
    end
    end
end

for ha=1:4
    for i=(DE1/4)/4+1:(DE1/2)/4
    lg=length(c{ind_DE1_2((ha-1)*48+i),2});
    %cons_ch=0.5;
    %cons_ch=0.5;
    dene=min(floor(min(gene_levels(ind_DE1_2((ha-1)*48+i),:))),lg);
    cons=1/dene;
    %cons=(1-cons_ch)/(lg-2);
    cons_ch=cons*2;
    test=0;
    %while test==0
        %
        sigma1=unifrnd(1,3);
        %sigma2=unifrnd(1,5);
        sigma2=unifrnd(2,5);
        %c1=unifrnd(0,10);
        c1=unifrnd(3,8);
        c2=unifrnd(c1-sigma1,c1+sigma1); 
        %alpha1=unifrnd(1,5);
        %alpha2=unifrnd(0,1);
        alpha1=unifrnd(0.6,0.9);
        alpha2=1-alpha1;
        fx_gene=3*(1+((alpha1*exp(-((t-c1)./sigma1).^2))+(alpha2*exp(-((t-c2)./sigma2).^2)))); 
        %
        %r=unifrnd(0.1,0.9);
        %fx_gene=0.1./(0.1+exp(-(r*t)));
        minmax=min(fx_gene)+max(fx_gene);
        fx_gene1=(cons_ch/minmax)*fx_gene;
        fx_gene2=repmat(cons_ch,10,1)-fx_gene1;
        
        mm1=max(fx_gene1)/min(fx_gene1);
        mm2=max(fx_gene2)/min(fx_gene2);
        AA_1=[AA_1;mm1];
        AA_2=[AA_2;mm2];
      %  if (mm1>=1.2) & (mm2>=1.2) & (mm1<=2.2) & (mm2<=2.2)
      %      test=1;
      %  end
    %end    
        %for h=1:2
        T(c{ind_DE1_2((ha-1)*48+i),2}(1),:)=gene_levels(ind_DE1_2((ha-1)*48+i),:).*fx_gene1';
        T(c{ind_DE1_2((ha-1)*48+i),2}(2),:)=gene_levels(ind_DE1_2((ha-1)*48+i),:).*fx_gene2';
        trans_type(c{ind_DE1_2((ha-1)*48+i),2}(1))=1;
        trans_type(c{ind_DE1_2((ha-1)*48+i),2}(2))=2;
        %end
      %  mm1=max(T(c{ind_DE1_2(i),2}(1),:))/min(T(c{ind_DE1_2(i),2}(1),:));
      %  mm2=max(T(c{ind_DE1_2(i),2}(2),:))/min(T(c{ind_DE1_2(i),2}(2),:));
      %  if (mm1>=2) & (mm2>=2)
      %      test=1;
      %  end
    %end           
    if lg>2 & dene>2
        for h=3:dene
            T(c{ind_DE1_2((ha-1)*48+i),2}(h),:)=gene_levels(ind_DE1_2((ha-1)*48+i),:)*cons;
            trans_type(c{ind_DE1_2((ha-1)*48+i),2}(h))=0;
        end
        if  dene<lg
            for h=dene+1:lg
                T(c{ind_DE1_2((ha-1)*48+i),2}(h),:)=repmat(0,1,10);
                trans_type(c{ind_DE1_2((ha-1)*48+i),2}(h))=4; % not expressed
            end
        end
    end
    end
end

%TT=round(T); % get integer absolute counts


LEN=repmat(len_,1,10);
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
normV=normT+((0.05)^2)*(normT.^2); % 0.01->0.1->0.5  %0.05>0.1>0.2
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

save('simulated_counts_005.mat')



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
