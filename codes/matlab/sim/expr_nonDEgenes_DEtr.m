function[T, trans_type]=expr_nonDEgenes_DEtr(c,t,n_nonDEgene_DEtr,ind_nonDEgene_DEtr,gene_levels,trans_type,T)
%% Constant genes but changing transcripts
AA1=[];
AA2=[];
for i=1:n_nonDEgene_DEtr
    lg=length(c{ind_nonDEgene_DEtr(i),2});
    mm=min(floor(gene_levels(ind_nonDEgene_DEtr(i),1)),lg);
    cons=1/mm;
    cons_ch=cons*2;
        
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
    scale=3*(1+((alpha1*exp(-((t-c1)./sigma1).^2))+(alpha2*exp(-((t-c2)./sigma2).^2)))); 
    %scale=10*(1+((alpha1*exp(-((t-c1)./sigma1).^2))+(alpha2*exp(-((t-c2)./sigma2).^2)))); 
    %scale=30*(1+((alpha1*exp(-((t-c1)./sigma1).^2))+(alpha2*exp(-((t-c2)./sigma2).^2)))); 
    %scale=100*(1+((alpha1*exp(-((t-c1)./sigma1).^2))+(alpha2*exp(-((t-c2)./sigma2).^2))));
    %
    %r=unifrnd(0.1,0.9);
    %scale=0.1./(0.1+exp(-(r*t)));
    minmax=min(scale)+max(scale);
    ratio1=(cons_ch/minmax)*scale;
    
    ratio2=repmat(cons_ch,10,1)-ratio1;
    %for h=1:2
    T(c{ind_nonDEgene_DEtr(i),2}(1),:)=gene_levels(ind_nonDEgene_DEtr(i),:).*ratio1';
    T(c{ind_nonDEgene_DEtr(i),2}(2),:)=gene_levels(ind_nonDEgene_DEtr(i),:).*ratio2';
    trans_type(c{ind_nonDEgene_DEtr(i),2}(1))=1;
    trans_type(c{ind_nonDEgene_DEtr(i),2}(2))=2;
    %end
    mm1=max(T(c{ind_nonDEgene_DEtr(i),2}(1),:))/min(T(c{ind_nonDEgene_DEtr(i),2}(1),:));
    mm2=max(T(c{ind_nonDEgene_DEtr(i),2}(2),:))/min(T(c{ind_nonDEgene_DEtr(i),2}(2),:));
    AA1=[AA1;max(ratio1)/min(ratio1)];
    AA2=[AA2;max(ratio2)/min(ratio2)];
    %  if (mm1>=1.2)  & (mm1<=2.2) & (mm2>=1.2) & (mm2<=2.2)
      %      test=1;
      %  end
    %end           
    if lg>2 & mm>2 
        for h=3:mm
            T(c{ind_nonDEgene_DEtr(i),2}(h),:)=gene_levels(ind_nonDEgene_DEtr(i),:)*cons;
            trans_type(c{ind_nonDEgene_DEtr(i),2}(h))=0;
        end
        if  mm<lg
            for h=mm+1:lg
                T(c{ind_nonDEgene_DEtr(i),2}(h),:)=repmat(0,1,10);
                trans_type(c{ind_nonDEgene_DEtr(i),2}(h))=4; % not expressed
            end
        end
    end
end
end