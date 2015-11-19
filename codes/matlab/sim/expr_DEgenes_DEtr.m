function[T, trans_type]=expr_DEgenes_DEtr(c,n_DEgene,gene_levels,t,trans_type,T,ind_DE1_1,ind_DE1_2)
AA_1=[];
AA_2=[];
%% Gene level changing
for ha=1:4
    %for i=1:(n_DEgene/4)/4
    for i=1:(n_DEgene/2)/4
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
    %for i=1:(n_DEgene/4)/4
    for i=1:(n_DEgene/2)/4
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

end

