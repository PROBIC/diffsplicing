function[gene_levels,gene_type]=expr_DEgenes(gene_levels,gene_type,ind_DEgene,n_DEgene,ind_DE1_1,ind_DE1_2)

%% genes with changing expression levels

AA=[];
for i=1:n_DEgene/2
    
    i
    test=0;
    while test==0
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
    gene_type(ind_DE1_1(i))=type_gene1;
    gene_type(ind_DE1_2(i))=type_gene2;
    
end
end


