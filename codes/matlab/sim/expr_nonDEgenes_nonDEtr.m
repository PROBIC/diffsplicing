function[T, trans_type]=expr_nonDEgenes_nonDEtr(c,ind_nonDEgene_nonDEtr,gene_levels,trans_type,T)
%% All constant

for i=1:length(ind_nonDEgene_nonDEtr)
    
    lg=length(c{ind_nonDEgene_nonDEtr(i),2});
    %ratios=drchrnd(lg,1);
    mm=min(floor(gene_levels(ind_CO1_gene(i),1)),lg); % min of the gene level or number of transcripts
    
    ratios=repmat((1/mm),1,mm);
  
    for h=1:mm       
        T(c{ind_nonDEgene_nonDEtr(i),2}(h),:)=gene_levels(ind_nonDEgene_nonDEtr(i),:)*ratios(h);
        trans_type(c{ind_nonDEgene_nonDEtr(i),2}(h))=0;       
    end
    
    % if there are more transcripts than the gene expression (such that each expressed transcript has
    % an expression of at least 1), the excess trancripts are assumed
    % unexpressed:
    
    if mm<lg  
        for h=mm+1:lg
            T(c{ind_nonDEgene_nonDEtr(i),2}(h),:)=repmat(0,1,10);
            trans_type(c{ind_nonDEgene_nonDEtr(i),2}(h))=4; % not expressed
        end
    end
        
end

end