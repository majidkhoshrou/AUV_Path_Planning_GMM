function Individual_Components = SeparateComponents(gmm)

estmu = gmm.mu ;
estcov = gmm.cov ;
prior = gmm.prior ;
Nk = gmm.Nk ;
num_clusters = gmm.K ;

Individual_Components = struct('mu',[],'cov',[],'prior',[],'Nk',[]) ;
Individual_Components = repmat(Individual_Components,1,num_clusters) ;

for i = 1 : num_clusters
    Individual_Components(i) = struct( 'mu',estmu(:,i),'cov',estcov(:,:,i),'prior', prior(i),'Nk',Nk(i) ) ;
end 