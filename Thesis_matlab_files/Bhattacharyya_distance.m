function gmm_dist = Bhattacharyya_distance(gmm1,gmm2)
% gmm1 , gmm2

BhGMM = zeros([gmm1.K,gmm2.K]);
for i = 1 : gmm1.K
    for j = 1 : gmm2.K
        term1 = (1/8)*(gmm1.mu(:,i)-gmm2.mu(:,j))'*((.5*(gmm1.cov(:,:,i)+gmm2.cov(:,:,j)))^-1)*(gmm1.mu(:,i)-gmm2.mu(:,j)) ;
        term2 = .5 * log(det((.5*(gmm1.cov(:,:,i)+gmm2.cov(:,:,j))))/sqrt(det(gmm1.cov(:,:,i))*det(gmm2.cov(:,:,j)))) ;
        BhGMM(i,j) = gmm1.prior(i)*gmm2.prior(j)*(term1 + term2) ;
    end
end

gmm_dist = sum(sum(BhGMM)) ;

