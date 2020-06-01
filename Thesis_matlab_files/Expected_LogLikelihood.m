
function expctd_loglikehd = Expected_LogLikelihood(gmm)

gmm2 = gmm ;
BhGMM = zeros([gmm.K,gmm2.K]);
dimens = size(gmm.mu,1) ;
dimensover2 = dimens/2 ;
for i = 1 : gmm.K
    for j = 1 : gmm2.K
          C = ( gmm.cov(:,:,i)^-1 + gmm2.cov(:,:,j)^-1 )^-1 ;
          muu = C * ((gmm.cov(:,:,i)^-1)*(gmm.mu(:,i)) + (gmm2.cov(:,:,j)^-1)*(gmm2.mu(:,j))) ;
          K = gmm.mu(:,i)' * (gmm.cov(:,:,i)^-1) * (gmm.mu(:,i)) + ...
              gmm2.mu(:,j)' * (gmm2.cov(:,:,j)^-1) * (gmm2.mu(:,j)) - ...
              muu' * C^-1 * muu ;
              nominator = exp(-K/2) ;
              denominator = ((2*pi)^dimensover2) * sqrt(det(gmm.cov(:,:,i)*gmm2.cov(:,:,j)*C)) ;
%               BhGMM(i,j) = gmm2.prior(j)*(nominator / denominator) ; 
              BhGMM(i,j) = nominator / denominator ;
    end
end
BhGMM = BhGMM.* repmat(gmm.prior,gmm.K,1) ;

expctd_loglikehd = sum(log2(sum(BhGMM,2)).*gmm.Nk') ;

