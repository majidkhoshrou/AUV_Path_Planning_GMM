
function BestHypothesis = DefineHypotheses1(NewPoint, time)

persistent MyHypotheses1

if isempty(MyHypotheses1)
    dimens = length(NewPoint) ;
    init_cov = 5 ; 
   MyHypotheses1  = struct('mu',NewPoint, 'cov',eye(dimens)* init_cov, 'prior',1, 'Nk',1, 'loglike',[], 'dl',[], 'time',time,...
       'birth',time, 'K',1) ;  
   indic = mvnpdf(NewPoint, NewPoint, eye(dimens)* init_cov) ;
   loglikelihood = sum(log2(indic)) ;
   MyHypotheses1(1).loglike = loglikelihood ;
   npars = (dimens + dimens*(dimens+1)/2); nparsover2 = npars / 2 ; 
   updated_pp = 1 ;  N_k = 1 ; k = 1 ;
   MyHypotheses1(1).dl = - loglikelihood + (nparsover2 * sum(log2(updated_pp))) + ...
                            (nparsover2 + 0.5)*k*log2(sum(N_k)) ; %%% +k*(nparsover2+0.5)*(1-log2(12));  
else
    MyHypotheses1 = OnlineAlgorithm_3(MyHypotheses1, NewPoint, time);
end
BestHypothesis = MyHypotheses1(1) ;


if ~mod(time,100)
%     clc
    disp('AUV 1')
    disp(BestHypothesis)
    disp(BestHypothesis.mu)
    disp(BestHypothesis.cov)
    disp('------------------------------------')
end