
function BestHypothesis = DefineHypotheses2(NewPoint, time)

persistent MyHypotheses2

if isempty(MyHypotheses2)
    dimens = length(NewPoint) ;
    init_cov = 5 ; 
   MyHypotheses2  = struct('mu',NewPoint, 'cov',eye(dimens)* init_cov, 'prior',1, 'Nk',1, 'loglike',[], 'dl',[], 'time',time,...
       'birth',time, 'K',1) ; % ,'tRaCe',[],'N0rM',[] 
   indic = mvnpdf( NewPoint, NewPoint, eye(dimens)* init_cov) ;
   loglikelihood = sum(log2(indic)) ;
   MyHypotheses2(1).loglike = loglikelihood ;
   npars = (dimens + dimens*(dimens+1)/2); nparsover2 = npars / 2; 
   updated_pp = 1 ;  N_k = 1 ; k = 1 ;
   MyHypotheses2(1).dl = - loglikelihood + (nparsover2 * sum(log2(updated_pp))) + ...
                            (nparsover2 + 0.5)*k*log2(sum(N_k)) ; %%%+k*(nparsover2+0.5)*(1-log2(12));  
else
    MyHypotheses2 = OnlineAlgorithm_3(MyHypotheses2, NewPoint, time);
end

BestHypothesis = MyHypotheses2(1) ;

if ~mod(time,100)    
    disp('AUV 2')
    disp(BestHypothesis)
    disp(BestHypothesis.mu)
    disp(BestHypothesis.cov)
    disp('------------------------------------')
end