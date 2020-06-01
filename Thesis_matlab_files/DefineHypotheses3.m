
function BestHypothesis = DefineHypotheses3(NewPoint, time)

persistent MyHypotheses3

if isempty(MyHypotheses3)
    dimens = length(NewPoint) ;
    init_cov = 5 ; 
   MyHypotheses3  = struct('mu',NewPoint, 'cov',eye(dimens)* init_cov, 'prior',1, 'Nk',1, 'loglike',[], 'dl',[], 'time',time,...
       'birth',time, 'K',1) ; 
   indic = mvnpdf( NewPoint, NewPoint, eye(dimens)* init_cov ) ;
   loglikelihood = sum(log2(indic)) ;
   MyHypotheses3(1).loglike = loglikelihood ;
   npars = (dimens + dimens*(dimens+1)/2); nparsover2 = npars / 2 ; 
   updated_pp = 1 ;  N_k = 1 ; k = 1 ;
   MyHypotheses3(1).dl = - loglikelihood + (nparsover2 * sum(log2(updated_pp))) + ...
                            (nparsover2 + 0.5)*k*log2(sum(N_k)) ; %%% +k*(nparsover2+0.5)*(1-log2(12));  
else
    MyHypotheses3 = OnlineAlgorithm_3(MyHypotheses3, NewPoint, time);
end

BestHypothesis = MyHypotheses3(1) ;

if ~mod(time,100)
    disp('AUV 3')
    disp(BestHypothesis)
    disp(BestHypothesis.mu)
    disp(BestHypothesis.cov)
    disp('------------------------------------')
end