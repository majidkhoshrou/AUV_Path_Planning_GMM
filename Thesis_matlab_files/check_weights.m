
function Priors = check_weights(Priors,varargin)
%
aux = Priors;
if nargin < 2
    prior_threshold = 1e-15 ;
else 
    prior_threshold = varargin{1} ;
end
prior_threshold = min( prior_threshold , 1/length(Priors) ) ;% In case the number of components were too many !
if sum( Priors < prior_threshold )    
    Priors( Priors < prior_threshold ) = prior_threshold ;   
    subtract_term = sum(aux - Priors) * Priors(aux > prior_threshold) / sum(Priors(aux > prior_threshold));
    Priors(aux > prior_threshold) = Priors(aux > prior_threshold) + subtract_term ;  
end
Priors = Priors / sum(Priors) ;
Priors(end) = 1 - sum(Priors(1 : end-1)) ;
