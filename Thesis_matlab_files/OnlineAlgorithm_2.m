function MyHypotheses = OnlineAlgorithm_2(Hypothesis, NewPoint, time)

% init_cov = var(NewPoint) /10 ;
init_cov = 8 ;
% merge_thresgold = 1e-2 ;
pp_extinction = 1 ;    pp_thresh = 0.001 ;
max_Hypothesis_K_record = 10 ;
Nk_record = Hypothesis(1).Nk ;
alpha_cte = 800 ;
dimens = length(NewPoint) ;
npars = (dimens + dimens*(dimens+1)/2); nparsover2 = npars / 2; 
min_cov_threshold = eye(dimens) * 1e-6 ;
            j = 1 ;
            while j <= length(Hypothesis) && ~isempty(Hypothesis(j).mu) % In case we did preallocation ... !
                clear updated_mu updated_cov updated_pp ;
                updated_mu = Hypothesis(j).mu ;
                updated_cov = Hypothesis(j).cov ;
                updated_pp = Hypothesis(j).prior ;
                k = length( updated_pp ) ;
                semi_indicator = zeros(1,k) ;
                clear N_k ;
                N_k = Hypothesis(j).Nk ; 
                Hypothesis(j).time = time ;          
                    for i = 1 : k  % Calculate responsibilities
                    semi_indicator(i) = mvnpdf( NewPoint, updated_mu(:,i), updated_cov(:,:,i) ) ;
                    end
                    indicator = semi_indicator.* updated_pp ;
                    indicator( indicator < 1e-84 ) = 1e-84 ;                     
                    normindicator = indicator./ sum(indicator) ; % responsibility
                      alphaa = 1 / ( alpha_cte + time - Hypothesis(j).birth );  % IMPORTANT !!!   
                      alphaa = max( alphaa , 1e-7 ) ;
                         for comp = 1 : k   % Update ...
                              diff = ( NewPoint - updated_mu(:,comp) );
                              updated_mu(:,comp) = updated_mu(:,comp) + alphaa * (normindicator(comp) / updated_pp(comp)) * diff ;                              
                              updated_cov(:,:,comp) = updated_cov(:,:,comp) + alphaa * ( normindicator(comp) / updated_pp(comp) ) * ((diff * diff') - updated_cov(:,:,comp)) ;
% % %                               updated_cov(:,:,comp) = updated_cov(:,:,comp) + min_cov_threshold ; % changed here!!! I AM GOOD, THANKS !!!
                              if det(updated_cov(:,:,comp)) <= 0
                                  disp('.... Determinant was NEGATIVE!')
                                  updated_cov(:,:,comp) = eye(dimens) * 5  ; % Don't put diag(diag(updated_cov(:,:,comp))), leads to some weird results (loglikelihood positive and dl negative and stuffs ...)
                              end                                
                              updated_pp(comp) = updated_pp(comp) + alphaa * ( normindicator(comp) - updated_pp(comp) );
                              % In some datasets DIAGONAL works better!
% % %                               updated_cov(:,:,comp) = diag(diag(updated_cov(:,:,comp))) ;   % + eye(dimens) * realmin
                         end
                              N_k = N_k + normindicator ;
                              updated_pp = check_weights(updated_pp);                             
                              N_k = sum(N_k) * updated_pp;
                     clear normindicator ; clear indicator ; clear semi_indicator ;                  
..... Again calculate indicator to update loglike, it really helps us in building a stronger case
                   semi_indicator = zeros(1,k);
                   for ii = 1 : k 
                        semi_indicator(ii) = mvnpdf( NewPoint , updated_mu(:,ii), updated_cov(:,:,ii));  
                   end 
                    indicator = semi_indicator.* updated_pp ;
                    indicator(indicator < 1e-84) = 1e-84 ;
%                    new_loglike = Hypothesis(j).loglike + ( log2( sum(indicator) ) ) ;  % realmin + sum(indicator)
                        Hypothesis(j).mu = updated_mu;
                    Hypothesis(j).cov = updated_cov;
                    Hypothesis(j).Nk = N_k;
                    Hypothesis(j).prior = updated_pp;
                    Hypothesis(j).K = k ;
                   
                    new_loglike = Expected_LogLikelihood(Hypothesis(j)) ;
                   dl_after_update = - new_loglike + ( nparsover2 * sum(log2(updated_pp)) ) + ...
                               ( nparsover2 + 0.5 ) * k * log2(sum(N_k));%%%+k*(nparsover2+0.5)*(1-log2(12));
                    
                    
                    Hypothesis(j).loglike = new_loglike;
                    Hypothesis(j).dl = dl_after_update; 
                    
                    clear N_k ;
                    clear updated_cov ;  clear updated_mu ;  clear updated_pp ;
              j = j + 1 ;
            end 
              j = j - 1 ;             
%   Add hypothesis : (New point as a new component) + (updated Hypothesis(1))
    clear mn Sigma ;
         mn = Hypothesis(1).mu ;
         Sigma = Hypothesis(1).cov ;
         Hypothesis(j+1).mu = [ mn , NewPoint ] ;
         kk = Hypothesis(1).K ; 
         Sigma(:,:,kk+1) = eye(dimens) * init_cov ;
         Hypothesis(j+1).cov = Sigma ;
         Hypothesis(j+1).Nk = [ Nk_record , 1 ] ;  
         clear Nk_record ; clear Sigma mn ;
         Hypothesis(j+1).prior = Hypothesis(j+1).Nk / sum( Hypothesis(j+1).Nk ) ; 
         updated_pp = Hypothesis(j+1).prior ;
                  k = length( updated_pp ) ;
         clear N_k ;
         N_k = Hypothesis(j+1).Nk ;
                  Hypothesis(j+1).time = time ;
         Hypothesis(j+1).birth = time ; 
         Hypothesis(j+1).K = k ; 
         
         Hypothesis(j+1).loglike = Expected_LogLikelihood(Hypothesis(j+1));
%          Hypothesis(j+1).loglike = Hypothesis(1).loglike + Variational_distance(Hypothesis(1),Hypothesis(j+1))  ; 
         
     

%          N_k = sum(N_k) * updated_pp ;
%          Hypothesis(j+1).Nk = N_k ;
         dl_new_H = - Hypothesis(j+1).loglike + (nparsover2 * sum(log2(updated_pp))) + ...
                    (nparsover2 + 0.5) * k * log2( sum(N_k) ) ;%%%+k*(nparsover2+0.5)*(1-log2(12));
         Hypothesis(j+1).dl = dl_new_H ; 
        
         clear dl_new_H ;
         clear updated_pp ; 
         clear N_k ;
...........................................................................
    
...........................................................................
        merge_threshold = 0.1 * (time/(time+500)) ; % CTD data [0 , 5]
    if 1 % To this point length(Hypothesis) = j+1
    %   Add hypothesis : Post-merge GMM of the best_Hypothesis    
     [ Merged_Cluster , merge_flag ] = Post_Merge_Pre_Merge_5( Hypothesis(1), nparsover2, merge_threshold, time ) ;     
       if merge_flag
          L = length(Hypothesis) ;
          s = sprintf(' Just Merged ... the length is %d', length(Merged_Cluster)) ;
          disp(s)
         Hypothesis(L+1 : L+length(Merged_Cluster)) = Merged_Cluster ;
       end
       clear Merged_Cluster ;   
    end
    
% ..................................................................... 
%   Add hypothesis : Prior Extinction in Hypothesis(1)
%   Be careful with it ! can cause overfitting if pp_thresh is too big !!!
    if  pp_extinction
        L = length(Hypothesis) ;
        if sum( Hypothesis(1).prior > pp_thresh ) && sum( Hypothesis(1).prior < pp_thresh )  
            Hypothesis(L+1).mu = Hypothesis(1).mu( : , Hypothesis(1).prior > pp_thresh );
            Hypothesis(L+1).cov = Hypothesis(1).cov( : , : , Hypothesis(1).prior > pp_thresh ); 
            aux_NK = sum( Hypothesis(1).Nk ) ;
            clear updated_pp ;
            updated_pp = Hypothesis(1).prior( Hypothesis(1).prior > pp_thresh ) ;
            a = (1 - sum(updated_pp)) * updated_pp / sum( updated_pp ); % Add sth, based on their weights! 
            updated_pp = updated_pp + a ;
            clear a ;
            updated_pp = updated_pp / sum( updated_pp ) ;
%             updated_pp = check_weights( updated_pp ) ;
            Hypothesis(L+1).Nk =  aux_NK * updated_pp;
            Hypothesis(L+1).prior = updated_pp ;
            k = length( updated_pp ) ;
            Hypothesis(L+1).K = k ;
            Hypothesis(L+1).loglike = Expected_LogLikelihood(Hypothesis(L+1)) ;            
            clear N_k ;
            N_k = Hypothesis(L+1).Nk ;
            
            dl_Comp = - Hypothesis(L+1).loglike + ( nparsover2 * sum(log2( updated_pp ))) + ...
                            (nparsover2 + 0.5) * k * log2(sum( N_k ));%%%+k*(nparsover2+0.5)*(1-log2(12));
            Hypothesis(L+1).dl = dl_Comp ;
            Hypothesis(L+1).time = time ;
            Hypothesis(L+1).birth = time ;   
%             Hypothesis(L+1).loglike = Hypothesis(1).loglike - Variational_distance(Hypothesis(1),Hypothesis(L+1)) ;

            disp('JUST DID COMPONENT EXTINCTION !')
            clear N_k ;
            clear updated_pp ;
            clear dl_Comp ;
        end
    end
    
    ....................................
      K_each_Hypothesis = arrayfun( @(x) (x.K), Hypothesis ) ;
      unique_K_each_Hypothesis = unique( K_each_Hypothesis ) ;
      find_min_dl = Inf( 1 , length( unique_K_each_Hypothesis ) );  
      Hypothesis_K_record = {1} ;    % Preallocation
      for i1 = 1 : length( unique_K_each_Hypothesis )
          aux_K_record = arrayfun( @(x) (x.K == unique_K_each_Hypothesis(i1)), Hypothesis ) ;
          Hypothesis_K_record{i1} = Hypothesis(aux_K_record) ;
          len_of_this = min( length( nonzeros(aux_K_record) ) , max_Hypothesis_K_record ) ;
          Hypothesis_K_record{i1} = Hypothesis_K_record{i1}( 1 : len_of_this ) ;
          find_min_dl(i1) = Hypothesis_K_record{i1}(1).dl ;
          clear aux_K_record ;
      end
      clear Hypothesis ;
      Hypothesis1 = horzcat( Hypothesis_K_record{:} ) ;
      dl_history = arrayfun( @(x) (x.dl), Hypothesis1 ) ;
      [~ , idx_dl] = sort( dl_history ) ;
      Hypothesis = Hypothesis1(idx_dl(1 : min( 60 , length(idx_dl) ))) ;
      clear Hypothesis1 ;
      clear idx_dl ;
      ...
         MyHypotheses = Hypothesis ;
         