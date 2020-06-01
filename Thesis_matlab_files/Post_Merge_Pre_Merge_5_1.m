function [Merged_Cluster ,merge_flag ] = Post_Merge_Pre_Merge_5_1( Cluster, nparsover2, merge_thresgold, time )

num_clusters = length( Cluster.prior ) ;
   if num_clusters == 1 % In case we have only one Cluster !
        Merged_Cluster = Cluster ;
        merge_flag = 0 ;
%         B_ij = 0 ;
        return ;
   end
estmu = Cluster.mu ;
estcov = Cluster.cov ;
prior = Cluster.prior ;
Nk = Cluster.Nk ;
old_loglike = Cluster.loglike ;
dimens = size(estmu,1) ;
.................              Preallocation    
Pre_individual_Cluster = struct('mu',[],'cov',[],'prior',[],'Nk',[]) ;
Pre_individual_Cluster = repmat(Pre_individual_Cluster,1,num_clusters) ;
..............................................
Merged_Cluster = struct('mu',[],'cov',[],'prior',[],'Nk',[],'loglike',[],'dl',[],'time',[],'birth',[],'K',[]) ;
..............................................
for i = 1 : num_clusters
    Pre_individual_Cluster(i) = struct( 'mu',estmu(:,i),'cov',estcov(:,:,i),'prior', prior(i),'Nk',Nk(i) ) ;
end                
            B_ij = Inf( num_clusters - 1 , num_clusters );
             for g = 1 : num_clusters
                 for gg = g + 1 : num_clusters
                    B_ij(g,gg) = My_KL( Pre_individual_Cluster(g),Pre_individual_Cluster(gg) );
                 end
             end
              IND = find( min(min(B_ij)) < merge_thresgold ) ;   %  IND = find(B_ij == min(nonzeros(B_ij))); 
                if isempty( IND )
                    Merged_Cluster = Cluster ;
                    merge_flag = 0 ; 
                    return ;
                end
         for id = 1 : length(IND)
                 [I,J] = ind2sub(size(B_ij),IND(id)) ;                            
            aux_Cluster = Pre_individual_Cluster ;
            aux_Nk = Pre_individual_Cluster(I).Nk + Pre_individual_Cluster(J).Nk ;
            record = 1 : num_clusters ;
            record = record( record ~= I & record ~= J ) ;
            Post_individual_Cluster = aux_Cluster( record ) ;
            [~, w_ij, mu_ij, p_ij] = My_KL( Pre_individual_Cluster(I), Pre_individual_Cluster(J) ) ;   % Now we merge two clusters ...
            kk = length(record);
            Post_individual_Cluster(kk+1)= struct( 'mu',mu_ij, 'cov',p_ij, 'prior',w_ij, 'Nk',aux_Nk ) ;
            num_clusters = length( Post_individual_Cluster ) ;
            updated_cov = zeros(dimens, dimens, num_clusters) ;
                for k = 1 : num_clusters
                        updated_cov( :,:,k ) = Post_individual_Cluster(k).cov ;                   
                end
                updated_mu = horzcat(Post_individual_Cluster.mu) ;     %  updated_mu = arrayfun( @(x) (x.mu), Post_individual_Cluster )           
                updated_pp = horzcat(Post_individual_Cluster.prior) ;
                N_k = horzcat(Post_individual_Cluster.Nk) ;
                updated_pp = check_weights(updated_pp);
                N_k = sum( N_k ) * updated_pp ;
                k = length(updated_pp) ;        
                dl = - (old_loglike - B_ij( I,J ) ) + (nparsover2*sum(log(updated_pp))) + ...
                                (nparsover2 + 0.5) * k * log(sum(N_k));%%%+k*(nparsover2+0.5)*(1-log(12));
               Merged_Cluster(id).mu = updated_mu ;
               Merged_Cluster(id).cov = updated_cov ;
               Merged_Cluster(id).prior = updated_pp ;
               Merged_Cluster(id).Nk = N_k ;
               Merged_Cluster(id).loglike = old_loglike - B_ij( I , J ) ;
               Merged_Cluster(id).dl = dl;
               Merged_Cluster(id).time = time ;
               Merged_Cluster(id).birth = time ;
               Merged_Cluster(id).K = length( updated_pp ) ;
               clear updated_mu updated_cov updated_pp N_k dl
               clear Post_individual_Cluster; clear aux_Cluster ;
         end        
           merge_flag = 1 ;
    