
function Variational_D = Variational_distance(gmm1, gmm2)

 Individual_Components_gmm1 = SeparateComponents(gmm1) ;
 Individual_Components_gmm2 = SeparateComponents(gmm2) ;

  pair_KL_gmm1 = zeros( gmm1.K  , gmm1.K ) ;
  pair_KL_gmm1_2 = zeros(gmm1.K, gmm2.K) ;
  
             for g = 1 : gmm1.K
                 for gg =  1 : gmm1.K
                    pair_KL_gmm1(g,gg) = My_KL( Individual_Components_gmm1(g),Individual_Components_gmm1(gg) );
                 end
                 for jj =  1 : gmm2.K
                    pair_KL_gmm1_2(g,jj) = My_KL( Individual_Components_gmm1(g),Individual_Components_gmm2(jj) );
                 end
             end
          
             % formula 20 KLdiv.pdf paper !
             
          nominator = sum( (repmat(gmm1.prior,gmm1.K,1)).* exp(-pair_KL_gmm1),2) ;
          denominator = sum( (repmat(gmm2.prior,gmm1.K,1)).* exp(-pair_KL_gmm1_2),2) ;
          
     Variational_D = sum( (gmm1.prior)'.* (log(nominator./denominator))) ;
      
      
      
      
      
      
             
             
          