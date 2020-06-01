function [B_ij,w_ij,mu_ij,p_ij] = My_KL(Cluster_i,Cluster_j)
% definitions are based on Kullback_Leibler Approach to Gaussian Mixture
% Reduction By Andrew R. Runnalls

w_i = Cluster_i.prior;
mu_i = Cluster_i.mu;
p_i = Cluster_i.cov;

w_j = Cluster_j.prior;
mu_j = Cluster_j.mu;
p_j = Cluster_j.cov;

% Eq. (2) to (4) in the paper
w_ij = w_i + w_j;
mu_ij = (w_i/w_ij)* mu_i + (w_j/w_ij)* mu_j;
p_ij = (w_i/w_ij)*p_i + (w_j/w_ij)*p_j + (w_i/w_ij)*(w_j/w_ij)*(mu_i-mu_j)*(mu_i-mu_j)';

% Eq. (21) in the paper
B_ij = 0.5*((w_i+w_j)*log(det(p_ij))-w_i*log(det(p_i))-w_j*log(det(p_j))) ;
% B_ij = abs(B_ij) ;