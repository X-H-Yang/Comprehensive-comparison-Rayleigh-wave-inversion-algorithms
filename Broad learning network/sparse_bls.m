function wk = sparse_bls(A,b,lam,itrs)
%
% note: built upon previously released codes (https://broadlearning.ai/)
%
% Chen and Liu 2017: Chen, C. P., & Liu, Z. (2017). Broad learning system: 
%                    An effective and efficient incremental learning system
%                    without the need for deep architecture. IEEE 
%                    transactions on neural networks and learning systems, 
%                    29(1), 10-24.
%
% Chen et al. 2018: Chen, C. P., Liu, Z., & Feng, S. (2018). Universal 
%                   approximation capability of broad learning system and 
%                   its structural variations. IEEE transactions on neural 
%                   networks and learning systems, 30(4), 1191-1204.
%
AA = (A') * A;
m = size(A,2);
n = size(b,2);
x = zeros(m,n);
wk = x; 
ok=x;uk=x;
L1=eye(m)/(AA+eye(m));
L2=L1*A'*b;

for i = 1:itrs,
    tempc=ok-uk;
  ck =  L2+L1*tempc;
 ok=shrinkage(ck+uk, lam);
 uk=uk+(ck-ok);
 wk=ok;
end
end
function z = shrinkage(x, kappa)
    z = max( x - kappa,0 ) - max( -x - kappa ,0);
end


