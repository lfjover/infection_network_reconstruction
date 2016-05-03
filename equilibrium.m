function x = equilibrium(M,r,a,K,phi,beta,m)
% equilibrium(M,r,a,K,phi,beta,m)

[nH,nV] = size(M);
cMatrix = [ a.*repmat(r,1,nH)/K, M.*phi; (M.*phi.*beta)', zeros(nV)];
x = cMatrix\[r;m];

end