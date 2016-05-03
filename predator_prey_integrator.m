function [t,x] = predator_prey_integrator(M,r,a,K,phi,beta,m,x0,tfinal)
%integrates the mutualist system presented in ode_mutualism
%output: time series x = [P;A]
%M: interaction matrix
    [nH,nV] = size(M);    %number of Plants and animals

    options = odeset('RelTol',10^-8); 
    [t,x] = ode45(@ode_predator_prey,[0 tfinal],x0,options);

%%
    function xDot = ode_predator_prey(t,x)
    %nestedmodelmatrix
            H= x(1:nH);
            V = x(nH+1:end);

            Hdot = H.*(r.*(1 - a*H./K) - (M.*phi)*V);
            Vdot = V.*( (M.*phi.*beta)'*H - m);

            xDot = [Hdot;Vdot];
    end
end