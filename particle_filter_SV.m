function [w_T,part_T] = particle_filter_SV(y,thetahat,N,Nth) 

%Model parameterization to be used in particle filtering
%y_t = exp(x_t/2)*eps_t
%x_t = alpha+beta*x_{t-1}+gamma*eta_t

%Model parameterization when estimating thetahat (SV model)
%y_t = mu+eps_y_t, eps_y_t ~ N(0,e^h_t)
%h_t = mu_h+phi_h*(h_{t-1}-mu_h)+eps_h_t, eps_h_t ~ N(0,omega_h^2)
%thetahat = [mu mu_h phi_h omega2_h]

obs = y-thetahat(1);
alpha = thetahat(2)-thetahat(3)*thetahat(2);
beta = thetahat(3);
gamma = sqrt(thetahat(4));

T = length(obs);
 
part = zeros(T,N); 
w = part; 
Neff = zeros(T,1); 
 
part(1,:) = alpha/(1-beta)+gamma/sqrt(1-beta^2)*randn(1,N); 
wei = (1./exp(part(1,:)/2)).*exp(-.5*(obs(1).^2)./exp(part(1,:))); 
w(1,:) = wei/sum(wei); 
Neff(1) = 1/sum(w(1,:).^2); 
 
for t = 2:T 
    if Neff(t-1) <= Nth 
        resampidx = randsample(N,N,true,w(t-1,:)); 
        part(t,:) = alpha+beta*part(t-1,resampidx)+gamma*randn(1,N); 
        wei = (1./exp(part(t,:)/2)).*exp(-.5*(obs(t).^2)./exp(part(t,:))); 
    else 
        part(t,:) = alpha+beta*part(t-1,:)+gamma*randn(1,N); 
        wei = w(t-1,:).*(1./exp(part(t,:)/2)).*exp(-.5*(obs(t).^2)./exp(part(t,:))); 
    end 
    w(t,:) = wei/sum(wei); 
    Neff(t) = 1/sum(w(t,:).^2); 
end

w_T = w(T,:);
part_T = part(T,:);

end
