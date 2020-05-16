function [w_T,part_T] = particle_filter_SV_L(y,thetahat,N,Nth) 

%Model parameterization to be used in particle filtering
%y_t = exp(x_t/2)*eps_t
%x_t = alpha+beta*x_{t-1}+gamma*eta_t

%Model parameterization when estimating thetahat (SV-L model)
%y_t = mu+eps_y_t
%h_t = mu_h+phi_h*(h_{t-1}-mu_h)+eps_h_t
%(eps_y_t,eps_h_t) follows a bivariate normal distribution
%    mean (0,0)
%    covariance matrix  (e^h_t                   rho*e^(h_t/2)*omega_h^2)
%                       (rho*e^(h_t/2)*omega_h^2               omega_h^2)
%thetahat = [mu rho mu_h phi_h omega2_h]

%The model above is equivalent to
%y_t = mu+eps_t*e^(h_t/2)
%h_t = mu_h*(1-phi_h)+phi_h*h_{t-1}+omega_h*eta_t
%(eps_t,eta_t) follows a bivariate normal distribution
%    mean (0,0)
%    covariance matrix  (1   rho)
%                       (rho   1)
%thetahat = [mu rho mu_h phi_h omega2_h]
%We can write eta_t = rho*eps_t + sqrt(1-rho^2)*N(0,1)

obs = y-thetahat(1);
alpha = thetahat(3)-thetahat(4)*thetahat(3);
beta = thetahat(4);
gamma = sqrt(thetahat(5));
rho = thetahat(2);

T = length(obs);
 
part = zeros(T,N); 
w = part; 
Neff = zeros(T,1); 
 
part(1,:) = alpha/(1-beta)+gamma/sqrt(1-beta^2)*(sqrt(1-(rho^2))*randn(1,N)); 
wei = (1./exp(part(1,:)/2)).*exp(-.5*(obs(1).^2)./exp(part(1,:))); 
w(1,:) = wei/sum(wei); 
Neff(1) = 1/sum(w(1,:).^2); 
 
for t = 2:T 
    if Neff(t-1) <= Nth 
        resampidx = randsample(N,N,true,w(t-1,:)); 
        part(t,:) = alpha+beta*part(t-1,resampidx)+gamma* ...
                         (rho*obs(t-1)*exp(-0.5*part(t-1,resampidx))+ ...
                          sqrt(1-(rho^2))*randn(1,N)); 
        wei = (1./exp(part(t,:)/2)).*exp(-.5*(obs(t).^2)./exp(part(t,:))); 
    else 
        part(t,:) = alpha+beta*part(t-1,:)+gamma* ...
                         (rho*obs(t-1)*exp(-0.5*part(t-1,:))+ ...
                          sqrt(1-(rho^2))*randn(1,N)); 
        wei = w(t-1,:).*(1./exp(part(t,:)/2)).*exp(-.5*(obs(t).^2)./exp(part(t,:))); 
    end 
    w(t,:) = wei/sum(wei); 
    Neff(t) = 1/sum(w(t,:).^2); 
end

w_T = w(T,:);
part_T = part(T,:);

end
