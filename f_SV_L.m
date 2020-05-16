function q_VAR = f_SV_L(y_end,thetahat,w_T,part_T,N)

%thetahat = [mu rho muh phih omega2_h]
alpha = thetahat(3)-thetahat(4)*thetahat(3);
beta = thetahat(4);
gamma = sqrt(thetahat(5));
rho = thetahat(2);

part = [part_T;zeros(20,N)];
yhat = [y_end*ones(1,N);zeros(20,N)];

for t = 2:21
    part(t,:) = alpha+beta*part(t-1,:)+gamma* ...
               (rho*(yhat(t-1,:)-thetahat(1)).*exp(-0.5*part(t-1,:))+ ...
                sqrt(1-(rho^2))*randn(1,N));
    resampidx = randsample(N,N,true,w_T); 
    part(t,:) = part(t,resampidx);
    yhat(t,:) = normrnd(thetahat(1)*ones(1,N),exp(part(t,:)/2));
end

q_VAR_T1 = quantile(yhat(2,:),[0.01 0.025 0.05 0.99 0.975 0.95]);
q_VAR_T5 = quantile(yhat(6,:),[0.01 0.025 0.05 0.99 0.975 0.95]);
q_VAR_T20 = quantile(yhat(21,:),[0.01 0.025 0.05 0.99 0.975 0.95]);

q_VAR = {q_VAR_T1 q_VAR_T5 q_VAR_T20};

end