function q_VAR = f_SV(thetahat,w_T,part_T,N)

%thetahat = [mu muh phih omega2_h]
alpha = thetahat(2)-thetahat(3)*thetahat(2);
beta = thetahat(3);
gamma = sqrt(thetahat(4));

part = [part_T;zeros(20,N)];

for t = 2:21
    part(t,:) = alpha+beta*part(t-1,:)+gamma*randn(1,N);
end

resampidx = randsample(N,N,true,w_T); 
part(2,:) = part(2,resampidx); 
yhat_T1 = normrnd(thetahat(1)*ones(1,N),exp(part(2,:)/2));
q_VAR_T1 = quantile(yhat_T1,[0.01 0.025 0.05 0.99 0.975 0.95]);

resampidx = randsample(N,N,true,w_T); 
part(6,:) = part(6,resampidx); 
yhat_T5 = normrnd(thetahat(1)*ones(1,N),exp(part(6,:)/2));
q_VAR_T5 = quantile(yhat_T5,[0.01 0.025 0.05 0.99 0.975 0.95]);

resampidx = randsample(N,N,true,w_T); 
part(21,:) = part(21,resampidx); 
yhat_T20 = normrnd(thetahat(1)*ones(1,N),exp(part(21,:)/2));
q_VAR_T20 = quantile(yhat_T20,[0.01 0.025 0.05 0.99 0.975 0.95]);

q_VAR = {q_VAR_T1 q_VAR_T5 q_VAR_T20};

end