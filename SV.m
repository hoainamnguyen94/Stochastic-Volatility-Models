function thetahat = SV(y,T,nloop,burnin)

phih0 = .97; Vphih = .1^2;
Vmu = 10;
muh0 = 1; Vmuh = 10;
nuh = 5; Sh = .2^2*(nuh-1);

muh = log(var(y)); phih = .98; 
omegah2 = .2^2;
h = muh + sqrt(omegah2)*randn(T,1);
exph = exp(h);

store_theta = zeros(nloop - burnin,4);
store_h = zeros(nloop - burnin,T);

Hphi = speye(T) - sparse(2:T,1:(T-1),phih*ones(1,T-1),T,T);
newnuh = T/2 + nuh;
counth = 0; countphi = 0;

randn('seed',sum(clock*100)); rand('seed',sum(clock*1000));  %#ok<RAND>
disp('Starting SV.... ');
disp(' ' );

start_time = clock;

for loop = 1:nloop 
    
    %% sample mu    
    invDmu = 1/Vmu + sum(1./exph);
    muhat = invDmu\sum(y./exph);   
    mu = muhat + chol(invDmu,'lower')'\randn;  
    
    %% sample h    
    HiSH = Hphi'*sparse(1:T,1:T,[(1-phih^2)/omegah2; 1/omegah2*ones(T-1,1)])*Hphi;
    deltah = Hphi\[muh; muh*(1-phih)*ones(T-1,1)];
    HiSHdeltah = HiSH*deltah;
    s2 = (y-mu).^2;
    errh = 1; ht = h;
    while errh> 10^(-3)
        expht = exp(ht);
        sinvexpht = s2./expht;        
        fh = -.5 + .5*sinvexpht;
        Gh = .5*sinvexpht;
        Kh = HiSH + sparse(1:T,1:T,Gh);
        newht = Kh\(fh+Gh.*ht+HiSHdeltah);
        errh = max(abs(newht-ht));
        ht = newht;          
    end 
    cholHh = chol(Kh,'lower');
    %AR-step: 
    lph = @(x) -.5*(x-deltah)'*HiSH*(x-deltah) -.5*sum(x) -.5*exp(-x)'*s2;
    logc = lph(ht) + log(3);    
    flag = 0;
    while flag == 0
        hc = ht + cholHh'\randn(T,1);        
        alpARc =  lph(hc) + .5*(hc-ht)'*Kh*(hc-ht) - logc;
        if alpARc > log(rand)
            flag = 1;
        end
    end        
    %MH-step    
    alpAR = lph(h) + .5*(h-ht)'*Kh*(h-ht) - logc;
    if alpAR < 0
        alpMH = 1;
    elseif alpARc < 0
        alpMH = - alpAR;
    else
        alpMH = alpARc - alpAR;
    end    
    if alpMH > log(rand) || loop == 1
        h = hc;
        exph = exp(h);
        counth = counth + 1;
    end 
    
    %% sample omegah2
    errh = [(h(1)-muh)*sqrt(1-phih^2);h(2:end)-phih*h(1:end-1)-muh*(1-phih)];    
    newSh = Sh+sum(errh.^2)/2;    
    omegah2 = 1/gamrnd(newnuh,1./newSh);  
    
    %% sample phih
    Xphi = h(1:end-1) - muh;
    yphi = h(2:end) - muh;
    Dphi = 1/(1/Vphih + Xphi'*Xphi/omegah2);
    phihat = Dphi*(phih0/Vphih + Xphi'*yphi/omegah2);
    phic = phihat + sqrt(Dphi)*randn;
    g = @(x) -.5*log(omegah2./(1-x.^2))-.5*(1-x.^2)/omegah2*(h(1)-muh)^2;
    if abs(phic)<.9999
        alpMH = exp(g(phic)-g(phih));
        if alpMH>rand
            phih = phic;
            countphi = countphi+1;
            Hphi = speye(T) - sparse(2:T,1:(T-1),phih*ones(1,T-1),T,T);
        end
    end    
    
    %% sample muh    
    Dmuh = 1/(1/Vmuh + ((T-1)*(1-phih)^2 + (1-phih^2))/omegah2);
    muhhat = Dmuh*(muh0/Vmuh + (1-phih^2)/omegah2*h(1) + ...
             (1-phih)/omegah2*sum(h(2:end)-phih*h(1:end-1)));
    muh = muhhat + sqrt(Dmuh)*randn;   

    if loop>burnin
        i = loop-burnin;        
        store_h(i,:) = h';         
        store_theta(i,:) =  [mu muh phih omegah2]; 
    end    
    
    if mod(loop,5000)==0
        disp([ num2str( loop ) ' loops... '])
    end        
end

disp(['MCMC takes '  num2str(etime(clock,start_time)) ' seconds']);
disp(' ');

thetahat = mean(store_theta)';

end



