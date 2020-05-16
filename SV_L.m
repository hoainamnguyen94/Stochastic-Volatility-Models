function thetahat = SV_L(y,T,nloop,burnin)

phih0 = .97; Vphih = .1^2;
mu0 = 0; Vmu = 10;
muh0 = 1; Vmuh = 10;
nuh = 5; Sh = .2^2*(nuh-1);
rho0 = 0; Vrho = 1;
lrhopri = @(x) -.5*(x-rho0).^2/Vrho;

mu = mean(y); %#ok<*NASGU>
phih = .98;
omegah2 = .2;
rho = -.5;
muh = -10;
Hphi = speye(T+1) - phih*sparse(2:T+1,1:T,ones(1,T),T+1,T+1);
HiSH = Hphi'*spdiags([(1-phih^2)/omegah2; 1/omegah2*ones(T,1)],0,T+1,T+1)*Hphi;
deltah = Hphi\[muh;(1-phih)*muh*ones(T,1)]; 
h = deltah + chol(HiSH,'lower')'\randn(T+1,1)/10;
y_temp = y;
y_temp(y_temp==0) = 0.0001;
h = (h + log([y_temp;y_temp(end)].^2))/2;

store_theta = zeros(nloop - burnin,5);
store_h = zeros(nloop - burnin,T+1);
counth = 0;
countphih = 0;
countomegah2 = 0;

disp('Starting SV-L.... ');
disp(' ' );

rand('state', sum(100*clock)); randn('state', sum(200*clock)); %#ok<*RAND>
start_time = clock;

for loop = 1:nloop  
    
    %% sample mu
    eh = h(2:end)-phih*h(1:end-1)-(1-phih)*muh;
    iexph = exp(-h(1:T));
    Dmu = 1/(1/Vmu + 1/(1-rho^2)*sum(iexph));
    muhat = Dmu*(mu0/Vmu + ...
        1/(1-rho^2)*iexph'*(y-rho/sqrt(omegah2)*exp(h(1:T)/2).*eh));    
    mu = muhat + sqrt(Dmu)*randn;    
     
    %% sample h  
    HiSH = Hphi'*spdiags([(1-phih^2)/omegah2; 1/omegah2*ones(T,1)],0,T+1,T+1)*Hphi;
    deltah = Hphi\[muh;(1-phih)*muh*ones(T,1)]; 
    errh = 1; ht = h;
    c1 = 1/(1-rho^2);
    c2 = rho/sqrt(omegah2);
    c3 = c2^2;        
    while errh > 10^(-3)
        tmp1 = (y-mu)./exp(ht(1:T)/2); 
        tmp2 = tmp1.^2;
        eh = ht(2:end)-phih*ht(1:end-1)-(1-phih)*muh;
        fh1 = -.5 - c1/2*(-tmp2 -2*c3*phih*eh + c2*tmp1.*(eh+2*phih));        
        fh2 = c1*c2*(tmp1 - c2*eh);        
        fh = [fh1;0] + [0;fh2];    
        Gh1 = c1/2*(tmp2 + 2*c3*phih^2 - c2/2*tmp1.*(eh + 4*phih));
        Gh2 = c1*c3*ones(T,1);
        Gh3 = -c1*c2*(c2*phih - .5*tmp1);
        Gh = sparse(1:T+1,1:T+1,[Gh1;0]+[0;Gh2]) ...
            + sparse(2:T+1,1:T,Gh3,T+1,T+1) + sparse(1:T,2:T+1,Gh3,T+1,T+1);      
        Kh = HiSH + Gh;
        newht = Kh\(HiSH*deltah + fh + Gh*ht); 
        errh = max(abs(newht-ht));
        ht = newht;
    end    
    % AR-step
    logc = like_sv_l(y,ht,mu,rho,muh,phih,omegah2,0) + log(3);
    CinvDh = chol(Kh,'lower');
    flag = 0;
    while flag == 0
        hc = ht + CinvDh'\randn(T+1,1);      
        alpARc = like_sv_l(y,hc,mu,rho,muh,phih,omegah2,0) + ...
            + .5*(hc-ht)'*Kh*(hc-ht) - logc;
        if alpARc > log(rand)
            flag = 1;
        end
    end       
    % MH-step    
    alpAR = like_sv_l(y,h,mu,rho,muh,phih,omegah2,0) + ...
        .5*(h-ht)'*Kh*(h-ht) - logc;
    if alpAR < 0 
        alpMH = 1;
    elseif alpARc < 0
        alpMH = - alpAR;
    else
        alpMH = alpARc - alpAR;
    end        
    if alpMH > log(rand) || loop<10
        h = hc;
        counth = counth + 1;    
    end
    
    %% sample rho - griddy Gibbs
    eh = h(2:end)-phih*h(1:end-1)-(1-phih)*muh;
    tmprho = exp(-h(1:T)/2).*(y-mu);
    k1 = tmprho'*tmprho;
    k2 = tmprho'*eh;
    k3 = eh'*eh; 
    grho = @(x) lrhopri(x) -T/2*log(1-x.^2) + ...
        - (k1 - 2*x/sqrt(omegah2)*k2 + x.^2/omegah2*k3)./(2*(1-x.^2));    
    rhogrid = linspace(-1+rand/100,1-rand/100,300)';
    rhopdf = grho(rhogrid);
    rhopdf = exp(rhopdf-max(rhopdf));
    rhocdf = cumsum(rhopdf);
    rhocdf = rhocdf/rhocdf(end);    
    rho = rhogrid(find(rhocdf>rand,1));
 
    %% sample omegah2
    eh = h(2:end)-phih*h(1:end-1)-(1-phih)*muh;
    eh_ss = eh'*eh;
    newSh = Sh + ((1-phih^2)*(h(1)-muh)^2 + eh_ss)/2;
    omegah2c = 1/gamrnd(nuh+(T+1)/2, 1/newSh);
    k1 = .5*rho^2/(1-rho^2)*eh_ss;
    k2 = rho/(1-rho^2)*eh'*((y-mu)./exp(h(1:T)/2));
    goh2 = @(x) -k1/x + k2/sqrt(x);
    alpMH = goh2(omegah2c) - goh2(omegah2);
    if alpMH > log(rand)
        omegah2 = omegah2c;
        countomegah2 = countomegah2 + 1;
    end
    
    %% sample phih
    tmph = h-muh;
    omegah = sqrt(omegah2);
    Dphih = 1/(1/Vphih + sum(tmph(1:T).^2)/(omegah2*(1-rho^2)));
    phihhat = Dphih*(phih0/Vphih + tmph(1:T)'*tmph(2:T+1)/omegah2 ...
        - rho/((1-rho^2)*omegah)*tmph(1:T)'*((y-mu)./exp(h(1:T)/2)-rho/omegah*tmph(2:T+1)));    
    phihc = phihhat + sqrt(Dphih)*randn;    
    gphih = @(x) -.5*log(omegah2./(1-x.^2)) -.5*(1-x.^2)/omegah2*(h(1)-muh)^2;
    if abs(phihc)<.999
        alpMH = gphih(phihc)-gphih(phih);
        if alpMH > log(rand)
            phih = phihc;
            countphih = countphih + 1;
        end
    end 
    Hphi = speye(T+1) - phih*sparse(2:T+1,1:T,ones(1,T),T+1,T+1);
    
    %% sample muh  
    tmph = h(2:end)-phih*h(1:end-1);
    omegah = sqrt(omegah2);
    Dmuh = 1/(1/Vmuh + T*(1-phih)^2/omegah2 + (1-phih^2)/omegah2 ...
        + T*(1-phih)^2*rho^2/(omegah2*(1-rho^2)));
    muhhat = Dmuh*(muh0/Vmuh + (1-phih^2)/omegah2*h(1) + (1-phih)/omegah2*sum(tmph) ...
        - rho*(1-phih)/(omegah*(1-rho^2))*sum((y-mu)./exp(h(1:T)/2) - rho/omegah*tmph));
    muh = muhhat + sqrt(Dmuh)*randn;

    if ( mod( loop, 2000 ) ==0 )
        disp(  [ num2str( loop ) ' loops... ' ] )
    end     
    
    if loop>burnin
        i = loop-burnin;
        store_h(i,:) = h';         
        store_theta(i,:) =  [mu rho muh phih omegah2];         
    end    
end

disp( ['MCMC takes '  num2str( etime( clock, start_time) ) ' seconds' ] );
disp(' ' );

hhat = mean(exp(store_h/2))';  %% plot std dev
thetahat = mean(store_theta)';

end

