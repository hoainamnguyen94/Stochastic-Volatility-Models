%select a market
%1-JCI 2-KLSE 3-PCOMP 4-SET 5-STI 6-VNI

for market = 1:6

switch market
    case 1
        %read the data for a market
        y = xlsread('ASEAN_RET.xlsx','JCI');
        start_1 = 1214;
        filename = 'JCI_SV.csv';
        fprintf('running SV for JCI \n\n');
    case 2
        y = xlsread('ASEAN_RET.xlsx','KLSE');
        %in-sample period: 01/07/2006 - 30/06/2011
        %out-of-sample period: 01/07/2011 - 30/06/2017
        %start_1 represents the index of the last in-sample observation
        %start to forecast 1 day ahead from start_1
        start_1 = 1236;
        filename = 'KLSE_SV.csv';
        fprintf('running SV for KLSE \n\n');
    case 3
        y = xlsread('ASEAN_RET.xlsx','PCOMP');
        start_1 = 1222;
        %set the name of the file to which the result will be written
        filename = 'PCOMP_SV.csv';
        fprintf('running SV for PCOMP \n\n');
    case 4
        y = xlsread('ASEAN_RET.xlsx','SET');
        start_1 = 1218;
        filename = 'SET_SV.csv';
        fprintf('running SV for SET \n\n');
    case 5
        y = xlsread('ASEAN_RET.xlsx','STI');
        start_1 = 1257;
        filename = 'STI_SV.csv';
        fprintf('running SV for STI \n\n');
    case 6
        y = xlsread('ASEAN_RET.xlsx','VNI');
        start_1 = 1240;
        filename = 'VNI_SV.csv';
        fprintf('running SV for VNI \n\n');
end
y = 100*y; %multiply by 100
T_total = length(y); %T_total represents the total number of data points
start_5 = start_1-4; %start to forecast 5 days ahead from start_5
start_20 = start_1-19; %start to forecast 20 days ahead from start_20

%set the parameters to be used in the estimation
nloop = 11000;
burnin = 1000;
%set the parameters to be used in particle filtering
N = 1000;
Nth = 667;

%estimate the parameters using in-sample observations from 1 to start_5
%and forecast 5 days ahead
thetahat = SV(y(1:start_5),start_5,nloop,burnin);
%the time is incremented by 1 up to start_5+3 (not start_5+4 as usual)
%the reason is to have synchronization with forecasting 1 day ahead
%each time we forecast 5 days ahead
result_5_ini = zeros(4,6); %create a matrix to hold the estimated VaR
for i = 0:3
    [w_T,part_T] = particle_filter_SV(y(1:(start_5+i)),thetahat,N,Nth);
    q_VAR = f_SV(thetahat,w_T,part_T,N);
    result_5_ini(i+1,:) = q_VAR{2};
end

%we do the same thing for start_20
%again, we stop at start_20+18 (not start_20+19) to obtain synchronization
result_20_ini = zeros(19,6);
for i = 0:18
    if mod(i,5)==0
        train = y(1:(start_20+i));
        T_train = length(train);
        thetahat = SV(train,T_train,nloop,burnin);
    end
    [w_T,part_T] = particle_filter_SV(y(1:(start_20+i)),thetahat,N,Nth);
    q_VAR = f_SV(thetahat,w_T,part_T,N);
    result_20_ini(i+1,:) = q_VAR{3};
end

%I divide the loops into blocks of 5 (the model is estimated every 5 days)
%suppose that (s, s+1, s+2, s+3, s+4) is a block
%the parameters are estimated at time s 
                  %and we forecast VaR at s+1, s+5 and s+20
%at time s+1, we skip the estimation and forecast VaR at s+2, s+6 and s+21
%...
%therefore, the iterations within the same block ARE NOT INDEPENDENT
%but the blocks ARE INDEPENDENT
%so parallel processing can be used to carry out computations
                  %on multiple blocks at the same time
iter = fix((T_total-1-start_1)/5);

%create a matrix to hold the estimated VaR (1 day ahead)
result_1 = zeros(iter+1,30);
%create a matrix to hold the estimated VaR (5 days ahead)
result_5 = zeros(iter+1,30); 
%create a matrix to hold the estimated VaR (20 days ahead)
result_20 = zeros(iter+1,30);

parfor i = 0:(iter-1) %iterate over the blocks  
    
    T = start_1+5*i; %T represents the current time

    %include all the data points up to T
    train = y(1:T); %#ok<*PFBNS>
    T_train = length(train); 
    %estimate the parameters
    thetahat = SV(train,T_train,nloop,burnin);

    %forecast the VaR and obtain the quantiles
    %the current time is incremented by 1
    [w_T,part_T] = particle_filter_SV(y(1:T),thetahat,N,Nth);
    q_VAR_1 = f_SV(thetahat,w_T,part_T,N);
    [w_T,part_T] = particle_filter_SV(y(1:(T+1)),thetahat,N,Nth);
    q_VAR_2 = f_SV(thetahat,w_T,part_T,N);
    [w_T,part_T] = particle_filter_SV(y(1:(T+2)),thetahat,N,Nth);
    q_VAR_3 = f_SV(thetahat,w_T,part_T,N);
    [w_T,part_T] = particle_filter_SV(y(1:(T+3)),thetahat,N,Nth);
    q_VAR_4 = f_SV(thetahat,w_T,part_T,N);
    [w_T,part_T] = particle_filter_SV(y(1:(T+4)),thetahat,N,Nth);
    q_VAR_5 = f_SV(thetahat,w_T,part_T,N);
    
    %store the quantiles in the pre-allocated matrices
    %note that all the quantiles must be stored on the same row
    %otherwise the parfor loop will flag an error 
    %the matrices will need to be reshaped later before writing to a csv file
    result_1(i+1,:) = [q_VAR_1{1} q_VAR_2{1} q_VAR_3{1} q_VAR_4{1} ...
                       q_VAR_5{1}]
    result_5(i+1,:) = [q_VAR_1{2} q_VAR_2{2} q_VAR_3{2} q_VAR_4{2} ...
                       q_VAR_5{2}]
    result_20(i+1,:) = [q_VAR_1{3} q_VAR_2{3} q_VAR_3{3} q_VAR_4{3} ...
                        q_VAR_5{3}]
                    
end

%stop using parallel processing after the parfor loop above
delete(gcp('nocreate'));

%the last block i = iter may not contain 5 iterations
%although it is possible to incorporate this case to the parfor loop above
                  %by using conditional statements
%this may slow down the loop because the condition always holds true
                  %except for the last block
%so I decide to treat it separately
T = start_1+5*iter;
train = y(1:T);
T_train = length(train);
thetahat = SV(train,T_train,nloop,burnin);
for i = 0:(T_total-1-T)
    [w_T,part_T] = particle_filter_SV(y(1:(T+i)),thetahat,N,Nth);
    q_VAR = f_SV(thetahat,w_T,part_T,N);
    start_point = 6*i+1;
    end_point = 6*(i+1);
    result_1(iter+1,start_point:end_point) = q_VAR{1};
    result_5(iter+1,start_point:end_point) = q_VAR{2};
    result_20(iter+1,start_point:end_point) = q_VAR{3};
end

%reshape the result into a more convenient form
reshaped_mat_1 = reshaping(result_1,iter);
reshaped_mat_5 = reshaping(result_5,iter);
reshaped_mat_20 = reshaping(result_20,iter);
%combine result_5_ini and result_20_ini with reshaped_mat_5 and
%reshaped_mat_20 to produce the final results
len = T_total-start_1;
final_mat_1 = reshaped_mat_1(1:len,:);
%note that we have to remove some values at the end of reshaped_mat_5 and
%reshaped_mat_20 because the predicted quantiles are outside our sample
final_mat_5 = [result_5_ini;reshaped_mat_5(1:(len-4),:)];
final_mat_20 = [result_20_ini;reshaped_mat_20(1:(len-19),:)];

%create a vector to hold the time T of the last observation
time = transpose(start_1:(T_total-1));
%append the times to the reshaped matrix
final_mat = [time final_mat_1 final_mat_5 final_mat_20];

%create the headers
header = {'T','1.00%','2.50%','5.00%','99.00%','97.50%','95.00%', ...
              '1.00%','2.50%','5.00%','99.00%','97.50%','95.00%', ...
              '1.00%','2.50%','5.00%','99.00%','97.50%','95.00%'};
%write the final matrix with time and header to a csv file
csvwrite_with_headers(filename,final_mat,header);

end
