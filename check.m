%a quick check for the coverage
y = xlsread('ASEAN_RET.xlsx','VNI'); %change the worksheet
y = y*100; %multiply by 100
filename = 'VNI_SV_t.csv'; %change the filename
data = csvread(filename,1,0);
threshold = 1240; %change the threshold
%JCI - 1214 KLSE - 1236 PCOMP - 1222 SET - 1218 STI - 1257 VNI - 1240
num_pred = length(y)-threshold;
q_VAR_1d_1 = data(:,2);
coverage = sum(q_VAR_1d_1<y((end-num_pred+1):end))/num_pred;

