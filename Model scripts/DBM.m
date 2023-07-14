function [ER] = DBM(R,PET)   % In this model, lag of one day has considered. see the part of model discharge calculation. where i+'1' represents the lag
% lag=0 ;
%% all parameter values
W=zeros(1,1);
H=zeros(1,1);

len = size(R,1);
% Calculate W
for iii= 1:len
    for j=1:size(R,2)
        if R(iii,j)> (PET(iii,1)/size(R,2))
            W(iii,j) = R(iii,j)-(PET(iii,1)/size(R,2));
        else
            W(iii,j) = 0;
        end
    end
end
% Calculate H
for i= 1:len
    for j=1:size(R,2)
        if R(i,j)< (PET(i,1)/size(R,2))
            H(i,j) = (PET(i,1)/size(R,2))-(R(i,j));
        else
            H(i,j) = 0;
        end
    end
end

% Calculate FW
FW = zeros(len,1);
for i = 1:len
    if i<366
        FW(i,1)=  nan;
    else
        count = 0;
        for j = 1:365
            FW(i,1)= FW(i,1)+ sum(W(i-j+1,:))/(1+count*0.4);
            count = count + 1;
        end
    end
end
%---------------------------**Important**----------------------------------------------
% in some cases FW is coming as zero. so the 'phi' value which is ratio
% of 'FH' by 'FW'  comes as infinity. so zero of FW by 0.00001 to avoid
% the nan value.
FW(FW==0)=0.0000001 ;
%-----------------------------END of-----------------------------------
% Calculate FH
FH = zeros(len,1);
for i = 1:len
    if i<366
        FH(i,1)=  nan;
    else
        count = 0;
        for j = 1:365
            FH(i,1)= FH(i,1)+ sum(H(i-j+1,:))/(1+count*0.4);
            count = count + 1;
        end
    end
end
% Calculate FW/FW+FH
FW_ratio = FW./(FW+FH);
FH_ratio = FH./(FW+FH);

% % Calculate Discharge
% 
% model_disch = zeros(len,1);
% Calculate
phi = FH./FW;
IR = zeros(len,size(R,2));
model_disch(1:366,1) = NaN ;
for i = 366:len
    for j=1:size(R,2)
        ER(i,j) = W(i,j)* (1-((phi(i,1)*tanh(phi(i,1)^(-1)))*(1-exp(-phi(i,1))))^0.5);
    end
end


