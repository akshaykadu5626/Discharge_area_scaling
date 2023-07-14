clear all
clc
tic
load Ken_parameters.mat %Optimal parameters for Kentucky basin
load Data_kenall.mat
nos_basins= 8; %number of basins
parameters_kenall=opt_paramgh;
[basin_area name]=xlsread('BasinArea.xlsx');
rs= 24; %Model Temporal Resolution in hours
klr=1
for klr=1:nos_basins
    klr
    parameters=parameters_kenall(klr,:);
    basin_id=Data_kenall(klr,:);
    P=basin_id{1,1}; %1/10/1975 to 30/09/2002
    T=basin_id{1,2};
    PET=basin_id{1,3};
    Q_d=basin_id{1,4};
    VLU=basin_id{1,5};
    VLD=basin_id{1,6};
    VFD=basin_id{1,7};
    parTT = parameters(1);
    parCFMAX = parameters(2);
    parSFCF = parameters(3);
    parCFR = parameters(4);
    parCWH = parameters(5);
    vs = parameters(6);
    vd= parameters(7);
    alpha= parameters(8);
    parPCORR = 1;
    %% ER generation using DB model
    %% Creating packs with all nan values and initial values as 0.0001
    SNOWPACK = nan(size(P));
    SNOWPACK(1,:) = 0.0001;
    MELTWATER = nan(size(P));
    MELTWATER(1,:) = 0.0001;
    temperary=nan(size(P)) ;   % added by me
    temperary(1,:)=0.0001 ;
    %%  Snow estimation
    P(:,1) = P(:,1).*parPCORR;  % multiplying by ppt correction factor. here it is one
    SNOW=zeros(1,1) ;  % creaing empty snow matrix
    RAIN=zeros(1,1) ;  % creaitng empty rain matrix
    for mm=1:length(P)
        if T(mm,1)>=parTT
            SNOW(mm,1)=0 ;
        else
            SNOW(mm,1)=P(mm,1) ;
        end
    end
    SNOW(:,1) = SNOW(:,1).*parSFCF ; % multiplying with parSFCF value to snow
    %  for rain
    for mm=1:length(P)
        if T(mm,1)<parTT
            RAIN(mm,1)=0 ;
        else
            RAIN(mm,1)=P(mm,1) ;
        end
    end
    %% Boxes
    for t = 2:length(P)
        % Snow
        SNOWPACK(t,:) = SNOWPACK(t-1,:)+SNOW(t,:);
        melt = parCFMAX .* (T(t,:)-parTT);
        melt(melt<0) = 0;
        melt = min(melt,SNOWPACK(t,:));
        MELTWATER(t,:) = MELTWATER(t-1,:)+melt;
        SNOWPACK(t,:) = SNOWPACK(t,:)-melt;
        refreezing = parCFR .* parCFMAX .* (parTT-T(t,:));
        refreezing(refreezing<0) = 0;
        refreezing = min(refreezing,MELTWATER(t,:));
        SNOWPACK(t,:) = SNOWPACK(t,:)+refreezing;
        MELTWATER(t,:) = MELTWATER(t,:)-refreezing;
        tosoil = MELTWATER(t,:) - (parCWH .* SNOWPACK(t,:));
        tosoil(tosoil<0) = 0;
        temperary(t,:)=tosoil ;
        MELTWATER(t,:) = MELTWATER(t,:)-tosoil;
        % final rain
    end
    R=temperary+RAIN ;
    ER=DBM(R,PET) ; % DB Model
    IPSF=alpha*ER;
    IMSSF=(1-alpha)*ER;
    %%
    timeI=zeros(1,1);
    timeF=zeros(1,1);
    dimrain=size(ER,1);
    dimluld=size(VLU,1);
    for k=1:dimluld
        timeI(k)= VLD(k)/vs/rs;
        timeF(k)= (VLU(k)+VLD(k))/vs/rs+ VLU(k)/vd/rs;
    end
    timeImax=max(timeI);
    timeFmax=max(timeF);
    sizeBF=floor(timeFmax+2);
    sizeWF=floor(timeImax+2);
    CWF=zeros(sizeWF,1);
    PSFIUH=zeros(sizeBF,1);
    MSSFIUH= zeros(sizeBF,1);
    %% PSFIUH computation
    tm=0;
    for i=1:sizeWF      %calculating cumulative width function
        for k=1:dimluld
            if timeI(k)>tm
                temp2= floor(VFD(k)/2);
                temp2=temp2*2;   %temp2 is an integer type variable.
                if(temp2~= VFD(k))
                    CWF(i)= CWF(i)+1;    %If VFD[k] is even
                else
                    CWF(i)= CWF(i)+sqrt(2);       %If VFD[k] is odd
                end
            end
        end
        tm=tm+1;
    end
    
    for i=1:sizeWF-1 %calculating  WFIUH from CWF
        PSFIUH(i)= (CWF(i) - CWF(i+1))/CWF(1);
    end
    %% MSSF Inst. unit hydrograph computation
    tm=1;
    for i=1:sizeBF
        for k=1:dimluld
            if (timeI(k)<tm && timeF(k)>tm)
                temp2= floor(VFD(k)/2);
                temp2=temp2*2;   %temp2 is an integer type variable.
                if(temp2~= VFD(k))
                    MSSFIUH(i)= MSSFIUH(i)+1;    %If VFD[k] is even
                else
                    MSSFIUH(i)= MSSFIUH(i)+sqrt(2);       %If VFD[k] is odd
                end
            end
        end
        tm=tm+1;
    end
    sumBF= sum(MSSFIUH);
    MSSFIUH=MSSFIUH/sumBF;
    % Convolution of ER with PSFIUH  and MSSFIUH
    US=zeros(dimrain+sizeBF,1);
    US=conv(IPSF,PSFIUH);
    UBF=zeros(dimrain+sizeBF,1);
    UBF=conv(IMSSF,MSSFIUH);
    szUS= length(US);
    szUBF=length(UBF);
    US(szUS+1:szUBF)=0;
    %Total flow hydrograph
    U=US+UBF; %(1:size(US));
    total_length=length(Q_d) ;
    U_test=U(6941:total_length);% Modelled discharge for testing period 6699=1/10/1994
    Q_test=Q_d(6941:total_length); % Observed discharge for testing period 6699=1/10/1994
    Q_Kentuckymod{klr,1}=U_test;
    Q_Kentuckyobs{klr,1}=Q_test;
end
%%
for bn=1:nos_basins
    qallmod(:,bn)=Q_Kentuckymod{bn,1}*basin_area(bn,1)*10^3;
    qallobs(:,bn)=Q_Kentuckyobs{bn,1}*basin_area(bn,1)*10^3;
end
[row col]=find(isnan(qallobs)); % finding NaN values from the observed data
qallobs(row,:)=[]; % omitting NaN from observed data 
qallmod(row,:)=[];% removing the corresponding data points from the modelled data
[sortedValueobs] = sort(qallobs,1,'descend');
[sortedValuemod] = sort(qallmod,1,'descend');
%% Scaling exponent for ranked observed Q values
for ev=1:length(sortedValueobs)
    rQobs(:,1)= sortedValueobs(ev,:);
    scaling_expobs(ev,:)= polyfit(log(basin_area),log(rQobs),1);
    rQobs=[];
end
%% Scaling exponent for ranked modelled Q values
for ev=1:length(sortedValuemod)
    rQmod(:,1)= sortedValuemod(ev,:);
    scaling_expmod(ev,:)= polyfit(log(basin_area),log(rQmod),1);
    rQmod=[];
end

semilogx(1:length(sortedValueobs), scaling_expobs(:,1),'b','LineWidth',1.5)
hold on
semilogx(1:length(sortedValuemod), scaling_expmod(:,1),'r','LineWidth',1.5)
xlabel('Rank');
ylabel('Theta');







