clear all
clc
tic
load Data_kensynt4ev1.mat
[basin_area name]=xlsread('BasinAreaAll.xlsx');
rs=4; % Resolution in hours
klr=1
Q=cell(23,1);
for klr=1:23
    klr
    basin_id=Data_kensyn(klr,:);
    P=basin_id{1,1};
    ER=P;
    VLU=basin_id{1,5};
    VLD=basin_id{1,6};
    VFD=basin_id{1,7};
    dimluld=size(VLU,1);
     vs=3;
     vd=0.01;
     alpha= 0.3;
    %%
    timeI=zeros(1,1);
    timeF=zeros(1,1);
    dimrain=size(ER,1);
    for k=1:dimluld
        timeI(k)=VLD(k)/vs/rs;
        timeF(k)=  (VLU(k)+VLD(k))/vs/rs+ VLU(k)/vd/rs; %dividing by 4 for 4 hourly time scale
    end
    timeImax=max(timeI);
    timeFmax=max(timeF);
    sizeBF=ceil(timeFmax)+1;
    sizeWF=ceil(timeImax)+1;
    CWF=zeros(sizeWF,1);
    PSFIUH=zeros(sizeBF,1);
    MSSFIUH= zeros(sizeBF,1);
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
    
    tm=0; %baseflow unit hydrograph computation
    for i=1:floor(sizeBF)
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
    %Partitioning the effective rainfall for fast and slow routing
    IS=alpha*ER;
    IBF=(1-alpha)*ER;
    % Convolution
    US=zeros(dimrain+sizeBF,1);
    US=conv(IS,PSFIUH);
    UBF=zeros(dimrain+sizeBF,1);
    UBF=conv(IBF,MSSFIUH);
    szUS= length(US);
    szUBF=length(UBF);
    %Total flow hydrograph
    U=US+UBF; 
    Q{klr,1}= U*basin_area(klr,1);
end
store_recession=cell(1,1);
for klr=1:23
idmax(klr,1)=find(Q{klr,1}==max(Q{klr,1}(:)),1);
store_recession{1,klr}= Q{klr,1}(idmax(klr):end); %Storing the discharge series from peak (maximum peak) onwards
end
Recconsi= store_recession;
log_area=log(basin_area)';
nos=size(Recconsi,1);
%%
for bn=1:size(Recconsi,2)
    szQ(bn,1)=size(Recconsi{1,bn},1);
end
maxsz=max(szQ);
minsz=min(szQ);
for bn=1:size(Recconsi,2)
    qu=Recconsi{1,bn}; 
    ct= zeros(maxsz-szQ(bn),1);
    qall(1:maxsz,bn)=[qu;ct];
end
[sortedValue] = sort(qall,1,'descend');
 para_day=zeros(1,2);
%% Scaling exponents for recession days for ranked Q values
for ev=1:length(sortedValue)
    rQ(1,:)= sortedValue(ev,:);
    para_day(ev,:)= polyfit(log_area,log(rQ),1);
    rQ=[];
end
semilogx(1:length( para_day), para_day(:,1),'m','LineWidth',1.5)
toc






