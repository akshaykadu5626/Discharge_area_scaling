clear all
load SampleLULD.mat %loading upstream and downstream length of the delineated channel pixels
VLU=LU;
VLD=LD;
dimluld=size(VLU,1);
Ts= 24; %time-step
timeI=zeros(dimluld,1);
timeF=zeros(dimluld,1);
 vs=2; % Surface flow velocity in Km/hr
 vd=0.005;% Desaturation rate Km/hr
alpha=0.3;% Splitting parameter
pixelsize=30; %Defining pixel size of DEM (30 m resolution)
for k=1:dimluld
    VLU(k)= VLU(k)*pixelsize/1000; %Converting lengths into km units
    VLD(k)= VLD(k)*pixelsize/1000;
    timeI(k)= VLD(k)/vs/Ts; %Time req. for channel pixel to start contributing PSF
    timeF(k)= (VLU(k)+VLD(k))/vs/Ts+ VLU(k)/vd/Ts; % Time till channel link will be active
end
timeImax=max(timeI);
sizeWF=floor(timeImax+2);
CWF=zeros(sizeWF,1);
WFIUH=zeros(sizeWF,1);
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

for i=2:sizeWF-1 %calculating  WFIUH from CWF
    WFIUH(i)= (CWF(i) - CWF(i+1))/CWF(1);
end
plot(WFIUH,'r','LineWidth',1.5)
timeFmax=max(timeF);
sizeBF=floor(timeFmax+2);
BFIUH= zeros(sizeBF,1);
tm=1; %baseflow unit hydrograph computation

for i=1:sizeBF
    for k=1:dimluld
        if (timeI(k)<tm && timeF(k)>tm)
            temp2= floor(VFD(k)/2);
            temp2=temp2*2;   %temp2 is an integer type variable.
            if(temp2~= VFD(k))
                BFIUH(i)= BFIUH(i)+1;    %If VFD[k] is even
            else
                BFIUH(i)= BFIUH(i)+sqrt(2);       %If VFD[k] is odd
            end
        end
    end
    tm=tm+1;
end
sumBF= sum(BFIUH);
BFIUH=BFIUH/sumBF;
plot(BFIUH,'b','LineWidth',1.5)
%%
rain=readmatrix('rain.txt');% Importing the effective rainfall
%%
dimrain=size(rain,1);
US=zeros(dimrain+sizeBF,1);
for i=1: dimrain+sizeWF %Convoluting fast flow UH with ER
    j=i;
    k=1;
    US(i,1)=0;
    while (j>i-dimrain)
        if(j>=1 && j<= sizeWF)
            US(i)= US(i)+rain(k)*WFIUH(j);
        end
   j=j-1;
   k=k+1;
    end
end

UBF=zeros(dimrain+sizeBF,1);
for i=1: dimrain+sizeBF %Convoluting slow flow UH with ER
    j=i;
    k=1;
    UBF(i,1)=0;
    while (j>i-dimrain)
        if(j>=1 && j<= sizeBF)
            UBF(i)= UBF(i)+rain(k)*BFIUH(j);
        end
   j=j-1;
   k=k+1;
    end
end

U= alpha*US+ (1-alpha)*UBF;

plot(US,'r','LineWidth',1.5)
hold on
plot(UBF,'b','LineWidth',1.5)
hold on
plot(U,'k','LineWidth',1.5)
