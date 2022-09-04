clear;clc;close all;
%����Ŀǰ��1000�˲������A����
m1=0.42875;m2=0.15125;m3=0.57000;m4=0.51375;
%��λ ÿ����ÿƽ���ײ��������ǧ��  kg / m^2 * month  �����˽������������������  
%����5322��Ӧ����m1������6029��Ӧ����m2 ��  ����2233��Ӧ����m3 �� ����2281��Ӧ����m4

outputcrop = [m1     1.6*m1   0;
              m2     1.05*m2  0;
              m3     0.8*m3   0.6*m3;
              m4     m4       0];%��ֲͬ�ﲻͬ�ַ���Ӧ�Ĳ���Ч�� 

%1�ǵ����� 2��֬�� 3��̼ˮ������ 4��ˮ
nutritionlist=[0.1110,  0.0200,  0.7100, 0.1202;
               0.3430,  0.1750,  0.2670,  0.1200;
               0.0934,  0.0136,  0.6045, 0.1196;
               0.0955,  0.0400,  0.7097, 0.1135]; % ��λkg/1kg
nutritioncommand = [36.5, 29.2, 182.5,   0]; %һ����һ�������Ӫ�� ��λkg
%�ٶ�����һ����ÿ����Ҫ0.1kg�ĵ����� 0.08kg��֬�� 0.4kg��̼ˮ������ ˮ�Ļ����Ժ�����ˮ �˴�����Ҫ�ر�������л�ȡ

nutritionweight = [4000.0 ,9000.0, 4000.0 ,0]; %������,�ҽ�Ȩ�ػ�Ϊ����   ��λkcal/kg

totalarea = 100; %�����ƽ����

%7�ִ����ַ�(6������5532��6������2281Ч���г����ظ� ��˲����ǣ�
%S1��Ӧ0.8*m3*6+m4*6   S2��Ӧ1.05*m2*8+0.6*m3*4  S3��Ӧ1.6*m1*12  S4��Ӧ 12*m2
%S5��Ӧ12*m3  S6��Ӧ12*m4  S7��Ӧ1.6*m1*6+0.8*m3*3 

max_binarylength=7; %�����������
NP=100;%��Ⱥ�и�������
Iteration=1000;%�ܵ�������
gen=1;%��ǰ��������
Pc=0.7;%�������
Pm=0.1;%�������
Gap = 0.9;%ÿ�ζ��������ٷ�֮10���Ա����ٷ�֮10���Ӵ�

X=zeros(NP,7*max_binarylength);
for i = 1:NP  %��ʼ����Ⱥ
    matchcode = encode(max_binarylength,totalarea); %����ͬ������ַ����һ����ֲ���
    X(i,:) = matchcode;
end

while (gen<=Iteration) %��������
    Y = decode(X,NP,max_binarylength) %�Ѷ�������ת��ʮ����
    Z = Fitness(Y,NP,nutritionlist,nutritionweight,nutritioncommand,outputcrop) %������Ӧ��
    parentselector = Select(Z,NP,X,Gap,max_binarylength)%���̶��㷨 ����90����һ��
    kidselector = cross (parentselector,NP,Pc,Gap,max_binarylength) %ͨ�������Ľ���������10���ӱ�
    X = [parentselector;kidselector];%�ӱ���������������һ��
    V = Variation(Pm,X,max_binarylength,NP);%��һ�����������
    X = V;
    gen=gen+1;
end







%%
function Y = decode(X,NP,max_binarylength)   %��������תΪʮ����
Y=zeros(NP,7);sum=0;
for i = 1:NP
    for j = 1:7*max_binarylength
        if (mod(j,7)==0)
            Y(i,fix(j/7))=sum;
            sum=0;
        end
        sum = X(i,j)+ sum*2;
    end
end

end

function Z = Fitness(Y,NP,nutritionlist,nutritionweight,nutritioncommand,outputcrop)  
Z=zeros(NP,1);

for i=1:NP
    realnutrition= zeros(1,4); %Ŀǰ�ַ�ʵ�ʻ�õ�Ӫ��
    
    chase = 0;%Ŀ������������������ǰ���� ʵ�ֻ�ȡ���������
    %Ӫ�������
    realnutrition(i,1) =  Y(i,1)*outputcrop(3,2)*6*nutritionlist(3,1) + Y(i,1)*outputcrop(4,2)*6*nutritionlist(4,1) ...
                       +  Y(i,2)*outputcrop(2,2)*8*nutritionlist(2,1) + Y(i,2)*outputcrop(3,3)*4*nutritionlist(3,1) ...
                       +  Y(i,3)*outputcrop(1,2)*12*nutritionlist(1,1)  ...
                       +  Y(i,4)*outputcrop(2,1)*12*nutritionlist(2,1) ...
                       +  Y(i,5)*outputcrop(3,1)*12*nutritionlist(3,1) ...
                       +  Y(i,6)*outputcrop(4,1)*12*nutritionlist(4,1) ...
                       +  Y(i,7)*outputcrop(1,2)*6*nutritionlist(1,1) + Y(i,7)*outputcrop(3,2)*6*nutritionlist(3,1);
    
     realnutrition(i,2) =  Y(i,1)*outputcrop(3,2)*6*nutritionlist(3,2) + Y(i,1)*outputcrop(4,2)*6*nutritionlist(4,2) ...
                       +  Y(i,2)*outputcrop(2,2)*8*nutritionlist(2,2) + Y(i,2)*outputcrop(3,3)*4*nutritionlist(3,2) ...
                       +  Y(i,3)*outputcrop(1,2)*12*nutritionlist(1,2)  ...
                       +  Y(i,4)*outputcrop(2,2)*12*nutritionlist(2,2) ...
                       +  Y(i,5)*outputcrop(3,1)*12*nutritionlist(3,2) ...
                       +  Y(i,6)*outputcrop(4,1)*12*nutritionlist(4,2) ...
                       +  Y(i,7)*outputcrop(1,2)*6*nutritionlist(1,2) + Y(i,7)*outputcrop(3,2)*6*nutritionlist(3,2);                            
    
      realnutrition(i,3) =  Y(i,1)*outputcrop(3,2)*6*nutritionlist(3,3) + Y(i,1)*outputcrop(4,2)*6*nutritionlist(4,3) ...
                       +  Y(i,2)*outputcrop(2,2)*8*nutritionlist(2,3) + Y(i,2)*outputcrop(3,3)*4*nutritionlist(3,3) ...
                       +  Y(i,3)*outputcrop(1,2)*12*nutritionlist(1,3)  ...
                       +  Y(i,4)*outputcrop(2,1)*12*nutritionlist(2,3) ...
                       +  Y(i,5)*outputcrop(3,1)*12*nutritionlist(3,3) ...
                       +  Y(i,6)*outputcrop(4,1)*12*nutritionlist(4,3) ...
                       +  Y(i,7)*outputcrop(1,2)*6*nutritionlist(1,3) + Y(i,7)*outputcrop(3,2)*6*nutritionlist(3,3);
    
      realnutrition(i,4) =  Y(i,1)*outputcrop(3,2)*6*nutritionlist(3,4) + Y(i,1)*outputcrop(4,2)*6*nutritionlist(4,4) ...
                       +  Y(i,2)*outputcrop(2,2)*8*nutritionlist(2,4) + Y(i,2)*outputcrop(3,3)*4*nutritionlist(3,4) ...
                       +  Y(i,3)*outputcrop(1,2)*12*nutritionlist(1,4)  ...
                       +  Y(i,4)*outputcrop(2,1)*12*nutritionlist(2,4) ...
                       +  Y(i,5)*outputcrop(3,1)*12*nutritionlist(3,4) ...
                       +  Y(i,6)*outputcrop(4,1)*12*nutritionlist(4,4) ...
                       +  Y(i,7)*outputcrop(1,2)*6*nutritionlist(1,4) + Y(i,7)*outputcrop(3,2)*6*nutritionlist(3,4); 
      
       chase = realnutrition(i,1)*nutritionweight(1,1)+realnutrition(i,2)*nutritionweight(1,2)+realnutrition(i,3)*nutritionweight(1,3)+realnutrition(i,4)*nutritionweight(1,4);
       %���޷��������Ӫ������ʱ ���ʸ�ֵΪ0
       if  nutritioncommand(1,1) >realnutrition(i,1)
            chase = 0; 
       elseif nutritioncommand(1,2) >realnutrition(i,2)
            chase = 0; 
       elseif nutritioncommand(1,3) >realnutrition(i,3)
            chase = 0; 
       elseif nutritioncommand(1,4) >realnutrition(i,4)
            chase = 0; 
       end    
       
       Z(i,1) = chase;
end
end

function parentselector = Select(Z,NP,X,Gap,max_binarylength) %���̶�ѡ����һ��
parentselector = zeros(NP*Gap,7*max_binarylength);
sumfitness = sum (Z);
accP = cumsum(Z/sumfitness); %�ۻ�����
for n = 1:NP*Gap
    matrix = find (accP>rand) %�ҵ����������ĸ���
    if isempty(matrix)
        continue;
    end
    temp = X(matrix(1),:); %��һ������������ʴ�Ļ�����
    parentselector(n,:) = temp;
end
end

function  kidselector = cross (parentselector,NP,Pc,Gap,max_binarylength) %����
kidselector = zeros(10,7*max_binarylength);n=1;
father = zeros(1,7*max_binarylength);
mother = zeros(1,7*max_binarylength);
while n<=10
    father(1,:) = parentselector(ceil(rand*NP*Gap),:);  %�����������ĸ��  
    mother(1,:) = parentselector(ceil(rand*NP*Gap),:);
    crossLocation = ceil(7*max_binarylength*rand);
    if rand < Pc %���ݸ��ʽ��н���
        father(1,crossLocation:7*max_binarylength) = mother(1,crossLocation:7*max_binarylength);
        kidselector(n,:) = father(1,:);
        n=n+1;
    end
end
end

function V = Variation(Pm,X,max_binarylength,NP)  %����
V=zeros(NP,7*max_binarylength);
for n=1:NP
    if (rand < Pm)
        location = ceil(7*max_binarylength*rand);
        X(n,location) = 1-X(n,location); %�������һλ��0 1�任
    end
end
V = X;
end

function  matchcode = encode(max_binarylength,totalarea)     
    matchcode = []
    
    digit = randi([0,100]); %�������7����Ϊ100������ �ֱ�Ϊ��Ӧ��7����ֲ���䷨�����
    count=[];
    count = [count digit];
    i = 1;
    maximum = 100 - digit
    while (i<=5)
        digit = randi([0,maximum])
        i = i+1;
        count = [count digit];
        maximum = maximum - digit;
    end
    count = [count maximum];
    
    for n=7:-1:1 %���ж����Ƶı���
        num=count(1,n);
        for j=1:max_binarylength
            digit = mod(num,2);
            matchcode  = [matchcode digit];
            num = fix(num /2)
        end
    end
    
end