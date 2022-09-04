clear;clc;close all;
%假设目前有1000人参与调研A星球
m1=0.42875;m2=0.15125;m3=0.57000;m4=0.51375;
%单位 每个月每平方米产出作物的千克  kg / m^2 * month  （用了金克拉，产量翻翻翻！  
%作物5322对应的是m1；作物6029对应的是m2 ；  作物2233对应的是m3 ； 作物2281对应的是m4

outputcrop = [m1     1.6*m1   0;
              m2     1.05*m2  0;
              m3     0.8*m3   0.6*m3;
              m4     m4       0];%不同植物不同种法对应的产量效益 

%1是蛋白质 2是脂肪 3是碳水化合物 4是水
nutritionlist=[0.1110,  0.0200,  0.7100, 0.1202;
               0.3430,  0.1750,  0.2670,  0.1200;
               0.0934,  0.0136,  0.6045, 0.1196;
               0.0955,  0.0400,  0.7097, 0.1135]; % 单位kg/1kg
nutritioncommand = [36.5, 29.2, 182.5,   0]; %一个人一年所需的营养 单位kg
%假定的是一个人每天需要0.1kg的蛋白质 0.08kg的脂肪 0.4kg的碳水化合物 水的话可以喝饮用水 此处不需要特别从作物中获取

nutritionweight = [4000.0 ,9000.0, 4000.0 ,0]; %此题中,我将权重化为热量   单位kcal/kg

totalarea = 100; %总体的平方数

%7种搭配种法(6个月种5532，6个月种2281效益有出现重复 因此不考虑）
%S1对应0.8*m3*6+m4*6   S2对应1.05*m2*8+0.6*m3*4  S3对应1.6*m1*12  S4对应 12*m2
%S5对应12*m3  S6对应12*m4  S7对应1.6*m1*6+0.8*m3*3 

max_binarylength=7; %二进制最长长度
NP=100;%种群中个体数量
Iteration=1000;%总迭代次数
gen=1;%当前迭代次数
Pc=0.7;%交叉概率
Pm=0.1;%变异概率
Gap = 0.9;%每次丢弃父代百分之10，以便加入百分之10的子代

X=zeros(NP,7*max_binarylength);
for i = 1:NP  %初始化种群
    matchcode = encode(max_binarylength,totalarea); %给不同搭配的种法随机一个种植面积
    X(i,:) = matchcode;
end

while (gen<=Iteration) %迭代过程
    Y = decode(X,NP,max_binarylength) %把二进制码转成十进制
    Z = Fitness(Y,NP,nutritionlist,nutritionweight,nutritioncommand,outputcrop) %计算适应度
    parentselector = Select(Z,NP,X,Gap,max_binarylength)%轮盘赌算法 产生90个下一代
    kidselector = cross (parentselector,NP,Pc,Gap,max_binarylength) %通过父辈的交叉随机算出10个子辈
    X = [parentselector;kidselector];%子辈父辈重新整合在一起
    V = Variation(Pm,X,max_binarylength,NP);%下一代的随机变异
    X = V;
    gen=gen+1;
end







%%
function Y = decode(X,NP,max_binarylength)   %将二进制转为十进制
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
    realnutrition= zeros(1,4); %目前种法实际获得的营养
    
    chase = 0;%目标是在满足基本需求的前提下 实现获取热量的最大化
    %营养计算表
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
       %当无法满足基本营养需求时 概率赋值为0
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

function parentselector = Select(Z,NP,X,Gap,max_binarylength) %轮盘赌选择下一代
parentselector = zeros(NP*Gap,7*max_binarylength);
sumfitness = sum (Z);
accP = cumsum(Z/sumfitness); %累积概率
for n = 1:NP*Gap
    matrix = find (accP>rand) %找到比随机数大的概率
    if isempty(matrix)
        continue;
    end
    temp = X(matrix(1),:); %第一个比随机数概率大的基因型
    parentselector(n,:) = temp;
end
end

function  kidselector = cross (parentselector,NP,Pc,Gap,max_binarylength) %交叉
kidselector = zeros(10,7*max_binarylength);n=1;
father = zeros(1,7*max_binarylength);
mother = zeros(1,7*max_binarylength);
while n<=10
    father(1,:) = parentselector(ceil(rand*NP*Gap),:);  %随机产生父辈母辈  
    mother(1,:) = parentselector(ceil(rand*NP*Gap),:);
    crossLocation = ceil(7*max_binarylength*rand);
    if rand < Pc %根据概率进行交叉
        father(1,crossLocation:7*max_binarylength) = mother(1,crossLocation:7*max_binarylength);
        kidselector(n,:) = father(1,:);
        n=n+1;
    end
end
end

function V = Variation(Pm,X,max_binarylength,NP)  %变异
V=zeros(NP,7*max_binarylength);
for n=1:NP
    if (rand < Pm)
        location = ceil(7*max_binarylength*rand);
        X(n,location) = 1-X(n,location); %变异采用一位的0 1变换
    end
end
V = X;
end

function  matchcode = encode(max_binarylength,totalarea)     
    matchcode = []
    
    digit = randi([0,100]); %随机产生7个和为100的整数 分别为对应的7种种植搭配法的面积
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
    
    for n=7:-1:1 %进行二进制的编码
        num=count(1,n);
        for j=1:max_binarylength
            digit = mod(num,2);
            matchcode  = [matchcode digit];
            num = fix(num /2)
        end
    end
    
end