clear;clc;close all;
p=[4,5;
   5,5;
   6,8;
   7,7]; %四种作物
Height=18;Length=19;value=[4,5,7,7];
max_binarylength=5; %二进制最长长度

NP=100;%种群中个体数量
Iteration=1000;%总迭代次数
gen=1;%当前迭代次数
Pc=0.7;%交叉概率
Pm=0.1;%变异概率
Gap = 0.9;%每次丢弃父代百分之10，以便加入百分之10的子代

X=zeros(NP,4*max_binarylength);
for i = 1:NP  %初始化种群
    boxcode=encode(p,Height,Length);%搭建可以实现的整箱，并返回每种箱子的二进制码
    X(i,:)=boxcode;
end

while (gen<=Iteration) %迭代过程
    Y = decode(X,NP) %把二进制码转成十进制
    Z = Fitness(Y,NP,value,p) %计算适应度
    parentselector = Select(Z,NP,X,Gap)%轮盘赌算法 产生90个下一代
    kidselector = cross (parentselector,NP,Pc,Gap) %通过父辈的交叉随机算出10个子辈
    X = [parentselector;kidselector];%子辈父辈重新整合在一起
    V = Variation(Pm,X);%下一代的随机变异
    X = V;
    gen=gen+1;
end

%%
function Y = decode(X,NP)   %将二进制转为十进制
Y=zeros(NP,4);sum=0;
for i = 1:NP
    for j = 1:20
        if (mod(j,5)==0)
            Y(i,fix(j/5))=sum;
            sum=0;
        end
        sum = X(i,j)+ sum*2;
    end
end

end

function Z = Fitness(Y,NP,value,p)  %计算利润 
Z=zeros(NP,1);
total_area = 342
for i=1:NP
    profit=0;
    area =0;
    for j=1:4
        %if判断不可能的情况，如果不可能，利润直接降为0 
        profit = profit + Y(i,j)*value(1,j);%作物总利润
        if Y(i,1)>16 %经过数学计算 第一种矩形最多只能放16个 因此大于16的情况全部归0
            profit=0;
            break;
        end
        if Y(i,2)>9
            profit=0;
            break;
        end
        if Y(i,3)>6
            profit=0;
            break;
        end
        if Y(i,4)>4
            profit=0;
            break;
        end
        if Y(i,1)>16
            profit=0;
            break;
        end
        if Y(i,1)*20+Y(i,2)*25+Y(i,3)*48+Y(i,4)*49>342
            profit=0;
            break;
        end
    end

    Z(i,1) = profit;
end
end

function parentselector = Select(Z,NP,X,Gap) %轮盘赌选择下一代
parentselector = zeros(NP*Gap,20);
sumfitness = sum (Z);
accP = cumsum(Z/sumfitness); %累积概率
for n = 1:NP*Gap
    matrix = find (accP>rand) %找到比随机数大的概率
    if isempty(matrix)
        continue;
    end
    temp = X(matrix(1),:);
    parentselector(n,:) = temp;
end
end

function  kidselector = cross (parentselector,NP,Pc,Gap) %交叉
kidselector = zeros(10,20);n=1;
father = zeros(1,20);
mother = zeros(1,20);
while n<=10
    father(1,:) = parentselector(ceil(rand*NP*Gap),:);  %随机产生父辈母辈  
    mother(1,:) = parentselector(ceil(rand*NP*Gap),:);
    crossLocation = ceil(20*rand);
    if rand < Pc %根据概率进行交叉
        father(1,crossLocation:20) = mother(1,crossLocation:20);
        kidselector(n,:) = father(1,:);
        n=n+1;
    end
end
end

function V = Variation(Pm,X)  %变异
V=zeros(100,20);
for n=1:100
    if (rand < Pm)
        location = ceil(20*rand);
        X(n,location) = 1-X(n,location); %变异采用一位的0 1变换
        
    end
end
V = X;
end

function boxcode=encode(p,Height,Length) %重点代码 初始化种群 搭矩形部分
    boxcode = [];
    count =zeros(1,4);
    x_max=Height;y_max=Length;
    x=1;y=1;
    A = ones(Height,Length);  % 1的地方是空 0的地方是有矩形
    kinds = randi(4);  %随机选择一个矩形类型
    m = round(rand)+1; %随机一种横竖方式
    count(1,kinds) = count(1,kinds)+1; %产生有效矩形后要加1 方面后期统计各个矩形的个数
    A (x:x+p(kinds,m)-1,y:y+p(kinds,3-m)-1) = 1-A (x:x+p(kinds,m)-1,y:y+p(kinds,3-m)-1); %把被矩形占据的地方实现0 1变换
    next_y= y+p(kinds,3-m);%下一个y的位置
    x = x + p(kinds,m);%首先遍历x轴
    while  true
        trylineflag = 0;
        for trylinetimes = 1:8        
            kinds = randi(4);
            m = round(rand)+1;
            if x+p(kinds,m)-1 <=x_max && y+p(kinds,3-m)-1<=y_max 
                position = find(A(x:x+p(kinds,m)-1,y+p(kinds,3-m)-1)==0); 
                if isempty(position) %预计空间没有矩形存在 即能放得下该矩形
                    A (x:x+p(kinds,m)-1,y:y+p(kinds,3-m)-1) = 1-A (x:x+p(kinds,m)-1,y:y+p(kinds,3-m)-1);
                    x = x + p(kinds,m) ;
                    trylineflag = 1;
                    count(1,kinds) = count(1,kinds)+1;
                    break;
                end 
            end
        end
        if (trylineflag==0) %当x轴遍历完无法再装矩形后 跳出第一次循环
            break;
       end
    end
   
    y = next_y;x=1;nextyflag=0; %第二次循环
    while  true
        trylineflag = 0;
        nextyflag=0;
        for trylinetimes = 1:16        
            kinds = randi(4);
            m = round(rand)+1;
            if x+p(kinds,m)-1 <=x_max && y+p(kinds,3-m)-1<=y_max 
                position = find(A(x:x+p(kinds,m)-1,y+p(kinds,3-m)-1)==0);
                if isempty(position)
                    A (x:x+p(kinds,m)-1,y:y+p(kinds,3-m)-1) = 1-A (x:x+p(kinds,m)-1,y:y+p(kinds,3-m)-1);
                    x = x + p(kinds,m) ;
                    trylineflag = 1;
                    count(1,kinds) = count(1,kinds)+1;
                    if (nextyflag==0)
                        nextyflag=1;
                        next_y= y+p(kinds,3-m);
                    end
                    break;
                end 
            end
        end
        if (trylineflag==0)
            break;
        end
    end
    
    y = next_y;x=1;nextyflag=0; %第三次循环
    while  true
        trylineflag = 0;
        nextyflag=0;
        for trylinetimes = 1:16        
            kinds = randi(4);
            m = round(rand)+1;
            if x+p(kinds,m)-1 <=x_max && y+p(kinds,3-m)-1<=y_max 
                position = find(A(x:x+p(kinds,m)-1,y+p(kinds,3-m)-1)==0);
                if isempty(position)
                    A (x:x+p(kinds,m)-1,y:y+p(kinds,3-m)-1) = 1-A (x:x+p(kinds,m)-1,y:y+p(kinds,3-m)-1);
                    x = x + p(kinds,m) ;
                    trylineflag = 1;
                    count(1,kinds) = count(1,kinds)+1;
                    if (nextyflag==0)
                        nextyflag=1;
                        next_y= y+p(kinds,3-m);
                    end
                    break;
                end 
            end
        end
        if (trylineflag==0)
            break;
        end
    end
    
    %三次循环后，样本显示基本完成箱子装搭,因此不再继续，现在开始编码
    for i=4:-1:1
        num=count(1,i);
        for j=1:5
            digit = mod(num,2);
            boxcode  = [boxcode digit];
            num = fix(num /2)
        end
    end
    
end