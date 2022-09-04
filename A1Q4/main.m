clear;clc;close all;
p=[4,5;
   5,5;
   6,8;
   7,7]; %��������
Height=18;Length=19;value=[4,5,7,7];
max_binarylength=5; %�����������

NP=100;%��Ⱥ�и�������
Iteration=1000;%�ܵ�������
gen=1;%��ǰ��������
Pc=0.7;%�������
Pm=0.1;%�������
Gap = 0.9;%ÿ�ζ��������ٷ�֮10���Ա����ٷ�֮10���Ӵ�

X=zeros(NP,4*max_binarylength);
for i = 1:NP  %��ʼ����Ⱥ
    boxcode=encode(p,Height,Length);%�����ʵ�ֵ����䣬������ÿ�����ӵĶ�������
    X(i,:)=boxcode;
end

while (gen<=Iteration) %��������
    Y = decode(X,NP) %�Ѷ�������ת��ʮ����
    Z = Fitness(Y,NP,value,p) %������Ӧ��
    parentselector = Select(Z,NP,X,Gap)%���̶��㷨 ����90����һ��
    kidselector = cross (parentselector,NP,Pc,Gap) %ͨ�������Ľ���������10���ӱ�
    X = [parentselector;kidselector];%�ӱ���������������һ��
    V = Variation(Pm,X);%��һ�����������
    X = V;
    gen=gen+1;
end

%%
function Y = decode(X,NP)   %��������תΪʮ����
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

function Z = Fitness(Y,NP,value,p)  %�������� 
Z=zeros(NP,1);
total_area = 342
for i=1:NP
    profit=0;
    area =0;
    for j=1:4
        %if�жϲ����ܵ��������������ܣ�����ֱ�ӽ�Ϊ0 
        profit = profit + Y(i,j)*value(1,j);%����������
        if Y(i,1)>16 %������ѧ���� ��һ�־������ֻ�ܷ�16�� ��˴���16�����ȫ����0
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

function parentselector = Select(Z,NP,X,Gap) %���̶�ѡ����һ��
parentselector = zeros(NP*Gap,20);
sumfitness = sum (Z);
accP = cumsum(Z/sumfitness); %�ۻ�����
for n = 1:NP*Gap
    matrix = find (accP>rand) %�ҵ����������ĸ���
    if isempty(matrix)
        continue;
    end
    temp = X(matrix(1),:);
    parentselector(n,:) = temp;
end
end

function  kidselector = cross (parentselector,NP,Pc,Gap) %����
kidselector = zeros(10,20);n=1;
father = zeros(1,20);
mother = zeros(1,20);
while n<=10
    father(1,:) = parentselector(ceil(rand*NP*Gap),:);  %�����������ĸ��  
    mother(1,:) = parentselector(ceil(rand*NP*Gap),:);
    crossLocation = ceil(20*rand);
    if rand < Pc %���ݸ��ʽ��н���
        father(1,crossLocation:20) = mother(1,crossLocation:20);
        kidselector(n,:) = father(1,:);
        n=n+1;
    end
end
end

function V = Variation(Pm,X)  %����
V=zeros(100,20);
for n=1:100
    if (rand < Pm)
        location = ceil(20*rand);
        X(n,location) = 1-X(n,location); %�������һλ��0 1�任
        
    end
end
V = X;
end

function boxcode=encode(p,Height,Length) %�ص���� ��ʼ����Ⱥ ����β���
    boxcode = [];
    count =zeros(1,4);
    x_max=Height;y_max=Length;
    x=1;y=1;
    A = ones(Height,Length);  % 1�ĵط��ǿ� 0�ĵط����о���
    kinds = randi(4);  %���ѡ��һ����������
    m = round(rand)+1; %���һ�ֺ�����ʽ
    count(1,kinds) = count(1,kinds)+1; %������Ч���κ�Ҫ��1 �������ͳ�Ƹ������εĸ���
    A (x:x+p(kinds,m)-1,y:y+p(kinds,3-m)-1) = 1-A (x:x+p(kinds,m)-1,y:y+p(kinds,3-m)-1); %�ѱ�����ռ�ݵĵط�ʵ��0 1�任
    next_y= y+p(kinds,3-m);%��һ��y��λ��
    x = x + p(kinds,m);%���ȱ���x��
    while  true
        trylineflag = 0;
        for trylinetimes = 1:8        
            kinds = randi(4);
            m = round(rand)+1;
            if x+p(kinds,m)-1 <=x_max && y+p(kinds,3-m)-1<=y_max 
                position = find(A(x:x+p(kinds,m)-1,y+p(kinds,3-m)-1)==0); 
                if isempty(position) %Ԥ�ƿռ�û�о��δ��� ���ܷŵ��¸þ���
                    A (x:x+p(kinds,m)-1,y:y+p(kinds,3-m)-1) = 1-A (x:x+p(kinds,m)-1,y:y+p(kinds,3-m)-1);
                    x = x + p(kinds,m) ;
                    trylineflag = 1;
                    count(1,kinds) = count(1,kinds)+1;
                    break;
                end 
            end
        end
        if (trylineflag==0) %��x��������޷���װ���κ� ������һ��ѭ��
            break;
       end
    end
   
    y = next_y;x=1;nextyflag=0; %�ڶ���ѭ��
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
    
    y = next_y;x=1;nextyflag=0; %������ѭ��
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
    
    %����ѭ����������ʾ�����������װ��,��˲��ټ��������ڿ�ʼ����
    for i=4:-1:1
        num=count(1,i);
        for j=1:5
            digit = mod(num,2);
            boxcode  = [boxcode digit];
            num = fix(num /2)
        end
    end
    
end