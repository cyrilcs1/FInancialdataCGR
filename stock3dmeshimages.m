nvrtx = 4;
filename='NQapril2020-may2020cgrnum.xlsx';
filename1='YMapril2020-may2020cgrnum.xlsx';
coarsedate = xlsread(filename);%reads data 1 corresponds to a 20 to t 7 to g and 3 to c
coarsedate1=xlsread(filename1);

binbound= [-0.01;0;0.01]; %specify binboundaries

%%%%give the bit
b1=4;b2=5;b3=6;b4=7;b5=8;b6=9;
bit1=2^b1;bit2=2^b2;
bit3=2^b3;bit4=2^b4;
bit5=2^b5;bit6=2^b6;

% Define the vertices
vrtx =    [0 0;bit1 0;0 bit1;bit1 bit1];
vrtxfive = [0 0;bit2 0;0 bit2;bit2 bit2];
vrtxsix= [0 0;bit3 0;0 bit3;bit3 bit3];
vrtxsev =    [0 0;bit4 0;0 bit4;bit4 bit4];
vrtxeit = [0 0;bit5 0;0 bit5;bit5 bit5];
vrtxnine= [0 0;bit6 0;0 bit6;bit6 bit6];
sz=size(coarsedate);
niter = sz(1,1);  % Specify the number of iterations
coarse=zeros(niter-1,1);
sz1=size(coarsedate1);
niter1 = sz1(1,1);  % Specify the number of iterations
coarse1=zeros(niter1-1,1);

percentchange=zeros(niter-1,1);

for i=2:niter
    percentchange(i-1)=((coarsedate(i)-coarsedate(i-1))/coarsedate(i))*100;
end

percentchange1=zeros(niter1-1,1);

for i=2:niter1
    percentchange1(i-1)=((coarsedate1(i)-coarsedate1(i-1))/coarsedate1(i))*100;
end


count=0;
for i=1:niter-1
    if percentchange(i)~= 0
        count=count+1;
    end
end

percentchangenew=zeros(count,1);
i=1;j=1;
while i<niter
    if percentchange(i)~=0
        percentchangenew(j)=percentchange(i);
        i=i+1;
        j=j+1;
    else
        i=i+1;
    end
end


count1=0;
for i=1:niter1-1
    if percentchange1(i)~= 0
        count1=count1+1;
    end
end

percentchangenew1=zeros(count1,1);
i=1;j=1;
while i<niter1
    if percentchange1(i)~=0
        percentchangenew1(j)=percentchange1(i);
        i=i+1;
        j=j+1;
    else
        i=i+1;
    end
end

percentchange=percentchangenew;
percentchange1=percentchangenew1;

niter=count;
niter1=count1;

i=1;



while i<=niter
    if percentchange(i) <binbound(1)
        coarse(i)=1;  
       
        i=i+1;
       % i=i+1;
        
    
    elseif percentchange(i)<=binbound(2)
            coarse(i)=2;
            
           
           i=i+1;
                
            elseif percentchange(i)<=binbound(3)
                coarse(i)=3;
                
               
                i=i+1;
          else
                coarse(i)=4;
                i=i+1;
    end
            
end

i=1;


while i<=niter1
    if percentchange1(i) <binbound(1)
        coarse1(i)=1;  
       
        i=i+1;
       % i=i+1;
        
    
    elseif percentchange1(i)<=binbound(2)
            coarse1(i)=2;
            
           
           i=i+1;
                
            elseif percentchange1(i)<=binbound(3)
                coarse1(i)=3;
                
               
                i=i+1;
          else
                coarse1(i)=4;
                i=i+1;
    end
            
end

%n01-no4 specify the number of CATand G respectivley
no1=0;no2=0;no3=0;no4=0;
for i=1:niter
    if coarse(i)==1
        no1=no1+1;
    elseif coarse(i)==2
        no2=no2+1;
    elseif coarse(i)==3
        no3=no3+1;
    else
        no4=no4+1;
    end
end
percentno1=no1*100/niter;
percentno2=no2*100/niter;
percentno3=no3*100/niter;
percentno4=no4*100/niter;

% four bit cgr

point = zeros(niter,2);
point(1,1)=bit1/2;point(1,2)=bit1/2;
for i = 2:niter                                              % Generate the points
    vIdx = coarse(i-1,1);
    point(i,:) = vrtx(vIdx,:) - (vrtx(vIdx,:) - point(i-1,:))/2;
end

frbit= zeros(bit1,bit1);
for i=1:niter
    for k=1:bit1
        if (point(i,1)<=k && point(i,1)>(k-1))
            for j=1:bit1
                if( point(i,2)<=j && point(i,2)>(j-1))
                    frbit(k,j)=frbit(k,j)+1;
                end
            end
            
            
        end
    end
end
frbitp=frbit;%frbitp will be filled or notfilled ie 0 or 1.
frbitp(logical(frbitp))=1;

%five bit

pointfive = zeros(niter,2) ;
pointfive(1,1)=bit2/2;pointfive(1,2)=bit2/2;
for i = 2:niter                                              % Generate the points
    vIdx = coarse(i-1,1);
    pointfive(i,:) = vrtxfive(vIdx,:) - (vrtxfive(vIdx,:) - pointfive(i-1,:))/2;
end

fivebit= zeros(bit2,bit2);
for i=1:niter
    for k=1:bit2
        if (pointfive(i,1)<=k && pointfive(i,1)>(k-1))
            for j=1:bit2
                if( pointfive(i,2)<=j && pointfive(i,2)>(j-1))
                    fivebit(k,j)=fivebit(k,j)+1;
                end
            end
            
            
        end
    end
end
fivebitp=fivebit;% will be filled or notfilled ie 0 or 1.
fivebitp(logical(fivebitp))=1;

%sixbit
pointsix = zeros(niter,2) ;
pointsix(1,1)=bit3/2;pointsix(1,2)=bit3/2;
for i = 2:niter                                              % Generate the points
    vIdx = coarse(i-1,1);
    pointsix(i,:) = vrtxsix(vIdx,:) - (vrtxsix(vIdx,:) - pointsix(i-1,:))/2;
end

sixbit= zeros(bit3,bit3);
for i=1:niter
    for k=1:bit3
        if (pointsix(i,1)<=k && pointsix(i,1)>(k-1))
            for j=1:bit3
                if( pointsix(i,2)<=j && pointsix(i,2)>(j-1))
                    sixbit(k,j)=sixbit(k,j)+1;
                end
            end
            
            
        end
    end
end
sixbitp=sixbit;% will be filled or notfilled ie 0 or 1.
sixbitp(logical(sixbitp))=1;

%sevenbit
pointseven = zeros(niter,2) ;
pointseven(1,1)=bit4/2;pointseven(1,2)=bit4/2;
for i = 2:niter                                              % Generate the points
    vIdx = coarse(i-1,1);
    pointseven(i,:) = vrtxsev(vIdx,:) - (vrtxsev(vIdx,:) - pointseven(i-1,:))/2;
end

sevenbit= zeros(bit4,bit4);
for i=1:niter
    for k=1:bit4
        if (pointseven(i,1)<=k && pointseven(i,1)>(k-1))
            for j=1:bit4
                if( pointseven(i,2)<=j && pointseven(i,2)>(j-1))
                    sevenbit(k,j)=sevenbit(k,j)+1;
                end
            end
            
            
        end
    end
end
sevenbitp=sevenbit;% will be filled or notfilled ie 0 or 1.
sevenbitp(logical(sevenbitp))=1;

%eitbit
pointeit = zeros(niter,2) ;
pointeit(1,1)=bit5/2;pointeit(1,2)=bit5/2;
for i = 2:niter                                              % Generate the points
    vIdx = coarse(i-1,1);
    pointeit(i,:) = vrtxeit(vIdx,:) - (vrtxeit(vIdx,:) - pointeit(i-1,:))/2;
end

eitbit= zeros(bit5,bit5);
% k=1;
% j=1;

for i=1:niter
    for k=1:bit5
        if (pointeit(i,1)<=k && pointeit(i,1)>(k-1))
            for j=1:bit5
                if( pointeit(i,2)<=j && pointeit(i,2)>(j-1))
                    eitbit(k,j)=eitbit(k,j)+1;
                end
            end
            
            
        end
    end
end
eitbitp=eitbit; %same as frbitp
eitbitp(logical(eitbit))=1;

% nine bit
pointnine = zeros(niter,2) ;
pointnine(1,1)=bit6/2;pointnine(1,2)=bit6/2;
for i = 2:niter                                              % Generate the points
    vIdx = coarse(i-1,1);
    pointnine(i,:) = vrtxnine(vIdx,:) - (vrtxnine(vIdx,:) - pointnine(i-1,:))/2;
end

ninebit= zeros(bit6,bit6);
for i=1:niter
    for k=1:bit6
        if (pointnine(i,1)<=k && pointnine(i,1)>(k-1))
            for j=1:bit6
                if( pointnine(i,2)<=j && pointnine(i,2)>(j-1))
                    ninebit(k,j)=ninebit(k,j)+1;
                end
            end
            
            
        end
    end
end
ninebitp=ninebit;% will be filled or notfilled ie 0 or 1.
ninebitp(logical(ninebitp))=1;

%n01-no4 specify the number of CATand G respectivley
noo1=0;noo2=0;noo3=0;noo4=0;
for i=1:niter1
    if coarse1(i)==1
        noo1=noo1+1;
    elseif coarse1(i)==2
        noo2=noo2+1;
    elseif coarse1(i)==3
        noo3=noo3+1;
    else
        noo4=noo4+1;
    end
end

percent1noo1=noo1*100/niter1;
percent1noo2=noo2*100/niter1;
percent1noo3=noo3*100/niter1;
percent1noo4=noo4*100/niter1;



point1 = zeros(niter1,2);
point1(1,1)=bit1/2;point1(1,2)=bit1/2;
for i = 2:niter1                                              % Generate the points
    vIdx1 = coarse1(i-1,1);
    point1(i,:) = vrtx(vIdx1,:) - (vrtx(vIdx1,:) - point1(i-1,:))/2;
end

frbit1= zeros(bit1,bit1);
for i=1:niter1
    for k=1:bit1
        if (point1(i,1)<=k && point1(i,1)>(k-1))
            for j=1:bit1
                if( point1(i,2)<=j && point1(i,2)>(j-1))
                    frbit1(k,j)=frbit1(k,j)+1;
                end
            end
            
            
        end
    end
end
frbitp1=frbit1;%frbitp will be filled or notfilled ie 0 or 1.
frbitp1(logical(frbitp1))=1;

%five bit

pointfive1 = zeros(niter1,2) ;
pointfive1(1,1)=bit2/2;pointfive1(1,2)=bit2/2;
for i = 2:niter1                                              % Generate the points
    vIdx1 = coarse1(i-1,1);
    pointfive1(i,:) = vrtxfive(vIdx1,:) - (vrtxfive(vIdx1,:) - pointfive1(i-1,:))/2;
end

fivebit1= zeros(bit2,bit2);
for i=1:niter1
    for k=1:bit2
        if (pointfive1(i,1)<=k && pointfive1(i,1)>(k-1))
            for j=1:bit2
                if( pointfive1(i,2)<=j && pointfive1(i,2)>(j-1))
                    fivebit1(k,j)=fivebit1(k,j)+1;
                end
            end
            
            
        end
    end
end
fivebitp1=fivebit1;% will be filled or notfilled ie 0 or 1.
fivebitp1(logical(fivebitp1))=1;

%sixbit
pointsix1 = zeros(niter1,2) ;
pointsix1(1,1)=bit3/2;pointsix1(1,2)=bit3/2;
for i = 2:niter1                                              % Generate the points
    vIdx1 = coarse1(i-1,1);
    pointsix1(i,:) = vrtxsix(vIdx1,:) - (vrtxsix(vIdx1,:) - pointsix1(i-1,:))/2;
end

sixbit1= zeros(bit3,bit3);
for i=1:niter1
    for k=1:bit3
        if (pointsix1(i,1)<=k && pointsix1(i,1)>(k-1))
            for j=1:bit3
                if( pointsix1(i,2)<=j && pointsix1(i,2)>(j-1))
                    sixbit1(k,j)=sixbit1(k,j)+1;
                end
            end
            
            
        end
    end
end
sixbitp1=sixbit1;% will be filled or notfilled ie 0 or 1.
sixbitp1(logical(sixbitp1))=1;

%sevenbit
pointseven1 = zeros(niter1,2) ;
pointseven1(1,1)=bit4/2;pointseven1(1,2)=bit4/2;
for i = 2:niter1                                              % Generate the points
    vIdx1 = coarse1(i-1,1);
    pointseven1(i,:) = vrtxsev(vIdx1,:) - (vrtxsev(vIdx1,:) - pointseven1(i-1,:))/2;
end

sevenbit1= zeros(bit4,bit4);
for i=1:niter1
    for k=1:bit4
        if (pointseven1(i,1)<=k && pointseven1(i,1)>(k-1))
            for j=1:bit4
                if( pointseven1(i,2)<=j && pointseven1(i,2)>(j-1))
                    sevenbit1(k,j)=sevenbit1(k,j)+1;
                end
            end
            
            
        end
    end
end
sevenbitp1=sevenbit1;% will be filled or notfilled ie 0 or 1.
sevenbitp1(logical(sevenbitp1))=1;

%eitbit
pointeit1 = zeros(niter1,2) ;
pointeit1(1,1)=bit5/2;pointeit1(1,2)=bit5/2;
for i = 2:niter1                                              % Generate the points
    vIdx1 = coarse1(i-1,1);
    pointeit1(i,:) = vrtxeit(vIdx1,:) - (vrtxeit(vIdx1,:) - pointeit1(i-1,:))/2;
end

eitbit1= zeros(bit5,bit5);
% k=1;
% j=1;

for i=1:niter1
    for k=1:bit5
        if (pointeit1(i,1)<=k && pointeit1(i,1)>(k-1))
            for j=1:bit5
                if( pointeit1(i,2)<=j && pointeit1(i,2)>(j-1))
                    eitbit1(k,j)=eitbit1(k,j)+1;
                end
            end
            
            
        end
    end
end
eitbitp1=eitbit1; %same as frbitp
eitbitp1(logical(eitbit1))=1;

% nine bit
pointnine1 = zeros(niter1,2) ;
pointnine1(1,1)=bit6/2;pointnine1(1,2)=bit6/2;
for i = 2:niter1                                              % Generate the points
    vIdx1 = coarse1(i-1,1);
    pointnine1(i,:) = vrtxnine(vIdx1,:) - (vrtxnine(vIdx1,:) - pointnine1(i-1,:))/2;
end

ninebit1= zeros(bit6,bit6);
for i=1:niter1
    for k=1:bit6
        if (pointnine1(i,1)<=k && pointnine1(i,1)>(k-1))
            for j=1:bit6
                if( pointnine1(i,2)<=j && pointnine1(i,2)>(j-1))
                    ninebit1(k,j)=ninebit1(k,j)+1;
                end
            end
            
            
        end
    end
end
ninebitp1=ninebit1;% will be filled or notfilled ie 0 or 1.
ninebitp1(logical(ninebitp1))=1;

frbitpercent=frbit*100/niter; frbitpercent1=frbit1*100/niter1;
fivebitpercent=fivebit*100/niter; fivebitpercent1=fivebit1*100/niter1;
sixbitpercent=sixbit*100/niter; sixbitpercent1=sixbit1*100/niter1;
sevenbitpercent=sevenbit*100/niter; sevenbitpercent1=sevenbit1*100/niter1;
eitbitpercent=eitbit*100/niter; eitbitpercent1=eitbit1*100/niter1;
ninebitpercent=ninebit*100/niter; ninebitpercent1=ninebit1*100/niter1;

sfrbit=frbitpercent-frbitpercent1;
sfivebit=fivebitpercent-fivebitpercent1;
ssixbit=sixbitpercent-sixbitpercent1;
ssevenbit=sevenbitpercent-sevenbitpercent1;
seitbit=eitbitpercent-eitbitpercent1;
sninebit=ninebitpercent-ninebitpercent1;

frbitpp=frbitpercent;
frbitpp(logical(frbitpp))=1;
fivebitpp=fivebitpercent;
fivebitpp(logical(fivebitpp))=1;
sixbitpp=sixbitpercent;
sixbitpp(logical(sixbitpp))=1;
sevenbitpp=sevenbitpercent;
sevenbitpp(logical(sevenbitpp))=1;
eitbitpp=eitbitpercent;
eitbitpp(logical(eitbitpp))=1;
ninebitpp=ninebitpercent;
ninebitpp(logical(ninebitpp))=1;

frbitpp1=frbitpercent1;
frbitpp1(logical(frbitpp1))=1;
fivebitpp1=fivebitpercent1;
fivebitpp1(logical(fivebitpp1))=1;
sixbitpp1=sixbitpercent1;
sixbitpp1(logical(sixbitpp1))=1;
sevenbitpp1=sevenbitpercent1;
sevenbitpp1(logical(sevenbitpp1))=1;
eitbitpp1=eitbitpercent1;
eitbitpp1(logical(eitbitpp1))=1;
ninebitpp1=ninebitpercent1;
ninebitpp1(logical(ninebitpp1))=1;


sfrbitpp=frbitpp-frbitpp1;
sfivebitpp=fivebitpp-fivebitpp1;
ssixbitpp=sixbitpp-sixbitpp1;
ssevenbitpp=sevenbitpp-sevenbitpp1;
seitbitpp=eitbitpp-eitbitpp1;
sninebitpp=ninebitpp-ninebitpp1;


figure,
surfl(1:bit1,1:bit1,(frbitpercent));shading interp;colormap(pink);

figure,
surfl(1:bit1,1:bit1,(frbitpercent1));shading interp;colormap(pink);


