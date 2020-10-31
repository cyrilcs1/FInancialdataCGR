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
subplot('Position',[0.0241581259150808 0.0337941628264208 0.437774524158126 0.918586789554531]);
imagesc(1:bit1,1:bit1,rot90(frbitpercent));shading flat;colormap(flipud(gray(10)));colorbar('Position',[0.469375304863213 0.0276497695852535 0.0156173743289409 0.918586789554531]);
caxis([0 0.01]);

subplot('Position',[0.513038338484903 0.0337941628264209 0.440841602949951 0.918586789554531]);
imagesc(1:bit2,1:bit2,rot90(fivebitpercent));shading flat;colormap(flipud(gray(10)));colorbar('Position',[0.961688628599317 0.0337941628264209 0.0156173743289409 0.918586789554531]);
caxis ([0 0.001]);


figure,
subplot('Position',[0.0241581259150808 0.0337941628264208 0.437774524158126 0.918586789554531]);
imagesc(1:bit3,1:bit3,rot90(sixbitpercent));shading flat;colormap(flipud(gray(10)));colorbar('Position',[0.469375304863213 0.0276497695852535 0.0156173743289409 0.918586789554531]);
caxis([0 0.01]);

subplot('Position',[0.513038338484903 0.0337941628264209 0.440841602949951 0.918586789554531]);
imagesc(1:bit4,1:bit4,rot90(sevenbitpercent));shading flat;colormap(flipud(gray(10)));colorbar('Position',[0.961688628599317 0.0337941628264209 0.0156173743289409 0.918586789554531]);
caxis ([0 0.001]);

figure,
subplot('Position',[0.0241581259150808 0.0337941628264208 0.437774524158126 0.918586789554531]);
imagesc(1:bit5,1:bit5,rot90(eitbitpercent));shading flat;colormap(flipud(gray(10)));colorbar('Position',[0.469375304863213 0.0276497695852535 0.0156173743289409 0.918586789554531]);
caxis([0 0.01]);

subplot('Position',[0.513038338484903 0.0337941628264209 0.440841602949951 0.918586789554531]);
imagesc(1:bit6,1:bit6,rot90(ninebitpercent));shading flat;colormap(flipud(gray(10)));colorbar('Position',[0.961688628599317 0.0337941628264209 0.0156173743289409 0.918586789554531]);
caxis ([0 0.001]);


figure,
subplot('Position',[0.0241581259150808 0.0337941628264208 0.437774524158126 0.918586789554531]);
imagesc(1:bit1,1:bit1,rot90(frbitpercent1));shading flat;colormap(flipud(gray(10)));colorbar('Position',[0.469375304863213 0.0276497695852535 0.0156173743289409 0.918586789554531]);
caxis([0 0.01]);

subplot('Position',[0.513038338484903 0.0337941628264209 0.440841602949951 0.918586789554531]);
imagesc(1:bit2,1:bit2,rot90(fivebitpercent1));shading flat;colormap(flipud(gray(10)));colorbar('Position',[0.961688628599317 0.0337941628264209 0.0156173743289409 0.918586789554531]);
caxis ([0 0.001]);


figure,
subplot('Position',[0.0241581259150808 0.0337941628264208 0.437774524158126 0.918586789554531]);
imagesc(1:bit3,1:bit3,rot90(sixbitpercent1));shading flat;colormap(flipud(gray(10)));colorbar('Position',[0.469375304863213 0.0276497695852535 0.0156173743289409 0.918586789554531]);
caxis([0 0.01]);

subplot('Position',[0.513038338484903 0.0337941628264209 0.440841602949951 0.918586789554531]);
imagesc(1:bit4,1:bit4,rot90(sevenbitpercent1));shading flat;colormap(flipud(gray(10)));colorbar('Position',[0.961688628599317 0.0337941628264209 0.0156173743289409 0.918586789554531]);
caxis ([0 0.001]);

figure,
subplot('Position',[0.0241581259150808 0.0337941628264208 0.437774524158126 0.918586789554531]);
imagesc(1:bit5,1:bit5,rot90(eitbitpercent1));shading flat;colormap(flipud(gray(10)));colorbar('Position',[0.469375304863213 0.0276497695852535 0.0156173743289409 0.918586789554531]);
caxis([0 0.01]);

subplot('Position',[0.513038338484903 0.0337941628264209 0.440841602949951 0.918586789554531]);
imagesc(1:bit6,1:bit6,rot90(ninebitpercent1));shading flat;colormap(flipud(gray(10)));colorbar('Position',[0.961688628599317 0.0337941628264209 0.0156173743289409 0.918586789554531]);
caxis ([0 0.001]);

c=jet(20);
d=jet(20);

figure,
subplot('Position',[0.0241581259150808 0.0337941628264208 0.437774524158126 0.918586789554531]);
imagesc(1:bit1,1:bit1,rot90(sfrbit));shading flat;colormap(c);colorbar('Position',[0.469375304863213 0.0276497695852535 0.0156173743289409 0.918586789554531]);
caxis([-0.05 0.05]);

subplot('Position',[0.513038338484903 0.0337941628264209 0.440841602949951 0.918586789554531]);
imagesc(1:bit2,1:bit2,rot90(sfivebit));shading flat;colormap(d);colorbar('Position',[0.961688628599317 0.0337941628264209 0.0156173743289409 0.918586789554531]);
caxis ([-0.05 0.05]);

figure,
subplot('Position',[0.0241581259150808 0.0337941628264208 0.437774524158126 0.918586789554531]);
imagesc(1:bit3,1:bit3,rot90(ssixbit));shading flat;colormap(c);colorbar('Position',[0.469375304863213 0.0276497695852535 0.0156173743289409 0.918586789554531]);
caxis([-0.05 0.05]);

subplot('Position',[0.513038338484903 0.0337941628264209 0.440841602949951 0.918586789554531]);
imagesc(1:bit4,1:bit4,rot90(ssevenbit));shading flat;colormap(d);colorbar('Position',[0.961688628599317 0.0337941628264209 0.0156173743289409 0.918586789554531]);
caxis ([-0.05 0.05]);

figure,
subplot('Position',[0.0241581259150808 0.0337941628264208 0.437774524158126 0.918586789554531]);
imagesc(1:bit5,1:bit5,rot90(seitbit));shading flat;colormap(c);colorbar('Position',[0.469375304863213 0.0276497695852535 0.0156173743289409 0.918586789554531]);
caxis([-0.05 0.05]);

subplot('Position',[0.513038338484903 0.0337941628264209 0.440841602949951 0.918586789554531]);
imagesc(1:bit6,1:bit6,rot90(sninebit));shading flat;colormap(d);colorbar('Position',[0.961688628599317 0.0337941628264209 0.0156173743289409 0.918586789554531]);
caxis ([-0.05 0.05]);


%calculating dissimilarity index

posindexfrbit=0;
negindexfrbit=0;


for i= 1:bit1
    for j= 1:bit1
 if (sfrbit(i,j)<=0)
     negindexfrbit=negindexfrbit+sfrbit(i,j);
 else
posindexfrbit=posindexfrbit+sfrbit(i,j);
 end
    end
end

posindexfivebit=0;
negindexfivebit=0;


for i= 1:bit2
    for j= 1:bit2
 if (sfivebit(i,j)<=0)
     negindexfivebit=negindexfivebit+sfivebit(i,j);
 else
posindexfivebit=posindexfivebit+sfivebit(i,j);
 end
    end
end

posindexsixbit=0;
negindexsixbit=0;


for i= 1:bit3
    for j= 1:bit3
 if (ssixbit(i,j)<=0)
     negindexsixbit=negindexsixbit+ssixbit(i,j);
 else
posindexsixbit=posindexsixbit+ssixbit(i,j);
 end
    end
end

posindexsevenbit=0;
negindexsevenbit=0;


for i= 1:bit4
    for j= 1:bit4
 if (ssevenbit(i,j)<=0)
     negindexsevenbit=negindexsevenbit+ssevenbit(i,j);
 else
posindexsevenbit=posindexsevenbit+ssevenbit(i,j);
 end
    end
end

posindexeitbit=0;
negindexeitbit=0;


for i= 1:bit5
    for j= 1:bit5
 if (seitbit(i,j)<=0)
     negindexeitbit=negindexeitbit+seitbit(i,j);
 else
posindexeitbit=posindexeitbit+seitbit(i,j);
 end
    end
end

posindexninebit=0;
negindexninebit=0;


for i= 1:bit6
    for j= 1:bit6
 if (sninebit(i,j)<=0)
     negindexninebit=negindexninebit+sninebit(i,j);
 else
posindexninebit=posindexninebit+sninebit(i,j);
 end
    end
end
%%%% matrx with similarity index %%%%%%%
similarityindexmatr(1,1)=posindexfrbit;similarityindexmatr(2,1)=negindexfrbit;
similarityindexmatr(1,2)=posindexfivebit;similarityindexmatr(2,2)=negindexfivebit;
similarityindexmatr(1,3)=posindexsixbit;similarityindexmatr(2,3)=negindexsixbit;
similarityindexmatr(1,4)=posindexsevenbit;similarityindexmatr(2,4)=negindexsevenbit;
similarityindexmatr(1,5)=posindexeitbit;similarityindexmatr(2,5)=negindexeitbit;
similarityindexmatr(1,6)=posindexninebit;similarityindexmatr(2,6)=negindexninebit;


%calculation of fractal dimension%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% amtrix with name fractal dime matrix is created 
% % row1 ->fd of filename virus
% % row2-> fd of filename1virus
% % row3-> fd of subtracton
% % row4-> fd of subtraction empty fill
frzero=0;fivezero=0;
sixzero=0;sevenzero=0;
eitzero=0;ninezero=0;
frnonz=0;fivenonz=0;
sixnonz=0;sevennonz=0;
eitnonz=0;ninenonz=0;

for i=1:bit1
    for j=1:bit1
        if frbit(i,j)==0
            frzero= frzero+1;

        else
            frnonz=frnonz+1;
        end
    end
end


for i=1:bit2
    for j=1:bit2
        if fivebit(i,j)==0
            fivezero= fivezero+1;

        else
            fivenonz=fivenonz+1;
        end
    end
end


for i=1:bit3
    for j=1:bit3
        if sixbit(i,j)==0
            sixzero= sixzero+1;

        else
            sixnonz=sixnonz+1;
        end
    end
end


for i=1:bit4
    for j=1:bit4
        if sevenbit(i,j)==0
            sevenzero= sevenzero+1;

        else
            sevennonz=sevennonz+1;
        end
    end
end


for i=1:bit5
    for j=1:bit5
        if eitbit(i,j)==0
            eitzero= eitzero+1;

        else
            eitnonz=eitnonz+1;
        end
    end
end


for i=1:bit6
    for j=1:bit6
        if ninebit(i,j)==0
            ninezero= ninezero+1;

        else
            ninenonz=ninenonz+1;
        end
    end
end



frzero1=0;fivezero1=0;
sixzero1=0;sevenzero1=0;
eitzero1=0;ninezero1=0;
frnonz1=0;fivenonz1=0;
sixnonz1=0;sevennonz1=0;
eitnonz1=0;ninenonz1=0;

for i=1:bit1
    for j=1:bit1
        if frbit1(i,j)==0
            frzero1= frzero1+1;

        else
            frnonz1=frnonz1+1;
        end
    end
end


for i=1:bit2
    for j=1:bit2
        if fivebit1(i,j)==0
            fivezero1= fivezero1+1;

        else
            fivenonz1=fivenonz1+1;
        end
    end
end


for i=1:bit3
    for j=1:bit3
        if sixbit1(i,j)==0
            sixzero1= sixzero1+1;

        else
            sixnonz1=sixnonz1+1;
        end
    end
end


for i=1:bit4
    for j=1:bit4
        if sevenbit1(i,j)==0
            sevenzero1= sevenzero1+1;

        else
            sevennonz1=sevennonz1+1;
        end
    end
end


for i=1:bit5
    for j=1:bit5
        if eitbit1(i,j)==0
            eitzero1= eitzero1+1;

        else
            eitnonz1=eitnonz1+1;
        end
    end
end


for i=1:bit6
    for j=1:bit6
        if ninebit1(i,j)==0
            ninezero1= ninezero1+1;

        else
            ninenonz1=ninenonz1+1;
        end
    end
end

frzerosub=0;fivezerosub=0;
sixzerosub=0;sevenzerosub=0;
eitzerosub=0;ninezerosub=0;
frnonzsub=0;fivenonzsub=0;
sixnonzsub=0;sevennonzsub=0;
eitnonzsub=0;ninenonzsub=0;

for i=1:bit1
    for j=1:bit1
        if sfrbit(i,j)==0
            frzerosub= frzerosub+1;

        else
            frnonzsub=frnonzsub+1;
        end
    end
end


for i=1:bit2
    for j=1:bit2
        if sfivebit(i,j)==0
            fivezerosub= fivezerosub+1;

        else
            fivenonzsub=fivenonzsub+1;
        end
    end
end


for i=1:bit3
    for j=1:bit3
        if ssixbit(i,j)==0
            sixzerosub= sixzerosub+1;

        else
            sixnonzsub=sixnonzsub+1;
        end
    end
end


for i=1:bit4
    for j=1:bit4
        if ssevenbit(i,j)==0
            sevenzerosub= sevenzerosub+1;

        else
            sevennonzsub=sevennonzsub+1;
        end
    end
end


for i=1:bit5
    for j=1:bit5
        if seitbit(i,j)==0
            eitzerosub= eitzerosub+1;

        else
            eitnonzsub=eitnonzsub+1;
        end
    end
end


for i=1:bit6
    for j=1:bit6
        if sninebit(i,j)==0
            ninezerosub= ninezerosub+1;

        else
            ninenonzsub=ninenonzsub+1;
        end
    end
end

fractaldimfrbit=log(frnonz)/log(bit1); %fractaldimension of virus given in filename bit1
fractaldimfivebit=log(fivenonz)/log(bit2);
fractaldimsixbit=log(sixnonz)/log(bit3);
fractaldimsevenbit=log(sevennonz)/log(bit4);
fractaldimeitbit=log(eitnonz)/log(bit5);%same bit2
fractaldimninebit=log(ninenonz)/log(bit6);


fractdimematrix(1,1)=fractaldimfrbit;fractdimematrix(1,2)=fractaldimfivebit;
fractdimematrix(1,3)=fractaldimsixbit;fractdimematrix(1,4)=fractaldimsevenbit;
fractdimematrix(1,5)=fractaldimeitbit;fractdimematrix(1,6)=fractaldimninebit;

fractaldimfrbit1=log(frnonz1)/log(bit1); %fractaldimension of virus given in filename1 bit1
fractaldimfivebit1=log(fivenonz1)/log(bit2);
fractaldimsixbit1=log(sixnonz1)/log(bit3);
fractaldimsevenbit1=log(sevennonz1)/log(bit4);
fractaldimeitbit1=log(eitnonz1)/log(bit5);
fractaldimninebit1=log(ninenonz1)/log(bit6);


fractdimematrix(2,1)=fractaldimfrbit1;fractdimematrix(2,2)=fractaldimfivebit1;
fractdimematrix(2,3)=fractaldimsixbit1;fractdimematrix(2,4)=fractaldimsevenbit1;
fractdimematrix(2,5)=fractaldimeitbit1;fractdimematrix(2,6)=fractaldimninebit1;

fractaldimfrbitsub=log(frnonzsub)/log(bit1); %fractaldimension of virus given in filename bit1
fractaldimfivebitsub=log(fivenonzsub)/log(bit2);
fractaldimsixbitsub=log(sixnonzsub)/log(bit3);
fractaldimsevenbitsub=log(sevennonzsub)/log(bit4);
fractaldimeitbitsub=log(eitnonzsub)/log(bit5);%same bit2
fractaldimninebitsub=log(ninenonzsub)/log(bit6);


fractdimematrix(3,1)=fractaldimfrbitsub;fractdimematrix(3,2)=fractaldimfivebitsub;
fractdimematrix(3,3)=fractaldimsixbitsub;fractdimematrix(3,4)=fractaldimsevenbitsub;
fractdimematrix(3,5)=fractaldimeitbitsub;fractdimematrix(3,6)=fractaldimninebitsub;


%fractal dimension considering matrices as empty fill

frzerosubpp=0;fivezerosubpp=0;
sixzerosubpp=0;sevenzerosubpp=0;
eitzerosubpp=0;ninezerosubpp=0;
frnonzsubpp=0;fivenonzsubpp=0;
sixnonzsubpp=0;sevennonzsubpp=0;
eitnonzsubpp=0;ninenonzsubpp=0;


for i=1:bit1
    for j=1:bit1
        if sfrbitpp(i,j)==0
            frzerosubpp= frzerosubpp+1;

        else
            frnonzsubpp=frnonzsubpp+1;
        end
    end
end


for i=1:bit2
    for j=1:bit2
        if sfivebitpp(i,j)==0
            fivezerosubpp= fivezerosubpp+1;

        else
            fivenonzsubpp=fivenonzsubpp+1;
        end
    end
end


for i=1:bit3
    for j=1:bit3
        if ssixbitpp(i,j)==0
            sixzerosubpp= sixzerosubpp+1;

        else
            sixnonzsubpp=sixnonzsubpp+1;
        end
    end
end


for i=1:bit4
    for j=1:bit4
        if ssevenbitpp(i,j)==0
            sevenzerosubpp= sevenzerosubpp+1;

        else
            sevennonzsubpp=sevennonzsubpp+1;
        end
    end
end


for i=1:bit5
    for j=1:bit5
        if seitbitpp(i,j)==0
            eitzerosubpp= eitzerosubpp+1;

        else
            eitnonzsubpp=eitnonzsubpp+1;
        end
    end
end


for i=1:bit6
    for j=1:bit6
        if sninebitpp(i,j)==0
            ninezerosubpp= ninezerosubpp+1;

        else
            ninenonzsubpp=ninenonzsubpp+1;
        end
    end
end

fractaldimfrbitsubpp=log(frnonzsubpp)/log(bit1); %fractaldimension of virus given in filename bit1
fractaldimfivebitsubpp=log(fivenonzsubpp)/log(bit2);
fractaldimsixbitsubpp=log(sixnonzsubpp)/log(bit3);
fractaldimsevenbitsubpp=log(sevennonzsubpp)/log(bit4);
fractaldimeitbitsubpp=log(eitnonzsubpp)/log(bit5);%same bit2
fractaldimninebitsubpp=log(ninenonzsubpp)/log(bit6);


fractdimematrix(4,1)=fractaldimfrbitsubpp;fractdimematrix(4,2)=fractaldimfivebitsubpp;
fractdimematrix(4,3)=fractaldimsixbitsubpp;fractdimematrix(4,4)=fractaldimsevenbitsubpp;
fractdimematrix(4,5)=fractaldimeitbitsubpp;fractdimematrix(4,6)=fractaldimninebitsubpp;


frbitn=rot90(frbit);
fivebitn=rot90(fivebit);
sixbitn=rot90(sixbit);
sevenbitn=rot90(sevenbit);
eitbitn=rot90(eitbit);
ninebitn=rot90(ninebit);


m=1;n=1;
frbitadz=zeros(frzero,2);frbitadnz=zeros(frnonz,2);






for i=1:bit1

    for j=1:bit1
        if frbitn(i,j)==0
            frbitadz(m,1)=i;
            frbitadz(m,2)=j;
            m=m+1;
        else
            frbitadnz(n,1)=i;
            frbitadnz(n,2)=j;
            n=n+1;
        end
    end 
end
m=1;n=1;
fivebitadz=zeros(fivezero,2);fivebitadnz=zeros(fivenonz,2);
for i=1:bit2

    for j=1:bit2
        if fivebitn(i,j)==0
            fivebitadz(m,1)=i;
            fivebitadz(m,2)=j;
            m=m+1;
        else
            fivebitadnz(n,1)=i;
            fivebitadnz(n,2)=j;
            n=n+1;
        end
    end 
end

m=1;n=1;
sixbitadz=zeros(sixzero,2);sixbitadnz=zeros(sixnonz,2);

for i=1:bit3

    for j=1:bit3
        if sixbitn(i,j)==0
            sixbitadz(m,1)=i;
            sixbitadz(m,2)=j;
            m=m+1;
        else
            sixbitadnz(n,1)=i;
            sixbitadnz(n,2)=j;
            n=n+1;
        end
    end 
end

m=1;n=1;
sevenbitadz=zeros(sevenzero,2);sevenbitadnz=zeros(sevennonz,2);
for i=1:bit4

    for j=1:bit4
        if sevenbitn(i,j)==0
            sevenbitadz(m,1)=i;
            sevenbitadz(m,2)=j;
            m=m+1;
        else
            sevenbitadnz(n,1)=i;
            sevenbitadnz(n,2)=j;
            n=n+1;
        end
    end 
end

m=1;n=1;
eitbitadz=zeros(eitzero,2);eitbitadnz=zeros(eitnonz,2);
for i=1:bit5

    for j=1:bit5
        if eitbitn(i,j)==0
            eitbitadz(m,1)=i;
            eitbitadz(m,2)=j;
            m=m+1;
        else
            eitbitadnz(n,1)=i;
            eitbitadnz(n,2)=j;
            n=n+1;
        end
    end
end

m=1;n=1;
ninebitadz=zeros(ninezero,2);ninebitadnz=zeros(ninenonz,2);

for i=1:bit6

    for j=1:bit6
        if ninebitn(i,j)==0
            ninebitadz(m,1)=i;
            ninebitadz(m,2)=j;
            m=m+1;
        else
            ninebitadnz(n,1)=i;
            ninebitadnz(n,2)=j;
            n=n+1;
        end
    end
end

frzbxad=zeros(frzero,b1);

for i=1:frzero
 
            frzbxad(i,:)=box_co(frbitadz(i,1),frbitadz(i,2),bit1);
end

frnzbxad=zeros(frnonz,b1);

for i=1:frnonz
  
            frnzbxad(i,:)=box_co(frbitadnz(i,1),frbitadnz(i,2),bit1);
end

fivezbxad=zeros(fivezero,b2);

for i=1:fivezero
 
            fivezbxad(i,:)=box_co(fivebitadz(i,1),fivebitadz(i,2),bit2);
end

fivenzbxad=zeros(fivenonz,b2);

for i=1:fivenonz
  
            fivenzbxad(i,:)=box_co(fivebitadnz(i,1),fivebitadnz(i,2),bit2);
end

sixzbxad=zeros(sixzero,b3);

for i=1:sixzero
 
            sixzbxad(i,:)=box_co(sixbitadz(i,1),sixbitadz(i,2),bit3);
end

sixnzbxad=zeros(sixnonz,b3);

for i=1:sixnonz
  
            sixnzbxad(i,:)=box_co(sixbitadnz(i,1),sixbitadnz(i,2),bit3);
end

sevenzbxad=zeros(sevenzero,b4);

for i=1:sevenzero
 
            sevenzbxad(i,:)=box_co(sevenbitadz(i,1),sevenbitadz(i,2),bit4);
end

sevennzbxad=zeros(sevennonz,b4);

for i=1:sevennonz
  
            sevennzbxad(i,:)=box_co(sevenbitadnz(i,1),sevenbitadnz(i,2),bit4);
end

eitzbxad=zeros(eitzero,b5);

for i=1:eitzero
    
            eitzbxad(i,:)=box_co(eitbitadz(i,1),eitbitadz(i,2),bit5);
end


eitnzbxad=zeros(eitnonz,b5);

for i=1:eitnonz
  
            eitnzbxad(i,:)=box_co(eitbitadnz(i,1),eitbitadnz(i,2),bit5);
end

ninezbxad=zeros(ninezero,b6);

for i=1:ninezero
 
            ninezbxad(i,:)=box_co(ninebitadz(i,1),ninebitadz(i,2),bit6);
end

ninenzbxad=zeros(ninenonz,b6);

for i=1:ninenonz
  
            ninenzbxad(i,:)=box_co(ninebitadnz(i,1),ninebitadnz(i,2),bit6);
end

frbnzv=zeros(frnonz,1);
for i=1:frnonz
    frbnzv(i)=frbitn(frbitadnz(i,1),frbitadnz(i,2));
end
frnzvpercent=frbnzv*100/niter;

fivebnzv=zeros(fivenonz,1);
for i=1:fivenonz
    fivebnzv(i)=fivebitn(fivebitadnz(i,1),fivebitadnz(i,2));
end
fivenzvpercent=fivebnzv*100/niter;

sixbnzv=zeros(sixnonz,1);
for i=1:sixnonz
    sixbnzv(i)=sixbitn(sixbitadnz(i,1),sixbitadnz(i,2));
end
sixnzvpercent=sixbnzv*100/niter;

sevenbnzv=zeros(sevennonz,1);
for i=1:sevennonz
    sevenbnzv(i)=sevenbitn(sevenbitadnz(i,1),sevenbitadnz(i,2));
end
sevennzvpercent=sevenbnzv*100/niter;

eitbnzv=zeros(eitnonz,1);
for i=1:eitnonz
    eitbnzv(i)=eitbitn(eitbitadnz(i,1),eitbitadnz(i,2));
end
eitnzvpercent=eitbnzv*100/niter;
ninebnzv=zeros(ninenonz,1);
for i=1:ninenonz
    ninebnzv(i)=ninebitn(ninebitadnz(i,1),ninebitadnz(i,2));
end
ninenzvpercent=ninebnzv*100/niter;



frbitn1=rot90(frbit1);
fivebitn1=rot90(fivebit1);
sixbitn1=rot90(sixbit1);
sevenbitn1=rot90(sevenbit1);
eitbitn1=rot90(eitbit1);
ninebitn1=rot90(ninebit1);


m=1;n=1;
frbitadz1=zeros(frzero1,2);frbitadnz1=zeros(frnonz1,2);






for i=1:bit1

    for j=1:bit1
        if frbitn1(i,j)==0
            frbitadz1(m,1)=i;
            frbitadz1(m,2)=j;
            m=m+1;
        else
            frbitadnz1(n,1)=i;
            frbitadnz1(n,2)=j;
            n=n+1;
        end
    end 
end
m=1;n=1;
fivebitadz1=zeros(fivezero1,2);fivebitadnz1=zeros(fivenonz1,2);
for i=1:bit2

    for j=1:bit2
        if fivebitn1(i,j)==0
            fivebitadz1(m,1)=i;
            fivebitadz1(m,2)=j;
            m=m+1;
        else
            fivebitadnz1(n,1)=i;
            fivebitadnz1(n,2)=j;
            n=n+1;
        end
    end 
end

m=1;n=1;
sixbitadz1=zeros(sixzero1,2);sixbitadnz1=zeros(sixnonz1,2);

for i=1:bit3

    for j=1:bit3
        if sixbitn1(i,j)==0
            sixbitadz1(m,1)=i;
            sixbitadz1(m,2)=j;
            m=m+1;
        else
            sixbitadnz1(n,1)=i;
            sixbitadnz1(n,2)=j;
            n=n+1;
        end
    end 
end

m=1;n=1;
sevenbitadz1=zeros(sevenzero1,2);sevenbitadnz1=zeros(sevennonz1,2);
for i=1:bit4

    for j=1:bit4
        if sevenbitn1(i,j)==0
            sevenbitadz1(m,1)=i;
            sevenbitadz1(m,2)=j;
            m=m+1;
        else
            sevenbitadnz1(n,1)=i;
            sevenbitadnz1(n,2)=j;
            n=n+1;
        end
    end 
end

m=1;n=1;
eitbitadz1=zeros(eitzero1,2);eitbitadnz1=zeros(eitnonz1,2);
for i=1:bit5

    for j=1:bit5
        if eitbitn1(i,j)==0
            eitbitadz1(m,1)=i;
            eitbitadz1(m,2)=j;
            m=m+1;
        else
            eitbitadnz1(n,1)=i;
            eitbitadnz1(n,2)=j;
            n=n+1;
        end
    end
end

m=1;n=1;
ninebitadz1=zeros(ninezero1,2);ninebitadnz1=zeros(ninenonz1,2);

for i=1:bit6

    for j=1:bit6
        if ninebitn1(i,j)==0
            ninebitadz1(m,1)=i;
            ninebitadz1(m,2)=j;
            m=m+1;
        else
            ninebitadnz1(n,1)=i;
            ninebitadnz1(n,2)=j;
            n=n+1;
        end
    end
end

frzbxad1=zeros(frzero1,b1);

for i=1:frzero1
 
            frzbxad1(i,:)=box_co(frbitadz1(i,1),frbitadz1(i,2),bit1);
end

frnzbxad1=zeros(frnonz1,b1);

for i=1:frnonz1
  
            frnzbxad1(i,:)=box_co(frbitadnz1(i,1),frbitadnz1(i,2),bit1);
end

fivezbxad1=zeros(fivezero1,b2);

for i=1:fivezero1
 
            fivezbxad1(i,:)=box_co(fivebitadz1(i,1),fivebitadz1(i,2),bit2);
end

fivenzbxad1=zeros(fivenonz1,b2);

for i=1:fivenonz1
  
            fivenzbxad1(i,:)=box_co(fivebitadnz1(i,1),fivebitadnz1(i,2),bit2);
end

sixzbxad1=zeros(sixzero1,b3);

for i=1:sixzero1
 
            sixzbxad1(i,:)=box_co(sixbitadz1(i,1),sixbitadz1(i,2),bit3);
end

sixnzbxad1=zeros(sixnonz1,b3);

for i=1:sixnonz1
  
            sixnzbxad1(i,:)=box_co(sixbitadnz1(i,1),sixbitadnz1(i,2),bit3);
end

sevenzbxad1=zeros(sevenzero1,b4);

for i=1:sevenzero1
 
            sevenzbxad1(i,:)=box_co(sevenbitadz1(i,1),sevenbitadz1(i,2),bit4);
end

sevennzbxad1=zeros(sevennonz1,b4);

for i=1:sevennonz1
  
            sevennzbxad1(i,:)=box_co(sevenbitadnz1(i,1),sevenbitadnz1(i,2),bit4);
end

eitzbxad1=zeros(eitzero1,b5);

for i=1:eitzero1
    
            eitzbxad1(i,:)=box_co(eitbitadz1(i,1),eitbitadz1(i,2),bit5);
end


eitnzbxad1=zeros(eitnonz1,b5);

for i=1:eitnonz1
  
            eitnzbxad1(i,:)=box_co(eitbitadnz1(i,1),eitbitadnz1(i,2),bit5);
end

ninezbxad1=zeros(ninezero1,b6);

for i=1:ninezero1
 
            ninezbxad1(i,:)=box_co(ninebitadz1(i,1),ninebitadz1(i,2),bit6);
end

ninenzbxad1=zeros(ninenonz1,b6);

for i=1:ninenonz1
  
            ninenzbxad1(i,:)=box_co(ninebitadnz1(i,1),ninebitadnz1(i,2),bit6);
end

frbnzv1=zeros(frnonz1,1);
for i=1:frnonz1
    frbnzv1(i)=frbitn1(frbitadnz1(i,1),frbitadnz1(i,2));
end
frnzvpercent1=frbnzv1*100/niter1;

fivebnzv1=zeros(fivenonz1,1);
for i=1:fivenonz1
    fivebnzv1(i)=fivebitn1(fivebitadnz1(i,1),fivebitadnz1(i,2));
end
fivenzvpercent1=fivebnzv1*100/niter1;

sixbnzv1=zeros(sixnonz1,1);
for i=1:sixnonz1
    sixbnzv1(i)=sixbitn1(sixbitadnz1(i,1),sixbitadnz1(i,2));
end
sixnzvpercent1=sixbnzv1*100/niter1;

sevenbnzv1=zeros(sevennonz1,1);
for i=1:sevennonz1
    sevenbnzv1(i)=sevenbitn1(sevenbitadnz1(i,1),sevenbitadnz1(i,2));
end
sevennzvpercent1=sevenbnzv1*100/niter1;

eitbnzv1=zeros(eitnonz1,1);
for i=1:eitnonz1
    eitbnzv1(i)=eitbitn1(eitbitadnz1(i,1),eitbitadnz1(i,2));
end

eitnzvpercent1=eitbnzv1*100/niter1;

ninebnzv1=zeros(ninenonz1,1);
for i=1:ninenonz1
    ninebnzv1(i)=ninebitn1(ninebitadnz1(i,1),ninebitadnz1(i,2));
end

ninenzvpercent1=ninebnzv1*100/niter1;

subfrbitn=rot90(sfrbit);
subfivebitn=rot90(sfivebit);
subsixbitn=rot90(ssixbit);
subsevenbitn=rot90(ssevenbit);
subeitbitn=rot90(seitbit);
subninebitn=rot90(sninebit);


m=1;n=1;
subfrbitadz=zeros(frzerosub,2);subfrbitadnz=zeros(frnonzsub,2);






for i=1:bit1

    for j=1:bit1
        if subfrbitn(i,j)==0
            subfrbitadz(m,1)=i;
            subfrbitadz(m,2)=j;
            m=m+1;
        else
            subfrbitadnz(n,1)=i;
            subfrbitadnz(n,2)=j;
            n=n+1;
        end
    end 
end
m=1;n=1;
subfivebitadz=zeros(fivezerosub,2);subfivebitadnz=zeros(fivenonzsub,2);
for i=1:bit2

    for j=1:bit2
        if subfivebitn(i,j)==0
            subfivebitadz(m,1)=i;
            subfivebitadz(m,2)=j;
            m=m+1;
        else
            subfivebitadnz(n,1)=i;
            subfivebitadnz(n,2)=j;
            n=n+1;
        end
    end 
end

m=1;n=1;
subsixbitadz=zeros(sixzerosub,2);subsixbitadnz=zeros(sixnonzsub,2);

for i=1:bit3

    for j=1:bit3
        if subsixbitn(i,j)==0
            subsixbitadz(m,1)=i;
            subsixbitadz(m,2)=j;
            m=m+1;
        else
            subsixbitadnz(n,1)=i;
            subsixbitadnz(n,2)=j;
            n=n+1;
        end
    end 
end

m=1;n=1;
subsevenbitadz=zeros(sevenzerosub,2);subsevenbitadnz=zeros(sevennonzsub,2);
for i=1:bit4

    for j=1:bit4
        if subsevenbitn(i,j)==0
            subsevenbitadz(m,1)=i;
            subsevenbitadz(m,2)=j;
            m=m+1;
        else
            subsevenbitadnz(n,1)=i;
            subsevenbitadnz(n,2)=j;
            n=n+1;
        end
    end 
end

m=1;n=1;
subeitbitadz=zeros(eitzerosub,2);subeitbitadnz=zeros(eitnonzsub,2);
for i=1:bit5

    for j=1:bit5
        if subeitbitn(i,j)==0
            subeitbitadz(m,1)=i;
            subeitbitadz(m,2)=j;
            m=m+1;
        else
            subeitbitadnz(n,1)=i;
            subeitbitadnz(n,2)=j;
            n=n+1;
        end
    end
end

m=1;n=1;
subninebitadz=zeros(ninezerosub,2);subninebitadnz=zeros(ninenonzsub,2);

for i=1:bit6

    for j=1:bit6
        if subninebitn(i,j)==0
            subninebitadz(m,1)=i;
            subninebitadz(m,2)=j;
            m=m+1;
        else
            subninebitadnz(n,1)=i;
            subninebitadnz(n,2)=j;
            n=n+1;
        end
    end
end

subfrzbxad=zeros(frzerosub,b1);

for i=1:frzerosub
 
            subfrzbxad(i,:)=box_co(subfrbitadz(i,1),subfrbitadz(i,2),bit1);
end

subfrnzbxad=zeros(frnonzsub,b1);

for i=1:frnonzsub
  
            subfrnzbxad(i,:)=box_co(subfrbitadnz(i,1),subfrbitadnz(i,2),bit1);
end

subfivezbxad=zeros(fivezerosub,b2);

for i=1:fivezerosub
 
            subfivezbxad(i,:)=box_co(subfivebitadz(i,1),subfivebitadz(i,2),bit2);
end

subfivenzbxad=zeros(fivenonzsub,b2);

for i=1:fivenonzsub
  
            subfivenzbxad(i,:)=box_co(subfivebitadnz(i,1),subfivebitadnz(i,2),bit2);
end

subsixzbxad=zeros(sixzerosub,b3);

for i=1:sixzerosub
 
            subsixzbxad(i,:)=box_co(subsixbitadz(i,1),subsixbitadz(i,2),bit3);
end

subsixnzbxad=zeros(sixnonzsub,b3);

for i=1:sixnonzsub
  
            subsixnzbxad(i,:)=box_co(subsixbitadnz(i,1),subsixbitadnz(i,2),bit3);
end

subsevenzbxad=zeros(sevenzerosub,b4);

for i=1:sevenzerosub
 
            subsevenzbxad(i,:)=box_co(subsevenbitadz(i,1),subsevenbitadz(i,2),bit4);
end

subsevennzbxad=zeros(sevennonzsub,b4);

for i=1:sevennonzsub
  
            subsevennzbxad(i,:)=box_co(subsevenbitadnz(i,1),subsevenbitadnz(i,2),bit4);
end

subeitzbxad=zeros(eitzerosub,b5);

for i=1:eitzerosub
    
            subeitzbxad(i,:)=box_co(subeitbitadz(i,1),subeitbitadz(i,2),bit5);
end


subeitnzbxad=zeros(eitnonzsub,b5);

for i=1:eitnonzsub
  
            subeitnzbxad(i,:)=box_co(subeitbitadnz(i,1),subeitbitadnz(i,2),bit5);
end

subninezbxad=zeros(ninezerosub,b6);

for i=1:ninezerosub
 
            subninezbxad(i,:)=box_co(subninebitadz(i,1),subninebitadz(i,2),bit6);
end

subninenzbxad=zeros(ninenonzsub,b6);

for i=1:ninenonzsub
  
            subninenzbxad(i,:)=box_co(subninebitadnz(i,1),subninebitadnz(i,2),bit6);
end

subfrbnzv=zeros(frnonzsub,1);
for i=1:frnonzsub
    subfrbnzv(i)=subfrbitn(subfrbitadnz(i,1),subfrbitadnz(i,2));
end
subfivebnzv=zeros(fivenonzsub,1);
for i=1:fivenonzsub
    subfivebnzv(i)=subfivebitn(subfivebitadnz(i,1),subfivebitadnz(i,2));
end

subsixbnzv=zeros(sixnonzsub,1);
for i=1:sixnonzsub
    subsixbnzv(i)=subsixbitn(subsixbitadnz(i,1),subsixbitadnz(i,2));
end


subsevenbnzv=zeros(sevennonzsub,1);
for i=1:sevennonzsub
    subsevenbnzv(i)=subsevenbitn(subsevenbitadnz(i,1),subsevenbitadnz(i,2));
end

subeitbnzv=zeros(eitnonzsub,1);
for i=1:eitnonzsub
    subeitbnzv(i)=subeitbitn(subeitbitadnz(i,1),subeitbitadnz(i,2));
end

subninebnzv=zeros(ninenonzsub,1);
for i=1:ninenonzsub
    subninebnzv(i)=subninebitn(subninebitadnz(i,1),subninebitadnz(i,2));
end


%%% non zero box adresses in ascending order col 1 adreess.. col2 numb of
%%% points and clo 3 percent of number of points%%%%%

letters='ABCD';

m=1;
frbitad=zeros(bit1^2,2);

for i=1:bit1

    for j=1:bit1
            frbitad(m,1)=i;
            frbitad(m,2)=j;
            m=m+1;
    
    end 
end

frbitbxad=zeros(bit1^2,b1);
for i=1:bit1^2
 
            frbitbxad(i,:)=box_co(frbitad(i,1),frbitad(i,2),bit1);
end

frbvlaue=zeros(bit1^2,1);
for i=1:bit1^2
    frbvlaue(i)=frbitn(frbitad(i,1),frbitad(i,2));
end
frbvaluepercent=frbvlaue*100/niter;


frbitbxadltr=letters(frbitbxad);
combinedfrbitbxad=cell(bit1^2,3);
for i= 1: bit1^2
    combinedfrbitbxad(i,1)={frbitbxadltr(i,:)};
   combinedfrbitbxad(i,2)={frbvlaue(i,:)};
   combinedfrbitbxad(i,3)={frbvaluepercent(i,:)};
end
combinedfrbitbxad = sortrows(combinedfrbitbxad,2);

m=1;
fivebitad=zeros(bit2^2,2);

for i=1:bit2

    for j=1:bit2
            fivebitad(m,1)=i;
            fivebitad(m,2)=j;
            m=m+1;
    
    end 
end

fivebitbxad=zeros(bit2^2,b2);
for i=1:bit2^2
 
            fivebitbxad(i,:)=box_co(fivebitad(i,1),fivebitad(i,2),bit2);
end

fivebvlaue=zeros(bit2^2,1);
for i=1:bit2^2
    fivebvlaue(i)=fivebitn(fivebitad(i,1),fivebitad(i,2));
end
fivebvaluepercent=fivebvlaue*100/niter;


fivebitbxadltr=letters(fivebitbxad);
combinedfivebitbxad=cell(bit2^2,3);
for i= 1:bit2^2
    combinedfivebitbxad(i,1)={fivebitbxadltr(i,:)};
   combinedfivebitbxad(i,2)={fivebvlaue(i,:)};
   combinedfivebitbxad(i,3)={fivebvaluepercent(i,:)};
end
combinedfivebitbxad = sortrows(combinedfivebitbxad,2);

m=1;
sixbitad=zeros(bit3^2,2);

for i=1:bit3

    for j=1:bit3
            sixbitad(m,1)=i;
            sixbitad(m,2)=j;
            m=m+1;
    
    end 
end

sixbitbxad=zeros(bit3^2,b3);
for i=1:bit3^2
 
            sixbitbxad(i,:)=box_co(sixbitad(i,1),sixbitad(i,2),bit3);
end

sixbvlaue=zeros(bit3^2,1);
for i=1:bit3^2
    sixbvlaue(i)=sixbitn(sixbitad(i,1),sixbitad(i,2));
end
sixbvaluepercent=sixbvlaue*100/niter;


sixbitbxadltr=letters(sixbitbxad);
combinedsixbitbxad=cell(bit3^2,3);
for i= 1:bit3^2
    combinedsixbitbxad(i,1)={sixbitbxadltr(i,:)};
   combinedsixbitbxad(i,2)={sixbvlaue(i,:)};
   combinedsixbitbxad(i,3)={sixbvaluepercent(i,:)};
end
combinedsixbitbxad = sortrows(combinedsixbitbxad,2);

m=1;
sevenbitad=zeros(bit4^2,2);

for i=1:bit4

    for j=1:bit4
            sevenbitad(m,1)=i;
            sevenbitad(m,2)=j;
            m=m+1;
    
    end 
end

sevenbitbxad=zeros(bit4^2,b4);
for i=1:bit4^2
 
            sevenbitbxad(i,:)=box_co(sevenbitad(i,1),sevenbitad(i,2),bit4);
end

sevenbvlaue=zeros(bit4^2,1);
for i=1:bit4^2
    sevenbvlaue(i)=sevenbitn(sevenbitad(i,1),sevenbitad(i,2));
end
sevenbvaluepercent=sevenbvlaue*100/niter;


sevenbitbxadltr=letters(sevenbitbxad);
combinedsevenbitbxad=cell(bit4^2,3);
for i= 1:bit4^2
    combinedsevenbitbxad(i,1)={sevenbitbxadltr(i,:)};
   combinedsevenbitbxad(i,2)={sevenbvlaue(i,:)};
   combinedsevenbitbxad(i,3)={sevenbvaluepercent(i,:)};
end
combinedsevenbitbxad = sortrows(combinedsevenbitbxad,2);

m=1;
eitbitad=zeros(bit5^2,2);

for i=1:bit5

    for j=1:bit5
            eitbitad(m,1)=i;
            eitbitad(m,2)=j;
            m=m+1;
    
    end 
end

eitbitbxad=zeros(bit5^2,b5);
for i=1:bit5^2
 
            eitbitbxad(i,:)=box_co(eitbitad(i,1),eitbitad(i,2),bit5);
end

eitbvlaue=zeros(bit5^2,1);
for i=1:bit5^2
    eitbvlaue(i)=eitbitn(eitbitad(i,1),eitbitad(i,2));
end
eitbvaluepercent=eitbvlaue*100/niter;


eitbitbxadltr=letters(eitbitbxad);
combinedeitbitbxad=cell(bit5^2,3);
for i= 1:bit5^2
    combinedeitbitbxad(i,1)={eitbitbxadltr(i,:)};
   combinedeitbitbxad(i,2)={eitbvlaue(i,:)};
   combinedeitbitbxad(i,3)={eitbvaluepercent(i,:)};
end
combinedeitbitbxad = sortrows(combinedeitbitbxad,2);
m=1;
ninebitad=zeros(bit6^2,2);

for i=1:bit6

    for j=1:bit6
            ninebitad(m,1)=i;
            ninebitad(m,2)=j;
            m=m+1;
    
    end 
end

ninebitbxad=zeros(bit6^2,b6);
for i=1:bit6^2
 
            ninebitbxad(i,:)=box_co(ninebitad(i,1),ninebitad(i,2),bit6);
end

ninebvlaue=zeros(bit6^2,1);
for i=1:bit6^2
    ninebvlaue(i)=ninebitn(ninebitad(i,1),ninebitad(i,2));
end
ninebvaluepercent=ninebvlaue*100/niter;


ninebitbxadltr=letters(ninebitbxad);
combinedninebitbxad=cell(bit6^2,3);
for i= 1:bit6^2
    combinedninebitbxad(i,1)={ninebitbxadltr(i,:)};
   combinedninebitbxad(i,2)={ninebvlaue(i,:)};
   combinedninebitbxad(i,3)={ninebvaluepercent(i,:)};
end
combinedninebitbxad = sortrows(combinedninebitbxad,2);



m=1;
frbitad1=zeros(bit1^2,2);

for i=1:bit1

    for j=1:bit1
            frbitad1(m,1)=i;
            frbitad1(m,2)=j;
            m=m+1;
    
    end 
end

frbitbxad1=zeros(bit1^2,b1);
for i=1:bit1^2
 
            frbitbxad1(i,:)=box_co(frbitad1(i,1),frbitad1(i,2),bit1);
end

frbvlaue1=zeros(bit1^2,1);
for i=1:bit1^2
    frbvlaue1(i)=frbitn1(frbitad1(i,1),frbitad1(i,2));
end
frbvaluepercent1=frbvlaue1*100/niter1;


frbitbxadltr1=letters(frbitbxad1);
combinedfrbitbxad1=cell(bit1^2,3);
for i= 1: bit1^2
    combinedfrbitbxad1(i,1)={frbitbxadltr1(i,:)};
   combinedfrbitbxad1(i,2)={frbvlaue1(i,:)};
   combinedfrbitbxad1(i,3)={frbvaluepercent1(i,:)};
end
combinedfrbitbxad1 = sortrows(combinedfrbitbxad1,2);

m=1;
fivebitad1=zeros(bit2^2,2);

for i=1:bit2

    for j=1:bit2
            fivebitad1(m,1)=i;
            fivebitad1(m,2)=j;
            m=m+1;
    
    end 
end

fivebitbxad1=zeros(bit2^2,b2);
for i=1:bit2^2
 
            fivebitbxad1(i,:)=box_co(fivebitad1(i,1),fivebitad1(i,2),bit2);
end

fivebvlaue1=zeros(bit2^2,1);
for i=1:bit2^2
    fivebvlaue1(i)=fivebitn1(fivebitad1(i,1),fivebitad1(i,2));
end
fivebvaluepercent1=fivebvlaue1*100/niter1;


fivebitbxadltr1=letters(fivebitbxad1);
combinedfivebitbxad1=cell(bit2^2,3);
for i= 1:bit2^2
    combinedfivebitbxad1(i,1)={fivebitbxadltr1(i,:)};
   combinedfivebitbxad1(i,2)={fivebvlaue1(i,:)};
   combinedfivebitbxad1(i,3)={fivebvaluepercent1(i,:)};
end
combinedfivebitbxad1 = sortrows(combinedfivebitbxad1,2);

m=1;
sixbitad1=zeros(bit3^2,2);

for i=1:bit3

    for j=1:bit3
            sixbitad1(m,1)=i;
            sixbitad1(m,2)=j;
            m=m+1;
    
    end 
end

sixbitbxad1=zeros(bit3^2,b3);
for i=1:bit3^2
 
            sixbitbxad1(i,:)=box_co(sixbitad1(i,1),sixbitad1(i,2),bit3);
end

sixbvlaue1=zeros(bit3^2,1);
for i=1:bit3^2
    sixbvlaue1(i)=sixbitn1(sixbitad1(i,1),sixbitad1(i,2));
end
sixbvaluepercent1=sixbvlaue1*100/niter1;


sixbitbxadltr1=letters(sixbitbxad1);
combinedsixbitbxad1=cell(bit3^2,3);
for i= 1:bit3^2
    combinedsixbitbxad1(i,1)={sixbitbxadltr1(i,:)};
   combinedsixbitbxad1(i,2)={sixbvlaue1(i,:)};
   combinedsixbitbxad1(i,3)={sixbvaluepercent1(i,:)};
end
combinedsixbitbxad1 = sortrows(combinedsixbitbxad1,2);

m=1;
sevenbitad1=zeros(bit4^2,2);

for i=1:bit4

    for j=1:bit4
            sevenbitad1(m,1)=i;
            sevenbitad1(m,2)=j;
            m=m+1;
    
    end 
end

sevenbitbxad1=zeros(bit4^2,b4);
for i=1:bit4^2
 
            sevenbitbxad1(i,:)=box_co(sevenbitad1(i,1),sevenbitad1(i,2),bit4);
end

sevenbvlaue1=zeros(bit4^2,1);
for i=1:bit4^2
    sevenbvlaue1(i)=sevenbitn1(sevenbitad1(i,1),sevenbitad1(i,2));
end
sevenbvaluepercent1=sevenbvlaue1*100/niter1;


sevenbitbxadltr1=letters(sevenbitbxad1);
combinedsevenbitbxad1=cell(bit4^2,3);
for i= 1:bit4^2
    combinedsevenbitbxad1(i,1)={sevenbitbxadltr1(i,:)};
   combinedsevenbitbxad1(i,2)={sevenbvlaue1(i,:)};
   combinedsevenbitbxad1(i,3)={sevenbvaluepercent1(i,:)};
end
combinedsevenbitbxad1 = sortrows(combinedsevenbitbxad1,2);

m=1;
eitbitad1=zeros(bit5^2,2);

for i=1:bit5

    for j=1:bit5
            eitbitad1(m,1)=i;
            eitbitad1(m,2)=j;
            m=m+1;
    
    end 
end

eitbitbxad1=zeros(bit5^2,b5);
for i=1:bit5^2
 
            eitbitbxad1(i,:)=box_co(eitbitad1(i,1),eitbitad1(i,2),bit5);
end

eitbvlaue1=zeros(bit5^2,1);
for i=1:bit5^2
    eitbvlaue1(i)=eitbitn1(eitbitad1(i,1),eitbitad1(i,2));
end
eitbvaluepercent1=eitbvlaue1*100/niter1;


eitbitbxadltr1=letters(eitbitbxad1);
combinedeitbitbxad1=cell(bit5^2,3);
for i= 1:bit5^2
    combinedeitbitbxad1(i,1)={eitbitbxadltr1(i,:)};
   combinedeitbitbxad1(i,2)={eitbvlaue1(i,:)};
   combinedeitbitbxad1(i,3)={eitbvaluepercent1(i,:)};
end
combinedeitbitbxad1 = sortrows(combinedeitbitbxad1,2);
m=1;
ninebitad1=zeros(bit6^2,2);

for i=1:bit6

    for j=1:bit6
            ninebitad1(m,1)=i;
            ninebitad1(m,2)=j;
            m=m+1;
    
    end 
end

ninebitbxad1=zeros(bit6^2,b6);
for i=1:bit6^2
 
            ninebitbxad1(i,:)=box_co(ninebitad1(i,1),ninebitad1(i,2),bit6);
end

ninebvlaue1=zeros(bit6^2,1);
for i=1:bit6^2
    ninebvlaue1(i)=ninebitn1(ninebitad1(i,1),ninebitad1(i,2));
end
ninebvaluepercent1=ninebvlaue1*100/niter1;


ninebitbxadltr1=letters(ninebitbxad1);
combinedninebitbxad1=cell(bit6^2,3);
for i= 1:bit6^2
    combinedninebitbxad1(i,1)={ninebitbxadltr1(i,:)};
   combinedninebitbxad1(i,2)={ninebvlaue1(i,:)};
   combinedninebitbxad1(i,3)={ninebvaluepercent1(i,:)};
end
combinedninebitbxad1 = sortrows(combinedninebitbxad1,2);


m=1;
frbitadsub=zeros(bit1^2,2);

for i=1:bit1

    for j=1:bit1
            frbitadsub(m,1)=i;
            frbitadsub(m,2)=j;
            m=m+1;
    
    end 
end

frbitbxadsub=zeros(bit1^2,b1);
for i=1:bit1^2
 
            frbitbxadsub(i,:)=box_co(frbitadsub(i,1),frbitadsub(i,2),bit1);
end

frbvlauesub=zeros(bit1^2,1);
for i=1:bit1^2
    frbvlauesub(i)=subfrbitn(frbitadsub(i,1),frbitadsub(i,2));
end



frbitbxadltrsub=letters(frbitbxadsub);
combinedfrbitbxadsub=cell(bit1^2,2);
for i= 1: bit1^2
    combinedfrbitbxadsub(i,1)={frbitbxadltrsub(i,:)};
   combinedfrbitbxadsub(i,2)={frbvlauesub(i,:)};
   
end
combinedfrbitbxadsub = sortrows(combinedfrbitbxadsub,2);

m=1;
fivebitadsub=zeros(bit2^2,2);

for i=1:bit2

    for j=1:bit2
            fivebitadsub(m,1)=i;
            fivebitadsub(m,2)=j;
            m=m+1;
    
    end 
end

fivebitbxadsub=zeros(bit2^2,b2);
for i=1:bit2^2
 
            fivebitbxadsub(i,:)=box_co(fivebitadsub(i,1),fivebitadsub(i,2),bit2);
end

fivebvlauesub=zeros(bit2^2,1);
for i=1:bit2^2
    fivebvlauesub(i)=subfivebitn(fivebitadsub(i,1),fivebitadsub(i,2));
end



fivebitbxadltrsub=letters(fivebitbxadsub);
combinedfivebitbxadsub=cell(bit2^2,2);
for i= 1:bit2^2
    combinedfivebitbxadsub(i,1)={fivebitbxadltrsub(i,:)};
   combinedfivebitbxadsub(i,2)={fivebvlauesub(i,:)};
   
end
combinedfivebitbxadsub = sortrows(combinedfivebitbxadsub,2);

m=1;
sixbitadsub=zeros(bit3^2,2);

for i=1:bit3

    for j=1:bit3
            sixbitadsub(m,1)=i;
            sixbitadsub(m,2)=j;
            m=m+1;
    
    end 
end

sixbitbxadsub=zeros(bit3^2,b3);
for i=1:bit3^2
 
            sixbitbxadsub(i,:)=box_co(sixbitadsub(i,1),sixbitadsub(i,2),bit3);
end

sixbvlauesub=zeros(bit3^2,1);
for i=1:bit3^2
    sixbvlauesub(i)=subsixbitn(sixbitadsub(i,1),sixbitadsub(i,2));
end



sixbitbxadltrsub=letters(sixbitbxadsub);
combinedsixbitbxadsub=cell(bit3^2,2);
for i= 1:bit3^2
    combinedsixbitbxadsub(i,1)={sixbitbxadltrsub(i,:)};
   combinedsixbitbxadsub(i,2)={sixbvlauesub(i,:)};
   
end
combinedsixbitbxadsub = sortrows(combinedsixbitbxadsub,2);

m=1;
sevenbitadsub=zeros(bit4^2,2);

for i=1:bit4

    for j=1:bit4
            sevenbitadsub(m,1)=i;
            sevenbitadsub(m,2)=j;
            m=m+1;
    
    end 
end

sevenbitbxadsub=zeros(bit4^2,b4);
for i=1:bit4^2
 
            sevenbitbxadsub(i,:)=box_co(sevenbitadsub(i,1),sevenbitadsub(i,2),bit4);
end

sevenbvlauesub=zeros(bit4^2,1);
for i=1:bit4^2
    sevenbvlauesub(i)=subsevenbitn(sevenbitadsub(i,1),sevenbitadsub(i,2));
end


sevenbitbxadltrsub=letters(sevenbitbxadsub);
combinedsevenbitbxadsub=cell(bit4^2,2);
for i= 1:bit4^2
    combinedsevenbitbxadsub(i,1)={sevenbitbxadltrsub(i,:)};
   combinedsevenbitbxadsub(i,2)={sevenbvlauesub(i,:)};
  
end
combinedsevenbitbxadsub = sortrows(combinedsevenbitbxadsub,2);

m=1;
eitbitadsub=zeros(bit5^2,2);

for i=1:bit5

    for j=1:bit5
            eitbitadsub(m,1)=i;
            eitbitadsub(m,2)=j;
            m=m+1;
    
    end 
end

eitbitbxadsub=zeros(bit5^2,b5);
for i=1:bit5^2
 
            eitbitbxadsub(i,:)=box_co(eitbitadsub(i,1),eitbitadsub(i,2),bit5);
end

eitbvlauesub=zeros(bit5^2,1);
for i=1:bit5^2
    eitbvlauesub(i)=subeitbitn(eitbitadsub(i,1),eitbitadsub(i,2));
end



eitbitbxadltrsub=letters(eitbitbxadsub);
combinedeitbitbxadsub=cell(bit5^2,2);
for i= 1:bit5^2
    combinedeitbitbxadsub(i,1)={eitbitbxadltrsub(i,:)};
   combinedeitbitbxadsub(i,2)={eitbvlauesub(i,:)};
  
end
combinedeitbitbxadsub = sortrows(combinedeitbitbxadsub,2);
m=1;
ninebitadsub=zeros(bit6^2,2);

for i=1:bit6

    for j=1:bit6
            ninebitadsub(m,1)=i;
            ninebitadsub(m,2)=j;
            m=m+1;
    
    end 
end

ninebitbxadsub=zeros(bit6^2,b6);
for i=1:bit6^2
 
            ninebitbxadsub(i,:)=box_co(ninebitadsub(i,1),ninebitadsub(i,2),bit6);
end

ninebvlauesub=zeros(bit6^2,1);
for i=1:bit6^2
    ninebvlauesub(i)=subninebitn(ninebitadsub(i,1),ninebitadsub(i,2));
end

ninebitbxadltrsub=letters(ninebitbxadsub);
combinedninebitbxadsub=cell(bit6^2,2);
for i= 1:bit6^2
    combinedninebitbxadsub(i,1)={ninebitbxadltrsub(i,:)};
   combinedninebitbxadsub(i,2)={ninebvlauesub(i,:)};
  
end
combinedninebitbxadsub = sortrows(combinedninebitbxadsub,2);


%matrix with percentage % first row for filename and 2nd row for filename 1
percentmatr(1,1)=percentno1(1,1);
percentmatr(1,2)=percentno2(1,1);
percentmatr(1,3)=percentno3(1,1);
percentmatr(1,4)=percentno4(1,1);

percentmatr(2,1)=percent1noo1(1,1);
percentmatr(2,2)=percent1noo2(1,1);
percentmatr(2,3)=percent1noo3(1,1);
percentmatr(2,4)=percent1noo4(1,1);


