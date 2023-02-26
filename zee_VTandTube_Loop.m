function [Z] = zee_VTandTube_Loop(f0s,dL,Area,Lstraw_loop,Astraw_loop)
%dL is tubelet length of vocal tract, single value assumed
%Area is vector of vocal tract tubelet areas
%Lstraw_loop is vector (or single number) of straw lengths (not related to dL or Area)
%Astraw_loop is vector (or single number) of straw areas(not related to dL or Area)
%f0s is vector (or single number) of fundamental frequencies to investigate

%output Z is an array of size length(Astraw_loop) x length(Lstraw_loop) x length(f0s)


%% some default paramaters
rho = 0.001225;%
c = 34300;
f = f0s; 
w = 2*pi*f;
k = w/c;

%% parameters for VT from Story 2000
Fw = 15;
FT = 200;
q=4;
alpha = sqrt(1j*w*q);
r = 408;
beta = (1j*w*(2*pi*FT).^2)./((1j*w+r).*1j.*w+(2*pi*Fw).^2) + alpha;
gamma = sqrt((r+1j*w)./(beta+1j*w));
sigma = gamma.*(beta+1j*w);

%% loop over tube segments of vocal tract here

Kn = [ones(1,1,size(f,2)),zeros(1,1,size(f,2));...
    zeros(1,1,size(f,2)),ones(1,1,size(f,2))];
Knn= Kn;

for ind = 1:length(Area)
    %sentacsKn+1
    A = reshape( cosh(sigma*dL/c), [1,1,size(f,2)] );
    B = reshape( -rho*c/Area(ind)*gamma.*sinh(sigma*dL/c), [1,1,size(f,2)] );
    C = reshape( -(Area(ind))./(rho*c*gamma).*sinh(sigma*dL/c),  [1,1,size(f,2)] );
    D = reshape( cosh(sigma*dL/c), [1,1,size(f,2)] );

    %Kn
    E = Kn(1,1,:);
    F = Kn(1,2,:);
    G = Kn(2,1,:);
    H = Kn(2,2,:);

    %matrix math (Kn+1)*(Kn), Story eqn 8
    Knn(1,1,:) = A.*E + B.*G;
    Knn(1,2,:) = A.*F + B.*H;
    Knn(2,1,:) = C.*E + D.*G;
    Knn(2,2,:) = C.*F + D.*H;

    %re-assign for next loop
    Kn = Knn;
end

%store VT to add straws to in next cell
KnVT = Kn;


%% add straws here

Fw = 15;
FT = 200;
q=4;
alpha = 0;%SET TO ZERO IN STRAW
r = 0; %%SET TO ZERO IN STRAW
beta = 0; %SET TO ZERO IN STRAW
gamma = sqrt((r+1j*w)./(beta+1j*w));
sigma = gamma.*(beta+1j*w);

%Impedance array
Z = zeros(length(Astraw_loop),length(Lstraw_loop),length(f0s)); %impedance array at select f0s

for ind=1:length(Astraw_loop)
    for jnd=1:length(Lstraw_loop)
        Astraw = Astraw_loop(ind);
        Lstraw = Lstraw_loop(jnd);
        
        A = reshape( cosh(sigma*Lstraw/c), [1,1,size(f,2)] );
        B = reshape( -rho*c/Astraw*gamma.*sinh(sigma*Lstraw/c), [1,1,size(f,2)] );
        C = reshape( -(Astraw)./(rho*c*gamma).*sinh(sigma*Lstraw/c),  [1,1,size(f,2)] );
        D = reshape( cosh(sigma*Lstraw/c), [1,1,size(f,2)] );

        %Kn
        E = KnVT(1,1,:);
        F = KnVT(1,2,:);
        G = KnVT(2,1,:);
        H = KnVT(2,2,:);

        %matrix math (Kn+1)*(KnVT), Story eqn 8
        Knn(1,1,:) = A.*E + B.*G;
        Knn(1,2,:) = A.*F + B.*H;
        Knn(2,1,:) = C.*E + D.*G;
        Knn(2,2,:) = C.*F + D.*H;

        %re-assign
        Kn = Knn;

        %add radiation impedance from end of tube
        b = sqrt(Astraw/pi);
        Zm = (rho*c)/(Astraw);
        R = 128*Zm/(9*pi^2);
        L = (8*b*Zm)./(3*pi*c);
        ZL = (1j*w.*R.*L)./(R+1j*w*L);
        ZL = reshape(ZL,[1,1,size(f,2)]);
        
        %impedance for VT and this specific straw
        Zin = (Kn(2,2,:).*ZL - Kn(1,2,:))./(Kn(1,1,:)-Kn(2,1,:).*ZL);

        %combine into impedance array for all straws, 3rd dimen is freq
        Z(ind,jnd,:) = Zin;
       
    end
end


end

