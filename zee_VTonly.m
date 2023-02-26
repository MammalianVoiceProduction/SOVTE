function [Zin] = zee_VTonly(f,dL,Area)
%uses story2000 paper, typo in his paper, referenced Sondhi paper is better

%% some default paramaters
rho = 0.001225;%;
c = 34300;
% f = 1:1:2500;
w = 2*pi*f;
k = w/c;

%% parameters for VT from Story
Fw = 15;
FT = 200;
q=4;
alpha = sqrt(1j*w*q);
r = 408;
beta = (1j*w*(2*pi*FT).^2)./((1j*w+r).*1j.*w+(2*pi*Fw).^2) + alpha;
gamma = sqrt((r+1j*w)./(beta+1j*w));
sigma = gamma.*(beta+1j*w);
%% loop over tube segments

Kn = [ones(1,1,size(f,2)),zeros(1,1,size(f,2));...
    zeros(1,1,size(f,2)),ones(1,1,size(f,2))];
Knn= Kn;

for ind = 1:length(Area)
    %Kn+1
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


%% add radiation impedance from end of mouth (only if there is no tube)
b = sqrt(Area(end)/pi);
Zm = (rho*c)/(Area(end));
R = 128*Zm/(9*pi^2);
L = (8*b*Zm)./(3*pi*c);
ZL = (1j*w.*R.*L)./(R+1j*w*L);
ZL = reshape(ZL,[1,1,size(f,2)]);

Zin = (Kn(2,2,:).*ZL - Kn(1,2,:))./(Kn(1,1,:)-Kn(2,1,:).*ZL);
Zin = reshape(Zin,[1,size(f,2)]);




end

