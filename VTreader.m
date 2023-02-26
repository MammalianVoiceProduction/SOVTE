function [Ax,Li,Lvt] = VTreader(Sex,Age,voweltype)
%Input sex as string 'M' or 'F'
%Input age as numeral 0,4,8,12,18 (which is adult)

%load male or female
if Sex == 'M'
    data = readmatrix('Story2018data - Male.csv');
elseif Sex == 'F'
    data = readmatrix('Story2018data - Female.csv');
else
    disp('input string M or F')
end

%load ages
if Age==18
    cols = 2;
elseif Age==12
    cols = 5;
elseif Age==8
    cols = 8;
elseif Age==4
    cols = 11;
elseif Age==0
    cols = 14;
else
    disp('Age should be numeral 0, 4, 8, 12, or 18')
end

%vowel coefficients from table
VowelData = readmatrix('Story2018data - VowelCoeffs.csv');
q1 = VowelData(voweltype,1);
q2 = VowelData(voweltype,2);

%build function, eqn1 story2018
N =     data(3:end-2,1);
Omega = data(3:end-2,cols);
Phi1 =  data(3:end-2,cols+1);
Phi2 =  data(3:end-2,cols+2);
Li =    data(end-1,cols);
Lvt=    data(end,cols);

Ax = pi/4*(Omega + q1*Phi1 + q2*Phi2).^2; %cm, see figure 1a, eqn 1







