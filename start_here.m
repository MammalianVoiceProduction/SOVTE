%% all the codes use cm, dyne, sec unit system.


%% bring in area function

%this function reads in Story 2018 vocal tract shapes
%see voweltypes.png in this folder

%this reads in Male, 18yo, vowel 10 (which is /u/)
%output is area function vector (Ax), element length (Li), and total length
%of vocal tract (Lvt)

[Ax,Li,Lvt] = VTreader('M',18,10);

%axial length of vocal tract
x = (0:length(Ax)-1)*Li;

figure(1)
plot(x,Ax);
xlabel('axial distance (cm)')
ylabel('vocal tract area (cm^2)')

%% Zin of just the vocal tract

f = 1:1000; %freq, Hz

%input frequency, Li, Ax. Output is input impedance of supraglottal tract
[Zvt] = zee_VTonly(f,Li,Ax);

figure(2)
subplot(2,1,1)
semilogy(f,abs(Zvt))
xlabel('f (Hz)'), ylabel('|Zin|')
subplot(2,1,2)
plot(f,rad2deg(angle(Zvt)))
xlabel('f (Hz)'), ylabel('angle(Zin)')

%% Zin with some straws

%now input vectors of straw lengths and areas to try out. Output Zstraw is
%output Zstraw is an array of size length(Astraw) x length(Lstraw) x length(f)

Lstraws = [0:5:30]; %cm, straw length
Astraws = pi*[0.4].^2; %cm^2, can be a vector but then is a lot to plot

[Zstraws] = zee_VTandTube_Loop(f,Li,Ax,Lstraws,Astraws);

%first diameter
figure(2)
subplot(2,1,1)
semilogy(f,squeeze(abs(Zstraws))) %have to use "squeeze" because "f" corresponds to third dimension of Zstraws
xlabel('f (Hz)'), ylabel('|Zin|')
subplot(2,1,2)
plot(f,squeeze(rad2deg(angle(Zstraws))))
xlabel('f (Hz)'), ylabel('angle(Zin)')

%%







