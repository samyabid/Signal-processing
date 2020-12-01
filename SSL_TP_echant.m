%% Codes du TP 4 SSL sur echantillonnage
% PG : 2019 / 2020

%% function in GeneSinus.m
% function [K,s] = GeneSinus(A,nu0,phi0,D,Fe)
% 
% M = D*Fe ;
% K = 0:M-1 ;
% t = K./Fe ;
% s = A*sin(2*pi*nu0*t+phi0)

%% Main program

nu0 = 7 ;
A = 10 ;
phi0 = 45 ;
D = 25*(1/nu0) ;
Fe = 19*nu0 ;

phi0_rad = phi0/180 * pi ;
[K,s] = GeneSinus(A,nu0,phi0_rad,D,Fe) ;
M = length(K) ;

figure(1), clf
subplot(211)
stem(K,s)
axis([0 M-1 -Inf Inf])
xlabel(['indices k = 0,1,...',num2str(M-1)])
title('s[k]')

subplot(212)
plot(K./Fe,s,K./Fe,s,'.','markersize',18)
axis([0 (M-1)/Fe -Inf +Inf])
xlabel(['temps : 0 \leq kT_e < ,',num2str((M)/Fe),' s.'])
title('s(t=kT_e)')

%% Effet de la frequence d'echantillonnage

nu0 = 800 ;
A = 1 ;
phi0 = 0 ;
D = 10 ;
phi0_rad = phi0/180*pi ;

figure(2), clf
j=1 ;
for Fe = [1e4 3e3 1200 600] ;
    [K,s] = GeneSinus(A,nu0,phi0_rad,D,Fe) ;
    t = K./Fe ;
%    [S,f] = TransFourierPower(s,t) ;
    [S,f] = TransFourier2(s,t,Fe) ;
    subplot(4,2,j)
    plot(t,s,t,s,'.','markersize',18)
    axis([0 10/nu0 -Inf Inf])
    title(['s(kT_e) avec F_e = ',num2str(Fe),'Hz']), 
    xlabel('kT_e')
    subplot(4,2,j+1) ;
    plot(f,abs(S),f,abs(S),'.'), grid
    axis([-Fe Fe -Inf Inf])
    xlabel(['Frequence ',num2str(-Fe),' \leq f \leq ',num2str(Fe)])
    title('|S(f)|')
    j=j+2;
end

