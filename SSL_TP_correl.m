%% TP 2 SSL - Correlation
% PG : 2019/2020

%%
clear, close all

%% Synthese des signaux 
% x : signal emis
% y : signal-echo reçu

f0 = 7 ;
Fs = 1000 ;

% sans bruit : 
    T = 5/f0 ;
% avec bruit (sigma = 10) : 
    T = 25/f0 ;  
% avec Doppler sans bruit: 
    T = 3/f0


D = 3*T ;
t0 = 1.5*T ;
t = 0:1/Fs:D-1/Fs ;


x = cos(2*pi*f0*t).*( (t>=0) & (t<T) ) ;
y = cos(2*pi*f0*(t-t0)).*( (t>= t0) & (t<T+t0) ) ;

%% Echo sans bruit 

[Ryx,lag] = xcorr(y,x) ; % attention a l'ordre des signaux!
Ryx = Ryx./Fs ;
lag = lag./Fs ;
figure(1), clf
subplot(211)
plot(t,x,t,y) ;
axis([0 max(t) -Inf Inf])
legend({'signal emis','echo recu'})
xlabel('temps (sec.)')
print -depsc corr1

figure(2), clf
subplot(212)
plot(lag,Ryx)
axis([-D D -inf inf])
title('Fonction d''inter-correlation \gamma_{y,x}(\tau)')
xlabel('retard \tau (sec.)')
imax = find(abs(Ryx) == max(abs(Ryx))) ;
t0_est = lag(imax) ; 
disp(['Retard = ',num2str(t0),' - Retard estime = ',num2str(t0_est),...
    ' (',num2str(abs(t0-t0_est)/t0*100),' %)'])
print -depsc corr2

%% Echo avec bruit

sigma = 10 ;
n = sigma*randn(1,length(t)) ;
z = y+n ;
[Rzx,lag] = xcorr(z,x,'coeff') ; % attention a l'ordre des signaux!
lag = lag./Fs ;
figure(1), clf
subplot(211)
plot(t,z,'r',t,x) ;
axis([0 max(t) -Inf Inf])
legend({'echo recu','signal emis'})
xlabel('temps (sec.)')
print -depsc corr3

figure(2), clf
subplot(212)
plot(lag,Rzx)
axis([-D D -inf inf])
title('Fonction d''inter-correlation \gamma_{z,x}(\tau)')
xlabel('retard \tau (sec.)')
print -depsc corr4

imax = find(abs(Rzx) == max(abs(Rzx))) ;
t0_est = lag(imax) ;
disp(['Retard = ',num2str(t0),' - Retard estime = ',num2str(t0_est),...
    ' (',num2str(abs(t0-t0_est)/t0*100),' %)'])


%% Echo avec effet Doppler

df = 0.13*f0 ;
w = cos(2*pi*(f0+df)*(t-t0)).*( (t>= t0) & (t<T+t0) ) ;
[Rwx,lag] = xcorr(w,x,'coeff') ; % attention a l'ordre des signaux!
lag = lag./Fs ;
figure(1), clf
subplot(211)
plot(t,w,'r',t,x) ;
axis([0 max(t) -Inf Inf])
legend({'echo recu','signal emis'})
xlabel('temps (sec.)')

print -depsc corr5

figure(2), clf
subplot(212)
plot(lag,Rwx)
axis([-D D -inf inf])
title('Fonction d''inter-correlation \gamma_{w,x}(\tau)')
xlabel('retard \tau (sec.)')
print -depsc corr6

imax = find(abs(Rwx) == max(abs(Rwx))) ;
t0_est = lag(imax) ;
disp(['Retard = ',num2str(t0),' - Retard estime = ',num2str(t0_est),...
    ' (',num2str(abs(t0-t0_est)/t0*100),' %)'])

figure(3)
[X,f] = TransFourier(x,t) ;
[W,f] = TransFourier(w,t) ;
subplot(211)
plot(f,abs(X).^2,f,abs(W).^2)
axis([0 3*f0 -Inf Inf])
xlabel('frequence (Hz)')
legend({'|X(f)|^{2}','|W(f)|^{2}'})
subplot(212)
plot(f,abs(X.*conj(W)))
axis([0 3*f0 -Inf Inf])
xlabel('frequence (Hz)')
title('Densite interspectrale \Gamma_{W,X}(f)')
print -depsc corr7
