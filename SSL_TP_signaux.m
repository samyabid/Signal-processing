%% Code Matalb TP1 SSL - Transformation élémentaire des Signaux
% PG - 2019/2020

%%
clear, close all

%% Construction vecteur temps
N = 2000 ;
D = 5 ;
dt = (2*D)/N        % Intervalle entre deux points consecutifs
t = -D:dt:D-dt ;    % construction du vecteur temps
% variante possible : 
% t = linspace(-D,D,N+1);   % construction d'un vecteur de valeurs entre 
                            % -D et +D (inclus) regulierement espacées sur
                            % N+1 points
% t = t(1:N) ;              % Il faut exclure le dernier point car
                            % l'intervalle demandé est ouvert à droite

%% Construction porte de largeur T centree en zero

T = 2 ;                     % Largeur de la porte
s = (t>=-T/2 & t<=T/2) ;    % Construction en utilisant un operateur booleen
                            % s(n) = 1 si t(n) >= -T/2 et t(n) <= T/2
                            % s(n) = 0 sinon
plot(t,s,'linewidth',2)     % Affichage du signal porte en fonction du 
                            % vecteur temps en secondes
title('Porte_{T}(t)')       % Toujours indiquer le titre des Figures
xlabel('temps (secondes)')  % et les noms (et unités) des variables sur les axes 

%% Transformee de Fourier de la porte

subplot(211)                    % Partage de la figure en sous-figures
[S,f] = TransFourier(s,t) ;     % Calcul de la transformée de Fourier de la porte
plot(f,real(S),...              % Affichage SUPERPOSE de la partie réelle de la TF
    f,imag(S),...               % de la partie imaginaire de la TF
    f(f>0),(pi*f(f>0)).^(-1),...% de l'enveloppe théorique de la TF pour les f>0 uniquement
    '--','linewidth',2)         % UTILISER PLOT.M AVEC PLUSIEURS ARGUMENTS (PLUTOT QUE LA
                                % COMMANDE 'HOLD ON / HOLD OFF'
legend({'partie reelle',...     % Affichage d'une légende pour discriminer les 
    'partie imaginaire',...     % différentes courbes affichées sur la même figure
    'enveloppe'})
axis([-10 10 -Inf Inf])         % Limitation des axes (X uniquement ici) pour zoomer sur 
                                % la partie intéressante des courbes
grid                            % Superposition d'une grille pour faciliter la lecture 
                                % des valeurs
title('Transformée de Fourier de la Porte_T(t)')
xlabel('Frequence (Hz)')

subplot(212)
plot(f,abs(S).^2,'linewidth',2)
axis([-10 10 -Inf Inf])
grid
title('\Gamma_p(f) : Densité spectrale d''energie de la Porte_T(t)')
xlabel('Frequence (Hz)')

print -dpng TP1-Fig1            % Enregistrement de la Figure dans un fichier au format pbng 
                                % (>> help print pour voir les autres formats possibles)
                                
% print -dpng /Users/pgoncalv/CPE/TP/3ETI-SSL/TP1_Signaux/TP1-Fig1 

%% Transformee de Fourier de la porte comprimee d'un facteur a

a = 2 ; 
aT = T/a ;                          % La fonction P(t) est changée d'échelle 
                                    % d'un facteur inverse 1/a
                                    % On calcule donc les nouvelles bornes
s_dil = t>=-aT/2 & t<=aT/2 ;        % La porte P(at) est synthétisée à partir 
                                    % du vecteur temps 't' inchangé
subplot(211)
plot(t,s,'--b',t,s_dil,'r','linewidth',2)
axis([-2 2 -Inf Inf])
legend({'P_T(t)','P_T(at)'})
title('Changement d''echelle de la porte')
xlabel('temps (sec.)')

[S_dil,f] = TransFourier(s_dil,t) ;
subplot(223)
plot(f,abs(S).^2,'b',f,abs(S_dil).^2,'r','linewidth',2)
axis([-4 4 -Inf Inf]), grid
title('DSE - echelle lineaire'), xlabel('frequence (Hz)')
legend({'P_t(t)','P_T(at)'})

subplot(224)
semilogy(f,abs(S).^2,'b',f,abs(S_dil).^2,'r','linewidth',1)  % Affichage logarthmique
axis([-4 4 10e-5 4]), grid
title('DSE - echelle logarithmique'), xlabel('frequence (Hz)')
legend({'P_t(t)','P_T(at)'})

print -dpng TP1-Fig3

%% Transformee de Fourier de la porte translatee

t0 = 2 ;
s_shift = t>=-T/2 + t0 & t<= T/2 + t0 ;         % Meme construction que pour la Porte
subplot(211)
plot(t,s,'b--',t,s_shift,'r','linewidth',2)
legend({'P_T(t)','P_T(t-t_0)'})
title('Translation de la porte Porte_T(t)')
xlabel('temps (sec.)')
[S_shift,f] = TransFourier(s_shift,t) ;
subplot(223)
plot(f,real(S),f,real(S_shift),'linewidth',2)
axis([-3 3 -Inf Inf])
legend({'P_T(t)','P_T(t-t_0)'})
title('TF (partie reelle)')
xlabel('frequence (Hz)')
subplot(224)
plot(f,imag(S),f,imag(S_shift),'linewidth',2)
axis([-3 3 -Inf Inf])
legend({'P_T(t)','P_T(t-t_0)'})
title('TF (partie imag.)')
xlabel('frequence (Hz)')

print -dpng /Users/pgoncalv/CPE/TP/3ETI-SSL/TP1_Signaux/TP1-Fig2

%% Multiplication de la porte par un signal sinusoidal

f0 = 20 ; A0 = 3 ; 
s2 = s.*(A0*cos(2*pi*f0*t)) ; 
clf, subplot(211)
plot(t,s,t,s2,'linewidth',1)
axis([-2 2 -Inf Inf])
legend({'P_T(t)','P_T(t) modulé'})
title('P_T(t) x [A cos(2\pi f_0 t)]'), xlabel('temps (sec.)')
[S2,f] = TransFourier(s2,t) ;
subplot(212)
plot(f,real(S2),f,imag(S2))
axis([-40 40 -Inf Inf])
title('Transformee de Fourier de \{ P_T(t) x [A cos(2\pi f_0 t)] \}'),
xlabel('frequence (Hz)')

print -dpng /Users/pgoncalv/CPE/TP/3ETI-SSL/TP1_Signaux/TP1-Fig4

%% Calcul de l'energie 

% a partir du signal en temps
Est = dt*sum(abs(s2).^2) ;          % Integration par la méthode des rectangles
                                    % hauteur = |s2(n)|^2
                                    % largeur = dt : largeur intervalle temporel
                                    
disp(['Energie calculee sur le signal en temps = ',num2str(Est),' (J)'])

% a partir du signal en frequence

df = f(2)-f(1) ;                    % calcul largeur intervalle frequentiel
Esf = df*sum(abs(S2).^2) ;          % Integration par la methode des rectangles

disp(['Energie calculee sur le signal en frequence = ',num2str(Esf),' (J)'] )







