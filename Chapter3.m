%Lab Orbital mechanics (lezione 23/11/2022)
%% Esercizio 1
%Pulizia
clc;
clear;
close all;

%Dati
r1 = [-21800 ; 37900 ; 0]; %initial position vector [km] (istante t1)
r2 = [ 27300 ; 27700 ; 0]; %final position vector [km] (istante t2)
Tof=54400; %15 h, 6 min, 40 s --> in secondi;sarebbe ToF (time of flight
              %(transfer time)= deltaT=t2-t1)
mu_E=astroConstants(13);
J2=astroConstants(9);
Re=astroConstants(23);

%gli input sono r1, r2, ToF, muE e altri parametri (posti pari a 0 (ho
%visto dalla funzione "Call_lamertMR.m"))
[a,p,e,ERROR,v1,v2,TPAR,THETA] = lambertMR( r1, r2, Tof, mu_E, 0, 0, 0, 0);
% Gli Output sono:
%	a        semi-major axis of the transfer orbit
% 	p        semi-latus rectum of the transfer orbit
%  	e        eccentricity of the transfer orbit
%	ERROR	Error flag
%                   0:	No error
%	v1       vector containing the initial velocity vector in Cartesian
%               coordinates --> v(t1)
%	v2		 vector containing the final velocity vector in Cartesian
%               coordinates --> v(t2)
%	TPAR 	 parabolic flight time between r1 and r2
%	THETA	 transfer angle [radians]

%Ora, grazie al Lambert solver, ho ottenuto v1 e v2 --> considero o (r1,v1)
% o (r2,v2) e, ricopiando la funzione Ode2bp, effettuo la propagazione
% dell'orbita partendo da una delle due (come condizioni iniziali)

T=2*pi*sqrt(a^3/mu_E); %periodo dell'orbita
options=odeset('RelTol',1e-13,'AbsTol',1e-14);
tspan=linspace(0,T,1000);

n = input('Enter a number: ');
%Caso 1 per le condizioni iniziali (r1,v1)
%Caso 2 per le condizioni iniziali (r2,v2)
switch n
    case 1 %(propagation forward in time; ma è la stessa orbita)
        y0=[r1;v1']; %vettore delle condizioni iniziali(deve essere colonna!)
                       %r1 è già colonna, mentre v1 è da trasporre
    case 2 %(propagation backward in time, ma è la stessa orbita)
        y0=[r2;v2']; %vettore delle condizioni iniziali(deve essere colonna!)
                       %r2 è già colonna, mentre v2 è da trasporre
end
[t y]=ode45(@(t,y)perturbed_ode_2bp( t, y, mu_E, J2, Re),tspan,y0,options);
figure(1);
earth_sphere();
hold on;
plot3(y(:,1),y(:,2),y(:,3),'-b','LineWidth',2); %plotto rx, ry, rz
scatter3(-21800, 37900, 0, 100, 'magenta', 'filled'); %plotto point P1
scatter3(27300, 27700, 0, 100, 'yellow', 'filled'); %plotto point P2
legend('','Orbita di trasferimento','P1','P2');
%Nel legend il primo '' è vuoto per evitare di considerare "Earth_sphere"
grid on;
hold off;
xlabel('x [km]');
ylabel('y [km]');
zlabel('z [km]');
title(['Solution for Lambert problem']);

%% Esercizio 2
%Pulizia
clc;
clear;
close all;

%Dati richiamati spesso
mu_E=astroConstants(13);
J2=astroConstants(9);
Re=astroConstants(23);

%Initial state (keplerian elements)
a1=12500;
e1=0;
i1=0;
Omega1=0;
w1=0;
f1=degtorad(120);

%Final state (keplerian elements)
a2=9500;
e2=0.3;
i2=0;
Omega2=0;
w2=0;
f2=degtorad(250);

%Trasformo in coordinate cartesiane col file 'kep2car' (gli angoli devono
%essere tutti in radianti --> ottengo r in km e v in km/s)
%f sarebbe "theta" (true anomaly)
%Initial state (cartesian elements)
%Le seguenti sono relative alle orbite 1 e 2
[r1,v1i]=kep2car(a1,e1,i1,Omega1,w1,f1, mu_E); %r e v (orbita 1) del punto P1
[r2,v2f]=kep2car(a2,e2,i2,Omega2,w2,f2, mu_E); %r e v (orbita 2) del punto P2
%v1i = velocità riferita al punto P1 sull'orbita iniziale (i --> orbita 1)
%v2f = velocità riferita al punto P2 sull'orbita finale (f --> orbita 2)

%Utilizzo Lambert solver per l'orbita di trasferimento invece (transfer
%arc) --> grazie ad esso posso dire che P1 e P2 appartengono ad una stessa
%orbita --> lo spacecraft passa in P1 all'istante t1 ed in P2 all'istante
%t2 --> passo r1 ed r2 al Lambert solver (vettori posizione di P1 e P2
%rispetto al fuoco dell'orbita di trasferimento) ed ottengo v1 e v2
%(diversi da quelli definiti prima (1 rispetto all'orbita 1 e l'altro
%rispetto all'orbita 2) --> ora invece ho v1 e v2 che fanno riferimento alla
%stessa orbita, ma a due istanti differenti)
tof = 3300; %Time of flight (tempo di trasferimento --> deltat=t2-t1)
[a,p,e,ERROR,v1t,v2t,TPAR,THETA] = lambertMR( r1, r2, tof, mu_E, 0, 0, 0, 0);
%da tale funzione mi servono gli OUTPUT v1t e v2t (rispettivamente le
%velocità relative al punto P1 ed al punto P2 SULL'orbita di trasferimento)

%Ora è necessario calcolare le differenze di velocità: 
% Injection manoeuvre (from initial orbit to transfer arc): deltav1=v1t-v1i
deltav1=v1t'-v1i; %devo trasporre v1t per ottenere un vettore colonna
%Arrival manoeuvre (from transfer arc to final orbit): deltav2=v2f-v2t
deltav2=v2f-v2t'; %devo trasporre v2f per ottenere un vettore colonna

%The total cost of the mission is:
deltavtot=norm(deltav1)+norm(deltav2); %si trova [km/s]

%Ora propagare la transfer orbit da t1 a t2 --> dunque parto dalle
%condizioni iniziali in t1 (r1,v1t) e propago per un deltat pari a t2-t1,
%ovvero pari a tof (=3300s) [utilizzo l'ODE (soluzione perturbata)]
options=odeset('RelTol',1e-13,'AbsTol',1e-14);
y0=[r1;v1t']; %deve essere colonna
tspan=linspace(0,tof,100); %come se t1 fosse pari a 0 e t2 pari a tof
[t yt]=ode45(@(t,y)perturbed_ode_2bp( t, y, mu_E, J2, Re),tspan,y0,options);

%Plot the initial and final orbits and the transfer arc
%Propago prima le orbite 1 e 2 per un periodo T (pari a quello delle
%orbite) e poi plotto
%Orbita iniziale (1)
T1=2*pi*sqrt(a1^3/mu_E); %periodo dell'orbita 1 (iniziale)
tspan1=linspace(0,T1,1000);
y01=[r1;v1i];
[t yi]=ode45(@(t,y)perturbed_ode_2bp( t, y, mu_E, J2, Re),tspan1,y01,options);
%Orbita finale (2)
T2=2*pi*sqrt(a2^3/mu_E); %periodo dell'orbita 2 (finale)
tspan2=linspace(0,T2,1000);
y02=[r2;v2f];
[t yf]=ode45(@(t,y)perturbed_ode_2bp( t, y, mu_E, J2, Re),tspan2,y02,options);
%Ricapitolando: yt è per l'orbita di trasferimento, yi è per l'orbita
%iniziale (1), yf è per l'orbita finale (2)

%Plotto anche la restante parte dell'orbita di trasferimento con linea
%tratteggiata
%Propago tutta l'orbita di trasferimento con condizioni iniziali pari a 
%(r1,v1t) o (r2,v2t)
[at,~,~,~,~,~] = car2kep(r1,v1t',mu_E); %ricavo 'at' da car2kep 
                                      %("at" semiasse maggiore dell'orbita)
Tt=2*pi*sqrt(at^3/mu_E);
tspant=linspace(0,Tt,1000);
[t yttot]=ode45(@(t,y)perturbed_ode_2bp( t, y, mu_E, J2, Re),tspant,y0,options);

%Ora plotto il tutto
figure(2);
earth_sphere();
hold on;
plot3(yi(:,1),yi(:,2),yi(:,3),'-m','LineWidth',2);
plot3(yf(:,1),yf(:,2),yf(:,3),'-g','LineWidth',2);
plot3(yt(:,1),yt(:,2),yt(:,3),'-b','LineWidth',2);
plot3(yttot(211:end,1),yttot(211:end,2),yttot(211:end,3),'--b','LineWidth',2);
plot3(r1(1), r1(2), r1(3),'oc','MarkerFaceColor','c','MarkerSize',10); %plotto point P1
plot3(r2(1), r2(2), r2(3),'oy','MarkerFaceColor','y','MarkerSize',10); %plotto point P2
%'MarkerFaceColor' è per colorare il cerchio internamente
grid on;
legend('','Orbita iniziale','Orbita finale','Orbita di trasferimento','unused transfer arc','P1','P2');
%In plot3 mettere "LineWidth", altrimenti legend non legge quella determinata curva 
%Nel legend il primo '' è vuoto per evitare di considerare "Earth_sphere"
xlabel('x [km]');
ylabel('y [km]');
zlabel('z [km]');
title(['The three orbits']);
%% Esercizio 3

%Definizione ulteriori costanti
mu_S=astroConstants(4); %costante gravitazionale del Sole
n=10; %è semplicemente un indice che utilizzo in linspace
%Definisco le finestre di partenza e di arrivo
jDPinf=date2mjd2000([2003 04 01 12 00 00]); %julian date di partenza (limite inferiore)
jDPsup=date2mjd2000([2003 08 01 12 00 00]); %julian date di partenza (limite superiore)
jDAinf=date2mjd2000([2003 09 01 12 00 00]); %julian date di arrivo (limite inferiore)
jDAsup=date2mjd2000([2004 03 01 12 00 00]); %julian date di arrivo (limite superiore)
finestraPartenza=linspace(jDPinf,jDPsup,n); %finestra di partenza 
finestraArrivo=linspace(jDAinf,jDAsup,n); %finestra di arrivo
for i=1:n
    for j=1:n
        [v_E,v_M,deltaVTOT(i,j)]=Transferdeltavtot(mu_S,finestraPartenza(i),finestraArrivo(j),3,4) %vedere la funzione "Transferdellavtot" per gli indici
        %queste v_E e v_M sono relative a t1finale e t2finale (perchè
        %vengono sovrascritte ad ogni iterazione)
    end
end