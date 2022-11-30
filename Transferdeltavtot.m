function [v_E,v_M,deltaVtot]=Transferdeltavtot(mu,t1,t2,iplanet1,iplanet2) %t1 e t2 devono essere in julian days
                                                                     %iplanet1
                                                                     %e
                                                                     %iplanet2
                                                                     %sono gli indicatori del pianeta
                                                                     %della funzione "uplanet"
                                                                     %mu è
                                                                     %relativo
                                                                     %al
                                                                     %corpo
                                                                     %celeste
                                                                     %attorno
                                                                     %a cui
                                                                     %ruotano
                                                                     %sia
                                                                     %planet1
                                                                     %che
                                                                     %planet2
    [kep_E,~]=uplanet(t1,iplanet1); %al julian day k-esimo (per pianeta iniziale --> nello script Terra) 
                                            %ho i 6 elementi kep_E (kep_E (alla fine) è una matrice)
    [kep_M,~]=uplanet(t2,iplanet2);  %al julian day k-esimo (per pianeta finale --> nello script Marte) 
                                           %ho i 6 elementi kep_E (kep_M (alla fine) è una matrice)
    %Ovviamente devo avere la Terra per la finestra di partenza e Marte per
    %la finestra di arrivo (lo s/c va dalla Terra a Marte)
    [r_E,v_E]=kep2car(kep_E(1),kep_E(2),kep_E(3),kep_E(4),kep_E(5),kep_E(6),mu); %A quell'istante di tempo ho quella determinata posizione e velocità della Terra
                                   %mu_S perchè la Terra ruota attorno al
                                   %Sole (stavolta m1=Earth) [corrispondono
                                   %anche a r e v dello spacecraft al tempo
                                   %iniziale(di partenza) --> quando
                                   %avviene il launch lo s/c si trova sulla
                                   %Terra e segue con essa la sua orbita]
    [r_M,v_M]=kep2car(kep_M(1),kep_M(2),kep_M(3),kep_M(4),kep_M(5),kep_M(6),mu); %A quell'istante di tempo ho quella determinata posizione e velocità della Terra
                                   %mu_S perchè la Terra ruota attorno al
                                   %Sole (stavolta m1=Earth) [corrispondono
                                   %anche a r e v dello spacecraft al tempo
                                   %iniziale(di partenza) --> quando
                                   %avviene il launch lo s/c si trova sulla
                                   %Terra e segue con essa la sua orbita]
    %Il kep2car restituisce r e v come vettori colonna
    %Ora devo propagare lungo l'arco di trasferimento --> per un tempo pari
    %al time of flight (dunque ToF=t2-t1)
    ToF=(t2-t1)*24*3600; %deve essere in secondi
    [A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(r_E,r_M,ToF,mu,0,0,0,0);
    %Ricordo che Lambert mi restituisce VI e VF in riga
    VI=VI'; %Ora è colonna
    VF=VF'; %Ora è colonna
    deltaVi=VI-v_E;
    deltaVf=v_M-VF;
    deltaVtot=norm(deltaVi)+norm(deltaVf);
end