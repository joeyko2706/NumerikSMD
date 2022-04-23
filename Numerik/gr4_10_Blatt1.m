%%  Numerische Mathematik fuer Physik und Ingenieurwissenschaften SS2022 
%
%   Prof. Dr. J. Stoeckler
%   M.Sc. M. Weimann
%
%   Programmierblatt 1

close all;  
clear all; 

%% Definition der Funktionen
f = @funktion;


% Polynom vierten Grades
n = 4;
x4 = linspace(-5,5, n+1); %... Knoten mittels linspace
y4 = f(x4);%... Daten
p4 = polyfit(x4, y4, n);%... Interpolationspolynom mittels polyfit		

% Polynom zehnten Grades
n = 10;%... Grad
x10 = linspace(-5,5, n+1);%...
y10 = f(x10);%...
p10 = polyfit(x10, y10, n)%...		

%% Auswertung
xx = linspace(-5,5,101); % diskrete Definitionsmenge
yyf = f(xx);%... Auswertung von f an den Stellen xx
yy4 = polyval(p4, xx);%... Auswertung von p4 an den Stellen xx mittels polyval	
yy10 = polyval(p10, xx);%...
r4 = abs(f(xx) - yy4);%... Interpolationsfehler zu p4
r10 = abs(f(xx) - yy10); %...

%% Ausgabe
%... 
subplot(3,1,1); %Drei Plots in einen packen
title('Funktion f(x) mit 4 Knoten')
hold on
plot(xx, yyf, "blue"); %Plot von der Funktion
plot(x4, y4, "magenta"); %Knotenpunkte
plot(xx, yy4, "green"); %Interpolationspolynom
legend('Funktion', "Knotenpunkte", 'Interpolationspolynom');
hold on

subplot(3,1,2); %Drei Plots in einen packen
title('Funktion f(x) mit 10 Knoten')
hold on
plot(xx, yyf, "blue"); %Plot von der Funktion
plot(x10, y10, "magenta"); %Knotenpunkte
plot(xx, yy10, "green"); %Interpolationspolynom
legend('Funktion', "Knotenpunkte", 'Interpolationspolynom');
hold on

subplot(3,1,3); %Drei Plots in einen packen
title('Fehlerfunktion')
hold on
plot(xx, r4, "blue"); %Interpolationspolynom
plot(xx, r10, "magenta");
legend('Fehlerpolynom fuer 4 Knoten', "Fehlerpolynom fuer 10 Knoten");
hold on

function a = funktion(x)
a = 10./(1+x.^2);
end