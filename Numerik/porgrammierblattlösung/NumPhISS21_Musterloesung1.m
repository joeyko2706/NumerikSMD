%%  Numerische Mathematik fuer Physik und Ingenieurwissenschaften SS2021 
%
%   Prof. Dr. J. Stoeckler
%   M.Sc. M. Weimann
%   Dipl.-Math. M. Bangert
%
%   Programmierblatt 1
%   Abgabe bis zum 21.04.2021
%  
%   Student*in 1: Vorname, Nachname, Matrikelnummer
%   Student*in 2: Vorname, Nachname, Matrikelnummer
%
%   Programmversion: z.B. Matlab R2021a oder Octave 6.2.0
%% 

close all; % close plots 
clear all; % delets all data in the Workspace

%% Definition der Funktionen
f = @(x) 10./(1+x.^2);

% Polynom vierten Grades
n = 4;
x4 = linspace(-5,5,n+1);
y4 = f(x4);
p4 = polyfit(x4, y4, n);	

% Polynom zehnten Grades
n = 10;
x10 = linspace(-5,5,n+1);
y10 = f(x10);
p10 = polyfit(x10, y10, n);	

%% Auswertung
xx = linspace(-5,5,101);
yyf = f(xx); 
yy4 = polyval(p4,xx);		
yy10 = polyval(p10,xx);
r4 = abs(yyf - yy4);
r10 = abs(yyf - yy10);

%% Ausgabe
subplot(2,2,1);
hold on
plot(xx,yyf,'b');			% Plotten von f
plot(xx,yy4, 'g');			% Plotten von p4
plot(x4, y4, 'r*','Markersize',10);

subplot(2,2,2);
plot(xx, yyf,'b');			% Plotten von f
hold on;
plot(xx, yy10, 'm');        % Plotten von p10
plot(x10, y10, 'r*','Markersize',10);

subplot(2,2,3);
plot(xx ,r4, 'g');			% Plotten von Fehlerfunktion r4
hold on;
plot(xx, r10, 'm');		    % Plotten von Fehlerfunktion r10