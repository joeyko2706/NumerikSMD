%Programmierblatt 02
%%

%Felix Symma und Joel Koch

clear all
close all

%% Definition der Funktion f als functin handle
f = @(x) sin(5*pi*x);
a = 0;
b = 1;

%% Berechnung der absoluten Fehler für N = 1,...,4
F_abs = zeros(4,4); % Zeilen: N=1,..,4; Spalten: Schemata (Mittelpunkt, ...)

for N = 1:4
    
    I = zeros(1,4);
    
    for k = 1:N
        y_start = a + (k-1) * (b-a)/N;
        y_end = a + k * (b-a)/N;
        I = I + [m(f,y_start,y_end), t(f,y_start,y_end), s(f,y_start,y_end), g(f,y_start,y_end)];
    end
    
    F_abs(N,:) = abs(2/(5*pi) - I);
    
end

%% Plot des Balkendiagramms analog zum Bild auf dem Blatt
figure('name','Programmierblatt 02')
bar(categorical({'N=1','N=2','N=3','N=4'}),F_abs)
legend('Mittelpunkt','Trapez','Simpson','2-pkt. Gauß')
title('Absolute Fehler der summierten Quadrationsformeln zu \int_0^1 sin(5\pi x) dx')

%% nichtsummierte Quadraturformeln als Funktionen
function I = m(f,c,d)
I = (d-c) * f(0.5 * (c+d));
end

function I = t(f,c,d)
I = 0.5 * (d-c) * (f(c) + f(d));
end

function I = s(f,c,d)
I = (d-c)/6 * (f(c) + 4 * f(0.5 * (c+d)) + f(d));
end

function I = g(f,c,d)

x0 = 0.5 * (d+c) - sqrt(3) * (d-c)/6;
x1 = 0.5 * (d+c) + sqrt(3) * (d-c)/6;
w0 = 0.5 * (d-c);
w1 = 0.5 * (d-c);

I = w0 * f(x0) + w1 * f(x1);
end

%% Interpretation

%Die Funktion wird nicht ausreichend abgetastet, wenn das N klein ist. Das
%bedeutet, dass es zu wenig Knoten gibt, um die Funktion genau zu
%approximieren. Ist das N das doppelte der Frequenz (2*2,5=5), werden die
%Ergebnisse besser.
