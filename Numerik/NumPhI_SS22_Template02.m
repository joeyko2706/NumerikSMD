

clear all
close all

%% Definition der Funktion f als functin handle 
f = @(x) sin(5*pi*x);
a = 0;
b = 1;

%% Berechnung der absoluten Fehler für N = 1,...,4 
F_abs = zeros(4,4); % Zeilen: N=1,..,4; Spalten: Schemata (Mittelpunkt, ...)

for N = 1:4
    
    I = [..., ..., ..., ...];% Integrale abhaengig von N
    
    F_abs(N,:) = %...
    
end

%% Plot des Balkendiagramms 
%...

%% summierte Quadraturformeln als Funktionen 
function I = sum_m(f,a,b,N)
y_k = %linspace...
z_k = %...
I = %...
end

function I = sum_t(f,a,b,N)
y_k = %...
I = %...
end

function I = sum_s(f,a,b,N)
y_k = %...
z_k = %...
I = %...
end

function I = sum_g(f,a,b,N)
y_k = %...
% Koordinatentransformation beachten!
I = 0;
for k = 1:N
    x0 = %...
    x1 = %...
    w0 = %...
    w1 = %...
    
    I = %...
end 
end

%% Interpretation 
%{
...
%}