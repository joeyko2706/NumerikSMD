%%  Numerische Mathematik fuer Physik und Ingenieurwissenschaften SS2021
%
%   Prof. Dr. J. Stoeckler
%   M.Sc. M. Weimann
%   Dipl.-Math. M. Bangert
%
%   Musterloesung 5
%   Abgabe bis zum 17.06.2021
%
%   Student*in 1: Vorname, Nachname, Matrikelnummer
%   Student*in 2: Vorname, Nachname, Matrikelnummer
%
%   Programmversion: z.B. Matlab R2021a oder Octave 6.2.0
%
close all
clear all
clc
%
%% 5a
% Definition von f
f = @(x) 1./(10 + x.*x);       % f als function handle

% Erstelle Figure
figure('name','Programmierblatt 5')
hold on

% Gauss-Approximation für verschiedne n
for n = 1:4                     % Dimension des Polynomraumes
    
    M = hilb(n+1);                % Hilbert-matrix
    b = zeros(n+1,1);             % Initialisierung von b als Null-Spaltenvektor
    s = ones(n+1,1);              % Stoervektor
    eps = 1e-5;                 % Faktor zum Stoervektor
    
    ret = zeros(1,2);           % Initialisierung der Ausgabematrix
    
    
    % ungestoert
    for i = 0:n
        b(i+1) = integral(@(x) (f(x) .* x.^i), 0,1);
    end
    
    c = invhilb(n+1) * b;
    nb = norm(b);
    nc = norm(c);
    
    g = @(x) 0;
    for i = 0:n
        g = @(x) g(x) + c(i+1) * x.^i;
    end
    
    ret(1,2) = sqrt(integral(@(x) abs(f(x) - g(x)).^2,0,1));
    
    % gestoert
    for l = 1:100
        
        k_l = l * eps;
        b_l = b + k_l * s;
        c_l = invhilb(n+1)*b_l;
        
        % Setze Ergebnis-Vektor g zusammen
        g_l = @(x) 0;
        for i = 0:n
            g_l = @(x) g_l(x) + c_l(i+1) * x.^i;
        end
        
        % Setze Ausgabematrix
        ret(end+1,1) = norm(b-b_l)/nb;
        ret(end,2) = sqrt(integral(@(x) abs(f(x) - g_l(x)).^2,0,1));
    end
    
    % plot für jedes n
     plot(ret(:,1),ret(:,2))
end

xlim([0,0.01])
legend('n=1','n=2','n=3','n=4')
title('Approximationsfehler in Abhaengigkeit der Stoerung von b')
xlabel('||b-b\_l||_2 / ||b||_2')
ylabel('(\int |f-g\_l|^2)^{1/2}')


%% Aufgabe 5b
%{
Die Gauß-Approximation verbessert sich fuer ungestoerte b mit steigendem
Polynomgrad (y-Achsen-Abschnitte). Allerdings verschlechtert die Stoerung die 
Approximationen fuer hohe Polynomgrade "schneller" (Steigungen der
Graphen), sodass sich bereits bei sehr kleinen Stoerungen die Guete der Approximationen
umdreht. Dann sind die Gauß-Approximationen mit kleineren Polynomgraden
besser, weil die Stoerungen wegen der kleineren Konditionszahl einen geringeren Einfluss haben.
%}

%% Aufgabe 5c
%{
Die Konditionszahl der Hilbert-Matrix waechst schnell mit der Raumdimension
n. Kleine Schwankungen in b koennen also fuer groesser werdende n grosse 
Variationen in c ausloesen. Die Variationen in c uebertragen sich auf die
Gauss-Approximation und somit auch auf den Approximationsfehler.
%}
for n = 1:4                    
    disp(['cond(H) fuer n = ' num2str(n) ':  ' num2str(cond(hilb(n+1)))]);
end