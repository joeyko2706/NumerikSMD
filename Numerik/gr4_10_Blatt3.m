%%  Numerische Mathematik fuer Physik und Ingenieurwissenschaften SS2022
%
%   Prof. Dr. J. Stoeckler
%   M.Sc. M. Weimann
%
%   Programmierblatt 3
%

close all
clear all


eps = 1e-5;                 % Faktor zum Stoervektor

%% 3a
% Definition von f
f = @(x) 1./(10 + x.*x);       % f als function handle

% Erstelle Figure
figure('name','Programmierblatt 5')
hold on

% Gauss-Approximationen für verschiedne n
for n = 1:4                     % Dimension des Polynomraumes
    
    M = hilb(n+1);                % Hilbert-matrix
    b = zeros(n+1,1);             % Initialisierung von b als Null-Spaltenvektor
    s = ones(n+1,1);              % Stoervektor
    
    ret = zeros(1,2);           % Initialisierung der Ausgabematrix
    
    for i = 0:n
        b(i+1) = integral(@(x) (f(x) .* x.^i), 0,1);    %... Berechnung des ungestoerten b
    end
    
    c = invhilb(n+1) * b;   %... Koeffizientenvektor

    g = @(x) 0;%... Gauss-Approximation als function handle
    for i = 0:n
        g = @(x) g(x) + c(i+1) * x.^i;
    end


    ret(1,2) = sqrt(integral(@(x) abs(f(x) - g(x)).^2,0,1)); %... Abstand von f und g (ungestoert)
    
    % gestoerte Faelle
    for j = 1:100
        
        k_j = j * eps;
        b_j = b + k_j * s;          %... gestoertes b
        c_j = invhilb(n+1) * b_j;     %... Koeffizientenvektor für gestoertes b
        
        %... Gauss-Approximation fuer gestoertes b
        g_j = @(x) 0;
        for i = 0:n
            g_j = @(x) g_j(x) + c_j(i+1) * x.^i;
        end


        % Setze Ausgabematrix
        ret(end+1,1) = norm(b-b_j)/norm(b);                             %... x-Achse des Plots
        ret(end,2) = sqrt(integral(@(x) abs(f(x) - g_j(x)).^2,0,1));    %... y-Achse des Plots
    end
   
     plot(ret(:,1),ret(:,2)) % plotte für jedes n
end

% Beschriftung des Plots
%...
xlim([0, 0.01])
ylim([0, 3*10^-3])
title("Approximationsfehler in Abhaengigkeit der Stoerung von b")
legend("n=1", "n=2", "n=3", "n=4")
xlabel('||b-b\_l||_2 / ||b||_2')
ylabel('(\int |f-g\_l|^2)^{1/2}')

%% Aufgabe 3b
%{
Der Plot sagt aus, dass die angewandten Gauß-Approximationen fuer
ungestoerte Faelle sich verbessern, wenn der Polynomgrad steigt.
Allerdings verschlechtert die Stoerung die Approxiamtion fuer hoehere
Polynomgrade schneller. Deshalb sind die Approximationen fuer geringe
Polynomgrade "besser".
%}