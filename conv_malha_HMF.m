%% Converg�ncia de malha para reator de leito fixo

clc; 
clear all;
close all;

% Vari�veis operacionais - par�metros da malha discretizada
% W = 100 kgs de massa de catalisador

M = 500;            % vaz�o m�ssica (kg/h) -> cond.otimizada (WHSV = 5 1/h)
ps = 1000;          % densidade fluido - �gua (kg/m�) 
d = 0.3;            % di�metro do reator (m)
A = (d^2)*pi/4;     % �rea da se��o transversal (m�)
T = 367.6;          % Temperatura (K)-> cond.otimizada
L = 1.5;            % Comprimento do reator (m)
pc = 2200;          % Densidade do catalisador (kg/m�)
Chmf0 = 0.1;        % Concentra��o inicial de HMF (mol/L)
Coxi = 0;           % Concentra��o inicial de produtos de oxida��o (mol/L)
por = 0.59;         % porosidade do leito

% Tempo de simula��o

tmin = 0;
tmax = 2;

% Par�metros iniciais para converg�ncia de malha

N1 = 10;      % n�mero de pontos
tol = 10^-1;  % toler�ncia do erro
erro = 1;

while erro >= tol | N1 <= 1000
    
    p1 = [M,ps,A,T,L,pc,Chmf0,Coxi,por,N1];
    X1 = zeros(N1,1);
    Y1 = zeros(N1,1);
    % discretiza��o do reator
    zspan = linspace(0,L,N1);
    % Resolu��o da EDO discretizada
    [t,X1] = ode45(@(t,X) modelo(t,X,p1),[tmin tmax],[X1 Y1]);
    xspan = X1(end,1:N1);
    

    % Novo n�mero de pontos da malha 
    N2 = N1*(2);
    
    % discretiza��o do reator - novo tamanho da malha
    znovo = linspace(0,L,N2);
    
    % Interpola��o dos pontos para compara��o da malha
    xant = interp1(zspan,xspan,znovo);
    
    
    p2 = [M,ps,A,T,L,pc,Chmf0,Coxi,por,N2];
    X2 = zeros(N2,1);
    Y2 = zeros(N2,1);
    
    % Resolu��o da EDO discretizada
    [t,X2] = ode45(@(t,X) modelo(t,X,p2),[tmin tmax], [X2 Y2]);
    xnovo = X2(end,1:N2);
    
    % M�trica para avaliar a converg�ncia da malha
    erro = sum(abs(xnovo - xant));
    
    % Perfil axial da concentra��o - Visualiza��o gr�fica
    
    figure(1)
    plot(znovo,xant,'-.','Linewidth',2,'DisplayName', strcat(' N = '," ", num2str(N1))); hold on; grid on
    xlabel('Espa�o discretizado (m)', 'Fontsize',12,'Fontweight','bold');
    ylabel('C_{EE} (mol L^{-1})', 'Fontsize',12, 'Fontweight','bold');
    title('\fontsize{16}Converg�ncia da malha');
    legend('Location', 'northeast');
    ylim([0.03, Chmf0 + Chmf0/10])
    
    % Recome�o do loop computacional
    N1 = N2;  
    
end

%% Fun��o discretizada por volumes finitos - PFR unidimensional sem dispers�o

function dxdt = modelo(t,X,p)

% Par�metros operacionais

M = p(1);    % vaz�o m�ssica (kg/s)
ps = p(2);   % densidade da solu��o (kg/m�)
A = p(3);    % �rea da se��o transversal (m�)
T = p(4);    % Temperatura (K) 
L = p(5);    % Comprimento do reator (m)
pc = p(6);   % densidade do catalisador (kg/m�)
C10 = p(7);  % concentra��o inicial de HMF (mol/L)
C20 = p(8);  % concentra��o inicial dos produtos oxidados (mol/L)
por = p(9);  % porosidade do leito
N = p(10);   % N�mero de pontos para discretiza��o

% Defini��o do passo 
delta = L /(N-1);

Chmf = X(1:N);    % Concentra��o de HMF
Coxi = X(N+1:2*N);    % Concentra��o dos produtos de oxida��o

% C�lculo das constantes cin�ticas

% Par�metros cin�ticos HMF 

E1 = 49.7*10^3;          % Energia de ativa��o HMF -> produtos (J/mol)
k01 = 4.73*10^4;         % fator pr�-exponencial (Kg/m�h)


% Par�metros cin�ticos produtos da oxida��o

E2 = 51.4*10^3;          % Energia de ativa��o HMF -> produtos de oxida��o (J/mol)
k02 = 56.8*10^3;        % fator pr�-exponencial (Kg/m�h)

R = 8.31;                % Constante universal dos gases (% mol K / J)

k1 = k01*exp(-E1/ (R*T));
k2 = k02*exp(- E2/ (R*T));

% Discretiza��o da malha espacial em N pontos dentro do 
% z1 < z2 < z3...zi... < zn (volumes finitos)


dxdt = zeros(2*N,1);    % Armazenar os valores calculados para derivadas (C1 e C2)
                    
               
% Condi��o de contorno na posi��o inicial z = 0 (i=1)

dxdt(1,1) =  + ( M / (ps * A * delta) )*( C10 - Chmf(1) ) - k1*Chmf(1)* pc*(1-por);
dxdt(N+1,1)=  + ( M / (ps * A * delta) )*(C20 - Coxi(1)) + k2*Chmf(1)* pc*(1-por);
   
% Discretiza��o no interior do reator com volumes finitos finitas ( 2 < z < N)

for i = 2:N
    
    % C�lculo da derivada temporal do HMF
    
    dC1 = ( M / (ps * A * delta) )*(Chmf(i-1) - Chmf(i));    % termo advectivo 
    r1 = (k1 * Chmf(i) * pc)*(1-por);                        % termo cin�tico 
    dxdt(i,1) = + dC1 - r1;                                  % derivada temporal
    
    % C�lculo da derivada temporal dos produtos de oxida��o
    
    dC2 = ( M / (ps * A * delta) )*(Coxi(i-1) - Coxi(i));    % termo advectivo
    r2 = (k2 * Chmf(i)* pc)*(1-por);                         % termo cin�tico
    dxdt(N+i,1) = + dC2 + r2;                                % derivada temporal
    
end
    

end