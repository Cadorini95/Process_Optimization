%% Convergência de malha para reator de leito fixo

clc; 
clear all;
close all;

% Variáveis operacionais - parâmetros da malha discretizada
% W = 100 kgs de massa de catalisador

M = 500;            % vazão mássica (kg/h) -> cond.otimizada (WHSV = 5 1/h)
ps = 1000;          % densidade fluido - água (kg/m³) 
d = 0.3;            % diâmetro do reator (m)
A = (d^2)*pi/4;     % área da seção transversal (m²)
T = 367.6;          % Temperatura (K)-> cond.otimizada
L = 1.5;            % Comprimento do reator (m)
pc = 2200;          % Densidade do catalisador (kg/m³)
Chmf0 = 0.1;        % Concentração inicial de HMF (mol/L)
Coxi = 0;           % Concentração inicial de produtos de oxidação (mol/L)
por = 0.59;         % porosidade do leito

% Tempo de simulação

tmin = 0;
tmax = 2;

% Parâmetros iniciais para convergência de malha

N1 = 10;      % número de pontos
tol = 10^-1;  % tolerância do erro
erro = 1;

while erro >= tol | N1 <= 1000
    
    p1 = [M,ps,A,T,L,pc,Chmf0,Coxi,por,N1];
    X1 = zeros(N1,1);
    Y1 = zeros(N1,1);
    % discretização do reator
    zspan = linspace(0,L,N1);
    % Resolução da EDO discretizada
    [t,X1] = ode45(@(t,X) modelo(t,X,p1),[tmin tmax],[X1 Y1]);
    xspan = X1(end,1:N1);
    

    % Novo número de pontos da malha 
    N2 = N1*(2);
    
    % discretização do reator - novo tamanho da malha
    znovo = linspace(0,L,N2);
    
    % Interpolação dos pontos para comparação da malha
    xant = interp1(zspan,xspan,znovo);
    
    
    p2 = [M,ps,A,T,L,pc,Chmf0,Coxi,por,N2];
    X2 = zeros(N2,1);
    Y2 = zeros(N2,1);
    
    % Resolução da EDO discretizada
    [t,X2] = ode45(@(t,X) modelo(t,X,p2),[tmin tmax], [X2 Y2]);
    xnovo = X2(end,1:N2);
    
    % Métrica para avaliar a convergência da malha
    erro = sum(abs(xnovo - xant));
    
    % Perfil axial da concentração - Visualização gráfica
    
    figure(1)
    plot(znovo,xant,'-.','Linewidth',2,'DisplayName', strcat(' N = '," ", num2str(N1))); hold on; grid on
    xlabel('Espaço discretizado (m)', 'Fontsize',12,'Fontweight','bold');
    ylabel('C_{EE} (mol L^{-1})', 'Fontsize',12, 'Fontweight','bold');
    title('\fontsize{16}Convergência da malha');
    legend('Location', 'northeast');
    ylim([0.03, Chmf0 + Chmf0/10])
    
    % Recomeço do loop computacional
    N1 = N2;  
    
end

%% Função discretizada por volumes finitos - PFR unidimensional sem dispersão

function dxdt = modelo(t,X,p)

% Parãmetros operacionais

M = p(1);    % vazão mássica (kg/s)
ps = p(2);   % densidade da solução (kg/m³)
A = p(3);    % área da seção transversal (m²)
T = p(4);    % Temperatura (K) 
L = p(5);    % Comprimento do reator (m)
pc = p(6);   % densidade do catalisador (kg/m³)
C10 = p(7);  % concentração inicial de HMF (mol/L)
C20 = p(8);  % concentração inicial dos produtos oxidados (mol/L)
por = p(9);  % porosidade do leito
N = p(10);   % Número de pontos para discretização

% Definição do passo 
delta = L /(N-1);

Chmf = X(1:N);    % Concentração de HMF
Coxi = X(N+1:2*N);    % Concentração dos produtos de oxidação

% Cálculo das constantes cinéticas

% Parâmetros cinéticos HMF 

E1 = 49.7*10^3;          % Energia de ativação HMF -> produtos (J/mol)
k01 = 4.73*10^4;         % fator pré-exponencial (Kg/m³h)


% Parâmetros cinéticos produtos da oxidação

E2 = 51.4*10^3;          % Energia de ativação HMF -> produtos de oxidação (J/mol)
k02 = 56.8*10^3;        % fator pré-exponencial (Kg/m³h)

R = 8.31;                % Constante universal dos gases (% mol K / J)

k1 = k01*exp(-E1/ (R*T));
k2 = k02*exp(- E2/ (R*T));

% Discretização da malha espacial em N pontos dentro do 
% z1 < z2 < z3...zi... < zn (volumes finitos)


dxdt = zeros(2*N,1);    % Armazenar os valores calculados para derivadas (C1 e C2)
                    
               
% Condição de contorno na posição inicial z = 0 (i=1)

dxdt(1,1) =  + ( M / (ps * A * delta) )*( C10 - Chmf(1) ) - k1*Chmf(1)* pc*(1-por);
dxdt(N+1,1)=  + ( M / (ps * A * delta) )*(C20 - Coxi(1)) + k2*Chmf(1)* pc*(1-por);
   
% Discretização no interior do reator com volumes finitos finitas ( 2 < z < N)

for i = 2:N
    
    % Cálculo da derivada temporal do HMF
    
    dC1 = ( M / (ps * A * delta) )*(Chmf(i-1) - Chmf(i));    % termo advectivo 
    r1 = (k1 * Chmf(i) * pc)*(1-por);                        % termo cinético 
    dxdt(i,1) = + dC1 - r1;                                  % derivada temporal
    
    % Cálculo da derivada temporal dos produtos de oxidação
    
    dC2 = ( M / (ps * A * delta) )*(Coxi(i-1) - Coxi(i));    % termo advectivo
    r2 = (k2 * Chmf(i)* pc)*(1-por);                         % termo cinético
    dxdt(N+i,1) = + dC2 + r2;                                % derivada temporal
    
end
    

end