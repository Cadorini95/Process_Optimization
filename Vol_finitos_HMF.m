%% Simulação dinâmica - Volumes Finitos

clc;
clear all;
close all; 

% Variáveis operacionais - parâmetros da malha discretizada
% W = 100 kgs de massa de catalisador

M = 500;            % vazão mássica (kg/h) -> cond.otimizada (WHSV = 5 1/h)
ps = 1000;          % densidade fluido - água (kg/m³) 
d = 0.3;            % diâmetro do reator (m)
A = (d^2)*pi/4;     % área da seção transversal (m²)
T = 367.6;          % Temperatura (K) -> cond.otimizada
L = 1.5;            % Comprimento do reator (m)
pc = 2200;          % Densidade do catalisador (kg/m³)
Chmf0 = 0.1;        % Concentração inicial de HMF (mol/L)
Coxi = 0;           % Concentração inicial de produtos de oxidação (mol/L)
por = 0.59;         % porosidade do leito
N = 300;            % Pontos da malha

p = [M,ps,A,T,L,pc,Chmf0,Coxi,por,N];

X0 = zeros(N,1);
Y0 = zeros(N,1);

% Tempo de simulação

tmin = 0;
tmax = 2;

% Resolução da EDO discretizada

[t,X] = ode45(@(t,X) modelo(t,X,p),[tmin tmax], [X0 Y0]);
zspan = linspace(0,L,N);   % discretização do reator 

% Perfil axial da concentração (HMF e produtos da oxidação) - Visualização gráfica

figure(1)

plot(zspan,X(end,1:N),'-.r', 'Linewidth',2,'DisplayName', 'HMF' ); hold on; grid on;
plot(zspan,X(end,N+1:2*N),'-.k', 'Linewidth',2, 'DisplayName', 'Produtos da oxidação'); 
xlabel('Comprimento (m)','Fontsize',12,'Fontweight','bold')
ylabel('C_{EE} (mol L^{-1})', 'Fontsize',12, 'Fontweight','bold')
title('\fontsize{16}Perfil axial de concentração')
legend('Location', 'northeast');
ylim([0,Chmf0 + 0.01]);

%  Visualização gráfica do perfil dinâmico  -  HMF e produtos da oxidação

figure(2)

plot(t, X(:,N),'-.r', 'Linewidth',2,'DisplayName', 'HMF'); hold on; grid on;
plot(t,X(:,2*N), '-.k', 'Linewidth',2, 'DisplayName', 'Produtos da oxidação');
xlabel('Tempo (h)','Fontsize',12,'Fontweight','bold')
ylabel('C_{EE} (mol L^{-1})', 'Fontsize',12, 'Fontweight','bold')
title('\fontsize{16}Evolução dinâmica no reator')
legend('Location', 'northeast');
ylim([0, 0.07]);

% Criação de matriz vazia 

X1 = X(end,1:N); 
X2 = X(end,N+1:2*N);
M = zeros(length(N),3);

% Preenchimento da matriz

for i = 1:length(X1)
    
    M(i,1) = abs(Chmf0-X1(1,i))*100/Chmf0; % conversão de HMF 
    M(i,2) = X2(1,i)*100/Chmf0;            % conversão à produtos oxidados
    M(i,3) = (X1(1,i)+ X2(1,i))*100/Chmf0; % Balanço de carbono 
    
end

% Perfil axial da concentração (Xhmf, Xoxi e BC) - Visualização gráfica

figure(3)

plot(zspan,M(:,1),'-.r', 'Linewidth',2,'DisplayName', 'X_{HMF}' ); hold on;
plot(zspan,M(:,2),'-.k', 'Linewidth',2, 'DisplayName', 'X_{OXI}'); hold on;
plot(zspan,M(:,3),'-.g', 'Linewidth',2, 'DisplayName', 'Balanço de carbono');
xlabel('Comprimento (m)','Fontsize',12,'Fontweight','bold')
ylabel('Conversão (%)','Fontsize',12, 'Fontweight','bold')
title('\fontsize{16}Perfil axial da conversão')
legend('Location', 'northeast');
ylim([0,120]);
grid on;

% Criação da matriz vazia 

W = zeros(length(t),3);
Chmf = X(:,N);
Coxi = X(:,2*N);

% Preenchimento da matriz 

for i = 1:length(t)
 
    W(i,1) = abs(Chmf0 - Chmf(i,1))/Chmf0*100;
    W(i,2) = Coxi(i,1)*100/Chmf0;
    W(i,3) = (Chmf(i,1) + Coxi(i,1))*100/Chmf0;
    
end

%  Visualização gráfica do perfil dinâmico  -  Xhmf, Xoxi e BC

figure(4)

plot(t,W(:,1),'-.r', 'Linewidth',2,'DisplayName','X_{HMF}'); hold on;
plot(t,W(:,2),'-.k', 'Linewidth',2, 'DisplayName', 'X_{oxi}'); hold on;
plot(t,W(:,3),'-.g', 'Linewidth',2, 'DisplayName', 'Balanço de carbono');
xlabel('Tempo (h)','Fontsize',12,'Fontweight','bold')
ylabel('Conversão (%)', 'Fontsize',12, 'Fontweight','bold')
title('\fontsize{16}Evolução dinâmica no reator')
legend('Location', 'northeast');
ylim([0,120]);
grid on;

%% Função discretizada por volumes finitos - PFR unidimensional sem dispersão

function dxdt = modelo(t,X,p)

% Parâmetros operacionais

M = p(1);    % vazão mássica (kg/h)
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

Chmf = X(1:N);        % Concentração de HMF
Coxi = X(N+1:2*N);    % Concentração dos produtos de oxidação

% Cálculo das constantes cinéticas

% Parâmetros cinéticos HMF 

E1 = 49.7*10^3;          % Energia de ativação HMF -> produtos (J/mol)
k01 = 4.73*10^4;         % fator pré-exponencial (Kg/m³h)

% Parâmetros cinéticos produtos da oxidação

E2 = 51.4*10^3;          % Energia de ativação HMF -> produtos de oxidação (J/mol)
k02 = 56.8*10^3;         % fator pré-exponencial (Kg/m³h)

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