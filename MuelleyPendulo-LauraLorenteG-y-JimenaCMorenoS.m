clear; 
clc;
close all;

%% 1. Parámetros del sistema.

m1 = 2; % Masa bloque
m2 = 1; % Masa del péndulo
l0 = 1; % Longitud natural del muelle
k = 10; % Constante elástica
l = 2; % Longitud del péndulo
g = 9.8; % Gravedad

%% 2. Condiciones iniciales

x0 = l0; % Posición inicial del bloque 
dotx0 = 0; % Velocidad inicial del bloque
theta0 = 90; % Ángulo inicial del péndulo
dottheta0 = 0; % Velocidad angular inicial del péndulo
t0 = 0; % Tiempo inicial
tf = 10; % Tiempo final

%% Salida

fprintf('Péndulo articulado sobre una masa con movimiento horizontal \n')

fprintf('--------------------------------------------------\n');
fprintf('Valores: \n')
fprintf('--------------------------------------------------\n');

fprintf('Masa del bloque en kg: %f\n', m1);
fprintf('Masa del péndulo en kg: %f\n', m2);
fprintf('Longitud natural del muelle en m: %f\n', l0);
fprintf('Constante elástica en N/m: %f\n', k);
fprintf('Longitud del péndulo en m: %f\n', l);

fprintf('--------------------------------------------------\n');
fprintf('Condiciones iniciales: \n')
fprintf('--------------------------------------------------\n');

fprintf('Posición inicial del bloque en m: %f\n', x0);
fprintf('Velocidad inicial del bloque en m/s: %f\n', dotx0);
fprintf('Ángulo inicial del péndulo en grados: %f\n', theta0);
fprintf('Velocidad inicial angular del péndulo en rad/s: %f\n', dottheta0);

fprintf('--------------------------------------------------\n');
fprintf('Ecuaciones de movimiento: \n')
fprintf('--------------------------------------------------\n');

fprintf('Ecuación de movimiento para theta: \n');
fprintf('m2*l^2* ddot(theta) + m2* ddot(x)*l*cos(theta) + m2*g*l*sen(theta)=0 \n');
fprintf('Ecuación de movimiento para x: \n');
fprintf('ddot(x) *(m1 + m2) + m2*l* ddot(theta)*cos(theta) - m2*l* dot(theta)^2*sen(theta) + kx=0 \n');

%% Pasamos a radianes el ángulo

theta0 = deg2rad(theta0);

%% Definimos el intervalo de tiempo

tinterval = [t0 tf];

%% 3. Definimos el vector columna de condiciones iniciales

X0 = [x0 ; dotx0 ; theta0 ; dottheta0];

%% 4. Resolvemos el sistema 

options = odeset('RelTol', 1e-7, 'AbsTol', 1e-9); % Vector opciones de tolerancia

derivadas = @(t, X) resolucion(t, X,m1,m2,l,g,k,l0); % Vector que contiene las primeras y segundass derivadas

tic
[t, xxpttp] = ode45(derivadas,tinterval,X0,options); % Obtenemos los valores de x, derivada primera de x, theta, derivada primera de theta
tiempo = toc; % Para calcular cuanto tiempo tarda en realizar las operaciones pertinentes

%% 5. Extraemos los resultados del vector xxpttp

xt = xxpttp(:,1);
xtp = xxpttp(:,2);
thetat = xxpttp(:,3);
thetatp = xxpttp(:,4);

xt_raw = xxpttp(:,1);
thetat_raw = xxpttp(:,3);
t_uniform = linspace(t0, tf, 600); % 600 fotogramas uniformes
xt = interp1(t, xt_raw, t_uniform, 'pchip');
thetat = interp1(t, thetat_raw, t_uniform, 'pchip');

function dxdt = resolucion(~, xin, M1, M2, L, G, K, L0)
    dxdt = zeros(4,1);
    x = xin(1);       
    xp = xin(2);      
    Theta = xin(3);   
    Thetap = xin(4);  
    
    sentheta = sin(Theta);
    costheta = cos(Theta);
    
    comdem = (M1 + M2 * (sentheta)^2);
    %% Fuerza de rebote contra la pared
    ancho_bloque = 1.0;                  
    borde_izq = x - (ancho_bloque / 2);  % Posición real del borde izquierdo del bloque
    f_pared = 0;

    if borde_izq < 0  
        k_pared = 100000;                % Constante alta
        f_pared = -k_pared * borde_izq;  % Fuerza de empuje
    end
    %% Ecuaciones de movimiento
    xpp = ((M2 * L * (Thetap^2) * sentheta) - (K * (x - L0)) + (M2 * G * costheta * sentheta) + f_pared) / comdem;
    thetapp = ((-G * sentheta * (M1 + M2)) - (M2 * L * (Thetap^2) * costheta * sentheta) + (K * (x - L0) * costheta)) / (L * comdem);
   
    dxdt(1) = xp;
    dxdt(2) = xpp;
    dxdt(3) = Thetap;
    dxdt(4) = thetapp;
end

%% 6. Elemento gráfico
figure ('Name','Péndulo articulado sobre una masa con movimiento horizontal')
n_puntos = length(xt); 

for i = 1:8:n_puntos 
    %% Posición central del bloque en este instante:
    x_bloque = xt(i);

    %% Dimensiones del bloque:
    ancho = 1;
    alto = 0.5;
    centro = alto / 2; % Altura del punto de enganche del péndulo (y = 0.25)

    %% Posición de la masa del péndulo (m2):
    % Se cuelga exactamente desde el centro del bloque (x_bloque)
    x2 = x_bloque + l * sin(thetat(i));
    y2 = centro - l * cos(thetat(i)); 

    %% Limpiamos la imagen anterior
    cla;
    hold on 

    %% Dibujar la pared fija en x = 0
    line([0, 0], [-2, 2], 'Color', [0.2 0.2 0.2], 'LineWidth', 4);
    for y_p = -2:0.1:2
        plot([0, -0.1], [y_p, y_p - 0.05], 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5);
    end

    % Línea horizontal del suelo
    longitud_suelo = max(xt) + l + 1.5;
    line([0, longitud_suelo], [0, 0], 'Color', [0.2 0.2 0.2], 'LineWidth', 2);
    for x_s = 0.2:0.3:longitud_suelo
        plot([x_s, x_s - 0.1], [0, -0.1], 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
    end

    %% Dibujar el muelle (desde la pared x=0 hasta el lado izquierdo del bloque)
    n_espiras = 10;
    xfmuelle = x_bloque - ancho / 2;
    xmuelle = linspace(0, xfmuelle, 40); 
    A = 0.30; 
    ymuelle = (alto / 2) + A * sin(linspace(0, n_espiras * 2 * pi, 40));
    plot(xmuelle, ymuelle, 'Color', [0 0.4 0.8], 'LineWidth', 1.2);

    %% Dibujar la caja (bloque m1)
    rectangle('Position', [x_bloque - ancho/2, 0, ancho, alto], ...
        'FaceColor', [0.8 0.8 0.8], ...
        'Curvature', 0.2);

    %% Dibujar la varilla rígida del péndulo (longitud 'l' constante)
    plot([x_bloque, x2], [centro, y2], 'k-', 'LineWidth', 2);

    %% Dibujar la masa del péndulo (m2)
    radiom2 = 10;
    plot(x2, y2, 'o', 'MarkerSize', radiom2, 'MarkerFaceColor', [1 0.6 0.6], 'MarkerEdgeColor', 'k');

    %% Configuración de ejes y proporción
    xlim([-1, max(xt) + l + 1.5]);
    ylim([-l - 1.5, alto + 1.5]);
    axis equal; 

    grid on;
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';

    drawnow;
    pause(0.01); 
end