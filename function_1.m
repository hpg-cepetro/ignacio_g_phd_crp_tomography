function [V_real,V_orig,V_real_interp,x_structure,z_structure,x_NIP_real,z_NIP_real,angulo_NIP_real,t_total_dado,x_dado,px_dado,x_matrix,z_matrix] = Constroi_camadas_sal_v2()

% #Copyright (C) 2021 Gustavo Barroso Dias Ignácio
% <ignaciogbd@ggaunicamp.com>
% <ignaciogbd@gmail.com>


%%%%%%%%% CRP Tomography V1 %%%%%%%%%%%
tic

%%%%%%%%%% Melhor Peso Inicial %%%%%%%%%%%%%%%%
peso_ini = 1/1000;
peso_ini_x = 1;
peso_ini_z = 1;

%%%%%%%%% Inicializa Variáveis %%%%%%%%
h=0.001;
%suporte_x = [0,0,0,0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5,3.75,4,4.25,4.5,4.75,5,5,5,5]; nx = length(suporte_x);
%suporte_z = [0,0,0,0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,4,4,4]; nz = length(suporte_z);
%suporte_x = [0,0,0,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5,5,5,5]; nx = length(suporte_x);
%suporte_z = [0,0,0,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4,4,4,4]; nz = length(suporte_z);
%suporte_x = [0,0,0,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5,5,5]; nx = length(suporte_x);
%suporte_z = [0,0,0,0,0.4,0.8,1.2,1.6,2,2.4,2.8,3.2,3.6,4,4,4,4]; nz = length(suporte_z);
suporte_x = [0,0,0,linspace(0,6,31),6,6,6]; nx = length(suporte_x);
suporte_z = [0,0,0,linspace(0,5,26),5,5,5]; nz = length(suporte_z);
espacamento = 0.01; tamanho_x = length(0:espacamento:suporte_x(nx-3)); tamanho_z = length(0:espacamento:suporte_z(nz-3));
V_real = zeros(tamanho_z,tamanho_x); Matriz_coef_real = zeros(nz-4,nx-4); %Inicializa componentes e par�metros das matrizes de velocidades

for j=1:tamanho_z
    for i=1:tamanho_x
        %if j*espacamento > 1.5 && j*espacamento < 3 && ((i*espacamento > 0.25 && i*espacamento < 1.75) || (i*espacamento > 2.25 && i*espacamento < 3.75) || (i*espacamento > 4.25 && i*espacamento < 5.75))
        if j*espacamento > 1.5 && j*espacamento < 3 && i*espacamento > 1.5 && i*espacamento < 4.5   
            V_real(j,i) = 3.5;
        else
            V_real(j,i) = 1.25 + 1.75*cos((-pi/2)+4*(3/2)*(j-1)/(tamanho_z-1)/2.3);
        end
        if ((mod(j-1,20) == 0 || j == tamanho_z) && (mod(i-1,20) == 0 || i == tamanho_x))
            Matriz_coef_real((j-1)/20+2,(i-1)/20+2) = V_real(j,i);
        end
    end
end
Matriz_coef_real(1,:) = Matriz_coef_real(2,:);
Matriz_coef_real(:,1) = Matriz_coef_real(:,2);
Matriz_coef_real(nz-4,:) = Matriz_coef_real(nz-5,:);
Matriz_coef_real(:,nx-4) = Matriz_coef_real(:,nx-5);

%Construção das B-Splines
Pre_Bx = struct('form',strings(2,1),'knots',zeros(5,1),'coefs',0,'number',0,'order',0,'dim',0);
Pre_Bz = struct('form',strings(2,1),'knots',zeros(5,1),'coefs',0,'number',0,'order',0,'dim',0);
Bx = repmat(Pre_Bx,nx-4,1); Bz = repmat(Pre_Bz,nz-4,1);
Vet_Bx = zeros(nx-4,tamanho_x); Vet_Bz = zeros(nz-4,tamanho_z);
Deriv_Beta_x = zeros(nx-4,tamanho_x); Deriv_Beta_z = zeros(nz-4,tamanho_z);
Deriv_Beta_x2 = zeros(nx-4,tamanho_x); Deriv_Beta_z2 = zeros(nz-4,tamanho_z);

cont=0;
for i=3:(nx-2)
    Bx(i-2) = spmak([suporte_x(i-2),suporte_x(i-1),suporte_x(i),suporte_x(i+1),suporte_x(i+2)],1);
    for j=0:espacamento:suporte_x(nx-3)
        cont = cont+1;
        Vet_Bx(i-2,cont) = fnval(Bx(i-2),j); %Matriz que armazena valores das bases B-splines ao longo do grid.
    end
    [Deriv_Beta_x(i-2,:)] = gradient(Vet_Bx(i-2,:),espacamento);
    [Deriv_Beta_x2(i-2,:)] = gradient(Deriv_Beta_x(i-2,:),espacamento);
    cont=0;
end
for i=3:(nz-2)
    Bz(i-2) = spmak([suporte_z(i-2),suporte_z(i-1),suporte_z(i),suporte_z(i+1),suporte_z(i+2)],1);
    for j=0:espacamento:suporte_z(nz-3)
        cont = cont + 1;
        Vet_Bz(i-2,cont) = fnval(Bz(i-2),j); %Matriz que armazena valores das bases B-splines ao longo do grid.
    end
    [Deriv_Beta_z(i-2,:)] = gradient(Vet_Bz(i-2,:),espacamento); Deriv_Beta_z(i-2,:) = -Deriv_Beta_z(i-2,:);
    [Deriv_Beta_z2(i-2,:)] = gradient(Deriv_Beta_z(i-2,:),espacamento); Deriv_Beta_z2(i-2,:) = -Deriv_Beta_z2(i-2,:);
    cont = 0;
end

%Interpolação do Modelo de Velocidade Via B-Splines
V_real_interp(:,:) = ((Matriz_coef_real(:,:)*Vet_Bx)'*Vet_Bz)'; %Velocidade Real
V_smoth = zeros(tamanho_z,tamanho_x); ini = 40;
for i=1:tamanho_x
   for j=1:tamanho_z
         V_smoth(j,i) = mean(mean(V_real_interp(max(j-ini,1):min(j+ini,501),max(i-ini,1):min(i+ini,601))));
         %V_smoth(j,i) = mean(mean(V_real(max(j-ini,1):min(j+ini,501),max(i-ini,1):min(i+ini,601))));
   end
end
V_orig = V_real;
V_real = V_smoth;
%V_real = V_real_interp;
%Calcular Derivadas
[Deriv_x_real,Deriv_z_real] = gradient(V_real(:,:),espacamento); Deriv_z_real = -Deriv_z_real;

%Constroir interfaces
pos = 1;
for i=1:600
       z_structure(pos) = 0.4; x_structure(pos) = i*espacamento; pos=pos+1;
       z_structure(pos) = 0.8; x_structure(pos) = i*espacamento; pos=pos+1;
       z_structure(pos) = 1.2 + 0.1*cos(i*espacamento*2*pi/6); x_structure(pos) = i*espacamento; pos=pos+1;
       z_structure(pos) = 1.6 + 0.15*cos(i*espacamento*2*pi/6); x_structure(pos) = i*espacamento; pos=pos+1;
       z_structure(pos) = 2.0 + 0.2*cos(i*espacamento*2*pi/6); x_structure(pos) = i*espacamento; pos=pos+1;
       z_structure(pos) = 2.5 - 0.2*cos(i*espacamento*2*pi/6); x_structure(pos) = i*espacamento; pos=pos+1;
       z_structure(pos) = 2.9 - 0.15*cos(i*espacamento*2*pi/6); x_structure(pos) = i*espacamento; pos=pos+1;
       z_structure(pos) = 3.3 - 0.1*cos(i*espacamento*2*pi/6); x_structure(pos) = i*espacamento; pos=pos+1;
       z_structure(pos) = 3.7; x_structure(pos) = i*espacamento; pos=pos+1;
       z_structure(pos) = 4.1; x_structure(pos) = i*espacamento; pos=pos+1;
       z_structure(pos) = 4.5; x_structure(pos) = i*espacamento; pos=pos+1;
       z_structure(pos) = 4.9; x_structure(pos) = i*espacamento; pos=pos+1;
end
raios_picks=20;
%Atribui pontos difratores
pos = 1; num_points = 11; delta=5.7/(num_points-1); 
for i=1:num_points
    for n = 1:raios_picks
       x_NIP_real(pos) = (i-1)*delta+0.15; z_NIP_real(pos) = 0.4; 
       if x_NIP_real(pos) > 0.6 && x_NIP_real(pos) < 5.4
           double_aperture = pi/3; 
       else
           double_aperture = pi/6;
       end
       angulo_NIP_real(1,pos) = ((n-(raios_picks+1)/2)/((raios_picks-1)/2))*double_aperture;
       pos = pos+1;
    end
end
num_points = 12; delta=5.7/(num_points-1); 
for i=1:num_points
    for n = 1:raios_picks
       x_NIP_real(pos) = (i-1)*delta+0.15; z_NIP_real(pos) = 0.8; 
       if x_NIP_real(pos) > 1 && x_NIP_real(pos) < 5
           double_aperture = pi/3; 
       elseif x_NIP_real(pos) > 0.5 && x_NIP_real(pos) < 5.5
           double_aperture = pi/6;
       else
           double_aperture = pi/12;
       end
       angulo_NIP_real(1,pos) = ((n-(raios_picks+1)/2)/((raios_picks-1)/2))*double_aperture;
       pos=pos+1;
    end
end
num_points = 13; delta=5.6/(num_points-1); 
for i=1:num_points
    for n = 1:raios_picks
       x_NIP_real(pos) = (i-1)*delta+0.2; z_NIP_real(pos) = 1.2 + 0.1*cos(x_NIP_real(pos)*2*pi/6); 
       if x_NIP_real(pos) > 1.8 && x_NIP_real(pos) < 4.2
           double_aperture = pi/3; 
       elseif x_NIP_real(pos) > 0.5 && x_NIP_real(pos) < 5.5
           double_aperture = pi/6;
       else
           double_aperture = pi/15;
       end
       angulo_NIP_real(1,pos) = ((n-(raios_picks+1)/2)/((raios_picks-1)/2))*double_aperture + atan(-0.1*(2*pi/6)*sin(x_NIP_real(pos)*2*pi/6));
       pos=pos+1;
    end
end
num_points = 14; delta=5.6/(num_points-1); 
for i=1:num_points
    for n = 1:raios_picks
        x_NIP_real(pos) = (i-1)*delta+0.2; z_NIP_real(pos) = 1.6 + 0.15*cos(x_NIP_real(pos)*2*pi/6); 
       if x_NIP_real(pos) > 2 && x_NIP_real(pos) < 4
           double_aperture = pi/3; 
       elseif x_NIP_real(pos) > 1 && x_NIP_real(pos) < 5
           double_aperture = pi/6;
       else
           double_aperture = pi/18;
       end
           angulo_NIP_real(1,pos) = ((n-(raios_picks+1)/2)/((raios_picks-1)/2))*double_aperture + atan(-0.15*(2*pi/6)*sin(x_NIP_real(pos)*2*pi/6));
           pos=pos+1;
        end
end
num_points = 15; delta=5.5/(num_points-1); 
for i=1:num_points
    for n = 1:raios_picks
        x_NIP_real(pos) = (i-1)*delta+0.25; z_NIP_real(pos) = 2.0 + 0.2*cos(x_NIP_real(pos)*2*pi/6); 
       if x_NIP_real(pos) > 2.5 && x_NIP_real(pos) < 3.5
           double_aperture = pi/3; 
       elseif x_NIP_real(pos) > 1.5 && x_NIP_real(pos) < 4.5
           double_aperture = pi/6;
       elseif x_NIP_real(pos) > 1 && x_NIP_real(pos) < 5
           double_aperture = pi/10;
       else
           double_aperture = pi/20;
       end
       angulo_NIP_real(1,pos) = ((n-(raios_picks+1)/2)/((raios_picks-1)/2))*double_aperture + atan(-0.2*(2*pi/6)*sin(x_NIP_real(pos)*2*pi/6));
       pos=pos+1;
    end
end
num_points = 16; delta=5.6/(num_points-1); 
for i=1:num_points
    for n = 1:raios_picks
        x_NIP_real(pos) = (i-1)*delta+0.2; z_NIP_real(pos) = 2.5 - 0.2*cos(x_NIP_real(pos)*2*pi/6); 
       if x_NIP_real(pos) > 2.5 && x_NIP_real(pos) < 3.5
           double_aperture = pi/4; 
       elseif x_NIP_real(pos) > 1.5 && x_NIP_real(pos) < 4.5
           double_aperture = pi/8;
       elseif x_NIP_real(pos) > 1 && x_NIP_real(pos) < 5
           double_aperture = pi/16;
       else
           double_aperture = pi/24;
       end
        angulo_NIP_real(1,pos) = ((n-(raios_picks+1)/2)/((raios_picks-1)/2))*double_aperture + atan(+0.2*(2*pi/6)*sin(x_NIP_real(pos)*2*pi/6));
        pos=pos+1;
    end
end
num_points = 17; delta=5.6/(num_points-1); 
for i=1:num_points
    for n = 1:raios_picks
        x_NIP_real(pos) = (i-1)*delta+0.2; z_NIP_real(pos) = 2.9 - 0.15*cos(x_NIP_real(pos)*2*pi/6); 
       if x_NIP_real(pos) > 2.5 && x_NIP_real(pos) < 3.5
           double_aperture = pi/5; 
       elseif x_NIP_real(pos) > 1.5 && x_NIP_real(pos) < 4.5
           double_aperture = pi/10;
       elseif x_NIP_real(pos) > 1 && x_NIP_real(pos) < 5
           double_aperture = pi/20;
       else
           double_aperture = pi/25;
       end
        angulo_NIP_real(1,pos) = ((n-(raios_picks+1)/2)/((raios_picks-1)/2))*double_aperture + atan(+0.15*(2*pi/6)*sin(x_NIP_real(pos)*2*pi/6));
        pos=pos+1;
    end
end
num_points = 18; delta=5.5/(num_points-1); 
for i=1:num_points
    for n = 1:raios_picks
        x_NIP_real(pos) = (i-1)*delta+0.25; z_NIP_real(pos) = 3.3 - 0.1*cos(x_NIP_real(pos)*2*pi/6); 
       if x_NIP_real(pos) > 2.5 && x_NIP_real(pos) < 3.5
           double_aperture = pi/6; 
       elseif x_NIP_real(pos) > 1.5 && x_NIP_real(pos) < 4.5
           double_aperture = pi/12;
       elseif x_NIP_real(pos) > 1 && x_NIP_real(pos) < 5
           double_aperture = pi/18;
       else
           double_aperture = pi/24;
       end
        angulo_NIP_real(1,pos) = ((n-(raios_picks+1)/2)/((raios_picks-1)/2))*double_aperture + atan(+0.1*(2*pi/6)*sin(x_NIP_real(pos)*2*pi/6));
        pos=pos+1;
    end
end
num_points = 19; delta=5.5/(num_points-1); 
for i=1:num_points
    for n = 1:raios_picks
        x_NIP_real(pos) = (i-1)*delta+0.25; z_NIP_real(pos) = 3.7; 
       if x_NIP_real(pos) > 2.5 && x_NIP_real(pos) < 3.5
           double_aperture = pi/7; 
       elseif x_NIP_real(pos) > 1.5 && x_NIP_real(pos) < 4.5
           double_aperture = pi/14;
       elseif x_NIP_real(pos) > 1 && x_NIP_real(pos) < 5
           double_aperture = pi/21;
       else
           double_aperture = pi/28;
       end
       angulo_NIP_real(1,pos) = ((n-(raios_picks+1)/2)/((raios_picks-1)/2))*double_aperture;
       pos=pos+1;
    end
end
num_points = 20; delta=5.4/(num_points-1); 
for i=1:num_points
    for n = 1:raios_picks
        x_NIP_real(pos) = (i-1)*delta+0.3; z_NIP_real(pos) = 4.1; 
       if x_NIP_real(pos) > 2.5 && x_NIP_real(pos) < 3.5
           double_aperture = pi/8; 
       elseif x_NIP_real(pos) > 1.5 && x_NIP_real(pos) < 4.5
           double_aperture = pi/16;
       elseif x_NIP_real(pos) > 1 && x_NIP_real(pos) < 5
           double_aperture = pi/24;
       else
           double_aperture = pi/32;
       end
       angulo_NIP_real(1,pos) = ((n-(raios_picks+1)/2)/((raios_picks-1)/2))*double_aperture;
       pos=pos+1;
    end
end
num_points = 21; delta=5.3/(num_points-1); 
for i=1:num_points
    for n = 1:raios_picks
        x_NIP_real(pos) = (i-1)*delta+0.35; z_NIP_real(pos) = 4.5; 
       if x_NIP_real(pos) > 2.5 && x_NIP_real(pos) < 3.5
           double_aperture = pi/8; 
       elseif x_NIP_real(pos) > 1.5 && x_NIP_real(pos) < 4.5
           double_aperture = pi/16;
       elseif x_NIP_real(pos) > 1 && x_NIP_real(pos) < 5
           double_aperture = pi/24;
       else
           double_aperture = pi/32;
       end
       angulo_NIP_real(1,pos) = ((n-(raios_picks+1)/2)/((raios_picks-1)/2))*double_aperture;
       pos=pos+1;
    end
end
num_points = 22; delta=5.2/(num_points-1); 
for i=1:num_points
    for n = 1:raios_picks
        x_NIP_real(pos) = (i-1)*delta+0.4; z_NIP_real(pos) = 4.9; 
       if x_NIP_real(pos) > 2.5 && x_NIP_real(pos) < 3.5
           double_aperture = pi/8; 
       elseif x_NIP_real(pos) > 1.5 && x_NIP_real(pos) < 4.5
           double_aperture = pi/16;
       elseif x_NIP_real(pos) > 1 && x_NIP_real(pos) < 5
           double_aperture = pi/24;
       else
           double_aperture = pi/32;
       end
       angulo_NIP_real(1,pos) = ((n-(raios_picks+1)/2)/((raios_picks-1)/2))*double_aperture;
       pos=pos+1;
    end
end

%Construir Data Space
picks = length(z_NIP_real); z_NIP_real = -z_NIP_real;
px_NIP_real = zeros(1,picks); pz_NIP_real = zeros(1,picks);
%P2_real = zeros(1,picks); Q1_real = zeros(1,picks);
for r=1:picks
    cord_x = max(1,min(round(x_NIP_real(1,r)/espacamento),tamanho_x));
    cord_z = max(1,min(round(-z_NIP_real(1,r)/espacamento),tamanho_z));
    px_NIP_real(1,r) = sin(angulo_NIP_real(1,r))/V_real(cord_z,cord_x);
    pz_NIP_real(1,r) = cos(angulo_NIP_real(1,r))/V_real(cord_z,cord_x);
end

%Componentes de Dado
x_matrix = zeros(4000,length(x_NIP_real)); z_matrix = zeros(4000,length(x_NIP_real)); 
x_matrix(1,:) = x_NIP_real; z_matrix(1,:) = z_NIP_real;
x_dado = x_NIP_real; px_dado = px_NIP_real;
z_dado = z_NIP_real; pz_dado = pz_NIP_real;
t_dado = zeros(1,picks); t_total_dado = zeros(1,picks/raios_picks);

%Tracamento de raios
for r=1:picks
    k = 1;
    while (z_matrix(k,r)<0)
        %Determinar Coordenada Atual
        cord_x = max(1,min(round(x_matrix(k,r)/espacamento),tamanho_x));
        cord_z = max(1,min(round(-z_matrix(k,r)/espacamento),tamanho_z));
        %Atualiza Coordenadas
        x_matrix(k+1,r) = x_matrix(k,r) + h*V_real(cord_z,cord_x)*px_dado(1,r);
        z_matrix(k+1,r) = z_matrix(k,r) + h*V_real(cord_z,cord_x)*pz_dado(1,r);
        %Atualiza vagarosidades
        px_dado(1,r) = px_dado(1,r) - h*(1/(V_real(cord_z,cord_x)^2))*Deriv_x_real(cord_z,cord_x);
        pz_dado(1,r) = pz_dado(1,r) - h*(1/(V_real(cord_z,cord_x)^2))*Deriv_z_real(cord_z,cord_x);
        %Atualiza tempo de trânsito
        t_dado(1,r) = t_dado(1,r) + h*(1/V_real(cord_z,cord_x));
        k = k+1;
        x_dado(1,r) = x_matrix(k,r); z_dado(1,r) = z_matrix(k,r);
    end
end

%Calcular tempo total de propagação de todos os raios que partem de um mesmo ponto em profundidade, da fonte até o receptor.
cont = 0;
for r=1:(picks/raios_picks)
    for rr=1:raios_picks
        t_total_dado(1,r) = t_total_dado(1,r) + t_dado(1,raios_picks*cont+rr);
    end
    cont = cont + 1;
end
Vetor_dado = [t_total_dado,px_dado,x_dado]; %A ser utilizado na avaliação da função objetivo.

end
