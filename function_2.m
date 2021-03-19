function [norma,FObj_sem_Reg,passo,V,x_modelo,z_modelo] = Tomografia_CRP_camadas_sal_4_oitava_passo_0002(x_dado,px_dado,t_total_dado,V_real)

% # Copyright (C) 2021 Gustavo Barroso Dias Ignácio
% <ignaciogbd@ggaunicamp.com>
% <ignaciogbd@gmail.com>

%%%%%%%%% CRP Tomography %%%%%%%%%%%
tic

%%%%%%%%%% Melhor Peso Inicial %%%%%%%%%%%%%%%%
peso_ini = 4/100000000;
peso_ini_x = 1;
peso_ini_z = 1;

%%%%%%%%% Inicializa Variáveis %%%%%%%%
qtd_it = 22; it = 1; max_num_it = 10; picks = length(x_dado); raios_picks = 20;
peso_reg=zeros(1,qtd_it); peso_x= peso_ini_x; peso_z=peso_ini_z; peso_v=0.0001; passo=zeros(1,qtd_it)+1; h=0.002; %Pode usar 0.005 sem prejuizo qualitativo. Com 0.001 fica um epsilon melhor
suporte_x = [0,0,0,linspace(0,6,11),6,6,6]; nx = length(suporte_x);
suporte_z = [0,0,0,linspace(0,5,11),5,5,5]; nz = length(suporte_z);
%suporte_x = [0,0,0,0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5,3.75,4,4.25,4.5,4.75,5,5,5,5]; nx = length(suporte_x);
%suporte_z = [0,0,0,0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,4,4,4]; nz = length(suporte_z);
espacamento = 0.01; tamanho_x = length(0:espacamento:suporte_x(nx-3)); tamanho_z = length(0:espacamento:suporte_z(nz-3));
V = zeros(tamanho_z,tamanho_x,qtd_it); Matriz_coef = zeros(nz-4,nx-4,qtd_it)+V_real(1,200); %V_norm = zeros(tamanho_z,tamanho_x,qtd_it);

% for j=1:tamanho_z
%     for i=1:tamanho_x   
%         if ((mod(j-1,20) == 0 || j == tamanho_z) && (mod(i-1,25) == 0 || i == tamanho_x))
%             Matriz_coef((j-1)/20+2,(i-1)/25+2,1) = 1.5 + 1.5*(j-1)/tamanho_z;% + ((2/5)*(i-1)/tamanho_x)^2;
%         end
%     end
% end
% Matriz_coef(1,:,1) = Matriz_coef(2,:,1);
% Matriz_coef(:,1,1) = Matriz_coef(:,2,1);
% Matriz_coef(nz-4,:,1) = Matriz_coef(nz-5,:,1);
% Matriz_coef(:,nx-4,1) = Matriz_coef(:,nx-5,1);

%Construção das B-Splines
Pre_Bx = struct('form',strings(2,1),'knots',zeros(5,1),'coefs',0,'number',0,'order',0,'dim',0);
Pre_Bz = struct('form',strings(2,1),'knots',zeros(5,1),'coefs',0,'number',0,'order',0,'dim',0);
Bx = repmat(Pre_Bx,nx-4,1); Bz = repmat(Pre_Bz,nz-4,1);
Vet_Bx = zeros(nx-4,tamanho_x); Vet_Bz = zeros(nz-4,tamanho_z);
Deriv_Beta_x = zeros(nx-4,tamanho_x); Deriv_Beta_z = zeros(nz-4,tamanho_z);
Deriv_Beta_x2 = zeros(nx-4,tamanho_x); Deriv_Beta_z2 = zeros(nz-4,tamanho_z);
Int_Prod_Beta_x0 = zeros(nx-4,nx-4); Int_Prod_Beta_z0 = zeros(nz-4,nz-4);
Int_Prod_Beta_x2 = zeros(nx-4,nx-4); Int_Prod_Beta_z2 = zeros(nz-4,nz-4);
D_xx = zeros((nx-4)*(nz-4),(nx-4)*(nz-4)); D_zz = zeros((nx-4)*(nz-4),(nx-4)*(nz-4)); D_int = zeros((nx-4)*(nz-4),(nx-4)*(nz-4));

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
%Construir Matrizes para Regulariza��o
for i=3:(nx-2)
    for j=i:(nx-2)
        Prod_Beta = Vet_Bx(i-2,:).*Vet_Bx(j-2,:);
        Int_Prod_Beta_x0(i-2,j-2) = trapz(espacamento,Prod_Beta); %Integral do produto entre as Beta Splines.
        Int_Prod_Beta_x0(j-2,i-2) = Int_Prod_Beta_x0(i-2,j-2); %Matriz � sim�trica.
        Prod_Beta_x2 = Deriv_Beta_x2(i-2,:).*Deriv_Beta_x2(j-2,:);
        Int_Prod_Beta_x2(i-2,j-2) = trapz(espacamento,Prod_Beta_x2); %Integral do produto entre as derivadas de segunda ordem das Beta Splines
        Int_Prod_Beta_x2(j-2,i-2) = Int_Prod_Beta_x2(i-2,j-2); %Matriz � sim�trica
    end
end
for i=3:(nz-2)
    for j=i:(nz-2)
        Prod_Beta = Vet_Bz(i-2,:).*Vet_Bz(j-2,:);
        Int_Prod_Beta_z0(i-2,j-2) = trapz(espacamento,Prod_Beta); %Integral do produto entre as Beta Splines.
        Int_Prod_Beta_z0(j-2,i-2) = Int_Prod_Beta_z0(i-2,j-2); %Matriz � sim�trica.
        Prod_Beta_z2 = Deriv_Beta_z2(i-2,:).*Deriv_Beta_z2(j-2,:);
        Int_Prod_Beta_z2(i-2,j-2) = trapz(espacamento,Prod_Beta_z2); %Integral do produto entre as derivadas de segunda ordem das Beta Splines
        Int_Prod_Beta_z2(j-2,i-2) = Int_Prod_Beta_z2(i-2,j-2); %Matriz � sim�trica
    end
end
for i=3:(nx-2)
    for j=3:(nz-2)
        for k=3:(nx-2)
            for l=3:(nz-2)
                D_xx((i-3)*(nz-4)+(j-2),(k-3)*(nz-4)+(l-2)) = Int_Prod_Beta_x2(i-2,k-2)*Int_Prod_Beta_z0(j-2,l-2);
                D_zz((i-3)*(nz-4)+(j-2),(k-3)*(nz-4)+(l-2)) = Int_Prod_Beta_x0(i-2,k-2)*Int_Prod_Beta_z2(j-2,l-2);
                D_int((i-3)*(nz-4)+(j-2),(k-3)*(nz-4)+(l-2)) = Int_Prod_Beta_x0(i-2,k-2)*Int_Prod_Beta_z0(j-2,l-2);
            end
        end
    end
end
D = peso_x*D_xx + peso_z*D_zz + peso_v*D_int;%Matriz a ser utilizada na regulariza��o.

%Interpolação do Modelo de Velocidade Via B-Splines
V(:,:,it) = ((Matriz_coef(:,:,it)*Vet_Bx)'*Vet_Bz)'; %Velocidade Inicial

%Calcular Derivadas
[Deriv_x,Deriv_z] = gradient(V(:,:,it),espacamento); Deriv_z = -Deriv_z;
[~,Deriv_z2] = gradient(Deriv_z(:,:),espacamento); Deriv_z2 = -Deriv_z2;
[Deriv_x2,~] = gradient(Deriv_x(:,:),espacamento);

%Traçamento Cinemático de Raio para inicializar as componentes de modelo.
% Sentido de propagação é da superfície para profundidade.

%Inicializa vetores e matrizes
x_modelo = zeros(qtd_it,picks); x_modelo(1,:) = x_dado(1,:);
z_modelo = zeros(qtd_it,picks); t_modelo = zeros(1,picks);
px_modelo = zeros(qtd_it,picks); px_modelo(1,:) = -px_dado(1,:); %Sinal negativo para inverter sentido da propagação.
pz_modelo = zeros(qtd_it,picks); pz_dado = zeros(1,picks);
for r=1:picks
    cord_x = max(1,min(round(x_modelo(1,r)/espacamento),tamanho_x));
    pz_dado(1,r) = sqrt((1/V(1,cord_x,1))^2 - (px_modelo(1,r))^2); 
end
pz_modelo(1,:) = -pz_dado(1,:); angulo_modelo = zeros(qtd_it,picks); %Sinal negativo para inverter sentido da propagação.


rr = 1; cont = 0;
for r=1:picks
    if cont == raios_picks
        cont = 0;
        rr = rr + 1;
    end
    cont = cont + 1;
    while (t_modelo(1,r) < t_total_dado(1,rr)/raios_picks)
        %Determinar Coordenada Atual
        cord_x = max(1,min(round(x_modelo(1,r)/espacamento),tamanho_x));
        cord_z = max(1,min(round(-z_modelo(1,r)/espacamento),tamanho_z));
        %Atualiza Coordenadas
        x_modelo(1,r) = x_modelo(1,r) + h*V(cord_z,cord_x,1)*px_modelo(1,r);
        z_modelo(1,r) = z_modelo(1,r) + h*V(cord_z,cord_x,1)*pz_modelo(1,r);
        %Atualiza vagarosidades
        px_modelo(1,r) = px_modelo(1,r) - h*(1/(V(cord_z,cord_x,1)^2))*Deriv_x(cord_z,cord_x);
        pz_modelo(1,r) = pz_modelo(1,r) - h*(1/(V(cord_z,cord_x,1)^2))*Deriv_z(cord_z,cord_x);
        %Atualiza tempo de trânsito
        t_modelo(1,r) = t_modelo(1,r) + h*(1/V(cord_z,cord_x,1));
    end
    px_modelo(1,r) = -px_modelo(1,r); pz_modelo(1,r) = -pz_modelo(1,r); %Alterar sinal para inverter sentido da propagação.
    angulo_modelo(1,r) = atan(px_modelo(1,r)/pz_modelo(1,r));
end

%Estipular um x_modelo e um z_modelo iniciais
cont = 0; x_aux = 0; z_aux = 0;
for r=1:picks
    x_aux = x_aux + x_modelo(1,r);
    z_aux = z_aux + z_modelo(1,r);
    cont = cont + 1;
    if cont == raios_picks
        x_modelo(1,(r-raios_picks+1):r) = x_aux/raios_picks;
        z_modelo(1,(r-raios_picks+1):r) = z_aux/raios_picks;
        cont = 0;
        x_aux = 0; z_aux = 0;
    end
end


% %%%%%%%%%%%% Início das Iterções Tomográficas %%%%%%%%%%%%%%%%%%
% 
%Inicializa Componentes de Dado Simuladas ao longo da Tomografia
x_simulado = zeros(qtd_it,picks); px_simulado = zeros(qtd_it,picks);
t_simulado = zeros(qtd_it,picks); t_total_simulado = zeros(qtd_it,picks/raios_picks);
z_simulado = zeros(qtd_it,picks); pz_simulado = zeros(qtd_it,picks); 
P1_simulado = zeros(qtd_it,picks); P2_simulado = zeros(qtd_it,picks); Q1_simulado = zeros(qtd_it,picks); Q2_simulado = zeros(qtd_it,picks);
%Inicializas derivadas de Frechet e parâmetros associados
cord_x_NIP = zeros(qtd_it,picks); cord_z_NIP = zeros(qtd_it,picks);
Delt_Delx = zeros(qtd_it,picks); Delt_Delz = zeros(qtd_it,picks); Delt_Delv = zeros(qtd_it,picks,nx-4,nz-4);
Delx_Delx = zeros(qtd_it,picks); Delx_Delz = zeros(qtd_it,picks); Delx_Delang = zeros(qtd_it,picks); Delx_Delv = zeros(qtd_it,picks,nx-4,nz-4);
Delp_Delx = zeros(qtd_it,picks); Delp_Delz = zeros(qtd_it,picks); Delp_Delang = zeros(qtd_it,picks); Delp_Delv = zeros(qtd_it,picks,nx-4,nz-4);
Int_Delq_Delv = zeros(qtd_it,picks,nx-4,nz-4); Int_Delp_Delv = zeros(qtd_it,picks,nx-4,nz-4); 
F = zeros(2*picks+picks/raios_picks,2*(picks/raios_picks)+picks+(nx-4)*(nz-4),qtd_it); Reg = zeros(1,qtd_it); Delta_modelo = zeros(qtd_it,2*(picks/raios_picks) + picks + (nx-4)*(nz-4));
Vetor_simulado = zeros(qtd_it,2*picks+picks/raios_picks); Vetor_delta = zeros(qtd_it,2*picks+picks/raios_picks); Funcao_Objetivo = zeros(1,qtd_it); FObj_sem_Reg = zeros(1,qtd_it);
sigma_t = 0.005; sigma_x = 0.005; sigma_px = 0.01;
%sigma_t = 0.01*raios_picks/2; sigma_x = 0.01; sigma_px = 0.002;
%sigma_t = 0.01; sigma_x = 0.01; sigma_px = 0.002;
vetor_coef_v = zeros(1,(nx-4)*(nz-4)); num_it = zeros(1,qtd_it)+1; perturbacao = zeros(1,qtd_it); perturbacao_2 = zeros(1,qtd_it);
%vetor_erro_t = max(t_total_dado.*0.01,0.001); vetor_erro_px = max(px_dado*0.1,0.001); vetor_erro_x = ones(1,picks)*0.001;
C=[eye(picks/raios_picks)*sigma_t,zeros(picks/raios_picks,picks),zeros(picks/raios_picks,picks);zeros(picks,picks/raios_picks),eye(picks)*sigma_px,zeros(picks,picks);zeros(picks,picks/raios_picks),zeros(picks,picks),eye(picks)*sigma_x]; 
%C=[diag(vetor_erro_t,0),zeros(picks/raios_picks,picks),zeros(picks/raios_picks,picks);zeros(picks,picks/raios_picks),diag(vetor_erro_px,0),zeros(picks,picks);zeros(picks,picks/raios_picks),zeros(picks,picks),diag(vetor_erro_x,0)]; 
%C=eye(2*picks+picks/raios_picks);
C_matrix = C^(-1); C_matrix_meio = C^(-1/2);
norma = zeros(1,qtd_it); Vetor_dado = [t_total_dado,px_dado,x_dado];
%norma_pos = zeros(1,qtd_it);
%t_total_delta = zeros(1,qtd_it); x_delta = zeros(1,qtd_it); px_delta = zeros(1,qtd_it);

%Traçamento dinâmico de raios para simular variáveis a partir das variáveis de modelo atuais
it = 1;
while (it <= qtd_it)
    controle = 0; %vari�vel que avalia se houve melhora na fun��o objetivo e utilizada como crit�rio de parada.
    while(controle ~= 2019)
        if(it>1)
            %Resolu��o do Problema de Otimiza��o Linearizado
            if((it>2) && (perturbacao(1,it)==0) && (perturbacao_2(1,it)==0) && (num_it(1,it) == 1))
                peso_reg(1,it) = min(sqrt(FObj_sem_Reg(1,it-1)/FObj_sem_Reg(1,it-2))*peso_reg(1,it-1),peso_reg(1,it-1));
                %peso_reg(1,it) = (peso_ini)*(Vetor_delta(it-1,:)*C_matrix*Vetor_delta(it-1,:)')/Reg(1,it-1);
                %peso_reg(1,it) = peso_reg(1,it-1);
            elseif((it == 2) && (perturbacao(1,it)==0) && (perturbacao_2(1,it)==0) && (num_it(1,it) == 1))
                peso_reg(1,it) = peso_reg(1,it-1);
            end
            if (num_it(1,it) == 1)%S� se calcula uma vez por itera��o
                B = chol(peso_reg(1,it)*D);
                F_ampliada = [C_matrix_meio*F(:,:,it-1);zeros((nx-4)*(nz-4),2*(picks/raios_picks)+picks),B];
                Vetor_delta_ampliado = [C_matrix_meio*Vetor_delta(it-1,:)';-B*vetor_coef_v(1,:)'];
                F_ampliada(isnan(F_ampliada)) = 0; F_ampliada(isinf(F_ampliada)) = 0;%Tratamento para NaN e Inf.
                [U_svd,S_svd,V_svd] = svd(F_ampliada,'econ'); Delta_modelo(it,:) = V_svd*(inv(S_svd))*U_svd'*Vetor_delta_ampliado;
            end
            %Atualiza Componentes de Modelo
            for rr = 1:picks/raios_picks
                x_modelo(it,((rr-1)*raios_picks+1):((rr-1)*raios_picks+raios_picks)) = x_modelo(it-1,((rr-1)*raios_picks+1):((rr-1)*raios_picks+raios_picks)) + passo(1,it)*Delta_modelo(it,rr);
                z_modelo(it,((rr-1)*raios_picks+1):((rr-1)*raios_picks+raios_picks)) = z_modelo(it-1,((rr-1)*raios_picks+1):((rr-1)*raios_picks+raios_picks)) + passo(1,it)*Delta_modelo(it,picks/raios_picks+rr);
            end
            angulo_modelo(it,:) = angulo_modelo(it-1,:) + passo(1,it)*Delta_modelo(it,(2*picks/raios_picks+1):(2*picks/raios_picks+picks));
            for r=1:picks
                if(z_modelo(it,r)>0)
                    z_modelo(it,r) = -espacamento;
                end
            end
            for i=3:(nx-2)
                Matriz_coef(:,i-2,it) = Matriz_coef(:,i-2,it-1) + passo(1,it)*Delta_modelo(it,((2*picks/raios_picks+picks+(nz-4)*(i-3)+1):2*picks/raios_picks+picks+(nz-4)*(i-2)))';
            end
            %Novo Modelo de Velocidade
            V(:,:,it) = ((Matriz_coef(:,:,it)*Vet_Bx)'*Vet_Bz)';
            [Deriv_x,Deriv_z] = gradient(V(:,:,it),espacamento); Deriv_z = -Deriv_z;
            [~,Deriv_z2] = gradient(Deriv_z(:,:),espacamento); Deriv_z2 = -Deriv_z2;
            [Deriv_x2,~] = gradient(Deriv_x(:,:),espacamento);
        end

        %Atribui Valores Iniciais
        x_simulado(it,:) = x_modelo(it,:); z_simulado(it,:) = z_modelo(it,:); t_simulado(it,:) = zeros(1,picks);
        P1_simulado(it,:) = zeros(1,picks); Q2_simulado(it,:) = zeros(1,picks);
        for r=1:picks
            cord_x_NIP(it,r) = max(1,min(round(x_simulado(it,r)/espacamento),tamanho_x));
            cord_z_NIP(it,r) = max(1,min(round(-z_simulado(it,r)/espacamento),tamanho_z));
            px_simulado(it,r) = sin(angulo_modelo(it,r))/V(cord_z_NIP(it,r),cord_x_NIP(it,r),it);
            pz_simulado(it,r) = cos(angulo_modelo(it,r))/V(cord_z_NIP(it,r),cord_x_NIP(it,r),it);
            P2_simulado(it,r) = 1/V(cord_z_NIP(it,r),cord_x_NIP(it,r),it); %Condição de fonte pontual.
            Q1_simulado(it,r) = 1; %Condição de onda plana.
        end
        
        %Inicializa ou reinicializa integrais e matriz de Frechet
        Delt_Delv(it,:,:,:) = zeros(1,picks,nx-4,nz-4); Int_Delq_Delv(it,:,:,:) = zeros(1,picks,nx-4,nz-4); Int_Delp_Delv(it,:,:,:) = zeros(1,picks,nx-4,nz-4); %Int_DelM_Delv(it,:,:,:) = zeros(1,picks,nx-4,nz-4); Int_DelM_Delq(it,:,:) = zeros(1,1,picks); Int_DelM_Delp(it,:,:) = zeros(1,1,picks);
        F(:,:,it) = zeros(2*picks+picks/raios_picks,2*(picks/raios_picks)+picks+(nx-4)*(nz-4),1);
        %Inicializa tempo total
        t_total_simulado(it,:) = zeros(1,picks/raios_picks);
        
        cont = 0; rr = 1;
        for r=1:picks
            aux = 0;
            while (z_simulado(it,r)<0 && aux < 1000000)
                aux = aux+1;
                %Determinar Coordenada Atual
                cord_x = max(1,min(round(x_simulado(it,r)/espacamento),tamanho_x));
                cord_z = max(1,min(round(-z_simulado(it,r)/espacamento),tamanho_z));
                %Atualiza Coordenadas
                x_simulado(it,r) = x_simulado(it,r) + h*V(cord_z,cord_x,it)*px_simulado(it,r); dx = h*V(cord_z,cord_x,it)*px_simulado(it,r); 
                z_simulado(it,r) = z_simulado(it,r) + h*V(cord_z,cord_x,it)*pz_simulado(it,r); dz = h*V(cord_z,cord_x,it)*pz_simulado(it,r);
                %Atualiza vagarosidades
                px_simulado(it,r) = px_simulado(it,r) - h*(1/(V(cord_z,cord_x,it)^2))*Deriv_x(cord_z,cord_x);
                pz_simulado(it,r) = pz_simulado(it,r) - h*(1/(V(cord_z,cord_x,it)^2))*Deriv_z(cord_z,cord_x);
                if pz_simulado(it,r) < 0 %controle para casos esp\FArios
                    z_simulado(it,r) = 0.001; %crit\E9rio de parada
                end
                %Atualiza tempo de trânsito
                t_simulado(it,r) = t_simulado(it,r) + h*(1/V(cord_z,cord_x,it));
                %Coordenadas de raio
                Q1_simulado(it,r) = Q1_simulado(it,r) + h*V(cord_z,cord_x,it)*P1_simulado(it,r);
                Q2_simulado(it,r) = Q2_simulado(it,r) + h*V(cord_z,cord_x,it)*P2_simulado(it,r);
                %Ângulo e Derivadas
                angulo = atan(px_simulado(it,r)/pz_simulado(it,r)); %Ângulo em relação ao eixo z (aponta para cima).
                Vqq = Deriv_x2(cord_z,cord_x)*(cos(angulo))^2 + Deriv_z2(cord_z,cord_x)*(sin(angulo))^2; %Derivadas parciais de segunda ordem da velocidade em relação às coordenadas centradas no raio.
                %Coordenadas de raio
                P2_simulado(it,r) = P2_simulado(it,r) - h*((1/V(cord_z,cord_x,it))^2)*Vqq*Q2_simulado(it,r);
                P1_simulado(it,r) = P1_simulado(it,r) - h*((1/V(cord_z,cord_x,it))^2)*Vqq*Q1_simulado(it,r);
                
                %if(num_it(1,it) == 1)%S� se calcula uma vez por itera��o
                    %Parâmetros para Derivadas de Frechet
                    delta = sqrt(dx^2 + dz^2);
                    Delv_Delq = cos(angulo)*Deriv_x(cord_z,cord_x) - sin(angulo)*Deriv_z(cord_z,cord_x);
                    %Integrais para Derivadas de Frechet
                    for i=3:(nx-2)
                        for k=3:(nz-2)
                            DelprodB_Delq = -sin(angulo)*Deriv_Beta_z(k-2,cord_z)*Vet_Bx(i-2,cord_x)+cos(angulo)*Deriv_Beta_x(i-2,cord_x)*Vet_Bz(k-2,cord_z);
                            Delt_Delv(it,r,i-2,k-2) = Delt_Delv(it,r,i-2,k-2) - (Vet_Bx(i-2,cord_x)*Vet_Bz(k-2,cord_z)/(V(cord_z,cord_x,it)^2))*delta; 
                            Int_Delq_Delv(it,r,i-2,k-2) = Int_Delq_Delv(it,r,i-2,k-2) - Q2_simulado(it,r)*(-(1/V(cord_z,cord_x,it)^2)*DelprodB_Delq + (1/(V(cord_z,cord_x,it)^3))*Delv_Delq*(Vet_Bx(i-2,cord_x)*Vet_Bz(k-2,cord_z)))*delta;
                            Int_Delp_Delv(it,r,i-2,k-2) = Int_Delp_Delv(it,r,i-2,k-2) + Q1_simulado(it,r)*(-(1/V(cord_z,cord_x,it)^2)*DelprodB_Delq + (1/(V(cord_z,cord_x,it)^3))*Delv_Delq*(Vet_Bx(i-2,cord_x)*Vet_Bz(k-2,cord_z)))*delta;
                        end
                    end
            end

            %%%%% Derivadas de Frechet %%%%%%
            %Tempo
            Delt_Delx(it,r) = -sin(angulo_modelo(it,r))/V(cord_z_NIP(it,r),cord_x_NIP(it,r),it); 
            Delt_Delz(it,r) = -cos(angulo_modelo(it,r))/V(cord_z_NIP(it,r),cord_x_NIP(it,r),it); 
            %Ponto de Emergência
            Delx_Delx(it,r) = Q1_simulado(it,r)*cos(angulo_modelo(it,r))/cos(angulo); 
            Delx_Delz(it,r) = -Q1_simulado(it,r)*sin(angulo_modelo(it,r))/cos(angulo); 
            Delx_Delang(it,r) = Q2_simulado(it,r)/(V(cord_z_NIP(it,r),cord_x_NIP(it,r),it)*cos(angulo)); 
            %Vagarosidade Horizontal
            Delp_Delx(it,r) = cos(angulo)*cos(angulo_modelo(it,r))*P1_simulado(it,r); 
            Delp_Delz(it,r) = -cos(angulo)*sin(angulo_modelo(it,r))*P1_simulado(it,r);
            Delp_Delang(it,r) = cos(angulo)*P2_simulado(it,r)/V(cord_z_NIP(it,r),cord_x_NIP(it,r),it); 
             
            %Derivadas em Relação aos nós de interpolacao do modelo de velocidade            
            if cont == raios_picks
                cont = 0;
                rr = rr + 1;
            end
            cont = cont + 1;
             
             for i=3:(nx-2)
                for k=3:(nz-2)
                    %F(r,3*picks+(i-3)*(nz-4)+(k-2),it) = Delt_Delv(it,r,i-2,k-2);
                    Delx_Delv(it,r,i-2,k-2) = (1/cos(angulo))*(Q1_simulado(it,r)*Int_Delq_Delv(it,r,i-2,k-2) + Q2_simulado(it,r)*Int_Delp_Delv(it,r,i-2,k-2)); %F(r+3*picks,3*picks+(i-3)*(nz-4)+(k-2),it)=Delx_Delv(it,r,i-2,k-2);
                    Delp_Delv(it,r,i-2,k-2) = cos(angulo)*(P1_simulado(it,r)*Int_Delq_Delv(it,r,i-2,k-2) + P2_simulado(it,r)*Int_Delp_Delv(it,r,i-2,k-2)); %F(r+2*picks,3*picks+(i-3)*(nz-4)+(k-2),it)=Delp_Delv(it,r,i-2,k-2);

                    %Derivadas Frechet - Derivadas em Relação à velocidade
                    F(rr,2*picks/raios_picks+picks+(i-3)*(nz-4)+(k-2),it) = F(rr,2*picks/raios_picks+picks+(i-3)*(nz-4)+(k-2),it) + Delt_Delv(it,r,i-2,k-2);
                    F(r+picks/raios_picks,2*picks/raios_picks+picks+(i-3)*(nz-4)+(k-2),it) = Delp_Delv(it,r,i-2,k-2);
                    F(r+picks/raios_picks+picks,2*picks/raios_picks+picks+(i-3)*(nz-4)+(k-2),it) = Delx_Delv(it,r,i-2,k-2);
                    
                end
            end

            %Matriz de Frechet

            %Derivdas do tempo total
            F(rr,rr,it) = F(rr,rr,it) + Delt_Delx(it,r);
            F(rr,rr+picks/raios_picks,it) = F(rr,rr+picks/raios_picks,it) + Delt_Delz(it,r);

            %Derivadas dos ângulos de emergência
            F(r+picks/raios_picks,rr,it) = Delp_Delx(it,r);
            F(r+picks/raios_picks,rr+picks/raios_picks,it) = Delp_Delz(it,r);
            F(r+picks/raios_picks,r+2*picks/raios_picks,it) = Delp_Delang(it,r);
            
            %Derivadas das posições de emergência
            F(r+picks/raios_picks+picks,rr,it) = Delx_Delx(it,r);
            F(r+picks/raios_picks+picks,rr+picks/raios_picks,it) = Delx_Delz(it,r);
            F(r+picks/raios_picks+picks,r+2*picks/raios_picks,it) = Delx_Delang(it,r);
            
            %Construir vetor de tempo total
            t_total_simulado(it,rr) = t_total_simulado(it,rr) + t_simulado(it,r);                 
        end 
        
        %%%Calcular Regularização
        for i=3:(nx-2)
            for j=3:(nz-2)
                vetor_coef_v(1,(i-3)*(nz-4)+(j-2)) = Matriz_coef(j-2,i-2,it);
            end
        end
        Reg(1,it) = vetor_coef_v*D*vetor_coef_v';
        %Calcular Função Objetivo
        Vetor_simulado(it,:) = [t_total_simulado(it,:),px_simulado(it,:),x_simulado(it,:)];
        Vetor_delta(it,:) = Vetor_dado - Vetor_simulado(it,:);
        Vetor_delta(isinf(Vetor_delta)) = 0; Vetor_delta(isnan(Vetor_delta)) = 0; %Tratamento para Nan e Inf.
        if(it==1 && (num_it(1,it) == 1))
            peso_reg(1,it) = (peso_ini)*(Vetor_delta(it,:)*C_matrix*Vetor_delta(it,:)')/Reg(1,it);
            %peso_reg(1,it) = peso_ini;
        end
        
        Funcao_Objetivo(1,it) = (1/2)*(Vetor_delta(it,:)*C_matrix*Vetor_delta(it,:)'+peso_reg(1,it)*Reg(1,it));
        FObj_sem_Reg(1,it) = Funcao_Objetivo(1,it) - (1/2)*(peso_reg(1,it)*Reg(1,it));
        %Avaliar se houve melhora na fun��o objetivo
        if(it>1)
            if (FObj_sem_Reg(1,it)>FObj_sem_Reg(1,it-1))
                passo(1,it) = passo(1,it)*(1/2); %diminui tamanho do passo caso houve piora.
                num_it(1,it) = num_it(1,it) + 1;
                %disp(num_it(1,it));
                if num_it(1,it)>max_num_it
                    num_it(1,it) = 1;
                    passo(1,it) = 1;
                    controle=2019;
                    FObj_sem_Reg(1,qtd_it) = FObj_sem_Reg(1,it-1); V(:,:,qtd_it) = V(:,:,it-1); 
                    x_modelo(qtd_it,:) = x_modelo(it-1,:); z_modelo(qtd_it,:) = z_modelo(it-1,:);
                    it = qtd_it;
                end
            else
                controle = 2019;
            end
        else
            controle = 2019;
        end
        if controle == 2019
            disp(it);
            %disp(passo); %disp(num_it);
            disp(FObj_sem_Reg);
%            V_norm(:,:,it) = V(:,:,it);
%             for j=1:tamanho_z
%                 for i=1:tamanho_x
%                     if V(j,i,it) < 1.5
%                         V_norm(j,i,it) = 1.5;
%                     end
%                     if V(j,i,it) > 4.15
%                         V_norm(j,i,it) = 4.15;
%                     end
%                 end
%             end
            norma(1,it) = norm(abs((V_real(1:460,20:580)-V(1:460,20:580,it))./V_real(1:460,20:580)),'fro');
            %norma_pos(1,it) = sqrt(sum((abs(x_NIP_real(1,:) - x_modelo(it,:))).^2) +  sum((abs(z_NIP_real(1,:) - z_modelo(it,:))).^2));
            %t_total_delta(1,it) = norm(t_total_dado - t_total_simulado(it,:));
            %x_delta(1,it) = norm(x_dado - x_simulado(it,:));
            %px_delta(1,it) = norm(px_dado - px_simulado(it,:));
            %norma_2 = norm(V_real(1:360,20:480)-V(1:360,20:480,it),'fro');
            disp(norma); disp(peso_reg(1,it)); %disp(norma_pos);
            %disp(t_total_delta); disp(x_delta); disp(px_delta);
            it = it + 1;
        end
    end
end
toc
end
