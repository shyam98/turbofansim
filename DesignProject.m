%Shyam Patel
%Design Project
%MAE112 Propulsion
%Professor Feng Liu
%Fall 2019

clear
close all
clc 



%Variables
n_d = 0.97;
n_c = 0.85;
n_b = 1.00;
n_t = 0.90;
n_n = 0.98;

gamma_d = 1.4;
gamma_c = 1.37;
gamma_b = 1.35;
gamma_t = 1.33;
gamma_n = 1.36;

t04 = 1400:50:1700;
R = 287;
q_r = 45000000;

mach_num = 0.80;
alt = 10668;
p_a = 23800;
t_a = 219;

pr_c = 16:2:40;
pr_b = 0.95;

%Fan conditions
pr_f = 1.5;
n_fc = 0.97;
n_fn = 0.99;
gamma_fc = 1.4;
gamma_fn = 1.4;
bypass_0 = 7;
%Max area intake
a_max = 3.14*1.0;

%Begin calculation for ST and TSFC

Cp = (gamma_b/(gamma_b-1))*R;
%Get airspeed
airspeed = sqrt(gamma_d*R*t_a);

u = airspeed * mach_num;
%Calculate inlet stag temperature
t_02 = t_a*(1+((gamma_d-1)/2)*(mach_num^2));
%Calculate inlet stag pressure
p_02 = p_a*(1 + n_d*(t_02/t_a - 1))^(gamma_d/(gamma_d-1));

%Calculate mass flow rate
delt_0 = p_02/p_a;
theta_0 = t_02/t_a;
m_in = a_max*231.8 * delt_0/sqrt(theta_0);

%Calculate compressor exit
p_03 = pr_c*p_02;
t_03 = t_02*(1 + (1/n_d)*(pr_c.^((gamma_c-1)/gamma_c)-1));

%Create zero index for ST and TSFC values
ST_f = zeros(length(t04),length(pr_c));
TSFC_f = ST_f;
    
    %Graph for TSFC and ST for zero bypass
    for tval = 1:1:length(t04)
        %Burner exit conditions for each t04 requirement
        f = ((t04(tval)./t_03) - 1)./(q_r./(Cp*t_03) - (t04(tval)./t_03));
        p_04 = p_03*pr_b;
        
        %Comparing different bpyass ratios
        %Fan exit
        p_08 = pr_f*p_02;
        t_08 = t_02*(1 + (1/n_fc)*((pr_f)^((gamma_fc-1)/gamma_fc)-1));

        %Fan nozzle exit
        ue_f = sqrt(2*n_fn*(gamma_fn/(gamma_fn-1))*R*t_08*(1-((p_a/p_08)^((gamma_fn-1)/gamma_fn))));

        %Turbine exit
        t_05 = t04(tval) - (t_03 - t_02) - (bypass_0*(t_08-t_02));
        p_05 = p_04.*(1 - (1/n_t)*(1 - t_05/t04(tval))).^(gamma_t/(gamma_t-1));

        %Core nozzle exit
        t_06 = t_05;
        p_06 = p_05;

        ue = sqrt(2*n_n*(gamma_n/(gamma_n-1))*R*t_06.*(1 - (p_a./p_06).^((gamma_n-1)/gamma_n)));

        %Specific Thrust and TSFC
        T = m_in*((1+f).*ue + bypass_0*ue_f - (1+bypass_0)*u);
        ST = ((1+f).*ue + bypass_0*ue_f - (1+bypass_0)*u);
        TSFC_bare = f./ST;
        TSFC = (1.04 + 0.01*(bypass_0-1))*TSFC_bare;

        ST_f(tval,:) = ST;
        TSFC_f(tval,:) = TSFC;




        
    end
%Plot graph to choose a t04 value
title('TSFC and ST for \beta = 0 vs. \pi_c')
yyaxis left
plt1 = semilogx(pr_c, ST_f/1000);
ylabel('Specific Thrust $\frac{kN \cdot s}{kg}$','Interpreter','latex')
legend( plt1(:), strcat(num2str(t04(:)),'k'))
yyaxis right
plt2 = semilogx(pr_c, TSFC_f*1000,'HandleVisibility','off');
ylabel('TSFC $\frac{kg}{kN \cdot s}$','Interpreter','latex')
xlabel('Compressor Pressure Ratio \pi_c')
hold on
p = semilogx(pr_c(9),TSFC_f(2,9)*1000,'o','MarkerFaceColor','red', 'DisplayName','Chosen T_0_4');
hold off;

%Chosen Parameters pr_c = 32 t04 = 1450k
t_04 = 1450;
pi_c = 32;
    
%Design for bypass and fan ratio
bypass = 1:1:7;
fpr = 1:0.1:2;
[x_bypass, y_fpr] = meshgrid(bypass, fpr);

%Begin calculation for ST and TSFC

%Get airspeed
u = airspeed * mach_num;

%Calculate compressor exit
p_03 = pi_c*p_02;
t_03 = t_02*(1 + (1/n_d)*(pi_c.^((gamma_c-1)/gamma_c)-1));

%Burner exit conditions for each t04 requirement
f = ((t_04./t_03) - 1)./(q_r./(Cp*t_03) - (t_04./t_03));
p_04 = p_03*pr_b;


for b = 1:1:length(bypass)
    for n = 1:1:length(fpr)
        %Comparing different bpyass ratios
        %Fan exit
        p_08 = fpr(n)*p_02;
        t_08 = t_02*(1 + (1/n_fc)*((fpr(n))^((gamma_fc-1)/gamma_fc)-1));

        %Fan nozzle exit
        ue_f = sqrt(2*n_fn*(gamma_fn/(gamma_fn-1))*R*t_08*(1-((p_a/p_08)^((gamma_fn-1)/gamma_fn))));

        %Turbine exit
        t_05 = t_04 - (t_03 - t_02) - (bypass(b)*(t_08-t_02));
        p_05 = p_04.*(1 - (1/n_t)*(1 - t_05/t_04)).^(gamma_t/(gamma_t-1));

        %Core nozzle exit
        t_06 = t_05;
        p_06 = p_05;
        ue = sqrt(2*n_n*(gamma_n/(gamma_n-1))*R*t_06.*(1 - (p_a./p_06).^((gamma_n-1)/gamma_n)));

        %Specific Thrust and TSFC
        T = m_in*((1+f)*ue + bypass(b)*ue_f - (1 + bypass(b))*u);
        ST = ((1+f).*ue + bypass(b)*ue_f - (1 + bypass(b))*u);
        TSFC_bare = f./ST;
        TSFC = (1.04 + 0.01*(bypass(b)-1))*TSFC_bare;   
        TSFC_tbf(n,b) = TSFC*1000;
        ST_tbf(n,b) = ST/1000;
        
        %Efficiencies
        n_th = (((1+f)*(ue^2/2)-(u^2/2))+(bypass(b)*(ue_f^2/2 - u^2/2)))/(f*q_r);
        n_th_f(n,b) = n_th;
        
        n_p = (T*u)/(((1+f)*(ue^2/2))+(bypass(b)*(ue_f^2/2-u^2/2)));
        n_p_f(n,b) = n_p/1000;
    end

end
figure('NumberTitle', 'off', 'Name', 'TSFC');
[mm nn] = contourf(x_bypass, y_fpr, TSFC_tbf);
title('Contour plot: TSFC vs \beta & \pi_f');
xlabel('Bypass Ratio (\beta)');
ylabel('Fan Pressure Ratio (\pi_f)');
nn.ShowText = 'on' ;
nn. LineColor ='r' ;
c = colorbar();
c.Label.Interpreter = 'latex';
c.Label.String = ('TSFC $\frac{kg}{kN \cdot s}$');

figure('NumberTitle', 'off', 'Name', 'Specific Thrust');
[mm nn] = contourf(x_bypass, y_fpr, ST_tbf);
title('Contour plot: ST vs \beta & \pi_f');
xlabel('Bypass Ratio (\beta)');
ylabel('Fan Pressure Ratio (\pi_f)');
nn.ShowText = 'on' ;
nn. LineColor ='r' ;
c = colorbar();
c.Label.Interpreter = 'latex';
c.Label.String = ('Specific Thrust $\frac{kN \cdot s}{kg}$');

figure('NumberTitle', 'off', 'Name', 'Thermal Efficiency');
[mm nn] = contourf(x_bypass, y_fpr, n_th_f);
title('Thermal Efficiency: \eta_t_h vs \beta & \pi_f');
xlabel('Bypass Ratio (\beta)');
ylabel('Fan Pressure Ratio (\pi_f)');
nn.ShowText = 'on' ;
nn. LineColor ='r' ;
c = colorbar();
c.Label.String = 'Thermal Efficiency \eta_t_h';

figure('NumberTitle', 'off', 'Name', 'Propulsion Efficiency');
[mm nn] = contourf(x_bypass, y_fpr, n_p_f);
title('Propulsion Efficiency: \eta_p vs \beta & \pi_f');
xlabel('Bypass Ratio (\beta)');
ylabel('Fan Pressure Ratio (\pi_f)');
nn.ShowText = 'on' ;
nn. LineColor ='r' ;
c = colorbar();
c.Label.String = 'Propulsion Efficiency \eta_p';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% **** Optimized Bypass and Fan Pressure Ratio **** %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_o = n_p_f.*n_th_f;
n_o_max = max(max(n_o));
[n_o_x,n_o_y]=find(n_o==n_o_max);

figure('NumberTitle', 'off', 'Name', 'Overall Efficiency');
[mm nn] = contourf(x_bypass, y_fpr, n_o);
title('Overall Efficiency: \eta_o vs \beta & \pi_f');
xlabel('Bypass Ratio (\beta)');
ylabel('Fan Pressure Ratio (\pi_f)');
nn.ShowText = 'on' ;
nn. LineColor ='r' ;
c = colorbar();
c.Label.String = 'Overall Efficiency \eta_o';

[bypass_final] = [bypass(n_o_y)]
[fpr_final] = fpr(n_o_x)
ST_final = ST_tbf(n_o_x,n_o_y)
TSFC_final = TSFC_tbf(n_o_x,n_o_y)


%Calculate compressor exit
p_03 = pr_c*p_02;
t_03 = t_02*(1 + (1/n_d)*(pr_c.^((gamma_c-1)/gamma_c)-1));
pr_c = 16:2:40;


    %Graph for TSFC and ST for final bypass and fan ratio
    for tval = 1:1:length(t04)
        %Burner exit conditions for each t04 requirement
        f = ((t04(tval)./t_03) - 1)./(q_r./(Cp*t_03) - (t04(tval)./t_03));
        p_04 = p_03*pr_b;
        
        %Comparing different bpyass ratios
        %Fan exit
        p_08 = fpr_final*p_02;
        t_08 = t_02*(1 + (1/n_fc)*((fpr_final)^((gamma_fc-1)/gamma_fc)-1));

        %Fan nozzle exit
        ue_f = sqrt(2*n_fn*(gamma_fn/(gamma_fn-1))*R*t_08*(1-((p_a/p_08)^((gamma_fn-1)/gamma_fn))));

        %Turbine exit
        t_05 = t04(tval) - (t_03 - t_02) - (bypass_final*(t_08-t_02));
        p_05 = p_04.*(1 - (1/n_t)*(1 - t_05/t04(tval))).^(gamma_t/(gamma_t-1));

        %Core nozzle exit
        t_06 = t_05;
        p_06 = p_05;

        ue = sqrt(2*n_n*(gamma_n/(gamma_n-1))*R*t_06.*(1 - (p_a./p_06).^((gamma_n-1)/gamma_n)));

        %Specific Thrust and TSFC
        T = m_in*((1+f).*ue + bypass_final*ue_f - (1+bypass_final)*u);
        ST = ((1+f).*ue + bypass_final*ue_f - (1+bypass_final)*u);
        TSFC_bare = f./ST;
        TSFC = (1.04 + 0.01*(bypass_final-1))*TSFC_bare;

        ST_f(tval,:) = ST/1000;
        TSFC_f(tval,:) = TSFC*1000;
    end
    
    [prc_x, t04_y] = meshgrid(pr_c, t04);
    o = CarpetPlot(prc_x,t04_y, ST_f, TSFC_f);
    label(o)

    