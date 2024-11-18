%% Mec E 541 Assignment 4 code

clc, close all, clear all
load Engine_Test_Data
%% Engine Geometry
D = 8.6/100;		    % Bore [m]
Lc = 14.55/100;	% Connecting rod length [m]
S = 8.6/100;   % Stroke [m]
Rc = S/2;	% Crank throw radius [m]
cr = 9.2;			% Compression ratio
vol_c = (pi*D^2/4*2*Rc)/(cr-1);	% Clearance volume [m^3]
Vss = (pi*D^2/4*2*Rc);	% Swept volume [m^3]
EVC = 338; % aTDC (22bTDC-Exhaust stroke)
IVC=-182;   % aTDC(2 bBDC-Intake stroke)
EVO = 144;    % 36 (bBDC-Expansion stroke)
IVO = 334.5; % aTDC (25.5 bTDC-Exhaust stroke)
Vivc= [vol_c + (pi*D^2)/4*(Lc + Rc - (Rc*cos(abs(IVC*pi/180)) + (Lc^2 - Rc^2*(sin(abs(IVC)*pi/180))^2)^(1/2)))];         %Volume [m^3]
Vevo= [vol_c + (pi*D^2)/4*(Lc + Rc - (Rc*cos(abs(EVO*pi/180)) + (Lc^2 - Rc^2*(sin(abs(EVO)*pi/180))^2)^(1/2)))];         %Volume [m^3]
Vivo= [vol_c + (pi*D^2)/4*(Lc + Rc - (Rc*cos(abs(IVO*pi/180)) + (Lc^2 - Rc^2*(sin(abs(IVO)*pi/180))^2)^(1/2)))];         %Volume [m^3]
Vevc= [vol_c + (pi*D^2)/4*(Lc + Rc - (Rc*cos(abs(EVC.*pi/180)) + (Lc^2 - Rc^2*(sin(abs(EVC)*pi/180)).^2).^(1/2)))];      %Volume [m^3]


%%  Engine Data
CAD_P(1,:) = Engine_Test_Data(1).CAD; % Crank Angle Degree after top dead center [aTDC]
P_cyc(1,:) = Engine_Test_Data(1).In_cyl_Pressure_kPa; % In-cylinder gas pressure [kPa]
mf_mg_cyc = Engine_Test_Data(1).mf_mg_cyc; %[mg]
LHV_fuel = 44560; % Lower heating value of fuel [kJ/kg]
n_mean_c1=1.3694; % Mean polytropic coefficient including both compression and exapnsion

%% Sorting P-trace engine data
 for i =1:100
     a(1) =1;
     b(1) =720;
     P_cyc_1(i,:) = P_cyc(1,a(i):b(i));   
     CAD_P_1(i,:) =CAD_P(1,a(i):b(i));    
     a(i+1)=a(i)+720;
     b(i+1)=b(i)+720;
 end

%% Plotting data
figure % Plotting pressure trace versus crank angle degrees
for i =1:100
figure(1), plot(CAD_P_1(1,:),P_cyc_1(i,:));hold on
xlabel('Crank Angle Degree [aTDC]'),ylabel('In-Cylinder Gas Pressure [kPa]'),title('Pressure vs. CAD for 100 engine cycles')
end
xlim([-180 540]); xticks([-180 -90 0 90 180 270 360 450 540]); 

%% Cylinder volume and dV
% Remember radians and degrees
for i_w = 1:720
    cos_term(i_w) = cos(CAD_P_1(1, i_w) * pi / 180);
    sin_term(i_w) = sin(CAD_P_1(1, i_w) * pi / 180);
    x(i_w) = Rc * (1 - cos_term(i_w)) + Lc - (sqrt((Lc ^ 2) - ((Rc ^ 2) * (sin_term(i_w)) ^ 2)));
    V_s(i_w) = (pi * D^2 * x(i_w) / 4) + vol_c;
    Vs(i_w) = (V_s(i_w) * (10 ^ 6));
end

% Takes the delta
dV = diff(V_s);

% Test if volume figures looks like the notes
figure(2), plot(CAD_P_1(1, 1:360), Vs(1:360))
xlabel("Crank Angle Degree [aTDC]"), ylabel("Cylinder Volume [cm^3]"), title("Cylinder Volume vs. CAD");

figure(3), plot(Vs(1:360), P_cyc_1(1, 1:360))
xlabel("Cylinder Volume [cm^3]"), ylabel("In-Cylinder Gas Pressure [kPa]"), title("Pressure (1st cycle) vs. Cylinder Volume");

%% Gross work, gross imep, gross indicated thermal efficiency
d_theta = diff(CAD_P_1(1, 1:720));
for i_p = 1:100
    work_cycle(i_p) = max(cumsum(P_cyc_1(i_p, 1:719) .* 10 ^ 3 .* dV(1:719)));
    imep(i_p) = (work_cycle(i_p) / Vss) ./ 1000;

    work_cycle_gross(i_p) = max(cumsum(P_cyc_1(i_p, 1:360) .* 10 ^ 3 .* dV(1:360)));
    imep_gross(i_p) = (work_cycle_gross(i_p) / Vss) ./ 1000;

    nth(i_p) = work_cycle(i_p) / (mf_mg_cyc * 0.001 * LHV_fuel) * 100;
end

figure(4), stem(imep);
xlabel("Cycle Number [-]"), ylabel("IMEP [kPa]"), title("Net IMEP per cycle");
ylim([450, 510]);

figure(5), stem(imep_gross);
xlabel("Cycle Number [-]"), ylabel("IMEP [kPa]"), title("Gross IMEP per cycle");
ylim([450, 510]);

fprintf("Results\n");
mu_imep = mean(imep);
fprintf('Mean net IMEP = %4.2f kPa.\n', mu_imep);
mu_imep_gross = mean(imep_gross);
fprintf('Mean gross IMEP = %4.2f kPa.\n', mu_imep);
sigma_imep = std(imep);
sigma_imep_gross = std(imep_gross);
cov_imep = (sigma_imep / mu_imep) * 100;
fprintf('COV_imep = %4.2f percent.\n', cov_imep);
cov_imep_gross = (sigma_imep_gross / mu_imep_gross) * 100;
fprintf('COV_imep_gross = %4.2f percent.\n', cov_imep_gross);

figure(6), stem(nth);
xlabel("Cycle Number [-]"), ylabel("nth [%]"), title("Thermal efficiency per cycle");
ylim([35, 42]);

mean_nth = mean(nth);
fprintf('Mean thermal efficiency = %4.2f percent.\n', mean_nth);

%% NHRR, CHR, MFB, CA50
% 2 stages heat release from 50 CAD bTDC to 100 CAD aTDC (130:280)
n_ratio = n_mean_c1 / (n_mean_c1 - 1);
n_ratio2 = 1 / (n_mean_c1 - 1);
for i_n =1:100
    d_P(i_n,:) = diff(P_cyc_1(i_n,:));
    NHRR_kJ = n_ratio .* P_cyc_1(i_n, 130:280) .* dV(130:280) ./ d_theta(130:280) + n_ratio2 .* V_s(130:280) .* d_P(i_n,130:280) ./ d_theta(130:280);
    NHRR = NHRR_kJ .* 1000;
    figure(7), plot(CAD_P_1(1,130:280),NHRR);hold on
    xlabel('Crank Angle Degree [aTDC]'),ylabel('NHRR [J/CAD]'),title('Net Heat Release Rate (NHRR) for 100 engine cycles')

    CHR = cumtrapz(1, NHRR);
    figure(8), plot(CAD_P_1(1,130:280),CHR);hold on
    xlabel('Crank Angle Degree [aTDC]'),ylabel('CHR [J]'),title('Cumulative Heat Release (NHRR) for 100 engine cycles')
    
    MFB = CHR / max(CHR) * 100;
    figure(9), plot(CAD_P_1(1,130:280),MFB);hold on
    xlabel('Crank Angle Degree [aTDC]'),ylabel('MFB [%]'),title('Fuel Mass Fraction Burnt for 100 engine cycles')
    ylim([0, 100]);

    MFB50 = abs(MFB - 50);
    CA = find(MFB50 == min(MFB50));
    CA50(i_n) = CAD_P_1(1, 130+CA);
end

mu_imep = mean(CA50);
fprintf('Mean CA50 = %4.2f.\n', mu_imep);
sigma_imep = std(CA50);
fprintf('Stdev CA50 = %4.2f.\n', cov_imep);

figure(10), stem(CA50);
xlabel('Cycle Number [-]'),ylabel('CA50 [CAD aTDC]'),title('Crank Angle for 50% fuel burnt (heat released) for 100 engine cycles')
ylim([0, 20]);