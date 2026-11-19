%% =====================================================================
%  Figure 4: normalised forward sensitivity indices of R0m and R0w
%  Manuscript: "A Fuzzy Mathematical Framework for Integrated Control
%  of Soil-Transmitted Helminth and Malaria Co-infection under Parameter
%  Uncertainty".
%  Standalone script: running it produces one figure (Fig04_Sensitivity.png).
%  Bars are drawn on a signed logarithmic axis so that small indices
%  (e.g. d_m) are visible; tick labels give the true index values.
%  Tested in MATLAB-compatible syntax (also runs in GNU Octave >= 6).
%% =====================================================================
clear; close all; clc;

%% ---------------- Model parameters (Table 1) ----------------
Ncen  = 421119;      % Magu District Council population, 2022 census (NBS 2022)
P.mu  = 1/(67.2*365); % natural mortality: life expectancy 67.2 years (AfDB 2025)
P.Lam = P.mu*Ncen;   % human recruitment rate (persons day^-1), = 17.17
P.muv = 0.1429;      % mosquito mortality rate (day^-1)
P.bvh = 0.8;         % vector-to-human effective transmission rate (day^-1)
P.bhv = 0.8;         % human-to-vector effective transmission rate (day^-1)
P.bsw = 4.6e-7;      % environmental STH transmission rate (larva^-1 day^-1)
P.am  = 0.038;       % malaria recovery rate (day^-1)
P.als = 1/28;        % STH recovery (treatment) rate (day^-1)
P.asm = 0.005;       % simultaneous clearance rate of co-infection (day^-1)
P.dm  = 8.66e-5;     % malaria-induced death rate (day^-1)
P.sm  = 0.005;       % waning rate of malaria immunity (day^-1)
P.thw = 0.09;        % egg shedding rate (eggs person^-1 day^-1, scaled)
P.pse = 0.033;       % egg mortality rate (day^-1)
P.phe = 0.1;         % egg-to-L3 maturation rate (day^-1)
P.pha = 0.033;       % L3 larval mortality rate (day^-1)
P.ka  = 1.5;         % co-infection susceptibility modifier (-)
P.eta = 1.2;         % co-infection infectiousness modifier (-)
P.eps = 0.8;         % co-infection treatment-efficacy factor (-)
P.Nv  = 5e4;         % constant adult female mosquito population
Nst   = P.Lam/P.mu;  % disease-free human population

%% ---------------- Triangular fuzzy numbers (Table 2) ----------------
TFN.bvh = [0.6    0.8    1.0   ];
TFN.bsw = [2.5e-7 4.6e-7 7.0e-7];
TFN.am  = [1/42   0.038  1/14  ];
TFN.als = [1/42   1/28   1/14  ];
alphas  = [0 0.25 0.5 0.75 1];
cutL = @(v,a) v(1) + a*(v(2)-v(1));     % lower end of alpha-cut
cutU = @(v,a) v(3) - a*(v(3)-v(2));     % upper end of alpha-cut

%% ---------------- Reproduction numbers (Sect. 2.5) ----------------
R0m = @(bvh,am,u1,km) sqrt(bvh*(1-u1)*P.bhv*(1-u1)*P.Nv*P.mu ./ ...
                          (P.Lam*P.muv*(P.mu+km*am+P.dm)));
R0w = @(bsw,als,u2,ks) sqrt(bsw*(1-u2)*P.thw*P.phe*P.Lam ./ ...
                          (P.mu*P.pha*(P.pse+P.phe)*(P.mu+ks*als)));

%% ---------------- Analytical sensitivity indices (Table 3) ----------------
k1 = P.mu+P.am+P.dm; k2 = P.mu+P.als; k5 = P.pse+P.phe;
namesM = {'\beta_{vh}','\beta_{hv}','N_v','\mu_h','\Lambda_h','\mu_v','\alpha_m','d_m'};
valsM  = [0.5, 0.5, 0.5, 0.5*(1-P.mu/k1), -0.5, -0.5, -0.5*P.am/k1, -0.5*P.dm/k1];
namesW = {'\beta_{sw}','\theta_w','\Lambda_h','\phi_e','\psi_e','\phi_a','\alpha_s','\mu_h'};
valsW  = [0.5, 0.5, 0.5, 0.5*P.pse/k5, -0.5*P.pse/k5, -0.5, -0.5*P.als/k2, -0.5*(1+P.mu/k2)];
fprintf('R0m indices:\n'); for i=1:numel(valsM), fprintf('  %-10s % .4f\n',namesM{i},valsM(i)); end
fprintf('R0w indices:\n'); for i=1:numel(valsW), fprintf('  %-10s % .4f\n',namesW{i},valsW(i)); end
%% ---------------- Plot on a signed logarithmic scale ----------------
% d_m has index -0.0011, too small to see on a linear axis. Bar heights use
% the transform T(v) = sign(v)*log10(1+|v|/c), which is linear near zero and
% logarithmic for larger |v|; tick labels show the true values.
c    = 1e-3;                                   % linear-to-log transition
T    = @(v) sign(v).*log10(1+abs(v)/c);
tkv  = [-0.5 -0.1 -0.01 -0.001 0 0.001 0.01 0.1 0.5];   % true tick values
tkl  = {'-0.5','-0.1','-0.01','-0.001','0','0.001','0.01','0.1','0.5'};
ymax = T(0.9);
figure('Color','w','Units','centimeters','Position',[2 2 26 10]);
dat = {valsM, valsW}; nm = {namesM, namesW};
ttl = {'(a) Sensitivity of \it{R}\rm_{0m}','(b) Sensitivity of \it{R}\rm_{0w}'};
for pnl = 1:2
    subplot(1,2,pnl); hold on; box on;
    v = dat{pnl};
    for k = 1:numel(tkv)                       % horizontal grid at true ticks
        plot([0 numel(v)+1],[T(tkv(k)) T(tkv(k))],'-','Color',[0.88 0.88 0.88]);
    end
    for i = 1:numel(v)
        if v(i) >= 0, col = [0.00 0.35 0.70]; 
        else,         col = [0.80 0.15 0.15];  end
        h = T(v(i));
        fill([i-0.32 i+0.32 i+0.32 i-0.32],[0 0 h h], col, 'EdgeColor', 'none');
    end
    plot([0 numel(v)+1],[0 0],'k-');
    set(gca,'XTick',1:numel(v),'XTickLabel',nm{pnl},'FontSize',10, ...
        'YTick',T(tkv),'YTickLabel',tkl,'Layer','top');
    xlim([0.3 numel(v)+0.7]); ylim([-ymax ymax]);
    ylabel('Sensitivity index \Upsilon_p (signed log scale)');
    title(ttl{pnl},'FontWeight','normal');
end
print(gcf,'Fig04_Sensitivity','-dpng','-r300'); fprintf('Figure saved as Fig04_Sensitivity.png\n');
