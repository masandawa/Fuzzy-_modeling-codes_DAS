%% =====================================================================
%  Figure 3: alpha-cut decomposition of a triangular fuzzy number
%  Manuscript: "A Fuzzy Mathematical Framework for Integrated Control
%  of Soil-Transmitted Helminth and Malaria Co-infection under Parameter
%  Uncertainty".
%  Standalone script: running it produces one figure (Fig03_AlphaCut.png).
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

%% ---------------- Alpha-cuts of the fuzzy transmission rate beta_vh ----------------
v = TFN.bvh; c = [0.00 0.35 0.70];
figure('Color','w','Units','centimeters','Position',[2 2 20 12]); hold on; box on; grid on;
fill([v(1) v(2) v(3)],[0 1 0],c,'FaceAlpha',0.12,'EdgeColor','none');
plot([v(1) v(2) v(3)],[0 1 0],'-','Color',c,'LineWidth',2.2);
for a = alphas
    lo = cutL(v,a); hi = cutU(v,a);
    plot([lo hi],[a a],'-','Color',[0.80 0.15 0.15],'LineWidth',3);
    plot([lo hi],[a a],'o','Color',[0.80 0.15 0.15],'MarkerFaceColor','w','MarkerSize',6,'LineWidth',1.4);
    text(v(3)+0.03, a, sprintf('[%.2f, %.2f]',lo,hi),'FontSize',10,'VerticalAlignment','middle');
end
xlim([v(1)-0.05, v(3)+0.17]); ylim([-0.05 1.1]);
set(gca,'FontSize',11,'YTick',alphas);
xlabel('\beta_{vh} (day^{-1})'); ylabel('Membership level \alpha');
print(gcf,'Fig03_AlphaCut','-dpng','-r300'); fprintf('Figure saved as Fig03_AlphaCut.png\n');
