%% =====================================================================
%  Figure 2: membership functions of the triangular fuzzy parameters (Table 2)
%  Manuscript: "A Fuzzy Mathematical Framework for Integrated Control
%  of Soil-Transmitted Helminth and Malaria Co-infection under Parameter
%  Uncertainty".
%  Standalone script: running it produces one figure (Fig02_Membership.png).
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

%% ---------------- Plot the four membership functions ----------------
names = {'\beta_{vh} (day^{-1})','\beta_{sw} (larva^{-1} day^{-1})','\alpha_m (day^{-1})','\alpha_s (day^{-1})'};
V = {TFN.bvh, TFN.bsw, TFN.am, TFN.als};
ttl = {'(a) \beta_{vh}','(b) \beta_{sw}','(c) \alpha_m','(d) \alpha_s'};
cols = [0.80 0.15 0.15; 0.10 0.55 0.25; 0.00 0.35 0.70; 0.90 0.50 0.05];
figure('Color','w','Units','centimeters','Position',[2 2 26 17]);
for k = 1:4
    subplot(2,2,k); hold on; box on; grid on;
    v = V{k}; c = cols(k,:);
    fill([v(1) v(2) v(3)],[0 1 0],c,'FaceAlpha',0.18,'EdgeColor','none');
    plot([v(1) v(2) v(3)],[0 1 0],'-','Color',c,'LineWidth',2.4);
    plot(v,[0 1 0],'o','Color',c,'MarkerFaceColor','w','MarkerSize',7,'LineWidth',1.5);
    lab = {'p^l','p^m','p^u'}; va = {'top','bottom','top'}; dy = [-0.06 0.05 -0.06];
    for j = 1:3
        yv = [0 1 0]; text(v(j), yv(j)+dy(j), sprintf('%s = %.3g',lab{j},v(j)), 'HorizontalAlignment','center','FontSize',10);
    end
    d = v(3)-v(1); xlim([v(1)-0.15*d, v(3)+0.15*d]); ylim([-0.15 1.15]);
    set(gca,'FontSize',11,'YTick',0:0.25:1);
    xlabel(names{k}); ylabel('Membership \mu(x)'); title(ttl{k},'FontWeight','normal');
end
print(gcf,'Fig02_Membership','-dpng','-r300'); fprintf('Figure saved as Fig02_Membership.png\n');
