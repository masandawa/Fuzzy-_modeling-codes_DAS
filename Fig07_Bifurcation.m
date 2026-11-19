%% =====================================================================
%  Figure 7: endemic prevalence versus reproduction numbers (forward bifurcation)
%  Manuscript: "A Fuzzy Mathematical Framework for Integrated Control
%  of Soil-Transmitted Helminth and Malaria Co-infection under Parameter
%  Uncertainty".
%  Standalone script: running it produces one figure (Fig07_Bifurcation.png).
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

%% ---------------- Model right-hand side (model system, Sect. 2.2) ----------------
% y = [Yss Ysi Ysr Yis Yii Yir Iv Fw Aw]'; u1 = vector control, u2 = sanitation,
% km, ks = treatment-intensification factors for malaria and STH.
lm  = @(y,Q,u1) Q.bvh*(1-u1)*y(7)/sum(y(1:6));   % malaria force of infection
lw  = @(y,Q,u2) Q.bsw*(1-u2)*y(9);               % STH force of infection
rhs = @(t,y,Q,u1,u2,km,ks) [ ...
  Q.Lam - lm(y,Q,u1)*y(1) - lw(y,Q,u2)*y(1) - Q.mu*y(1) + Q.sm*y(3) + ks*Q.als*y(4); ...
  lm(y,Q,u1)*y(1) - (Q.mu+km*Q.am+Q.dm)*y(2) - Q.ka*lw(y,Q,u2)*y(2) + Q.eps*ks*Q.als*y(5); ...
  km*Q.am*y(2) - (Q.mu+Q.sm)*y(3) - lw(y,Q,u2)*y(3) + Q.asm*y(5) + ks*Q.als*y(6); ...
  lw(y,Q,u2)*y(1) - (Q.mu+ks*Q.als)*y(4) - Q.ka*lm(y,Q,u1)*y(4) + Q.sm*y(6); ...
  Q.ka*lm(y,Q,u1)*y(4) + Q.ka*lw(y,Q,u2)*y(2) - (Q.mu+Q.eps*ks*Q.als+Q.eps*km*Q.am+Q.asm+Q.dm)*y(5); ...
  lw(y,Q,u2)*y(3) + Q.eps*km*Q.am*y(5) - (Q.mu+ks*Q.als+Q.sm)*y(6); ...
  Q.bhv*(1-u1)*(y(2)+Q.eta*y(5))*(Q.Nv-y(7))/sum(y(1:6)) - Q.muv*y(7); ...
  Q.thw*(y(4)+Q.eta*y(5)+y(6)) - (Q.pse+Q.phe)*y(8); ...
  Q.phe*y(8) - Q.pha*y(9) ];
opts = odeset('RelTol',1e-8,'AbsTol',1e-3,'InitialStep',1e-2);
y0intro = [0.90*Nst; 0.04*Nst; 0; 0.05*Nst; 0.01*Nst; 0; 0.02*P.Nv; 0; 0];
prevM = @(Y) (Y(:,2)+Y(:,5))./sum(Y(:,1:6),2);          % malaria prevalence
prevW = @(Y) (Y(:,4)+Y(:,5)+Y(:,6))./sum(Y(:,1:6),2);   % STH prevalence

%% ---------------- Endemic (coexistence) equilibrium, modal parameters ----------------
[~,Yl] = ode15s(@(t,y) rhs(t,y,P,0,0,1,1), [0 150*365], y0intro, opts);
ye = Yl(end,:)';
for it = 1:8                               % Newton refinement
    f0 = rhs(0,ye,P,0,0,1,1); Jn = zeros(9);
    for j = 1:9
        h = 1e-6*max(1,abs(ye(j))); e = zeros(9,1); e(j) = h;
        Jn(:,j) = (rhs(0,ye+e,P,0,0,1,1) - f0)/h;
    end
    ye = ye - Jn\f0;
end

%% ---------------- Analytical boundary equilibria  ----------------
D  = P.dm/P.mu; Cc = 1 + P.am/(P.mu+P.sm); B = P.bhv/P.muv;
fprintf('Bifurcation coefficient sign: D - C - B = %.4f (negative => forward)\n', D-Cc-B);
Rg = linspace(0.2,6,300);
xm = zeros(size(Rg));
for i = 1:numel(Rg)
    r = Rg(i)^2; rt = roots([r*D*Cc, -(r*(D-Cc)-B), 1-r]);
    rt = rt(abs(imag(rt))<1e-12 & real(rt)>0); if isempty(rt), xm(i)=0; else, xm(i)=max(real(rt)); end
end
yw = max(0, 1 - 1./Rg.^2);           % STH-only endemic prevalence

%% ---------------- Co-infection invasion threshold  ----------------
k1 = P.mu+P.am+P.dm; k2 = P.mu+P.als; k3 = P.mu+P.eps*(P.als+P.am)+P.asm+P.dm;
Rw = R0w(P.bsw,P.als,0,1); s = 1/Rw^2; ii = 1-s; lamw = k2*(Rw^2-1);
Dl = (k1+P.ka*lamw)*k3 - P.eps*P.als*P.ka*lamw;
Phi = k1*(s*(k3+P.eta*P.ka*lamw) + P.ka*ii*(P.eps*P.als+P.eta*(k1+P.ka*lamw)))/Dl;
fprintf('Phi = %.4f, malaria invasion threshold in STH-endemic setting R0m = %.4f\n', Phi, 1/sqrt(Phi));

%% ---------------- Numerical coexistence branches (continuation in time) ----------------
Rm_base = R0m(P.bvh,P.am,0,1); Rw_base = R0w(P.bsw,P.als,0,1);
Rgrid = [0.5:0.05:1.2, 1.3:0.1:2, 2.25:0.25:6];
cm = zeros(size(Rgrid)); yc = ye;
for i = numel(Rgrid):-1:1                    % malaria prevalence, STH at baseline
    Q = P; Q.bvh = P.bvh*(Rgrid(i)/Rm_base)^2;
    yc(7) = max(yc(7),1); yc(2) = max(yc(2),1);
    [~,Y] = ode15s(@(t,y) rhs(t,y,Q,0,0,1,1), [0 300*365], yc, opts);
    yq = Y(end,:)';
    for it = 1:10                          % Newton refinement of the equilibrium
        f0 = rhs(0,yq,Q,0,0,1,1); Jn = zeros(9);
        for jj = 1:9
            hh = 1e-6*max(1,abs(yq(jj))); ee = zeros(9,1); ee(jj) = hh;
            Jn(:,jj) = (rhs(0,yq+ee,Q,0,0,1,1) - f0)/hh;
        end
        yq = yq - Jn\f0;
    end
    if all(yq > -1e-6) && max(abs(rhs(0,yq,Q,0,0,1,1))) < 1e-6, yc = yq; else, yc = Y(end,:)'; end
    cm(i) = max(0,prevM(yc'));
end
cw = zeros(size(Rgrid)); yc = ye;
for i = numel(Rgrid):-1:1                    % STH prevalence, malaria at baseline
    Q = P; Q.bsw = P.bsw*(Rgrid(i)/Rw_base)^2;
    yc(9) = max(yc(9),1); yc(4) = max(yc(4),1);
    [~,Y] = ode15s(@(t,y) rhs(t,y,Q,0,0,1,1), [0 300*365], yc, opts);
    yq = Y(end,:)';
    for it = 1:10                          % Newton refinement of the equilibrium
        f0 = rhs(0,yq,Q,0,0,1,1); Jn = zeros(9);
        for jj = 1:9
            hh = 1e-6*max(1,abs(yq(jj))); ee = zeros(9,1); ee(jj) = hh;
            Jn(:,jj) = (rhs(0,yq+ee,Q,0,0,1,1) - f0)/hh;
        end
        yq = yq - Jn\f0;
    end
    if all(yq > -1e-6) && max(abs(rhs(0,yq,Q,0,0,1,1))) < 1e-6, yc = yq; else, yc = Y(end,:)'; end
    cw(i) = max(0,prevW(yc'));
end

figure('Color','w','Units','centimeters','Position',[2 2 26 10.5]);
subplot(1,2,1); hold on; box on; grid on;
plot([0.2 1],[0 0],'k-','LineWidth',2); plot([1 6],[0 0],'k--','LineWidth',1.2);
h1 = plot(Rg(Rg>=1), 100*xm(Rg>=1), '-', 'Color',[0.80 0.15 0.15],'LineWidth',2);
h2 = plot(Rgrid, 100*cm, 'o-', 'Color',[0.00 0.35 0.70],'MarkerSize',4,'MarkerFaceColor',[0.00 0.35 0.70]);
plot([1 1],[0 16],'k:'); plot(1/sqrt(Phi)*[1 1],[0 16],':','Color',[0 0.35 0.7]);
xlabel('\it{R}\rm_{0m}'); ylabel('Equilibrium malaria prevalence (%)'); ylim([0 16]); xlim([0.2 6]);
title('(a) Malaria','FontWeight','normal');
legend([h1 h2],{'Malaria-only equilibrium (analytical)','Co-endemic equilibrium (R_{0w} = 3.00)'},'Location','southeast','FontSize',8);
subplot(1,2,2); hold on; box on; grid on;
plot([0.2 1],[0 0],'k-','LineWidth',2); plot([1 6],[0 0],'k--','LineWidth',1.2);
h1 = plot(Rg(Rg>=1), 100*yw(Rg>=1), '-', 'Color',[0.80 0.15 0.15],'LineWidth',2);
h2 = plot(Rgrid, 100*cw, 'o-', 'Color',[0.10 0.55 0.25],'MarkerSize',4,'MarkerFaceColor',[0.10 0.55 0.25]);
plot([1 1],[0 100],'k:');
xlabel('\it{R}\rm_{0w}'); ylabel('Equilibrium STH prevalence (%)'); ylim([0 100]); xlim([0.2 6]);
title('(b) STH','FontWeight','normal');
legend([h1 h2],{'STH-only equilibrium: 1-1/R_{0w}^2','Co-endemic equilibrium (R_{0m} = 4.15)'},'Location','southeast','FontSize',8);
fprintf('Co-endemic malaria prevalence at R0m = 0.80: %.4f\n', interp1(Rgrid,cm,0.8));
print(gcf,'Fig07_Bifurcation','-dpng','-r300'); fprintf('Figure saved as Fig07_Bifurcation.png\n');
