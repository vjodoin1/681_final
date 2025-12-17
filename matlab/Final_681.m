%% ================================================================
%  Cancer Benefit Application Model — Decisions + Bifurcations
% ================================================================
function a681_final(G_t, m)

%% ========================================================================
outdir = sprintf('figs/G%d_m%d', G_t, m);
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

fig_counter = 0;  




%% ========================================================================
%Economic aspect
B     = 1.0;    
rho   = 0.03;   
L     = 2;      
Dlab  = 12;     

%Heterogeneity
eta_sd = 0.15;  
%m      = 2;     

%Social dynamics
alpha = 0.15;   
delta = 0.05;   
%G_t   = 1;      
N0    = 0.6;    
N_min = 0; 
N_max = 1;

k_mean = 0.30;

p_mean = 0.3*(8/(8+2)) + 0.4*(5/(5+3)) + 0.3*(2/(2+5));

present_value = @(p,L_,rho_) ((p./(1+rho_)).^L_) ./ (1 - (p./(1+rho_)));
finite_pv     = @(p,D_,rho_) (1 - (p./(1+rho_)).^D_) ./ (1 - (p./(1+rho_)));

PV_B    = present_value(p_mean, L, rho);
PV_stig = finite_pv(p_mean, Dlab, rho);
A   = B.*PV_B - k_mean;   
C   = PV_stig;            
eta = eta_sd;             

Phi = @(z) 0.5.*(1 + erf(z./sqrt(2)));
phi = @(z) (1./sqrt(2*pi)).*exp(-0.5.*z.^2);



%% ========================================================================
sigma_values = [0.30 0.60 0.90];
T = 40; 
N_sim = 1000;
colors = lines(numel(sigma_values));

figure('Units','inches','Position',[1 1 6 4]); hold on; box on; grid on;

for jj = 1:numel(sigma_values)
    sig = sigma_values(jj);
    color_base = colors(jj,:);

    a_runs = zeros(T, N_sim);

    for n = 1:N_sim
        N_series = zeros(T+1,1); 
        N_series(1) = N0;

        for t = 1:T
            N_t = N_series(t);

            eta_draw = eta * (1 + 0.05*randn);
            a_t = Phi((A - sig*C*(1 - N_t)^m) / eta_draw);

            N_next = N_t + alpha*(a_t - N_t) + delta*G_t*(1 - N_t);
            N_series(t+1) = max(N_min, min(N_max, N_next));

            a_runs(t,n) = a_t;
        end
    end

    a_runs_storage{jj} = a_runs;

    a_mean = mean(a_runs,2);
    a_p10  = prctile(a_runs,10,2);
    a_p90  = prctile(a_runs,90,2);

    fill([1:T, fliplr(1:T)], [a_p10', fliplr(a_p90')], ...
        color_base, 'FaceAlpha',0.10, 'EdgeColor','none');

    plot(1:T, a_mean, 'Color', color_base, 'LineWidth',2.5, ...
         'DisplayName', sprintf('\\sigma=%.2f', sig));
end

xlabel('Time (t)');
ylabel('Average application rate a_t');
title('Application decisions for different stigma levels');
legend('Location','southeast');

fig_counter = fig_counter + 1;
exportgraphics(gcf, fullfile(outdir, ...
    sprintf('fig%d.pdf', fig_counter)), ...
    'ContentType','vector');


%% ========================================================================
figure('Units','inches','Position',[1 1 6 4]); hold on; box on; grid on;

for jj = 1:numel(sigma_values)
    color_base = colors(jj,:);

    a_runs = a_runs_storage{jj};

    a_mean = mean(a_runs,2);
    a_p10  = prctile(a_runs,10,2);
    a_p90  = prctile(a_runs,90,2);

    fill([1:T, fliplr(1:T)], [a_p10', fliplr(a_p90')], ...
        color_base, 'FaceAlpha', 0.10, 'EdgeColor', 'none');

    plot(1:T, a_mean, 'Color', color_base, 'LineWidth', 2.5, ...
         'DisplayName', sprintf('\\sigma=%.2f', sigma_values(jj)));
end

xlim([0 10]);
xlabel('Time (t)');
ylabel('Average application rate a_t');
title('Early convergence (t ≤ 10)');
legend('Location','southeast');

fig_counter = fig_counter + 1;
exportgraphics(gcf, fullfile(outdir, sprintf('fig%d.pdf', fig_counter)), ...
    'ContentType','vector');

%% ========================================================================
sigmas_scan = linspace(0.2, 1.0, 80);
Ngrid = linspace(0,1,1000);

stable_pts   = [];
unstable_pts = [];

for sig = sigmas_scan
    aN = @(N) Phi((A - sig.*C.*(1 - N).^m)./eta);
    fN = @(N) aN(N) - N;

    rootsN = find_roots_unit(fN, Ngrid);

    for Nstar = rootsN(:).'
        z = (A - sig.*C.*(1 - Nstar).^m)./eta;
        aprime = phi(z) .* (sig.*C.*m.*(1 - Nstar).^(m-1)) ./ eta;

        if abs(aprime) < 1
            stable_pts   = [stable_pts;   sig, Nstar];
        else
            unstable_pts = [unstable_pts; sig, Nstar];
        end
    end
end

figure('Color','w','Units','inches','Position',[1 1 6 4]); hold on; grid on; box on;
if ~isempty(stable_pts)
    scatter(stable_pts(:,1), stable_pts(:,2), 25, 'b', 'filled', 'DisplayName','Stable');
end
if ~isempty(unstable_pts)
    scatter(unstable_pts(:,1), unstable_pts(:,2), 25, 'r', 'o', 'DisplayName','Unstable');
end
xlabel('\sigma_{stigma}');
ylabel('Fixed Point N^*');
title('Bifurcation Diagram');
legend('Location','best');

fig_counter = fig_counter + 1;
exportgraphics(gcf, fullfile(outdir, sprintf('fig%d.pdf', fig_counter)), ...
    'ContentType','vector');

%% ========================================================================
beta_over_alpha  = [0.25 0.5 0.75 1];
delta_over_alpha = [0 0.25 0.5 0.75];

y_vals = linspace(-1.0, 1.5, 75);

sigma_fix = 0.60;

figure('Units','normalized','Position',[0.06 0.06 0.88 0.86]);
for r = 1:numel(delta_over_alpha)
  for c = 1:numel(beta_over_alpha)
    subplot(numel(delta_over_alpha), numel(beta_over_alpha), ...
            (r-1)*numel(beta_over_alpha)+c);
    hold on; box on; grid on;

    beta  = beta_over_alpha(c)  .* alpha;
    delt  = delta_over_alpha(r) .* alpha;

    for yi = 1:3:numel(y_vals)
        y = y_vals(yi);
        Atemp = y .* eta;

        aN = @(N) Phi((Atemp - sigma_fix.*C.*(1 - N).^m + beta.*N)./eta);
        fN = @(N) aN(N) + (delt./alpha) - N;

        rootsN = find_roots_unit(fN, Ngrid);

        for Nstar = rootsN(:).'
            z = (Atemp - sigma_fix.*C.*(1 - Nstar).^m + beta.*Nstar)./eta;
            aprime = phi(z) .* ((sigma_fix.*C.*m.*(1 - Nstar).^(m-1) + beta) ./ eta);

            if abs(aprime) < 1
                plot(y, Nstar, 'kd', 'MarkerFaceColor','k', 'MarkerSize', 6);
            else
                plot(y, Nstar, 'kd', 'MarkerFaceColor','w', ...
                     'MarkerSize', 4, 'LineWidth', 1.2);
            end
        end
    end

    xlim([min(y_vals) max(y_vals)]); ylim([0 1]);
    if r==numel(delta_over_alpha), xlabel('y = A/\eta'); end
    if c==1, ylabel('N^*'); end
    title(sprintf('\\beta/\\alpha = %.2g, \\delta/\\alpha = %.2g', ...
    beta_over_alpha(c), delta_over_alpha(r)));

  end
end
sg = sgtitle('Equilibria Across Peer (\beta) and Authority (\delta) Forces', ...
    'Color', 'k', ...           % solid black title
    'FontWeight', 'bold', ...
    'FontSize', 16, ...
    'Interpreter','tex');

uistack(sg, 'top');  % ensure it is drawn above subplots

fig_counter = fig_counter + 1;
exportgraphics(gcf, fullfile(outdir, sprintf('fig%d.pdf', fig_counter)), ...
    'ContentType','vector');

%% ========================================================================
figure('Color','w','Units','inches','Position',[1 1 6.5 5]); 
hold on; grid on; box on;

sigma_fp = 0.60;

N_nc = linspace(0,1,200);
a_nc = arrayfun(@(N) Phi((A - sigma_fp*C*(1 - N)^m)/eta), N_nc);

plot(N_nc, a_nc, 'b-', 'LineWidth', 2.2, 'DisplayName','a(N) nullcline');
plot(N_nc, N_nc, 'r--', 'LineWidth', 2.2, 'DisplayName','N = a line');

N0_list = [0.15, 0.45, 0.75];
Tsim = 30;

for N0_init = N0_list
    N_path = zeros(Tsim+1,1);
    a_path = zeros(Tsim+1,1);

    N_path(1) = N0_init;

    for t = 1:Tsim
        N = N_path(t);
        a = Phi((A - sigma_fp*C*(1 - N)^m)/eta);
        a_path(t) = a;

        N_next = N + alpha*(a - N) + delta*G_t*(1 - N);
        N_path(t+1) = max(0, min(1, N_next));
    end

    plot(N_path(1:Tsim), a_path(1:Tsim), '-o', 'LineWidth',1.5, ...
        'DisplayName', sprintf('Trajectory N_0=%.2f', N0_init));
end

xlabel('Social Norm N_t');
ylabel('Application Probability a_t');
title('Phase Plane (N_t, a_t): Nullclines and Sample Trajectories');
legend('Location','southeast');

fig_counter = fig_counter + 1;
exportgraphics(gcf, fullfile(outdir, sprintf('fig%d.pdf', fig_counter)), ...
    'ContentType','vector');

%% ========================================================================
function rootsN = find_roots_unit(fN, Ngrid)
    vals = fN(Ngrid);
    vals(~isfinite(vals)) = NaN;
    epsv = 1e-10; 
    vals(abs(vals) < epsv) = 0;
    sgn = sign(vals); 
    sgn(~isfinite(sgn)) = 0;

    idx = find(sgn(1:end-1) .* sgn(2:end) < 0);

    rootsN = [];
    for k = idx(:).'
        a = Ngrid(k); 
        b = Ngrid(k+1);

        try
            r0 = fzero(fN, [a b]);
        catch
            mid = 0.5*(a+b);
            try
                r0 = fzero(fN, mid);
            catch
                continue;
            end
        end

        if isfinite(r0) && r0>=0 && r0<=1
            rootsN(end+1,1) = r0;
        end
    end
    rootsN = unique(round(rootsN, 8));
end
end