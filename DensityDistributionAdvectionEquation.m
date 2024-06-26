%% Solution of constant growth based on an advection equation for the density distribution
% Author = Plath, Timo
% E-mail: t.plath@utwente.nl
% Version = 1.0

% Solves a constant advection equation for the density distribution, plots
% its first 6 moments over time and compares it to their analytical solution.

% INPUT: tmax           a double determining the maximum simulation time
%         dt            a double determining the simulation timestep
%         g             a double determining the advection speed (growth
%                       rate)

% OUTPUT: Moments      a matrix containing the first 6 moments of each
%                      timestep

% Run for example: DensityDistributionAdvectionEquation(100,1,1)

function Moments = DensityDistributionAdvectionEquation(tmax,dt,g)

CDFv_L = ReadData('Data/InitialLactoseMCCPVP_CDF-v_L.csv',1e6,0,0.995);
v_L = convertCDFtoPDF(CDFv_L);
% initial distribution function
x = v_L(:,1);
y = v_L(:,2);

n = length(x);
N_delta = 3;
mMax = 2*N_delta;

%% Analytic
% Number of time steps
nt=round(tmax/dt)+1;
% time
time = (1:nt)*dt;

% Compute moments of the VPSD
M = ComputeMoments(v_L(:,1),v_L(:,2),mMax);
% we compute weights and nodes using the Wheeler algorithm
[xi,w] = Wheeler(M(1:2*N_delta),N_delta);

% Realize a moving-grid for the nodes and weights and the distribution n_L
for t = 1:nt
    % Moments
    m0(1,t) = dot(w,xi.^0);
    m1(1,t) = dot(w,xi.^1);
    m2(1,t) = dot(w,xi.^2);
    m3(1,t) = dot(w,xi.^3);
    m4(1,t) = dot(w,xi.^4);
    m5(1,t) = dot(w,xi.^5);

    for i = 1:N_delta
        xi(i,1) = xi(i,1) + g*dt;
    end
    for i = 1:n
        x(i,1) = x(i,1) + g*dt;
    end
    % plot moving distribution
    figure(1);
    plot (x,y)
    hold off
    title('Moving grid')
    xlabel('Granule Size [µm]')
    ylabel('Volume density [%]')
    grid on
    legend('Moving grid','initial','Location','north west')
    drawnow
end

Moments = [m0;m1;m2;m3;m4;m5];
Moments = Moments';

stepsize = 5;
% plot moments
figure(2)
axes2 = axes('Parent', figure(2));
plot(time(1:stepsize:tmax+1),m0(1:stepsize:end)/max(m0),'x', 'MarkerSize', 14, 'LineWidth',2)
hold on
plot(time(1:stepsize:tmax+1),m1(1:stepsize:end)/max(m1),'x-', 'MarkerSize', 14, 'LineWidth',2)
plot(time(1:stepsize:tmax+1),m2(1:stepsize:end)/max(m2),'x-', 'MarkerSize', 14, 'LineWidth',2)
plot(time(1:stepsize:tmax+1),m3(1:stepsize:end)/max(m3),'x-', 'MarkerSize', 14, 'LineWidth',2)
plot(time(1:stepsize:tmax+1),m4(1:stepsize:end)/max(m4),'x-', 'MarkerSize', 14, 'LineWidth',2)
plot(time(1:stepsize:tmax+1),m5(1:stepsize:end)/max(m5),'x-', 'MarkerSize', 14, 'LineWidth',2)
xlabel('time[s]')
ylabel('Normalized moments')
legend('m0','m1','m2','m3','m4','Analytical', 'location','best')

% analytical solution to constant growth
syms M_0Growth(t_)
M_0Growth(t_) = Moments(1,1);
M_1Growth(t_) = Moments(1,2) + g(1,1)*Moments(1,1)*t_;
M_2Growth(t_) = Moments(1,3) + 2*g(1,1)*Moments(1,2)*t_ + g(1,1)^2*Moments(1,1)*t_.^2;
M_3Growth(t_) = Moments(1,4) + 3*g(1,1)*Moments(1,3)*t_ + 3*g(1,1)^2*Moments(1,2)*t_.^2 + g(1,1)^3*Moments(1,1)*t_.^3;
M_4Growth(t_) = Moments(1,5) + 4*g(1,1)*Moments(1,4)*t_ + 6*g(1,1)^2*Moments(1,3)*t_.^2 + g(1,1)^3*4*Moments(1,2)*t_.^3 + g(1,1)^4*Moments(1,1)*t_.^4;
M_5Growth(t_) = Moments(1,6) + 5*g(1,1)*Moments(1,5)*t_ + 10*g(1,1)^2*Moments(1,4)*t_.^2 + 2*g(1,1)^3*5*Moments(1,3)*t_.^3 + g(1,1)^4*5*Moments(1,2)*t_.^4 + g(1,1)^5*Moments(1,1)*t_.^5;

% Reset the color order to repeat the defined colors
set(axes2, 'ColorOrderIndex', 1);
fplot(subs(M_0Growth)/max(subs(M_0Growth([0,time(1:nt)*dt]))),[0 tmax], 'LineWidth', 2);
hold on
fplot(subs(M_1Growth)/max(subs(M_1Growth([0,time(1:nt)*dt]))),[0 tmax], 'LineWidth', 2)
fplot(subs(M_2Growth)/max(subs(M_2Growth([0,time(1:nt)*dt]))),[0 tmax], 'LineWidth', 2)
fplot(subs(M_3Growth)/max(subs(M_3Growth([0,time(1:nt)*dt]))),[0 tmax], 'LineWidth', 2)
fplot(subs(M_4Growth)/max(subs(M_4Growth([0,time(1:nt)*dt]))),[0 tmax], 'LineWidth', 2)
fplot(subs(M_5Growth)/max(subs(M_5Growth([0,time(1:nt)*dt]))),[0 tmax], 'LineWidth', 2)
legend("$\widetilde{M}_{0}^{g\prime}$","$\widetilde{M}_{1,\mathrm{V}}^{g\prime}$","$\widetilde{M}_{2,\mathrm{V}}^{g\prime}$","$\widetilde{M}_{3,\mathrm{V}}^{g\prime}$","$\widetilde{M}_{4,\mathrm{V}}^{g\prime}$","$\widetilde{M}_{5,\mathrm{V}}^{g\prime}$",'location','best', 'FontSize', 14,'Interpreter','Latex')

% Compute error
for i = 1:length(Moments)
    E_V(i,1) = abs(Moments(i,1) - M_0Growth((i-1)*dt))/M_0Growth((i-1)*dt);
    E_V(i,2) = abs(Moments(i,2) - M_1Growth((i-1)*dt))/M_1Growth((i-1)*dt);
    E_V(i,3) = abs(Moments(i,3) - M_2Growth((i-1)*dt))/M_2Growth((i-1)*dt);
    E_V(i,4) = abs(Moments(i,4) - M_3Growth((i-1)*dt))/M_3Growth((i-1)*dt);
    E_V(i,5) = abs(Moments(i,5) - M_4Growth((i-1)*dt))/M_4Growth((i-1)*dt);
    E_V(i,6) = abs(Moments(i,6) - M_5Growth((i-1)*dt))/M_5Growth((i-1)*dt);
end

% plot error
figure(5)
axes5 = axes('Parent', figure(5));
for i = 1:size(E_V,2)
    plot(time(1:stepsize:length(E_V))*dt,E_V(1:stepsize:end,i), "-",'MarkerSize', 14, 'LineWidth',2)
    hold on
end
xlabel('$t^{\prime}$','Interpreter','Latex')
ylabel('$E$ [\%]','Interpreter','Latex')
set(gcf, 'Color', 'w');
set(axes5,'FontSize', 16)