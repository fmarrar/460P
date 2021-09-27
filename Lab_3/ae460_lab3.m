% ae460_lab3.m ... Code for Lab 3

% clean up
clear
close all

set(0,'defaultaxesfontname','times')
set(0,'defaultaxesfontsize',10)
set(0,'defaultFigurePosition',[360 198 450 350])

%% Initialize
file = 'data3.txt';
import = importdata(file);
% enumerate variables
enum = split(import.textdata{1});
% initialize variables in workspace
self = init(enum(2:end), import.data);
self.p_t = self.p_t*144;
R = 1716.554;
gamma = 1.4;
p_amb = 14.93*144;
T_amb = 75.2 + 459.67;
A = (0.1/12)^2*pi;

p_e = table;
rho_e = table;
V_e = table;

%% Jet Exit Analyses
% Bernoulli
p_e.Bernoulli = p_amb*ones(size(self.p_t));
rho_e.Bernoulli = p_e.Bernoulli/R/T_amb;
V_e.Bernoulli = sqrt(2*(self.p_t - p_amb)./rho_e.Bernoulli);

% full isentropic
p_e.isentropic = p_amb*ones(size(self.p_t));
M_e = sqrt(2/(gamma-1)*((self.p_t/p_amb).^((gamma-1)/gamma) - 1));
T_e = T_amb./(1 + (gamma-1)/2*M_e.^2);
V_e.isentropic = M_e.*sqrt(gamma*R*T_e);
rho_e.isentropic = p_e.isentropic/R./T_e;

% choked isentropic
M_1 = M_e;
M_1(M_1 > 1) = 1;
T_e = T_amb./(1 + (gamma-1)/2*M_1.^2);
V_e.choked = M_1.*sqrt(gamma*R*T_e);
p_e.choked = self.p_t./(1 + (gamma-1)/2*M_1.^2).^(gamma/(gamma-1));
rho_e.choked = p_e.choked/R./T_e;

Thrust = array2table((p_e{:,:}-p_amb+V_e{:,:}.^2.*rho_e{:,:})*A);

%% imp2metric
V_e = varfun(@(x) x*0.3048, V_e);
rho_e = varfun(@(x) x*515.37882, rho_e);
Thrust = varfun(@(x) x*4.4482216, Thrust);
% p_e = varfun(@(x) x*47.88, p_e);
V_e.Properties.VariableNames = {'Bernoulli','isentropic','choked'};
rho_e.Properties.VariableNames = {'Bernoulli','isentropic','choked'};
Thrust.Properties.VariableNames = {'Bernoulli','isentropic','choked'};
% p_e.Properties.VariableNames = {'Bernoulli','isentropic','choked'};
self.Load = self.Load*4.4482216;


%%
figure(1), clf
grid on
hold on
plot(self.p_t, M_e, 'LineStyle', '-', 'LineWidth', 0.75,...
    'Marker', '^', 'MarkerSize', 5,...
    'Color', 'blue', 'MarkerFaceColor', 'blue')

xlim([self.p_t(1), self.p_t(end)])
xlabel('Chamber {\itp_t} (Pa)', 'FontSize', 12)
ylabel('\itM_e', 'FontSize', 12)
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02,.02],...
    'XMinorTick', 'on', 'YMinorTick', 'on')
set(get(gca,'ylabel'), 'rotation', 0,...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right')
ax = gca;
ax.XAxis.LineWidth = 0.75;
ax.YAxis.LineWidth = 0.75;

exportgraphics(gcf,'Mach_isen.png','Resolution',300)

plot(self.p_t(M_1>0.9), M_1(M_1>0.9), 'LineStyle', '--',...
    'Color', 'blue', 'LineWidth', 0.75, 'Marker', '^', 'MarkerSize', 5)
hold off

exportgraphics(gcf,'Mach_choked.png','Resolution',300)

figure(2), clf
hold on
plot(self.p_t, V_e.Bernoulli, 'LineStyle', '-',...
    'LineWidth', 0.75, 'Marker', '^', 'MarkerSize', 5,...
    'Color', [0.1412, 0.5490, 0.0392],...
    'MarkerFaceColor', [0.1412, 0.5490, 0.0392])
plot(self.p_t, V_e.isentropic, 'LineStyle', '-',...
    'LineWidth', 0.75, 'Marker', 'o', 'MarkerSize', 5,...
    'Color', 'blue', 'MarkerFaceColor', 'blue')
plot(self.p_t, V_e.choked, 'LineStyle', '-',...
    'LineWidth', 0.75, 'Marker', 's', 'MarkerSize', 5,...
    'Color', [0.9490, 0.4431, 0],...
    'MarkerFaceColor', [0.9490, 0.4431, 0])
hold off
grid on
xlim([self.p_t(1), self.p_t(end)])
xlabel('Chamber {\itp_t} (Pa)', 'FontSize', 12)
ylabel('{\itV_e} (m/s)', 'FontSize', 12)
legend('Bernoulli', 'Isentropic', 'Isentropic + Choked' ,...
    'Location', 'NorthWest')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02,.02],...
    'XMinorTick', 'on', 'YMinorTick', 'on')
ax = gca;
ax.XAxis.LineWidth = 0.75;
ax.YAxis.LineWidth = 0.75;

exportgraphics(gcf,'Velocity.png','Resolution',300)

figure(3), clf
hold on
plot(self.p_t, self.Load, 'LineStyle', '-',...
    'LineWidth', 0.75, 'Marker', '^', 'MarkerSize', 5,...
    'Color', [0.1412, 0.5490, 0.0392],...
    'MarkerFaceColor', [0.1412, 0.5490, 0.0392])
plot(self.p_t, Thrust.Bernoulli, 'LineStyle', '--',...
    'LineWidth', 0.75, 'Marker', '^', 'MarkerSize', 5,...
    'Color', 'blue')
plot(self.p_t, Thrust.isentropic, 'LineStyle', '--',...
    'LineWidth', 0.75, 'Marker', 'o', 'MarkerSize', 5,...
    'Color', [0.9490, 0.4431, 0])
plot(self.p_t, Thrust.choked, 'LineStyle', '--',...
    'LineWidth', 0.75, 'Marker', 's', 'MarkerSize', 5,...
    'Color', [0.6902, 0.0431, 0.4118])
hold off
grid on
xlim([self.p_t(1), self.p_t(end)])
xlabel('Chamber {\itp_t} (Pa)', 'FontSize', 12)
ylabel('Thrust (N)', 'FontSize', 12)
legend('Experiment', 'Bernoulli', 'Isentropic', 'Isentropic + Choked' ,...
    'Location', 'NorthWest')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02,.02],...
    'XMinorTick', 'on', 'YMinorTick', 'on')
ax = gca;
ax.XAxis.LineWidth = 0.75;
ax.YAxis.LineWidth = 0.75;

exportgraphics(gcf,'Thrust.png','Resolution',300)

%% Functions
function self = init(name, value)
    for n = 1:numel(name)
        self.(name{n}) = value(:,n);
    end
end