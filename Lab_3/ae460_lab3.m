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

%%
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
M_e(M_e > 1) = 1;
T_e = T_amb./(1 + (gamma-1)/2*M_e.^2);
V_e.choked = M_e.*sqrt(gamma*R*T_e);
p_e.choked = self.p_t./(1 + (gamma-1)/2*M_e.^2).^(gamma/(gamma-1));
rho_e.choked = p_e.choked/R./T_e;

Thrust = array2table((p_e{:,:}-p_amb+V_e{:,:}.^2.*rho_e{:,:})*A);

psf2Pa = @(x) x*47.88;
degR2K = @(x) x*5/9;

V_e = varfun(@(x) x*0.3048, V_e);
rho_e = varfun(@(x) x*515.37882, rho_e);
Thrust = varfun(@(x) x*4.4482216, Thrust);
V_e.Properties.VariableNames = {'Bernoulli','isentropic','choked'};
rho_e.Properties.VariableNames = {'Bernoulli','isentropic','choked'};
Thrust.Properties.VariableNames = {'Bernoulli','isentropic','choked'};
self.Load = self.Load*4.4482216;


%%
% figure(1), clf
% hold on
% plot(self.p_t, self.Load, 'LineStyle', 'none', 'LineWidth', 0.75,...
%     'Marker', '^', 'MarkerSize', 5)
% hold off
% grid on
% set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02,.02],...
%     'XMinorTick', 'on', 'YMinorTick', 'on')
% ax = gca;
% ax.XAxis.LineWidth = 0.75;
% ax.YAxis.LineWidth = 0.75;



function self = init(name, value)
    for n = 1:numel(name)
        self.(name{n}) = value(:,n);
    end
end