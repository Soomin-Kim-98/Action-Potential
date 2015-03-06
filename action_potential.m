function V = action_potential(initial, I, dt)
% initial: initial voltage
% I: applied current
% dt: time step size (ms)
% Example function call: action_potential(0, [zeros(1, 50), -1000*ones(1,
% 50), zeros(1, 900)], 10^-3);

% action_potential(.0227, V, 10^-2)

V = zeros(length(I), 4);
V(1,:) = [initial n_stable(initial) m_stable(initial) h_stable(initial)]';
time = (1:length(I))*dt;


% =========================================================================
% BEGIN SIMULATION
% =========================================================================
tic;
for ii = 1:length(I)-1
    if ii > 50
        V(ii, 1) = -50;
    end
    k1 = dt * fox(V(ii,:), I(ii));
    k2 = dt * fox(V(ii,:) + k1/2, I(ii));
    k3 = dt * fox(V(ii,:) + k2/2, I(ii));
    k4 = dt * fox(V(ii,:) + k3,  I(ii));
    V(ii+1,:) = V(ii,:) + 1/6*(k1 + 2*k2 + 2*k3 + k4);
end
record_time = toc;
disp(['time per step = ' num2str(record_time/(length(I)-1))]);

% =========================================================================
% Plot Things
% =========================================================================
% Plot rate constants as a function of voltage and steady-state rate
% constants as a function of voltage
% figure, plot(-100:1:50, alphan(-100:1:50), -100:1:50, betan(-100:1:50));
% xlabel('V(mV)'); ylabel('Rate Constant (msec-1)');
% legend('alphan', 'betan');
% figure, plot(-100:1:50, alpham(-100:1:50), -100:1:50, betam(-100:1:50));
% xlabel('V(mV)'); ylabel('Rate Constant (msec-1)');
% legend('alpham', 'betam');
% figure, plot(-100:1:50, alphah(-100:1:50), -100:1:50, betah(-100:1:50));
% xlabel('V(mV)'); ylabel('Rate Constant (msec-1)');
% legend('alphah', 'betah');
% figure, plot(-100:1:50, n_stable(-100:1:50));
% xlabel('V(mV)'); ylabel('n');
% figure, plot(-100:1:50, m_stable(-100:1:50));
% xlabel('V(mV)'); ylabel('m');
% figure, plot(-100:1:50, h_stable(-100:1:50));
% xlabel('V(mV)'); ylabel('h');

% Plot Voltage vs time and Current vs time on different axes to show step response
% figure, [ax, h1, h2] = plotyy(time, -V(:,1), time, -I);
% set(get(ax(1), 'Ylabel'), 'String', 'Voltage (mV)');
% set(get(ax(2), 'Ylabel'), 'String', 'Current (\muA)');

% Plot voltage vs time
% figure, plot(time, -V(:,1));
% xlabel('Time (ms)');
% ylabel('Voltage (mV)');

% Plot m, n, h, and voltage
% figure, hold on, plot(1:length(I), V(:,2)), title('n, time');
% plot(1:length(I), V(:,3)), title('m, time');
% plot(1:length(I), V(:,4)), title('h, time');
% hold off;
C = 1; % uF
GNa = 120; % mSiemens
GK = 36;
GL = .3;
EL = 10.613;
ENa = 115;
EK = -12;
% Plot Potassium and Sodium Conductance on one axis and Voltage on another
% figure;
% [ax, h1, h2] = plotyy(time, [GK*V(:,2).^4 GNa*(V(:,3).^3).*V(:,4)], time, -V(:,1));
% set(get(ax(1), 'Ylabel'), 'String', 'Conductance (mSiemens)');
% set(get(ax(2), 'Ylabel'), 'String', 'Voltage');
% legend('Potassium','Sodium');
% xlabel('Time (ms)');

% Plot Conductances individually
% figure, plot(time, GK*V(:,2).^4);
% xlabel('Time (ms)');
% ylabel('Conductance (mSiemens)');
% title('Potassium Conductance');
% figure, plot(time, GNa*(V(:,3).^3).*V(:,4));
% xlabel('Time (ms)');
% ylabel('Conductance (mSiemens)');
% title('Sodium Conductance');

% Plot current terms
% figure, plot(time, GK*V(:,2).^4.*(V(:,1)-EK), time, GNa*(V(:,3).^3).*V(:,4).*(V(:,1)-ENa), time, GL*(V(:,1)-EL));
% title('Ionic currents');
% legend('Potassium', 'Sodium', 'Leakage');
% xlabel('Time (ms)');
% ylabel('Current (\muA)');

end
% =========================================================================
% FUNCTION DEFINITIONS
% =========================================================================
function dV = fox(V_array, I)
% Governing equation
% V_array:
V = V_array(1); % V_array(1) is voltage
n = V_array(2); % V_array(2) is n
m = V_array(3); % V_array(3) is m
h = V_array(4); % V_array(4) is h
% I: Applied current

% Nmax: appears in the definition of the variance of the noise
Nmax = 100;
% C = membrane Capacitance
C = 1; % uF
% G = conductance of a single ion channel
GNa = 120; % mSiemens
GK = 36;
GL = .3;
% E = Electric field due to ion concentration
EL = 10.613;
ENa = 115;
EK = -12;

% rates constants
an = alphan(V);
bn = betan(V);
am = alpham(V);
bm = betam(V);
ah = alphah(V);
bh = betah(V);

% actual function
dV(1) = 1/C*(I - GNa*m^3*h*(V-ENa) - GK*n^4*(V-EK) - GL*(V-EL));
sigmam = 2/Nmax*(am*bm)/(am+bm);
sigmah = 2/Nmax*(ah*bh)/(ah+bh);
sigman = 2/Nmax*(an*bn)/(an+bn);
dV(2) = an*(1-n)-bn*n + normrnd(0, sigman);
dV(3) = am*(1-m)-bm*m + normrnd(0, sigmam);
dV(4) = ah*(1-h)-bh*h + normrnd(0, sigmah);
end
function an = alphan(V)
an = .01*(V+10)./(exp((V+10)/10)-1);
an(V == -10) = .1;
end
function bn = betan(V)
bn = .125*exp(V/80);
end
function n = n_stable(V)
n = alphan(V)./(alphan(V)+betan(V));
end
function am = alpham(V)
am = .1*(V+25)./(exp((V+25)/10)-1);
am(V == -25) = 1;
end
function bm = betam(V)
bm = 4*exp(V/18);
end
function m = m_stable(V)
m = alpham(V)./(alpham(V)+betam(V));
end
function ah = alphah(V)
ah = .07*exp(V/20);
end
function bh = betah(V)
bh = 1./(exp((V+30)/10)+1);
end
function h = h_stable(V)
h = alphah(V)./(alphah(V)+betah(V));
end
function z = steady_state(V)
GNa = 120; % mSiemens
GK = 36;
GL = .3;
% E = Electric field due to ion concentration
EL = 10.613;
ENa = 115;
EK = -12;
z = - GL*(V-EL) - GK*n_stable(V).^4.*(V-EK) - GNa*m_stable(V).^3.*h_stable(V).*(V-ENa);
end