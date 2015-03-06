function [Na_state_matrix K_state_matrix] = markov_process(N_Na, N_K, initial, I, dt)
% N_Na: number of sodium gates.
% N_K: number of potassium gates. The ratio of densities reported by H&H 
% was 10:3 sodium:potassium.
% Initial: Initial voltage
% I: Current input
% dt: time step size (ms)
% Example function call: markov_process(600, 180, 0, [zeros(1, 250) -100000*ones(1, 250) zeros(1, 1500)], 10^-3);

% define voltage; set initial voltage
V = zeros(1, length(I)); V(1) = initial;
% variable used for plotting
time = (1:length(I))*dt;

% Initialize states to 1
Na_state = ones(N_Na,1);                % Sodium state matrix
K_state = ones(N_K,1);                  % Potassium state matrix
Na_state_matrix = zeros(N_Na, length(I));% Sodium state history
Na_state_matrix(:, 1) = Na_state;       % Initialize history
K_state_matrix = zeros(N_K, length(I));  % Potassium state history
K_state_matrix(:, 1) = K_state;         % Initialize

% =========================================================================
% BEGIN SIMULATION
% =========================================================================
tic; % start timer
for time_step = 1:length(I)-1 % for every time step
    % disp([num2str(time_step) '/' num2str(length(I))]);

    % calculate alphas and betas to be used in the transition matrices
    an = alphan(V(time_step));
    bn = betan(V(time_step));
    am = alpham(V(time_step));
    bm = betam(V(time_step));
    ah = alphah(V(time_step));
    bh = betah(V(time_step));
    
    % transition rate matrix
    Na_transition_rates = [0 3*am 0 0 ah 0 0 0; bm 0 2*am 0 0 ah 0 0; ...
        0 2*bm 0 am 0 0 ah 0; 0 0 3*bm 0 0 0 0 ah; ...
        bh 0 0 0 0 3*am 0 0; 0 bh 0 0 bm 0 2*am 0; ...
        0 0 bh 0 0 2*bm 0 am; 0 0 0 bh 0 0 3*bm 0];
    K_transition_rates = [0 4*an 0 0 0; bn 0 3*an 0 0; ...
        0 2*bn 0 2*an 0; 0 0 3*bn 0 an; ...
        0 0 0 4*bn 0];
    
    % transition probability matrix
    Na_next_state = diag(1./sum(Na_transition_rates, 2))*Na_transition_rates;   % divide each row by the sum of the row to create uniform distribution
    Na_next_state = cumsum(Na_next_state, 2);                                   % cumsum creates a CDF from the PDF
    K_next_state = diag(1./sum(K_transition_rates, 2))*K_transition_rates;
    K_next_state = cumsum(K_next_state, 2);
    
    for each_gate = 1:N_Na % for each gate
        % generate a random number from the exponential distribution with
        % parameter (sum(transition rates(current state))). 
        if (exprnd(sum(Na_transition_rates(Na_state(each_gate), :), 2)) < 1000*dt)
            % if that number is less than the time step, the state
            % switches. pick a new state.
            Na_state(each_gate) = find(rand < Na_next_state(Na_state(each_gate), :), 1, 'first');
        end
    end
    for each_gate = 1:N_K
        if (exprnd(sum(K_transition_rates(K_state(each_gate), :), 2)) < 1000*dt)
            K_state(each_gate) = find(rand < K_next_state(K_state(each_gate), :), 1, 'first');
        end
    end
    
    % perform runge-kutta step, using sum(K_state == 5) and sum(Na_state ==
    % 4) to count the number of open potassium and sodium channels.
    k1 = dt * transmembranepotential(V(time_step), I(time_step), sum(K_state == 5), sum(Na_state == 4));
    k2 = dt * transmembranepotential(V(time_step) + k1/2, I(time_step), sum(K_state == 5), sum(Na_state == 4));
    k3 = dt * transmembranepotential(V(time_step) + k2/2, I(time_step), sum(K_state == 5), sum(Na_state == 4));
    k4 = dt * transmembranepotential(V(time_step) + k3, I(time_step), sum(K_state == 5), sum(Na_state == 4));
    V(time_step + 1) = V(time_step) + 1/6*(k1 + 2*k2 + 2*k3 + k4);
    
    % record the states to the state matrix to plot later
    Na_state_matrix(:, time_step) = Na_state;
    K_state_matrix(:, time_step) = K_state;
end
t = toc/(length(I)-1); % divide total time by number of steps to get average time per step
disp(['time per step = ' num2str(t)]);

% =========================================================================
% Plot Things
% =========================================================================
% figure, [ax, h1, h2] = plotyy(time, -V, time, -I/1000);
% set(get(ax(1), 'Ylabel'), 'String', 'Voltage (mV)');
% % set(get(ax(2), 'Ylabel'), 'String', 'Current (\muA)');
% figure, plot(time, -V/50);
% xlabel('Time (ms)');
% ylabel('Voltage (mV)');
% figure, plot(time, sum(Na_state_matrix == 1), ...
%     time, sum(Na_state_matrix == 2), time, sum(Na_state_matrix == 3), ...
%     time, sum(Na_state_matrix == 4), time, sum(Na_state_matrix == 5), ...
%     time, sum(Na_state_matrix == 6), time, sum(Na_state_matrix == 7), ...
%     time, sum(Na_state_matrix == 8));
% xlabel('Time (ms)'); ylabel('Number of Sodium Channels');
% legend('1', '2', '3', '4', '5', '6', '7', '8'); 
% title('Potassium states');
% figure, plot(time, sum(K_state_matrix == 1), ...
%     time, sum(K_state_matrix == 2), time, sum(K_state_matrix == 3), ...
%     time, sum(K_state_matrix == 4), time, sum(K_state_matrix == 5));
% xlabel('Time (ms)'); ylabel('Number of Potassium Channels');
% legend('1', '2', '3', '4', '5'); 
% title('Sodium States');
% figure, plot(time, .2*(sum(Na_state_matrix == 4)), time, .2*(sum(K_state_matrix == 5)));
% xlabel('Time (ms)');
% ylabel('Conductance (mS)*');
% legend('Potassium Open States', 'Sodium Open States');

end
% =========================================================================
% FUNCTION DEFINITIONS
% =========================================================================

function dV = transmembranepotential(V, I, n4, m3h1)
% Governing equation.
% V: input voltage (mV)
% I: input current (uA)
% n4: number of open potassium channels
% m3h1: number of open sodium channels

% C = membrane Capacitance
C = 1; % uF/um
% G = conductance of a single channel
GNa = 1*m3h1; % pS*m3h1
GK = 1*n4;
GL = 10; %pS
% E = Electric field due to ion concentration
EL = 10.613; % mV
ENa = 115;
EK = -12;

% actual function
dV = 1/C*(I - GNa*(V-ENa) - GK*(V-EK) - GL*(V-EL));
end
function an = alphan(V)
an = .01*(V+10)./(exp((V+10)/10)-1);
an(V == -10) = .1; % the function is singular, but the limit exists
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