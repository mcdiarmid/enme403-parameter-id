tutorial_setup

% Forward simulation function
Temp_f = @(ts, t0, th, ph) exp(-th(1)*ts) .* (t0+cumtrapz(ts, ...
  exp(th(1)*ts).*((ph*th(2:end))+(th(1)*T_amb))'));

% Initializations
[~, n_data] = size(T_data); 
J = zeros(n_data, cols);
n_iterations = 30;
T0 = mean(T_data(1:2));

dt_HR = 0.025;
t_HR = 0:dt_HR:10;
n_HR = 10/dt_HR;
phi = zeros(size(t_HR'));

for c = 1:1:3
  for ii = 1:1:n_HR
    if t_impulse(c) <= t_HR(ii)
      phi(ii, c) = 1/dt_HR;
      break
    end
  end
end

theta = zeros(cols, n_iterations);
psi = zeros(rows, 1);
TT = 0:dt_data/dt_HR:n_HR;

% Initial guess
delta = [0.05, 100, 10, 10]';  % Pertubations
lambda = 1;

% Guess theta(1), calculate psi(1)
theta(:, 1) = [0 0 0 0]';  % Initial conditions/guess
T_fwd = Temp_f(t_HR, T0, theta(:, 1), phi);
psi(:, 1) = T_fwd(TT+1)' - T_data';

for it = 1:n_iterations-1
  % Calculate J: same for all Gauss-Newton, Levenberg, Levenberg-Marquardt
  for col = 1:4
    dtheta = zeros(4,1); dtheta(col, 1) = delta(col, 1);
    dpsidth_i =  (Temp_f(t_HR, T0, theta(:, it) + dtheta, phi) - T_fwd)/delta(col);
    J(:, col) = dpsidth_i(TT+1)';
  end
  % Next part is different for next iteration of theta!
  % Levenberg
  % theta(:, it+1) = theta(:, it) - (((J'*J)+lambda*eye(4))^-1) * (J'*psi(:, it));
  % T_fwd = Temp_f(t_HR, T0, theta(:, it+1), phi);
  % psi(:, it+1) = T_fwd(TT+1)' - T_data';
  
  % Gauss-Newton
  % theta(:, it+1) = theta(:, it) - ((J'*J)^-1)*(J'*psi(:, it));
  % T_fwd = Temp_f(t_HR, T0, theta(:, it+1), phi);
  % psi(:, it+1) = T_fwd(TT+1)' - T_data';
  
  % Levenberg-Marquardt
  while true
    % Calculate theta(ii+1)
    theta(:, it+1) = theta(:, it) - (((J'*J)+lambda*eye(4).*(J'*J))^-1) * (J'*psi(:, it));
    
    % Forward simulate for psi(ii+1)
    T_fwd = Temp_f(t_HR, T0, theta(:, it+1), phi);
    psi(:, it+1) = T_fwd(TT+1)' - T_data';
    
    % Levenberg-Marquardt
    % Update lambda, repeat if norm(psi(ii+1)) > norm(psi(ii))
    if norm(psi(:, it+1)) <= norm(psi(:, it))
      lambda = lambda / 5;
      break
    else
      lambda = lambda * 5;
    end
  end
    
  % Calculate the Paramaters to ID
  k = theta(1, it+1);
  Jx = theta(2:end, it+1)*V*Cp;
  fprintf('GN: k=%.2f, J1=%.0f, J2=%.0f, J3=%.0f\n', k, Jx(1), Jx(2), Jx(3));
end

fprintf('\n');

% Steepest Descent
n_iterationsSD = 20000;  % Takes a real fucking long time
thetaSD = zeros(cols, n_iterationsSD); thetaSD(:, 1) = [10 1 1 1]';
deltaSD = theta(:, end)/1000;
JSD = zeros(cols, 1);
alpha = 0.01;

for it = 1:n_iterationsSD-1
  T_fwdSD = Temp_f(t_HR, T0, thetaSD(:, it), phi);
  psi_thi = norm(T_fwdSD(TT+1)' - T_data');
  
  for ii = 1:4
    dtheta = zeros(4,1); dtheta(ii, 1) = deltaSD(ii, 1);
    T_fwdSD_i = Temp_f(t_HR, T0, thetaSD(:, it) + dtheta, phi);
    [T_fwdSD(TT+1)' T_fwdSD_i(TT+1)' T_data'];
    psi_thi_i =  norm(T_fwdSD_i(TT+1)' - T_data');
    dpsi_dth_i = (psi_thi_i - psi_thi)/deltaSD(ii);
    [dpsi_dth_i psi_thi_i psi_thi];
    JSD(ii, 1) = dpsi_dth_i;  % [4x1]
  end
  JSD;
  
  thetaSD(:, it+1) = thetaSD(:, it) - alpha*JSD;
  k = thetaSD(1, it+1);
  Jx = thetaSD(2:end, it+1)*V*Cp;
  fprintf('SD: k=%.2f, J1=%.0f, J2=%.0f, J3=%.0f\n', k, Jx(1), Jx(2), Jx(3));
end

figure(403);
yyaxis left;
plot([0 n_iterationsSD], [parent(2) parent(2)]/1000, 'k', ...
    [0 n_iterationsSD], [parent(3) parent(3)]/1000, 'k', ...
    [0 n_iterationsSD], [parent(4) parent(4)]/1000, 'k');
hold on
plot(1:n_iterationsSD, thetaSD(2, :)*V*Cp/1000, '.b', ...
  1:n_iterationsSD, thetaSD(3, :)*V*Cp/1000, '.b', ...
  1:n_iterationsSD, thetaSD(4, :)*V*Cp/1000, '.b'); 
ylim([0, 550]);
xlabel('Iterations');ylabel('J_x [kJ]')
hold all;

yyaxis right;
plot(1:n_iterationsSD, thetaSD(1, :), '.r');
plot([0 n_iterationsSD], [parent(1) parent(1)], 'r-');
xlabel('Iterations');ylabel('k [s^-^1]')
ylim([0, 4]);

figure(401);
yyaxis left;
plot([0 n_iterations], [parent(2) parent(2)]/1000, 'k', ...
    [0 n_iterations], [parent(3) parent(3)]/1000, 'k', ...
    [0 n_iterations], [parent(4) parent(4)]/1000, 'k');
hold on
plot(1:n_iterations, theta(2, :)*V*Cp/1000, '.-b', ...
  1:n_iterations, theta(3, :)*V*Cp/1000, '.b', ...
  1:n_iterations, theta(4, :)*V*Cp/1000, ':b'); 
ylim([0, 550]);
xlabel('Iterations');ylabel('J_x [kJ]')
hold all;

yyaxis right;
plot(1:n_iterations, theta(1, :), '.r');
plot([0 n_iterations], [parent(1) parent(1)], 'r-');
xlabel('Iterations');ylabel('k [s^-^1]')
ylim([0, 4]);

figure(402);
plot(t_HR, T_fwd, 'ob', t_HR, T_fwdSD, '--r'); ylim([20, 45]);
ylabel('Temp (^oC)'); xlabel('Time (s)');



