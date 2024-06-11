tutorial_setup

% Part a

[~, n_data] = size(T_data);

b = T_data(2:end)';
A = zeros(n_data-1, 5);
A(:, 1) = T_data(1:end-1)';
A(:, 5) = 1;

for c = 2:1:4
  for ii = 1:1:n_data
    if t_impulse(c-1) <= t_meas(ii)
      A(ii-1, c) = 1;
      break
    end
  end
end

theta = A\b;
spy(A)

T_arx1 = zeros(size(T_data));
T_arx1(1) = T_data(1);
for t = 2:n_data
  T_arx1(t) = theta(1) * T_arx1(t-1) + theta(2:4)' * A(t-1, 2:4)' + theta(5);
end

figure(302);
plot(t_meas, T_arx1, 'r', t_meas, T_data, '+k');
ylim([20, 45]); xlabel('Time [s]'); ylabel('Temperature [^oC]');
hold on;

% Part b
t_s = 0:0.25:3;
[~, rows] = size(t_s);
phi = A(:, 2) + 2*A(:, 3) + 1.5*A(:, 4);  % Combine inputs into one fn
A2 = zeros(rows, 3);
A2(:, 1) = T_data(1:rows);  
A2(:, 2) = A(1:rows, 2); A2(end, 2) = 2;
A2(:, 3) = 1;
b2 = T_data(2:rows+1)';

theta2 = A2\b2;  % Used for simulating, not param ID

T_arx2 = zeros(size(T_data));
T_arx2(1) = T_data(1);

for t = 2:n_data
  T_arx2(t) = theta2(1)*T_arx2(t-1) + theta2(2)*phi(t-1) + theta2(3);
end

plot(t_meas, T_arx2, 'b');

