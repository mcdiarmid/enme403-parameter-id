% Set-up
tutorial_setup

% Populate A and b matrix
b = T_data';
A(:, 1) = 1;
tri = 0:dt_data/2:1;
tri = flip(tri)';
[trilen, ~] = size(tri);

% Triangle of width 2s
for c = 2:1:4
    for ii = 1:1:rows
        if t_impulse(c-1) <= t_meas(ii)
            A(ii:ii+trilen-1, c) = tri;
            break
        end
    end
end

theta = A\b;
k = theta(1); J = theta(2:end) * V * Cp;
fprintf('k=%.2f, J1=%.0f, J2=%.0f, J3=%.0f\n', k, J(1), J(2), J(3));

% Populate A and b matrix
b = T_data';
A = zeros(rows, cols);
A(:, 1) = 1;

% Exponentially decaying: exp(-2*(t-t_impulse(x))
for c = 2:1:4
  for ii = 1:1:rows
    if t_impulse(c-1) <= t_meas(ii)
      A(ii:end, c) = exp(-2*(t_meas(ii:end) - t_meas(ii)));
      break
    end
  end
end

theta = A\b;
k = theta(1); J = theta(2:end) * V * Cp;
fprintf('k=%.2f, J1=%.0f, J2=%.0f, J3=%.0f\n', k, J(1), J(2), J(3));