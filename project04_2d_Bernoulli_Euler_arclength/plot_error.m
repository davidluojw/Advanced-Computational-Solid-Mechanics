% plot the error
figure;
log_h = [1/2,1/4,1/6,1/8,1/10,1/12,1/14];
log_h = log(log_h);


% =========================================================================
% Hermite
log_L2_error = [0.0230918627872429,0.00147847288916694,0.000293314667631010,...
    9.29468966628678e-05,3.80976187593281e-05,1.83796513878782e-05,9.92314184585987e-06];
log_L2_error = log(log_L2_error);
log_H1_error = [0.0289704282723901,0.00370697038891401,0.00110299801072476,...
    0.000466009941701609,0.000238758694684451,0.000138221315983554,8.70624353127395e-05];
log_H1_error = log(log_H1_error);
% =========================================================================



plot(log_h, log_L2_error,'b-', 'LineWidth', 2);
hold on;
plot(log_h, log_H1_error,'r-', 'LineWidth', 2);
xlabel("log-h");
ylabel("log-error");
legend('L2-error', 'H1-error','Location', 'Best', 'FontSize', 14, 'Box', 'on');

% rate of error
L2_rate_error = zeros(length(log_h)-1,1);
H1_rate_error = zeros(length(log_h)-1,1);
for ii = 1:length(log_h)-1
    L2_rate_error(ii) = (log_L2_error(ii+1) - log_L2_error(ii))/(log_h(ii+1)-log_h(ii));
    H1_rate_error(ii) = (log_H1_error(ii+1) - log_H1_error(ii))/(log_h(ii+1)-log_h(ii));
end