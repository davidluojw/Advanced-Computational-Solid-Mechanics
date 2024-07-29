% plot the error
figure;
log_h = [1/12,1/16,1/20,1/24,1/28];
log_h = log(log_h);


% =========================================================================
% deg=1
log_L2_error = [0.0979075908235058,0.0830596670159659,0.0778037277190330,...
    0.0756722952085487,0.0746991192539494];
log_L2_error = log(log_L2_error);
log_H1_error = [0.276047499591993,0.212887585430118,0.176257650160837,...
    0.152928304061990,0.137089232566351];
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