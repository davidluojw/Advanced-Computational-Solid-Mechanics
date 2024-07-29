% plot the error
figure;
log_h = [1/12,1/16,1/20,1/24,1/28];
log_h = log(log_h);



% =========================================================================
% deg=1
log_L2_error = [0.0941442113024744,0.0520062864820762,0.0329835779523977,...
    0.0227881966607843,0.0166895611764418];
log_L2_error = log(log_L2_error);
log_H1_error = [0.572190296751414,0.426447287526867,0.340145308777427,...
    0.282991609245218,0.242323222242737];
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