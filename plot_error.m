% plot the error
clear all;
clc
figure;
% log_h = [1/12,1/16,1/20,1/24,1/28];
% log_h = log(log_h);

N = [1, 2, 3, 4, 5];

% N = log10(N);

% =========================================================================
% deg=1
% log_L2_error = [0.0941442113024744,0.0520062864820762,0.0329835779523977,...
%     0.0227881966607843,0.0166895611764418];
% log_L2_error = log(log_L2_error);
% log_H1_error = [0.572190296751414,0.426447287526867,0.340145308777427,...
%     0.282991609245218,0.242323222242737];
% log_H1_error = log(log_H1_error);

% creep m = 1, n = 0
log_Rn = [8.784337e+01, 9.799771e-01, 2.334465e-02, 5.656858e-04, 1.378789e-05,...
    3.364457e-07];
log_Rn = log(log_Rn);

log_Rnp1 = [9.799771e-01, 2.334465e-02, 5.656858e-04, 1.378789e-05,...
    3.364457e-07,  1.638541e-08];
log_Rnp1 = log(log_Rnp1);

% creep m = 2, n = 0
% log_Rn = [8.062854e+01 , 5.471516e-02, 7.528160e-04, 1.725369e-05, 3.963512e-07];
% log_Rn = log(log_Rn);
% 
% log_Rnp1 = [5.471516e-02, 7.528160e-04, 1.725369e-05, 3.963512e-07, 2.806284e-08];
% log_Rnp1 = log(log_Rnp1);


% creep m = 3, n = 0
% log_Rn = [8.134012e+01, 1.801351e-01, 4.262534e-03, 9.321802e-05, 2.040403e-06];
% log_Rn = log(log_Rn);
% 
% log_Rnp1 = [1.801351e-01, 4.262534e-03, 9.321802e-05, 2.040403e-06, 4.818888e-08];
% log_Rnp1 = log(log_Rnp1);
% 
% log_Rnp1dRn = zeros(length(log_Rn),1);
% for ii = 1:length(log_Rn)
%     log_Rnp1dRn(ii) = log_Rnp1(ii) / log_Rn(ii);
% end


% local 
% log_Rn = [0.00422579, 2.19807e-07];
% log_Rn = log(log_Rn);
% 
% log_Rnp1 = [ 2.19807e-07, 6.25366e-15];
% log_Rnp1 = log(log_Rnp1);


% shear, m = 0, n = 0
% log_Rn = [3.045879e-03, 1.515003e-06, 2.751721e-09];
% log_Rn = log(log_Rn);
% 
% log_Rnp1 = [1.515003e-06, 2.751721e-09, 5.783669e-12];
% log_Rnp1 = log(log_Rnp1);


% shear, m = 1, n = 0
log_Rn = [3.084982e-03, 1.637719e-06, 3.045452e-09];
log_Rn = log(log_Rn);

log_Rnp1 = [1.637719e-06, 3.045452e-09, 6.434694e-12];
log_Rnp1 = log(log_Rnp1);

% old
% log_Rn = [9.776327e+03, 4.856300e+01, 3.073925e-02, 1.461154e-05];
% log_Rn = log(log_Rn);
% 
% log_Rnp1 = [4.856300e+01, 3.073925e-02, 1.461154e-05, 1.502559e-07];
% log_Rnp1 = log(log_Rnp1);


% old 2
% log_Rn = [9.776327e+03, 1.367437e+02, 8.827810e+00, 1.285547e+00, 2.322286e-01,...
%     4.872784e-02, 1.111891e-02, 2.668642e-03, 6.625706e-04, 1.677938e-04];
% % log_Rn = log(log_Rn);
% 
% log_Rnp1 = [1.367437e+02, 8.827810e+00, 1.285547e+00, 2.322286e-01, 4.872784e-02,...
%     1.111891e-02, 2.668642e-03, 6.625706e-04, 1.677938e-04, 4.312769e-05];
% % log_Rnp1 = log(log_Rnp1);
% 
% log_Rnp1dRn = zeros(length(log_Rn),1);
% for ii = 1:length(log_Rn)
%     log_Rnp1dRn(ii) = log_Rnp1(ii) / log_Rn(ii);
% end

% log_Rn = log(log_Rn);

% =========================================================================



plot(log_Rn, log_Rnp1,'b-', 'LineWidth', 2);
xlabel("log-Rn");
ylabel("log-Rnp1");
legend('log-Residual', 'Location', 'Best', 'FontSize', 14, 'Box', 'on');

% rate of error
R_rate_error = zeros(length(log_Rn)-1,1);
for ii = 1:length(log_Rn)-1
    R_rate_error(ii) = (log_Rnp1(ii+1) - log_Rnp1(ii))/(log_Rn(ii+1)-log_Rn(ii));
end






