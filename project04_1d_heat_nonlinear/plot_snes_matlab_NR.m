snes_uh = zeros(21, 10);

snes_uh(:,1) = [0;
                0;
                0;
                0;
                0;
                0;
                0;
                0;
                0;
                0;
                0;
                0;
                0;
                0;
                0;
                0;
                0;
                0;
                0;
                0;
                -0.54402111088937;
                ];

snes_uh(:, 2) = [0.0469949278790754;
                 0.57501104145373;
                 1.08057529775141;
                 1.3479061750176;
                 1.19767401011823;
                 0.715239713568341;
                 0.167463772260764;
                 -0.353022301177344;
                 -0.875284309963495;
                 -1.2541969406547;
                 -1.24263291106098;
                 -0.819175896179746;
                 -0.270823780105879;
                 0.234746862947046;
                 0.749431277873764;
                 1.20386942379745;
                 1.31813268224586;
                 0.96681092762493;
                 0.419903037643202;
                 -0.0847367177318517;
                 -0.54402111088937;
                ];
snes_uh(:, 3) = [0.0549059060536684;
                 0.537451767690144;
                 0.890500858245932;
                 1.06064965233445;
                 0.965876375311208;
                 0.64241967254036;
                 0.172530387965614;
                 -0.338155362023379;
                 -0.744861656584688;
                 -0.995226704871788;
                 -0.988503507832703;
                 -0.70622594523455;
                 -0.255920918055023;
                 0.242682224540777;
                 0.665482989973536;
                 0.969464032615956;
                 1.04267692962429;
                 0.815527096607978;
                 0.404975411850918;
                 -0.0838552694652632;
                 -0.54402111088937;
                ];

snes_uh(:, 4) = [0.0542837868156485;
                 0.535380042111042;
                 0.870394471632741;
                 1.01148186057102;
                 0.935349043317339;
                 0.63927999640147;
                 0.172864406484101;
                 -0.337254627751419;
                 -0.736363587021987;
                 -0.955996513147625;
                 -0.951783289435057;
                 -0.700441198608224;
                 -0.257546843211863;
                 0.240613683850835;
                 0.661258363722813;
                 0.937200735202064;
                 0.998349213609998;
                 0.803212909024789;
                 0.404923109812608;
                 -0.0838565839587557;
                 -0.54402111088937;
                ];

snes_uh(:, 5) = [0.0542778231408796;
                 0.53535955074199;
                 0.870198861138691;
                 1.01016393825826;
                 0.934901532684726;
                 0.63925394100857;
                 0.172881627945673;
                 -0.337236348547008;
                 -0.73632003446316;
                 -0.955182285809951;
                 -0.951120486420764;
                 -0.700405853487952;
                 -0.257574594331239;
                 0.24058379672858;
                 0.661229181098472;
                 0.936661482748501;
                 0.997340984616896;
                 0.803104079267787;
                 0.404923109095211;
                 -0.0838565839939879;
                 -0.54402111088937;
                ];

snes_uh(:, 6) = [0.0542778232026844;
                 0.535359549307029;
                 0.870198843193999;
                 1.01016293836254;
                 0.934901458548465;
                 0.639253961383673;
                 0.172881667801138;
                 -0.337236311408099;
                 -0.736320006848061;
                 -0.955181917449884;
                 -0.951120267588364;
                 -0.700405826984313;
                 -0.257574583745076;
                 0.240583807090849;
                 0.661229188452795;
                 0.936661344745857;
                 0.997340483036329;
                 0.803104036506538;
                 0.404923109095211;
                 -0.0838565839939878;
                 -0.54402111088937;
                ];

% plot the solution
figure;
X_h = omega_l: hh/pp :omega_r;
for i=1:counter
    Y_h = snes_uh(:,i);
    h_fem = plot(X_h, Y_h,'b-', 'LineWidth', 2);
    hold on;
end
X = omega_l:0.01:omega_r;
Y = exact(X);
h_exact = plot(X, Y,'r--', 'LineWidth', 2);

legend([h_fem, h_exact], {'FEM', 'EXACT'}, 'Location', 'Best', 'FontSize', 14, 'Box', 'on');
xlabel("X");
ylabel("Temperature");


snes_error = zeros(10, 1);

snes_error(1,:) = 2.900846405359e+00;
snes_error(2,:) = 3.371646644088e+00;
snes_error(3,:) = 5.530154587081e-01;
snes_error(4,:) = 1.698946870813e-02;
snes_error(5,:) = 1.410032504737e-05;
snes_error(6,:) = 8.732341124723e-12;


log_error = log(error_set);
figure;
plot(log_error(1:counter-1), log_error(2:counter),'b-', 'LineWidth', 2);

snes_log_error = log(snes_error);
hold on;
plot(snes_log_error(1:counter-1), snes_log_error(2:counter),'r--', 'LineWidth', 2);
xlabel("ln(r_{i})");
ylabel("ln(r_{i+1})");

legend('Matlab','PetscSNES','Location', 'Best', 'FontSize', 14, 'Box', 'on');

uh_error = snes_uh(:, 1:counter) - uh_set(:, 1:counter); 

uh_error_norm = zeros(counter, 1);
for i = 1:counter
    uh_error_norm(i) = norm(uh_error(:,i));
end

