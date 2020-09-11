%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%   Paper: Impact of UAV Wobbling on the Air-to-Ground                %%%
%%%          Wireless Channel                                           %%%
%%%   Authors: Morteza Banagar, Harpreet S. Dhillon, and                %%%
%%%            Andreas F. Molisch                                       %%%
%%%   Emails: mbanagar@vt.edu, hdhillon@vt.edu, molisch@usc.edu         %%%
%%%                                                                     %%%
%%%   This code is used to generate all the plots in this paper.        %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datetime('now')
%% Fig 1
T = 0.1;
tauVec = 0 : 0.0001 : T;
tauLen = length(tauVec);
aD = 40; % [cm]
aD = aD / 100;
Freq = 6e9;
Lambda = 3e8 / Freq;
ThetaMax = deg2rad(5);
sigma = 1;
K = 11.5;
N = 20;
NumIter = 1e3;
phi_0 = deg2rad(20);
phi_n = unifrnd(0, pi / 2, NumIter, 1);
cWin = 2 * pi ^ 2 * aD ^ 2 / Lambda ^ 2;
cSin = 2 / Lambda * aD * ThetaMax;
f_st = 5;
f_en = 25;

R_Wiener_1 = zeros(tauLen, 1);
i_tau = 0;
for tau = tauVec
    i_tau = i_tau + 1;
    R_Wiener_1(i_tau) = K * N * 1 / (2 * sigma) * exp(-cWin * cos(phi_0) ^ 2 * tau) * mean(exp(-abs(phi_n - phi_0) / sigma)) + N * 1 / (2 * sigma) * mean(exp(-abs(phi_n - phi_0) / sigma) .* exp(-cWin * cos(phi_n) .^ 2 * tau));
end

tVec = [0, 0.03, 0.07, 0.1];
tLen = length(tVec);
R_SineAF_1_All = zeros(tauLen, tLen);
i_t = 0;
for t = tVec
    disp(['t = ', num2str(t)])
    i_tau = 0;
    i_t = i_t + 1;
    for tau = tauVec
        i_tau = i_tau + 1;
        Int2 = zeros(NumIter, 1);
        for jj = 1 : NumIter
            Int2(jj) = integral(@(f) sinc(cSin * cos(phi_n(jj)) * (sin(2 * pi * f * (t + tau)) - sin(2 * pi * f * t))) / (f_en - f_st), f_st, f_en);
        end
        R_SineAF_1_All(i_tau, i_t) = K * N * 1 / (2 * sigma) * integral(@(f) sinc(cSin * cos(phi_0) * (sin(2 * pi * f * (t + tau)) - sin(2 * pi * f * t))) / (f_en - f_st), f_st, f_en) * mean(exp(-abs(phi_n - phi_0) / sigma)) + N * 1 / (2 * sigma) * mean(exp(-abs(phi_n - phi_0) / sigma) .* Int2);
    end
end
R_SineAF_1_All_Normalized = R_SineAF_1_All ./ repmat(R_SineAF_1_All(1, :), tauLen, 1);

figure(301)
hold on
plot(tauVec * 1000, R_Wiener_1 / R_Wiener_1(1), 'LineWidth', 2, 'LineStyle', '-')
plot(tauVec * 1000, R_SineAF_1_All_Normalized(:, 1), 'LineWidth', 2, 'LineStyle', '--')
plot(tauVec * 1000, R_SineAF_1_All_Normalized(:, 2 : end), 'LineWidth', 2, 'LineStyle', ':')
ylim([-0.05, 1])
hold off
xlabel('$\tau$ (ms)', 'Interpreter', 'latex')
ylabel('$R(0, \tau)$', 'Interpreter', 'latex')
legend({'Wiener Process'; 'Sinusoidal Process ($t = 0$ ms)'; 'Sinusoidal Process ($t = 30$ ms)'; 'Sinusoidal Process ($t = 70$ ms)'; 'Sinusoidal Process ($t = 100$ ms)'}, 'Location', 'northeast', 'FontSize', 12, 'Interpreter', 'latex')
grid on
box on
%% Fig 2
T = 0.1;
tauVec = 0 : 0.001 : T;
tauLen = length(tauVec);
aD = 40; % [cm]
aD = aD / 100;
Freqs = [2.4, 6, 30] * 1e9;
NumFreqs = length(Freqs);
Lambda = 3e8 ./ Freqs;
ThetaMax = deg2rad(5);
N = [20, 20, 10];
R_SineAF_2 = zeros(tauLen, NumFreqs);
for iLam = 1 : NumFreqs
    cSin = 2 / Lambda(iLam) * aD * ThetaMax;
    i_tau = 0;
    for tau = tauVec
        i_tau = i_tau + 1;
        Int2 = zeros(NumIter, 1);
        for jj = 1 : NumIter
            Int2(jj) = integral(@(f) sinc(cSin * cos(phi_n(jj)) * sin(2 * pi * f * tau)) / (f_en - f_st), f_st, f_en);
        end
        R_SineAF_2(i_tau, iLam) = K * N(iLam) * 1 / (2 * sigma) * integral(@(f) sinc(cSin * cos(phi_0) * sin(2 * pi * f * tau)) / (f_en - f_st), f_st, f_en) * mean(exp(-abs(phi_n - phi_0) / sigma)) + N(iLam) * 1 / (2 * sigma) * mean(exp(-abs(phi_n - phi_0) / sigma) .* Int2);
    end
end

figure(302)
plot(tauVec * 1000, R_SineAF_2 ./ repmat(R_SineAF_2(1, :), tauLen, 1), 'LineWidth', 2)
xlabel('$\tau$ (ms)', 'Interpreter', 'latex')
ylabel('$R(\tau)$', 'Interpreter', 'latex')
grid on
box on
%% Fig 3
T = 0.1;
tauVec = 0 : 0.001 : T;
tauLen = length(tauVec);
aD = 40; % [cm]
aD = aD / 100;
Freq = 2.4e9;
Lambda = 3e8 / Freq;
ThetaMax = deg2rad([5, 7, 10]);
NumTheta = length(ThetaMax);
N = 20;
R_SineAF_3 = zeros(tauLen, NumTheta);
for iTheta = 1 : NumTheta
    cSin = 2 / Lambda * aD * ThetaMax(iTheta);
    i_tau = 0;
    for tau = tauVec
        i_tau = i_tau + 1;
        Int2 = zeros(NumIter, 1);
        for jj = 1 : NumIter
            Int2(jj) = integral(@(f) sinc(cSin * cos(phi_n(jj)) * sin(2 * pi * f * tau)) / (f_en - f_st), f_st, f_en);
        end
        R_SineAF_3(i_tau, iTheta) = K * N * 1 / (2 * sigma) * integral(@(f) sinc(cSin * cos(phi_0) * sin(2 * pi * f * tau)) / (f_en - f_st), f_st, f_en) * mean(exp(-abs(phi_n - phi_0) / sigma)) + N * 1 / (2 * sigma) * mean(exp(-abs(phi_n - phi_0) / sigma) .* Int2);
    end
end

figure(303)
plot(tauVec * 1000, R_SineAF_3 ./ repmat(R_SineAF_3(1, :), tauLen, 1), 'LineWidth', 2)
xlabel('$\tau$ (ms)', 'Interpreter', 'latex')
ylabel('$R(\tau)$', 'Interpreter', 'latex')
grid on
box on

datetime('now')