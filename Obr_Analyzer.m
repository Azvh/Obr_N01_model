%% ================== Obr_Analyzer.m ============================
% Повний аналіз Ia, Ib, Ic з моделі Simulink
% Автор: А.В. (з оптимізацією від AI)
% Версія: стабільнащзжд 

clc; close all;

%% ---------- НАЛАШТУВАННЯ ---------------------------------------
model    = 'Obr_N02.slx';   % ваша модель
f0       = 50;              % частота мережі
Ncycles  = 5;               % кількість періодів у вікні
tstart   = 0.20;            % час початку вікна (режим from)
refMode  = 'mean';          % 'mean' або 'max'
removeDC = false;           % true → видалити DC перед RMS

fprintf('\nЗапуск %s ...\n', model);
simOut = sim(model, 'ReturnWorkspaceOutputs','on');

%% ---------- ВИТЯГ ДАНИХ ---------------------------------------
[tA, iA] = fetch_tx(simOut,'Ia');
[tB, iB] = fetch_tx(simOut,'Ib');
[tC, iC] = fetch_tx(simOut,'Ic');

%% ---------- Вибір останніх N циклів ---------------------------
[tAw, xAw] = select_window(tA,iA,Ncycles,f0);
[tBw, xBw] = select_window(tB,iB,Ncycles,f0);
[tCw, xCw] = select_window(tC,iC,Ncycles,f0);

if removeDC
    xAw = xAw - mean(xAw);
    xBw = xBw - mean(xBw);
    xCw = xCw - mean(xCw);
end

%% ---------- RMS ------------------------------------------------
IA = rms_trapz(tAw,xAw);
IB = rms_trapz(tBw,xBw);
IC = rms_trapz(tCw,xCw);

Ivec = [IA, IB, IC];
switch refMode
    case 'mean', Iref = mean(Ivec);
    case 'max',  Iref = max(Ivec);
end

kA = Iref/IA; 
kB = Iref/IB; 
kC = Iref/IC;

iAeq = kA * iA;
iBeq = kB * iB;
iCeq = kC * iC;

ikw  = iA.^2 + iB.^2 + iC.^2;
ikwq = iAeq.^2 + iBeq.^2 + iCeq.^2;

%% ---------- Вікно за tstart -----------------------------------
[twA, fwA] = select_window_from(tA,iA ,f0,tstart,Ncycles);
[twB, fwB] = select_window_from(tB,iB ,f0,tstart,Ncycles);
[twC, fwC] = select_window_from(tC,iC ,f0,tstart,Ncycles);

%% ---------- ВИДІЛЕННЯ ФУНДАМЕНТАЛИ (ФАЗОРИ) -------------------
w = 2*pi*f0;
T = Ncycles/f0;

% вікно Ганна
wA = hann(length(tAw)); 
wB = hann(length(tBw)); 
wC = hann(length(tCw));

% фінальні вибірки вікна
XAw = xAw .* wA;
XBw = xBw .* wB;
XCw = xCw .* wC;

% комплексні фазори
IAc = (2/T)*trapz(tAw, XAw.*exp(-1j*w*tAw));
IBc = (2/T)*trapz(tBw, XBw.*exp(-1j*w*tBw));
ICc = (2/T)*trapz(tCw, XCw.*exp(-1j*w*tCw));

phiA_deg = rad2deg(angle(IAc));
phiB_deg = rad2deg(angle(IBc));
phiC_deg = rad2deg(angle(ICc));

wrap = @(x) mod(x+180,360)-180;

phiAB = wrap(phiA_deg - phiB_deg);
phiBC = wrap(phiB_deg - phiC_deg);
phiCA = wrap(phiC_deg - phiA_deg);

%% ---------- ВІКНО tstart (друга серія фаз) --------------------
Aw = hann(length(twA)); 
Bw = hann(length(twB)); 
Cw = hann(length(twC));

IA2 = (2/T)*trapz(twA, fwA.*Aw.*exp(-1j*w*twA));
IB2 = (2/T)*trapz(twB, fwB.*Bw.*exp(-1j*w*twB));
IC2 = (2/T)*trapz(twC, fwC.*Cw.*exp(-1j*w*twC));

phiA2 = rad2deg(angle(IA2));
phiB2 = rad2deg(angle(IB2));
phiC2 = rad2deg(angle(IC2));

phiAB2 = wrap(phiA2 - phiB2);
phiBC2 = wrap(phiB2 - phiC2);
phiCA2 = wrap(phiC2 - phiA2);

%% ========= ВИВІД -------------------------------------------------

fprintf('\n=== ФАЗИ (останні %g циклів) ===\n',Ncycles);
fprintf('φA=%.2f°, φB=%.2f°, φC=%.2f°\n',phiA_deg,phiB_deg,phiC_deg);
fprintf('ΔAB=%.2f°, ΔBC=%.2f°, ΔCA=%.2f°\n',phiAB,phiBC,phiCA);

fprintf('\n=== ФАЗИ (tstart=%.3f, %g циклів) ===\n',tstart,Ncycles);
fprintf('φA=%.2f°, φB=%.2f°, φC=%.2f°\n',phiA2,phiB2,phiC2);
fprintf('ΔAB=%.2f°, ΔBC=%.2f°, ΔCA=%.2f°\n',phiAB2,phiBC2,phiCA2);

fprintf('\n=== RMS ===\n');
fprintf('Ia=%.4f A, Ib=%.4f A, Ic=%.4f A\n', IA,IB,IC);
fprintf('Після вирівнювання: Ia=%.4f A, Ib=%.4f A, Ic=%.4f A\n',...
    rms_trapz(tAw,kA*xAw), rms_trapz(tBw,kB*xBw), rms_trapz(tCw,kC*xCw));

%% ------------------ ГРАФІКИ -----------------------------------

figure('Name','Фазні струми');
subplot(3,1,1); plot(tA,iA,'b',tA,iAeq,'r'); grid on; title('Ia');
subplot(3,1,2); plot(tB,iB,'b',tB,iBeq,'r'); grid on; title('Ib');
subplot(3,1,3); plot(tC,iC,'b',tC,iCeq,'r'); grid on; title('Ic');

figure('Name','Сума квадратів');
plot(tA,ikw,'b',tA,ikwq,'r'); grid on;

figure('Name','Векторна діаграма');
hold on; axis equal; grid on;
quiver(0,0,real(IAc), imag(IAc),'r','LineWidth',2);
quiver(0,0,real(IBc), imag(IBc),'g','LineWidth',2);
quiver(0,0,real(ICc), imag(ICc),'b','LineWidth',2);
legend('Ia','Ib','Ic'); title('Фазори (50 Гц)');


figure('Name','tA,iA');
plot(tA,iA); grid on;
%% === FFT аналіз 3-ї та 5-ї гармонік у вікні tstart ===

H_A = fft_harmonics(twA, fwA, f0);
H_B = fft_harmonics(twB, fwB, f0);
H_C = fft_harmonics(twC, fwC, f0);

fprintf('\n=== FFT гармоніки (вікно tstart=%.3f, %g циклів) ===\n', tstart, Ncycles);

fprintf('\nФаза A:\n');
fprintf('H1 = %.3f A,  H3 = %.3f A (%.2f%%),  H5 = %.3f A (%.2f%%)\n', ...
    H_A.H1, H_A.H3, H_A.H3p, H_A.H5, H_A.H5p);

fprintf('\nФаза B:\n');
fprintf('H1 = %.3f A,  H3 = %.3f A (%.2f%%),  H5 = %.3f A (%.2f%%)\n', ...
    H_B.H1, H_B.H3, H_B.H3p, H_B.H5, H_B.H5p);

fprintf('\nФаза C:\n');
fprintf('H1 = %.3f A,  H3 = %.3f A (%.2f%%),  H5 = %.3f A (%.2f%%)\n', ...
    H_C.H1, H_C.H3, H_C.H3p, H_C.H5, H_C.H5p);
%% ===========================================================
figure('Name','FFT спектр (до вирівнювання)');
hold on; grid on;

plot(H_A.f, H_A.X, 'r');
plot(H_B.f, H_B.X, 'g');
plot(H_C.f, H_C.X, 'b');

xlim([0 500]);  % показати до 10-ї гармоніки
xlabel('Frequency, Hz');
ylabel('Amplitude');
title('FFT (до вирівнювання)');
legend('Ia','Ib','Ic');

%% ============ ЛОКАЛЬНІ ФУНКЦІЇ ==================================

function [t,x] = fetch_tx(simOut,name)
    if ismember(name,simOut.who)
        v = simOut.get(name);
        [t,x] = normalize_tx_from_var(v,simOut); return;
    end
    if ismember('logsout',simOut.who)
        try sig = simOut.logsout.get(name); catch, sig=[]; end
        if ~isempty(sig)
            [t,x]=normalize_tx_from_var(sig.Values,simOut); return;
        end
    end
    error('Сигнал %s не знайдено',name);
end

function [t,x] = normalize_tx_from_var(v,simOut)
    if isa(v,'timeseries')
        t=v.Time; x=v.Data; return;
    end
    if isnumeric(v)
        if size(v,2)>=2
            t=v(:,1); x=v(:,2); return;
        else
            tt=simOut.get('tout'); t=tt(:); x=v(:); return;
        end
    end
    if isstruct(v)
        t=v.time; x=v.signals.values; return;
    end
end

function [tw,xw] = select_window(t,x,Ncycles,f0)
    T=Ncycles/f0;
    t1=t(end)-T;
    m=(t>=t1)&(t<=t(end));
    if ~any(m), m=t>=t(end)-T; end
    tw=t(m); xw=x(m);
end

function [tw,xw] = select_window_from(t,x,f0,tStart,Ncycles)
    T=Ncycles/f0;
    [~,iStart]=min(abs(t-tStart));
    t0=t(iStart);
    tEnd=t0+T;
    m=(t>=t0)&(t<=tEnd);
    tw=t(m); xw=x(m);
end

function r=rms_trapz(t,x)
    r=sqrt(trapz(t,x.^2)/(t(end)-t(1)));
end

function Harm = fft_harmonics(t, x, f0)
% FFT-аналіз рівня 1, 3, 5 гармонік для нерівномірної сітки
% Вхід:
%   t – часова сітка (Nx1)
%   x – сигнал (Nx1)
%   f0 – фундаментальна частота
%
% Вихід:
%   Harm.H1, Harm.H3, Harm.H5 – амплітуди
%   Harm.H3p = % від фундаментали
%   Harm.H5p = % від фундаментали
%   Harm.f, Harm.X – спектр для графіка

    % --- інтерполяція на рівномірну сітку (потрібна для FFT)
    fs = 1/mean(diff(t));
    N  = length(t);

    % нова рівномірна сітка
    t_uni = linspace(t(1), t(end), N).';
    x_uni = interp1(t, x, t_uni, 'linear');

    % --- вікно Ганна для зменшення leakage
    w = hann(N);
    xw = x_uni .* w;

    % --- FFT
    X = fft(xw);
    Xmag = abs(X)/ (sum(w)/2);   % правильна нормалізація амплітуди
    f = (0:N-1).' * (fs/N);

    % --- знайти гармоніки
    [~, k1] = min(abs(f - f0));
    [~, k3] = min(abs(f - 3*f0));
    [~, k5] = min(abs(f - 5*f0));

    H1 = Xmag(k1);
    H3 = Xmag(k3);
    H5 = Xmag(k5);

    Harm.H1  = H1;
    Harm.H3  = H3;
    Harm.H5  = H5;

    Harm.H3p = 100 * H3/H1;
    Harm.H5p = 100 * H5/H1;

    Harm.f = f;
    Harm.X = Xmag;
end