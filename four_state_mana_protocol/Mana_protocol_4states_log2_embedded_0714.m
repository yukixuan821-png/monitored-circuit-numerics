%% Mana_protocol_4states_log2_embedded.m
% Corrected public-frame simulation for the monitored Clifford protocol.
%
% Goal:
%   Compare the exact Gross mana dynamics, using base-two logarithms,
%   for four different initial
%   ensembles at fixed d=3, N=5, theta_M=0.2:
%       (1) Haar-random state
%       (2) tensor product of qutrit T-type magic states
%       (3) exp(-i H_GUE t)|0>^{\otimes N}, with spectrum(H_GUE) in [-2,2]
%       (4) computational-basis product state |0>^{\otimes N}
%
% Protocol implemented:
%   Original one-cycle protocol:
%       U_C -> R_X^{(d)}(theta_M) on qudit 1
%       -> measure Z_1
%       -> U_C^\dagger.
%
%   Public-frame equivalent implementation:
%       sample the labelled Pauli pair
%           P_X = U_C^\dagger X_1 U_C,
%           P_Z = U_C^\dagger Z_1 U_C
%       induced by the same random Clifford frame, apply the qutrit
%       T-type small-angle kick along P_X, then measure P_Z by Born rule.
%
%   For a uniform projective Clifford draw, the phase-free labels
%   (P_X,P_Z) are uniformly distributed over ordered symplectic Pauli
%   pairs with <P_X,P_Z>=1.  The Pauli displacement part induces
%   independent uniform eigenvalue-label offsets for P_X and P_Z; these
%   offsets are sampled explicitly below.  This makes the public-frame
%   simulation match the labelled Clifford-frame protocol rather than the
%   weaker "fixed kick + independent random Pauli measurement" model.
%
% Output:
%   The script creates a local folder "figures" and saves:
%       .png  high-resolution image
%       .pdf  vector figure with embedded/vector text when supported
%       .eps  vector figure for LaTeX workflows
%       .fig  MATLAB figure
%       .mat  numerical data
%
% Notes:
%   - Mana convention: M(psi) = log_2 ||W_psi||_1.
%   - Exact Gross mana is expensive: the default T=2000, Nr=1000,
%     stride=1 reproduces the intended figure-quality experiment but can
%     take a long time.  For debugging, first try Nr=5, T=50.
%   - This file is self-contained; no external toolbox is required.

clear; clc;

%% ----------------- parameters -----------------
d       = 3;        % odd prime; this script uses the qutrit T-type gate
N       = 5;        % exact mana scales with d^(2N)
thetaM  = 0.2;      % magic-injection angle
Tsteps  = 2000;     % number of measurement cycles
Nr      = 1000;     % number of independent trajectories
stride  = 1;        % compute mana every 'stride' cycles
seed    = 1;        % RNG seed
tGUE    = 0.1;      % GUE evolution time
manaLogBase = 2;   % convention: M(psi)=log_2 ||W_psi||_1

assert(mod(d,2)==1 && isprime(d), 'Gross Wigner exact mana requires odd prime d.');
assert(d == 3, 'This script is written for the qutrit d=3 T-type gate.');
assert(mod(Tsteps, stride) == 0, 'Choose Tsteps divisible by stride.');

D     = d^N;
tgrid = 0:stride:Tsteps;
K     = numel(tgrid);

%% ----------------- caches -----------------
cache = qudit_cache(d);

%% ----------------- initial states -----------------
rng(seed, 'twister');

ket0N = zeros(D,1);
ket0N(1) = 1;

psi_init = cell(4,1);
init_names = { ...
    'Global Haar', ...
    '$|T\rangle^{\otimes N}$', ...
    '$e^{-iH_{\rm GUE}t}|0\rangle^{\otimes N}$', ...
    '$|0\rangle^{\otimes N}$'};

% (1) Global Haar-random state
psi_init{1} = haar_state(D);

% (2) Tensor product of qutrit T-type magic states:
%     |T> = (|0> + exp(2*pi*i/9)|1> + exp(-2*pi*i/9)|2>)/sqrt(3).
psi_init{2} = qutrit_T_product_state(N);

% (3) GUE-evolved |0>^{\otimes N}, with spectrum normalized to [-2,2].
psi_init{3} = gue_evolved_state(ket0N, D, tGUE);

% (4) Computational-basis product state |0>^{\otimes N}.
psi_init{4} = ket0N;

%% ----------------- run simulations -----------------
nInit    = numel(psi_init);
MmeanAll = zeros(nInit, K);
Mtraj1   = zeros(nInit, K);  % one representative trajectory for each initial state

fprintf('Public-frame monitored Clifford protocol\n');
fprintf('Mana convention: M(psi) = log_2 ||W_psi||_1\n');
fprintf('d=%d, N=%d, theta_M=%.6g, Tsteps=%d, Nr=%d, stride=%d\n', ...
    d, N, thetaM, Tsteps, Nr, stride);
fprintf('Each cycle: labelled random Pauli pair (P_X,P_Z), T-kick along P_X, Born measurement of P_Z.\n\n');

for s = 1:nInit
    fprintf('Initial ensemble (%d/%d): %s\n', s, nInit, regexprep(init_names{s}, '\$|\\', ''));

    % Use the same random stream for each initial-state block. This gives
    % comparable Clifford-frame randomness across the four initial states.
    rng(seed + 1000, 'twister');

    Msum = zeros(1, K);

    for r = 1:Nr
        psi = psi_init{s};

        ksave = 1;
        Mt = mana_gross_pure(psi, N, d);
        Msum(ksave) = Msum(ksave) + Mt;
        if r == 1
            Mtraj1(s, ksave) = Mt;
        end

        for t = 1:Tsteps
            % Sample a labelled public Clifford frame:
            %   P_X = omega^xOffset P(xLabel),
            %   P_Z = omega^zOffset P(zLabel),
            % with symplectic bracket <xLabel,zLabel>=1.
            [xA, xB, zA, zB, xOffset, zOffset] = sample_labelled_hyperbolic_pair(N, d);

            % Apply the qutrit T-type small-angle kick along P_X.
            psi = apply_qutrit_T_axis_kick(psi, xA, xB, xOffset, thetaM, N, d, cache);

            % Measure P_Z by Born rule and update the pure state.
            [psi, ~] = measure_pauli_spectral(psi, zA, zB, zOffset, N, d, cache);

            if mod(t, stride) == 0
                ksave = ksave + 1;
                Mt = mana_gross_pure(psi, N, d);
                Msum(ksave) = Msum(ksave) + Mt;
                if r == 1
                    Mtraj1(s, ksave) = Mt;
                end
            end
        end

        if mod(r, max(1, round(Nr/10))) == 0
            fprintf('  done %d/%d trajectories\n', r, Nr);
        end
    end

    MmeanAll(s,:) = Msum / Nr;
    fprintf('\n');
end

%% ----------------- plot -----------------
% Use a vector-friendly figure setup for manuscript insertion.
% The PDF/EPS outputs below use the painters renderer so that labels and
% curves remain vector objects whenever supported by the MATLAB backend.
fig = figure('Color','w', 'Units','inches', 'Position',[1 1 6.8 5.0]);
set(fig, 'Renderer', 'painters', 'InvertHardcopy', 'off', 'PaperPositionMode', 'auto');
ax = axes(fig);
hold(ax, 'on');

set(ax, 'XScale', 'log', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 13, ...
    'LineWidth', 0.9, ...
    'TickLabelInterpreter', 'latex');
cols = lines(nInit);

% Light representative trajectories.
for s = 1:nInit
    c_light = lighten_color(cols(s,:), 0.25);
    plot(ax, tgrid+1, Mtraj1(s,:), ...
        'Color', c_light, 'LineWidth', 1.0, 'HandleVisibility','off');
end

% Dark ensemble means.
for s = 1:nInit
    plot(ax, tgrid+1, MmeanAll(s,:), ...
        'Color', cols(s,:), 'LineWidth', 2.0, ...
        'DisplayName', init_names{s});
end

xlabel(ax, 'measurement steps $t+1$', 'Interpreter','latex');
ylabel(ax, 'exact mana $\mathcal{M}(\psi_t)=\log_{2}\|W_{\psi_t}\|_1$', 'Interpreter','latex');

title(ax, sprintf('$d=%d$, $N=%d$, $\\theta_M=%.3g$, $N_r=%d$, Monitored Quantum Circuits', ...
    d, N, thetaM, Nr), 'Interpreter','latex');

legend(ax, 'Location','best', 'Interpreter','latex', 'FontSize', 12);
grid(ax, 'on');
box(ax, 'on');
xlim(ax, [1, Tsteps+1]);

ymax = max([MmeanAll(:); Mtraj1(:)]);
if ~isfinite(ymax) || ymax <= 0
    ymax = 1;
end
ylim(ax, [0, 1.05*ymax]);

%% ----------------- save outputs -----------------
outdir = fullfile(pwd, 'figures');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

tag = sprintf('Mana_log2_protocol_d%d_N%d_theta_%s_Nr%d', ...
    d, N, strrep(sprintf('%.3g',thetaM),'.','p'), Nr);

pngFile = fullfile(outdir, [tag '.png']);
pdfFile = fullfile(outdir, [tag '.pdf']);
epsFile = fullfile(outdir, [tag '.eps']);
figFile = fullfile(outdir, [tag '.fig']);
matFile = fullfile(outdir, [tag '.mat']);

try
    exportgraphics(fig, pngFile, 'Resolution', 600, 'BackgroundColor', 'white');
catch
    print(fig, pngFile, '-dpng', '-r600');
end

try
    exportgraphics(fig, pdfFile, 'ContentType', 'vector', 'BackgroundColor', 'white');
catch
    print(fig, pdfFile, '-dpdf', '-painters');
end

try
    print(fig, epsFile, '-depsc', '-painters');
catch
    warning('Could not export EPS file. PDF and PNG outputs were still attempted.');
end

savefig(fig, figFile);
save(matFile, 'd', 'N', 'thetaM', 'Tsteps', 'Nr', 'stride', 'seed', 'tGUE', 'manaLogBase', ...
    'tgrid', 'init_names', 'MmeanAll', 'Mtraj1');

fprintf('Saved figure/data files in:\n  %s\n', outdir);
fprintf('PNG: %s\n', pngFile);
fprintf('PDF: %s\n', pdfFile);
fprintf('EPS: %s\n', epsFile);
fprintf('FIG: %s\n', figFile);
fprintf('MAT: %s\n', matFile);

%% ========================= local functions =========================

function cache = qudit_cache(d)
    omega = exp(2*pi*1i/d);

    % X|j> = |j+1>, Z|j> = omega^j |j>, j=0,...,d-1.
    X = zeros(d,d);
    for j = 0:d-1
        X(mod(j+1,d)+1, j+1) = 1;
    end
    Z = diag(omega.^(0:d-1));

    % Fourier transform F_{jk} = omega^{jk}/sqrt(d).
    j = (0:d-1).';
    k = 0:d-1;
    F = omega.^(j*k) / sqrt(d);

    Xp = cell(d,1);
    Zp = cell(d,1);
    Xp{1} = eye(d);
    Zp{1} = eye(d);
    for a = 1:d-1
        Xp{a+1} = Xp{a} * X;
        Zp{a+1} = Zp{a} * Z;
    end

    cache = struct();
    cache.d = d;
    cache.omega = omega;
    cache.X = X;
    cache.Z = Z;
    cache.F = F;
    cache.Xp = Xp;
    cache.Zp = Zp;
end

function psi = haar_state(D)
    psi = randn(D,1) + 1i*randn(D,1);
    psi = psi / norm(psi);
end

function psi = qutrit_T_product_state(N)
    % |T> = (|0> + exp(2*pi*i/9)|1> + exp(-2*pi*i/9)|2>)/sqrt(3).
    v = [1; exp(2*pi*1i/9); exp(-2*pi*1i/9)] / sqrt(3);
    psi = 1;
    for q = 1:N
        psi = kron(psi, v);
    end
    psi = psi(:);
    psi = psi / norm(psi);
end

function psi = gue_evolved_state(ket0N, D, t)
    % Draw a GUE-like Hermitian matrix and affine-normalize its eigenvalues
    % to lie in [-2,2], then compute exp(-i H t)|0...0>.
    X = randn(D,D) + 1i*randn(D,D);
    H = (X + X')/2;

    [V, lam] = eig(H, 'vector');
    lam = real(lam);

    lam_min = min(lam);
    lam_max = max(lam);
    width   = lam_max - lam_min;

    if width <= 0 || ~isfinite(width)
        psi = expm(-1i*H*t) * ket0N;
        psi = psi / norm(psi);
        return;
    end

    center = (lam_max + lam_min)/2;
    scale  = width/4;             % maps endpoints to ±2
    lamN   = (lam - center) / scale;

    psi = V * (exp(-1i * lamN * t) .* (V' * ket0N));
    psi = psi / norm(psi);
end

function [xA, xB, zA, zB, xOffset, zOffset] = sample_labelled_hyperbolic_pair(N, d)
    % Uniformly sample an ordered phase-free Pauli pair (x,z) satisfying
    % symplectic bracket <x,z> = 1, then sample independent eigenvalue
    % offsets in F_d.
    %
    % Label convention: P(a,b)=X^a Z^b on each qudit.
    % Symplectic bracket:
    %   <x,z> = a_x · b_z - a_z · b_x  mod d,
    % so that P_z P_x = omega^{<x,z>} P_x P_z.

    [xA, xB] = random_pauli_label_uniform_nonzero(N, d);

    while true
        zA = randi([0, d-1], 1, N);
        zB = randi([0, d-1], 1, N);
        if symplectic_bracket(xA, xB, zA, zB, d) == 1
            break;
        end
    end

    xOffset = randi([0, d-1]);
    zOffset = randi([0, d-1]);
end

function [aVec, bVec] = random_pauli_label_uniform_nonzero(N, d)
    % Uniform over (F_d)^{2N}\{0}, without rejection.
    M = d^(2*N);
    u = randi([1, M-1]);  % exclude the all-zero label

    digits = zeros(1, 2*N);
    for k = 1:(2*N)
        digits(k) = mod(u, d);
        u = floor(u / d);
    end

    aVec = digits(1:N);
    bVec = digits(N+1:end);
end

function val = symplectic_bracket(a1, b1, a2, b2, d)
    val = mod(sum(a1 .* b2 - a2 .* b1), d);
end

function psi = apply_qutrit_T_axis_kick(psi, aVec, bVec, phaseOffset, thetaM, N, d, cache)
    % Apply the qutrit T-type small-angle gate along the spectral
    % decomposition of the labelled Pauli P = omega^phaseOffset X^a Z^b.
    %
    % For d=3:
    %   D(theta) = diag(1, exp(i theta 2pi/9), exp(-i theta 2pi/9)).
    %
    % If Pi_m is the eigenspace projector of P with eigenvalue omega^m,
    % then the kick is sum_m phase_m Pi_m.

    assert(d == 3, 'apply_qutrit_T_axis_kick is written for d=3.');

    phaseVec = [1, exp(1i * thetaM * 2*pi/9), exp(-1i * thetaM * 2*pi/9)];

    components = pauli_spectral_components(psi, aVec, bVec, phaseOffset, N, d, cache);

    psi2 = zeros(size(psi));
    for m = 1:d
        psi2 = psi2 + phaseVec(m) * components(:,m);
    end

    % Numerically renormalize; the operator is unitary up to roundoff.
    psi = psi2 / norm(psi2);
end

function [psi, outcome] = measure_pauli_spectral(psi, aVec, bVec, phaseOffset, N, d, cache)
    % Projective measurement of P = omega^phaseOffset X^a Z^b.
    % Outcome m corresponds to eigenvalue omega^m, m=0,...,d-1.

    components = pauli_spectral_components(psi, aVec, bVec, phaseOffset, N, d, cache);
    probs = sum(abs(components).^2, 1);
    probs = real(probs);
    probs = max(probs, 0);

    psum = sum(probs);
    if psum <= 0 || ~isfinite(psum)
        error('Invalid Born probabilities encountered.');
    end
    probs = probs / psum;

    outcome = sample_discrete(probs);  % 0,...,d-1
    pb = probs(outcome+1);

    psi = components(:, outcome+1) / sqrt(max(pb, realmin));
    psi = psi / norm(psi);
end

function components = pauli_spectral_components(psi, aVec, bVec, phaseOffset, N, d, cache)
    % components(:,m+1) = Pi_m psi, where Pi_m is the spectral projector
    % of P = omega^phaseOffset X^a Z^b with eigenvalue omega^m.
    %
    % Pi_m psi = (1/d) sum_{s=0}^{d-1} omega^{-m s} P^s psi.
    % MATLAB fft along the second dimension implements this convention.

    D = numel(psi);
    Phi = complex(zeros(D, d));

    Phi(:,1) = psi;
    for s = 2:d
        Phi(:,s) = apply_pauli_tensor(Phi(:,s-1), aVec, bVec, N, d, cache);
        Phi(:,s) = cache.omega^phaseOffset * Phi(:,s);
    end

    components = fft(Phi, [], 2) / d;
end

function psi = apply_pauli_tensor(psi, aVec, bVec, N, d, cache)
    % Apply tensor_q X^{a_q} Z^{b_q}.
    for q = 1:N
        a = aVec(q);
        b = bVec(q);
        if a == 0 && b == 0
            continue;
        end
        U = cache.Xp{a+1} * cache.Zp{b+1};
        psi = apply1(psi, U, q, N, d);
    end
end

function psi = apply1(psi, U, q, N, d)
    % Apply a single-qudit matrix U to qudit q of an N-qudit vector.
    dims = d * ones(1,N);
    Psi  = reshape(psi, dims);

    rest  = [1:q-1, q+1:N];
    order = [q, rest];
    PsiP  = permute(Psi, order);

    PsiMat = reshape(PsiP, d, []);
    PsiMat = U * PsiMat;

    PsiP = reshape(PsiMat, [d, dims(rest)]);
    invOrder = invert_permutation(order);
    Psi = permute(PsiP, invOrder);

    psi = Psi(:);
end

function invOrder = invert_permutation(order)
    invOrder = zeros(size(order));
    for i = 1:numel(order)
        invOrder(order(i)) = i;
    end
end

function k0 = sample_discrete(p)
    r = rand();
    c = cumsum(p);
    k = find(r <= c, 1, 'first');
    if isempty(k)
        k = numel(p);  % guard against roundoff in cumsum
    end
    k0 = k - 1;
end

function c2 = lighten_color(c, amount)
    % amount in [0,1]: 0 -> white, 1 -> original color.
    amount = max(0, min(1, amount));
    c2 = 1 - (1 - c) * amount;
end

function M = mana_gross_pure(psi, N, d)
    % Exact Gross Wigner mana for a pure state in odd prime local dimension:
    %   M(psi) = log_2 ||W_psi||_1.
    persistent C
    if isempty(C) || C.d ~= d || C.N ~= N
        C = precompute_mana_cache(N, d);
    end

    dims = C.dims;
    D    = C.D;

    Psi = reshape(psi, dims);

    % G(q,x) = psi(q+x/2) conj(psi(q-x/2)), represented by flattened q
    % and flattened x.  Here x/2 is implemented modulo d.
    G = complex(zeros(D, D));
    for idx = 1:D
        s = C.shift(idx,:);
        PsiPlus  = circshift(Psi, -s);
        PsiMinus = circshift(Psi, +s);
        G(:, idx) = PsiPlus(:) .* conj(PsiMinus(:));
    end

    % Fourier transform over x variables. MATLAB ifft contributes the
    % normalization d^{-N}; the sign convention does not affect ||W||_1.
    W = reshape(G, [D, dims]);  % [q_flat, x1, ..., xN]
    for k = 1:N
        W = ifft(W, [], 1+k);
    end

    Wvec = W(:);
    M = log2(sum(abs(Wvec)));  % base-two mana, in bits
end

function C = precompute_mana_cache(N, d)
    assert(mod(d,2)==1 && isprime(d), 'Gross Wigner cache requires odd prime d.');

    inv2 = modInverse(2, d);
    dims = d * ones(1,N);
    D    = d^N;

    shift = zeros(D, N);
    subs  = cell(1,N);

    for idx = 1:D
        [subs{:}] = ind2sub(dims, idx);
        xvec = zeros(1,N);
        for k = 1:N
            xvec(k) = subs{k} - 1;
        end
        shift(idx,:) = mod(inv2 * xvec, d);
    end

    C = struct();
    C.N = N;
    C.d = d;
    C.inv2 = inv2;
    C.dims = dims;
    C.D = D;
    C.shift = shift;
end

function x = modInverse(a, m)
    a = int64(a);
    m = int64(m);
    [g,u,~] = gcd(a, m);
    if g ~= 1
        error('modInverse: a=%d has no inverse mod m=%d', a, m);
    end
    x = double(mod(u, m));
end
