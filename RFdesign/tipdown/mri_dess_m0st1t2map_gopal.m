  function [M0s_mom, M0s_med, M0s_reg, T1_reg, T2_reg, wf_reg, cost] = ...
      mri_dess_m0st1t2map(E1_init, E2_init, wf_init, flip_im,...
      yp_im, ym_im, mask, T1max, T1min, T2max, T2min, TR, TE, Rm, R1, R2,...
      n_outer, n_innerM, n_inner1, n_inner2, tolm, tol1, tol2, T2meth, disp)
% function [M0s_mom, M0s_med, M0s_reg, T1_reg, T2_reg, wf_reg, cost] = ...
%     mri_dess_m0st1t2map(E1_init, E2_init, wf_init, flip_im,...
%     yp_im, ym_im, mask, T1max, T1min, T2max, T2min, TR, TE, Rm, R1, R2,...
%     n_outer, n_innerM, n_inner1, n_inner2, tolm, tol1, tol2, T2meth, disp)
%
% Inputs: 
%   E1_init     [nx ny nz]          E1 initial guess
%   E2_init     [nx ny nz]          E2 initial guess
%   wf_init     [nx ny nz]          wf (B0 inhomogeneity) initial guess
%   flip_im     [nx ny nz nf]       Flip angle maps
%   yp_im       [nx ny nz nf]       DESS first echo data
%   ym_im       [nx ny nz nf]       DESS second echo data
%   mask        [nx ny nz]          Data mask to operate on
%   T1max       [1]                 Maximum T1 expected (lower --> faster)
%   T1min       [1]                 Minimum T1 expected (higher --> faster)
%   T2max       [1]                 Maximum T2 expected (lower --> faster)
%   T2min       [1]                 Minimum T2 expected (higher --> faster)
%   TR          [1]                 Repetition time
%   TE          [1]                 Echo time
%   Rm          {1} strum           M0s Regularizer Object
%   R1          {1} strum           E1 Regularizer Object
%   R2          {1} strum           E2 Regularizer Object
%   n_outer     [1]                 Number of outer iterations
%   n_innerM    [1]                 Max number of M0s inner iterations
%   n_inner1    [1]                 Max number of E1 inner iterations
%   n_inner2    [1]                 Max number of E2 inner iterations
%   tolm        [1]                 M0s stop tolerance
%   tol1        [1]                 E1 stop tolerance
%   tol2        [1]                 E2 stop tolerance
%   T2meth      [1] char            Method of T2 update
%   disp        [1]                 Toggle display on/off of outer iterates
%
% Outputs:
%   M0s_mom     [nx ny nz]          M0s initial guess
%   M0s_med     [nx ny nz]          M0s median-filtered initial guess
%   M0s_reg     [nx ny nz]          Output M0star map
%   T1_reg      [nx ny nz]          Output T1 map
%   T2_reg      [nx ny nz]          Output T2 map
%   wf_reg      [nx ny nz]          Output wf map
%   cost        [n_outer+1 1]       Cost function vs. outer iteration
% 
% Written By: Gopal Nataraj and Jeffrey A. Fessler; Copyright 2014

% Datafit cost function
datacost = @(y, f) (1/2) * norm(y - f).^2;

% Initial estimates
E1 = E1_init(mask);
E2 = E2_init(mask);
wf = wf_init(mask); 
np = length(E1);

% Thresholding values
E1max = exp(-TR/T1max);
E1min = exp(-TR/T1min);
E2max = exp(-TR/T2max);
E2min = exp(-TR/T2min);

% Vectorize DESS and flip angle data
nf = size(flip_im, 4);
yp = reshape(yp_im(repmat(mask, [1 1 1 nf])), [np nf]);
ym = reshape(ym_im(repmat(mask, [1 1 1 nf])), [np nf]);
flip = reshape(flip_im(repmat(mask, [1 1 1 nf])), [np nf]);

% Weighting matrix, W 
W = Gdiag(col(repmat(ones(1,nf), [2*length(E1) 1])));

for outer = 1:n_outer
    % Display Outer Iteration
    printm('Outer Iteration: %d of %d', outer, n_outer);
    
    %% M0s update
    % Update the system matrices
    tmp_p = cell(nf, 1);
    tmp_m = cell(nf, 1);
    for a = 1:nf
        tmp_p{a} = Gdiag(dfp_M0s(E1, E2, flip(:,a), TR, TE, wf));
        tmp_m{a} = Gdiag(dfm_M0s(E1, E2, flip(:,a), TR, TE, wf));
    end
    Am = [vertcat(tmp_p{:}); vertcat(tmp_m{:})];
    
    % On the first outer iteration, get a better initial M0s estimate
    if (outer == 1)
        M0s = ((Am'*W) * col([yp ym])) ./ ((Am'*W*Am) * ones(size(Am,2), 1));
        M0s_mom = embed(M0s, mask);
        if (disp)
            figure(1); im(abs(M0s_mom), 'cbar'); drawnow;
        end
        
        % M0s_med median filtering, for comparison only (might not be right!)
        M0s_med = medfilt2(real(M0s_mom)) + 1i*medfilt2(imag(M0s_mom));
        
        % Cost function initialization
        cost = zeros(3*n_outer+1, 1);
        cost_idx = 1;
        cost(cost_idx) = Rm.penal(Rm, M0s) + R1.penal(R1, E1)+ R2.penal(R2, E2);
        for a = 1:nf
            cost(cost_idx) = cost(cost_idx) ...
                + datacost(yp(:,a), fp(M0s, E1, E2, flip(:,a), TR, TE, wf)) ...
                + datacost(ym(:,a), fm(M0s, E1, E2, flip(:,a), TR, TE, wf));
        end
        cost_idx = cost_idx + 1;
    end
    
    % Update preconditioning matrix
    Dm = (Am' * Am) + Gdiag(Rm.denom(Rm, M0s));
    Pm = Gdiag(1 ./ (Dm * ones(size(Dm,1), 1)));
    
    % Preconditioned conjugate gradient (equiv. to WLS for beta_m = 0)
	M0s = pwls_pcg1(M0s, Am, W, col([yp ym]), Rm, 'niter', n_innerM,...
        'precon', Pm, 'stop_diff_tol', tolm, 'chat', 1); 
    
    % Cost function evaluation (after M0s update)
    cost(cost_idx) = Rm.penal(Rm, M0s) + R1.penal(R1, E1)+ R2.penal(R2, E2);
    for a = 1:nf
        cost(cost_idx) = cost(cost_idx) ...
            + datacost(yp(:,a), fp(M0s, E1, E2, flip(:,a), TR, TE, wf)) ...
            + datacost(ym(:,a), fm(M0s, E1, E2, flip(:,a), TR, TE, wf));
    end
    cost_idx = cost_idx + 1;
    
    %% E1 update
    % BBGM parameters
    alpha = 1; gam = 1e-4; eps = 1e-10; M = 10; sig1 = 0.1; sig2 = 0.5; 
    
    % Inner cost function initialization
    cost_E1 = zeros(n_inner1+1, 1);
    cost_E1(1) = R1.penal(R1, E1);
    for a = 1:nf
        cost_E1(1) = cost_E1(1) ...
            + datacost(yp(:,a), fp(M0s, E1, E2, flip(:,a), TR, TE, wf)) ...
            + datacost(ym(:,a), fm(M0s, E1, E2, flip(:,a), TR, TE, wf));
    end
    
    % Compute initial gradient
    g = R1.cgrad(R1, E1);
    for a = 1:nf
        g_p = real(conj(fp(M0s, E1, E2, flip(:,a), TR, TE, wf)...
            - yp(:,a)) .* dfp_E1(M0s, E1, E2, flip(:,a), TR, TE, wf));
        g_m = real(conj(fm(M0s, E1, E2, flip(:,a), TR, TE, wf)...
            - ym(:,a)) .* dfm_E1(M0s, E1, E2, flip(:,a), TR, TE, wf));
        g = g + g_p + g_m;
    end
    
    for inner1 = 1:n_inner1
        % Rescale alpha
        if ((alpha <= eps) || (alpha >= (1/eps)))
            alpha = max(1, 1/norm(g));
            alpha = min(10e5, alpha);
        end
        
        % Update lambda, step size
        lam = 1/alpha;
        
        % Non-monotonic line-search
        while (true)
            % Next candidate solution
            E1_next = E1 - lam * g;
           
        	% Compute and store cost
            cost_E1(inner1+1) = R1.penal(R1, E1_next);
            for a = 1:nf
                cost_E1(inner1+1) = cost_E1(inner1+1) ...
                    + datacost(yp(:,a), fp(M0s, E1_next, E2, flip(:,a), TR, TE, wf)) ...
                    + datacost(ym(:,a), fm(M0s, E1_next, E2, flip(:,a), TR, TE, wf));
            end
            
            % If line search criterion is satisfied, update the solution
            cost_rel = cost_E1(inner1-min(inner1-1,M-1) : inner1);
            if (cost_E1(inner1+1) <= max(cost_rel) - (gam * lam * (g' * g)))
                E1_prev = E1;
                E1 = E1_next;
                
                % Update gradient
                g_prev = g;
                g = R1.cgrad(R1, E1);
                for a = 1:nf
                    g_p = real(conj(fp(M0s, E1, E2, flip(:,a), TR, TE, wf)...
                        - yp(:,a)) .* dfp_E1(M0s, E1, E2, flip(:,a), TR, TE, wf));
                    g_m = real(conj(fm(M0s, E1, E2, flip(:,a), TR, TE, wf)...
                        - ym(:,a)) .* dfm_E1(M0s, E1, E2, flip(:,a), TR, TE, wf));
                    g = g + g_p + g_m;
                end
                
                % Update inverse step size
                alpha = -(g_prev' * (g-g_prev)) / (lam * (g_prev' * g_prev));
                
                % Project between [E1min, E1max]
                E1 = max(E1, E1min);
                E1 = min(E1, E1max);
                
                % Exit line search
                break;
            % Otherwise, rescale the step size and try again
            else
                % Choose sig, rescaling parameter, and rescale
                sig = mean([sig1 sig2]);
                lam = sig*lam;
                
                % Error catch: if lam < machine precision, exit early
                if (lam < 10^-32)
                    E1_prev = E1;
                    break;
                end
            end 
        end
        
        % Display progress
        if (rem(inner1,5)==0)
            printm('E1 BBGM: %d of %d', inner1, n_inner1); 
        end
        
        % Exit early if change in solution < tol1
        change = norm(E1 - E1_prev) / norm(E1); 
        if (change < tol1)
            printm('Exited at inner E1 iteration %u of %u', inner1, n_inner1);
            break;
        end
    end
    
    % Cost function evaluation (after E1 update)
    cost(cost_idx) = Rm.penal(Rm, M0s) + R1.penal(R1, E1)+ R2.penal(R2, E2);
    for a = 1:nf
        cost(cost_idx) = cost(cost_idx) ...
            + datacost(yp(:,a), fp(M0s, E1, E2, flip(:,a), TR, TE, wf)) ...
            + datacost(ym(:,a), fm(M0s, E1, E2, flip(:,a), TR, TE, wf));
    end
    cost_idx = cost_idx + 1;
    
    %% E2 update
    if (strcmp(T2meth,'BBGM'))
        % METHOD ONE: BBGM (FAST, BUT NONMONOTONIC)
        % BBGM parameters
        alpha = 1; gam = 1e-4; eps = 1e-10; M = 10; sig1 = 0.1; sig2 = 0.5; 

        % Inner cost function initialization
        cost_E2 = zeros(n_inner2+1, 1);
        cost_E2(1) = R2.penal(R2, E2);
        for a = 1:nf
            cost_E2(1) = cost_E2(1) ...
                + datacost(yp(:,a), fp(M0s, E1, E2, flip(:,a), TR, TE, wf)) ...
                + datacost(ym(:,a), fm(M0s, E1, E2, flip(:,a), TR, TE, wf));
        end

        % Compute initial gradient
        g = R2.cgrad(R2, E2);
        for a = 1:nf
            g_p = real(conj(fp(M0s, E1, E2, flip(:,a), TR, TE, wf)...
                - yp(:,a)) .* dfp_E2(M0s, E1, E2, flip(:,a), TR, TE, wf));
            g_m = real(conj(fm(M0s, E1, E2, flip(:,a), TR, TE, wf)...
                - ym(:,a)) .* dfm_E2(M0s, E1, E2, flip(:,a), TR, TE, wf));
            g = g + g_p + g_m;
        end

        for inner2 = 1:n_inner2
            % Rescale alpha
            if ((alpha <= eps) || (alpha >= (1/eps)))
                alpha = max(1, 1/norm(g));
                alpha = min(10e5, alpha);
            end

            % Update lambda, step size
            lam = 1/alpha;

            % Non-monotonic line-search
            while (true)
                % Next candidate solution
                E2_next = E2 - lam * g;

                % Compute and store cost
                cost_E2(inner2+1) = R2.penal(R2, E2_next);
                for a = 1:nf
                    cost_E2(inner2+1) = cost_E2(inner2+1) ...
                        + datacost(yp(:,a), fp(M0s, E1, E2_next, flip(:,a), TR, TE, wf)) ...
                        + datacost(ym(:,a), fm(M0s, E1, E2_next, flip(:,a), TR, TE, wf));
                end

                % If line search criterion is satisfied, update the solution
                cost_rel = cost_E2(inner2-min(inner2-1,M-1) : inner2);
                if (cost_E2(inner2+1) <= max(cost_rel) - (gam * lam * g' * g))
                    E2_prev = E2;
                    E2 = E2_next;

                    % Update gradient
                    g_prev = g;
                    g = R2.cgrad(R2, E2);
                    for a = 1:nf
                        g_p = real(conj(fp(M0s, E1, E2, flip(:,a), TR, TE, wf)...
                            - yp(:,a)) .* dfp_E2(M0s, E1, E2, flip(:,a), TR, TE, wf));
                        g_m = real(conj(fm(M0s, E1, E2, flip(:,a), TR, TE, wf)...
                            - ym(:,a)) .* dfm_E2(M0s, E1, E2, flip(:,a), TR, TE, wf));
                        g = g + g_p + g_m;
                    end

                    % Update inverse step size
                    alpha = -(g_prev' * (g-g_prev)) / (lam * (g_prev' * g_prev));

                    % Project between [E2min, E2max]
                    E2 = max(E2, E2min);
                    E2 = min(E2, E2max);
                    
                    % Exit line search
                    break;
                % Otherwise, rescale the step size and try again
                else
                    % Choose sig, rescaling parameter, and rescale
                    sig = mean([sig1 sig2]);
                    lam = sig*lam;

                    %Error catch: if lam < machine precision, exit early
                    if (lam < 10^-32)
                        E2_prev = E2;
                        break;
                    end
                end 
            end

            % Display progress
            if (rem(inner2,5)==0)
                printm('E2 BBGM: %d of %d', inner2, n_inner2); 
            end

            % Exit early if change in solution < tol2
            change = norm(E2 - E2_prev) / norm(E2); 
            if (change < tol2)
                printm('Exited at inner E2 iteration %u of %u', inner2, n_inner2);
                break;
            end
        end
    elseif (strcmp(T2meth,'PGD'))
        % METHOD TWO: PGD (SLOW, BUT MONOTONIC)
        % After a few iterations, tighten bound of [E2min, E2max] 
        if (outer > 10)
            if (min(E2) > E2min)
                fprintf('E2min changed from %0.4f to %0.4f at iter %u.\n',...
                    E2min, min(E2), outer);
                E2min = min(E2);
            end
            if (max(E2) < E2max)
                fprintf('E2max changed from %0.4f to %0.4f at iter %u.\n',...
                    E2max, max(E2), outer);
                E2max = max(E2);
            end
        end

        % Update preconditioner on outer iteration
        d2 = R2.denom(R2, E2);
        for a = 1:nf
            % Compute derivatives of cost function at endpoints 
            d2_pmax = real(conj(fp(M0s, E1, E2max, flip(:,a), TR, TE, wf)...
                - yp(:,a)) .* ddfp_E2(M0s, E1, E2max, flip(:,a), TR, TE, wf))...
                + abs(dfp_E2(M0s, E1, E2max, flip(:,a), TR, TE, wf) .^ 2);
            d2_pmin = real(conj(fp(M0s, E1, E2min, flip(:,a), TR, TE, wf)...
                - yp(:,a)) .* ddfp_E2(M0s, E1, E2min, flip(:,a), TR, TE, wf))...
                + abs(dfp_E2(M0s, E1, E2min, flip(:,a), TR, TE, wf) .^ 2);
            d2_mmax = real(conj(fm(M0s, E1, E2max, flip(:,a), TR, TE, wf)...
                - ym(:,a)) .* ddfm_E2(M0s, E1, E2max, flip(:,a), TR, TE, wf))...
                + abs(dfm_E2(M0s, E1, E2max, flip(:,a), TR, TE, wf) .^ 2);
            d2_mmin = real(conj(fm(M0s, E1, E2min, flip(:,a), TR, TE, wf)...
                - ym(:,a)) .* ddfm_E2(M0s, E1, E2min, flip(:,a), TR, TE, wf))...
                + abs(dfm_E2(M0s, E1, E2min, flip(:,a), TR, TE, wf) .^ 2);

            % Add the maximum curvature 
            d2 = d2 + max(d2_pmax, d2_pmin) + max(d2_mmax, d2_mmin);
        end
        P2 = Gdiag(1 ./ d2);

        % PGD with Nesterov Momentum for E2 update
        % Initialize Nesterov auxiliary variable and momentum parameter
        z2 = E2; 
        t = 1;

        for inner2 = 1:n_inner2
            % Compute the next momentum parameter
            t_next = (1 + sqrt(1 + 4*t^2))/2;

            % Compute gradient at z2
            grad = R2.cgrad(R2, z2);
            for a = 1:nf
                grad_p = real(conj(fp(M0s, E1, z2, flip(:,a), TR, TE, wf)...
                    - yp(:,a)) .* dfp_E2(M0s, E1, z2, flip(:,a), TR, TE, wf));
                grad_m = real(conj(fm(M0s, E1, z2, flip(:,a), TR, TE, wf)...
                    - ym(:,a)) .* dfm_E2(M0s, E1, z2, flip(:,a), TR, TE, wf));
                grad = grad + grad_p + grad_m;
            end

            % PGD over z2 is stored in E2
            E2_prev = E2;
            E2 = z2 - (P2 * grad);

            % Project between [E2min, E2max]
            E2 = max(E2, E2min);
            E2 = min(E2, E2max);

            % Momentum over E2 is stored in z2
            z2 = E2 + ((t-1)/t_next) * (E2 - E2_prev);

            % Update the momentum parameter
            t = t_next;

            % Display progress
            if (rem(inner2,5)==0)
                printm('E2 PGD: %d of %d', inner2, n_inner2); 
            end

            % Exit early if change in solution < tol2
            change = norm(E2 - E2_prev) / norm(E2); 
            if (change < tol2)
                printm('Exited at inner E2 iteration %u of %u', inner2, n_inner2);
                break;
            end
        end
    end
    
    % Cost function evaluation (after E2 update)
    cost(cost_idx) = Rm.penal(Rm, M0s) + R1.penal(R1, E1) + R2.penal(R2, E2);
    for a = 1:nf
        cost(cost_idx) = cost(cost_idx) ...
            + datacost(yp(:,a), fp(M0s, E1, E2, flip(:,a), TR, TE, wf)) ...
            + datacost(ym(:,a), fm(M0s, E1, E2, flip(:,a), TR, TE, wf));
    end
    cost_idx = cost_idx + 1;
    
    %% Observe images over outer iteration
    if (disp)
        figure(1); im(embed(abs(M0s), mask), 'cbar'); drawnow;
        figure(2); im(embed(-TR./log(E1), mask), [0 3000], 'cbar'); drawnow;
        figure(3); im(embed(-TR./log(E2), mask), [0 200], 'cbar'); drawnow;
    end
end

%% Extract M0s, T1, T2, wf estimates
M0s_reg = embed(M0s, mask);
T1_reg = embed(-TR./log(E1), mask);
T2_reg = embed(-TR./log(E2), mask);
wf_reg = embed(wf, mask);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DESS signal models
function sp = fp(M0s, E1, E2, a, ~, TE, wf)
sp = M0s.*tan(a.*(1.0./2.0)).*exp(pi.*5.0e-1i+TE.*wf.*1i).*...
    ((sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*1.0./(E1.*cos(a)-1.0).^...
    2-1.0)).*(E1-cos(a)))./(E1.*cos(a)-1.0)+1.0);
sp(isnan(sp)) = 0;
end

function sm = fm(M0s, E1, E2, a, TR, TE, wf) 
sm = -E2.^((TE.*-2.0)./TR).*M0s.*tan(a.*(1.0./2.0)).*exp(pi.*-5.0e-1i-TE.*...
    wf.*1i).*(sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*1.0./...
    (E1.*cos(a)-1.0).^2-1.0))-1.0);
sm(isnan(sm)) = 0;
end

%% DESS first derivatives w.r.t. M0s
function dsp_M0s = dfp_M0s(E1, E2, a, TR, TE, wf)
dsp_M0s = tan(a.*(1.0./2.0)).*exp(TE.*wf.*1i).*...
    ((sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*1.0./(E1.*cos(a)-1.0).^...
    2-1.0)).*(E1-cos(a)))./(E1.*cos(a)-1.0)+1.0).*1i;
dsp_M0s(isnan(dsp_M0s)) = 0;
end

function dsm_M0s = dfm_M0s(E1, E2, a, TR, TE, wf)
dsm_M0s = E2.^((TE.*-2.0)./TR).*tan(a.*(1.0./2.0)).*exp(TE.*wf.*-1i).*...
    (sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*1.0./...
    (E1.*cos(a)-1.0).^2-1.0))-1.0).*1i;
dsm_M0s(isnan(dsm_M0s)) = 0;
end

%% DESS first derivatives w.r.t. E1 
function dsp_E1 = dfp_E1(M0s, E1, E2, a, TR, TE, wf)
dsp_E1 = M0s.*tan(a.*(1.0./2.0)).*exp(TE.*wf.*1i).*(E2.^2-1.0).*1.0./...
    sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*1.0./(E1.*cos(a)-1.0).^...
    2-1.0)).*(cos(a).^2-1.0).*(E1.*cos(a)-1.0).^2.*1.0./(E1.^2.*cos(a).^...
    2-E2.^2.*cos(a).^2-E1.*cos(a).*2.0-E1.^2.*E2.^2+E1.*E2.^2.*...
    cos(a).*2.0+1.0).^2.*-1i;
dsp_E1(isnan(dsp_E1)) = 0;
end

function dsm_E1 = dfm_E1(M0s, E1, E2, a, TR, TE, wf)
dsm_E1 = E2.^((TE.*-2.0)./TR).*E2.^2.*M0s.*tan(a.*(1.0./2.0)).*...
    exp(TE.*wf.*-1i).*(E2.^2-1.0).*1.0./sqrt((E2.^2-1.0)./(E2.^2.*...
    (E1-cos(a)).^2.*1.0./(E1.*cos(a)-1.0).^2-1.0)).*(cos(a).^2-1.0).*...
    1.0./(E2.^2.*(E1-cos(a)).^2.*1.0./(E1.*cos(a)-1.0).^2-1.0).^2.*...
    (E1-cos(a)).*1.0./(E1.*cos(a)-1.0).^3.*-1i;
dsm_E1(isnan(dsm_E1)) = 0;
end

%% DESS first derivatives w.r.t. E2
function dsp_E2 = dfp_E2(M0s, E1, E2, a, ~, TE, wf) 
dsp_E2 = E2.*M0s.*tan(a.*(1.0./2.0)).*exp(TE.*wf.*1i).*(E1.^2-1.0).*1.0./...
    sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*1.0./(E1.*cos(a)-1.0).^2-...
    1.0)).*(cos(a).^2-1.0).*(E1-cos(a)).*(E1.*cos(a)-1.0).*1.0./(E1.^2.*...
    cos(a).^2-E2.^2.*cos(a).^2-E1.*cos(a).*2.0-E1.^2.*E2.^2+E1.*E2.^2.*...
    cos(a).*2.0+1.0).^2.*-1i;
dsp_E2(isnan(dsp_E2)) = 0;
end

function dsm_E2 = dfm_E2(M0s, E1, E2, a, TR, TE, wf)
dsm_E2 = (E2.^(-(TE.*2.0+TR)./TR).*M0s.*TE.*tan(a.*(1.0./2.0)).*exp(TE.*wf.*-1i).*...
    (sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*1.0./(E1.*cos(a)-1.0).^2-...
    1.0))-1.0).*-2.0i)./TR-E2.^(-(TE.*2.0-TR)./TR).*M0s.*tan(a.*(1.0./2.0)).*...
    exp(TE.*wf.*-1i).*(E1.^2-1.0).*1.0./sqrt((E2.^2-1.0)./(E2.^2.*...
    (E1-cos(a)).^2.*1.0./(E1.*cos(a)-1.0).^2-1.0)).*(cos(a).^2-1.0).*...
    (E1.*cos(a)-1.0).^2.*1.0./(E1.^2.*cos(a).^2-E2.^2.*cos(a).^2-E1.*...
    cos(a).*2.0-E1.^2.*E2.^2+E1.*E2.^2.*cos(a).*2.0+1.0).^2.*1i;
dsm_E2(isnan(dsm_E2)) = 0;
end

%% DESS second derivatives w.r.t. E2
function ddsp_E2 = ddfp_E2(M0s, E1, E2, a, ~, TE, wf)
ddsp_E2 = -M0s.*tan(a.*(1.0./2.0)).*exp(TE.*wf.*1i).*(E1.^2-1.0).*1.0./...
    ((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*1.0./(E1.*cos(a)-1.0).^2-1.0)).^...
    (3.0./2.0).*(cos(a).^2-1.0).*(E1-cos(a)).*(E1.*cos(a)-1.0).^3.*1.0./...
    (E1.^2.*cos(a).^2-E2.^2.*cos(a).^2-E1.*cos(a).*2.0-E1.^2.*E2.^2+E1.*...
    E2.^2.*cos(a).*2.0+1.0).^4.*(E1.^2.*cos(a).^2.*1i+E2.^2.*cos(a).^2.*...
    2.0i-E2.^4.*cos(a).^2.*3.0i-E1.*cos(a).*2.0i+E1.^2.*E2.^2.*2.0i-E1.^...
    2.*E2.^4.*3.0i-E1.*E2.^2.*cos(a).*4.0i+E1.*E2.^4.*cos(a).*6.0i+1i);
ddsp_E2(isnan(ddsp_E2)) = 0;
end

function ddsm_E2 = ddfm_E2(M0s, E1, E2, a, TR, TE, wf) 
ddsm_E2 = E2.^((TE.*-2.0-TR.*2.0)./TR).*M0s.*TE.*1.0./TR.^2.*tan(a.*(1.0./...
    2.0)).*exp(TE.*wf.*-1i).*(sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*...
    1.0./(E1.*cos(a)-1.0).^2-1.0))-1.0).*(TE.*2.0+TR).*2.0i-E2.^...
    ((TE.*-2.0+TR.*2.0)./TR).*M0s.*tan(a.*(1.0./2.0)).*exp(TE.*wf.*-1i).*...
    (E1.^2-1.0).^2.*1.0./((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*1.0./...
    (E1.*cos(a)-1.0).^2-1.0)).^(3.0./2.0).*(cos(a).^2-1.0).^2.*(E1.*...
    cos(a)-1.0).^4.*1.0./(E1.^2.*cos(a).^2-E2.^2.*cos(a).^2-E1.*cos(a).*...
    2.0-E1.^2.*E2.^2+E1.*E2.^2.*cos(a).*2.0+1.0).^4.*1i-E2.^...
    ((TE.*-2.0+TR.*2.0)./TR).*M0s.*tan(a.*(1.0./2.0)).*exp(TE.*wf.*-1i).*...
    (E1.^2-1.0).*1.0./sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*1.0./...
    (E1.*cos(a)-1.0).^2-1.0)).*(cos(a).^2-1.0).*(E1-cos(a)).^2.*(E1.*...
    cos(a)-1.0).^2.*1.0./(E1.^2.*cos(a).^2-E2.^2.*cos(a).^2-E1.*cos(a).*...
    2.0-E1.^2.*E2.^2+E1.*E2.^2.*cos(a).*2.0+1.0).^3.*4.0i+(E2.^...
    ((TE.*-2.0)./TR).*M0s.*TE.*tan(a.*(1.0./2.0)).*exp(TE.*wf.*-1i).*...
    (E1.^2-1.0).*1.0./sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*1.0./...
    (E1.*cos(a)-1.0).^2-1.0)).*(cos(a).^2-1.0).*(E1.*cos(a)-1.0).^2.*...
    1.0./(E1.^2.*cos(a).^2-E2.^2.*cos(a).^2-E1.*cos(a).*2.0-E1.^2.*...
    E2.^2+E1.*E2.^2.*cos(a).*2.0+1.0).^2.*2.0i)./TR+(E2.^((TE.*-2.0)./...
    TR).*M0s.*tan(a.*(1.0./2.0)).*exp(TE.*wf.*-1i).*(E1.^2-1.0).*1.0./...
    sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*1.0./(E1.*cos(a)-1.0).^...
    2-1.0)).*(cos(a).^2-1.0).*(TE.*2.0-TR).*(E1.*cos(a)-1.0).^2.*1.0./...
    (E1.^2.*cos(a).^2-E2.^2.*cos(a).^2-E1.*cos(a).*2.0-E1.^2.*E2.^2+E1.*...
    E2.^2.*cos(a).*2.0+1.0).^2.*1i)./TR;
ddsm_E2(isnan(ddsm_E2)) = 0;
end