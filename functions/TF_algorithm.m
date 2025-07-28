function [tfc, tfrsq, tfrsq_proposed, tcrtic, tfrtic, ttrtic] = TF_algorithm(x, fs, lowFreq, highFreq, chirpMax, alpha, tDS, chirpReso, h, Dh, DDh, DDDh)
% x: input signal
% lowFreq: minimal frequency to observe, smallest as 0
% highFreq: maximal frequency to observe, largest as 0.5
% alpha: resolution for frequency and chirp rate ; example 中取 alpha = 2/length(x)
% tDS: hopping in time (時間間隔?)
% h: window function, a column vector
% Dh: h'
% DDh: h''
% DDDh: h'''
thres_q0 = 1e-12;
[xrow,xcol] = size(x) ; % should be only one column.

% Time resolution:
dt = 1/fs; % time interval
t = [0:length(x)-1]*dt ; % t=[1 2 3 4 ... length(x)];
tau = (0:tDS:length(x)-1)*dt; % Outpus time 
ttrtic = tau; % Outpus time 

% Frequency resolution:
N = length([-0.5+alpha:alpha:0.5]) ; % N : Number of frequency bins
Q = floor(N/2); % Half points
df = 1/N/dt; % frequency interval
f = (0:Q)'*df; % half real frequency
Lidx = knnsearch(f, lowFreq, 'k', 1); % Low frequency index
Hidx = knnsearch(f, highFreq, 'k', 1); % High frequency index

% Chirp-rate resolution:
dc = 2*(fs/N)^2;  % chirp rate interval
crate = ([1:N]-ceil(N/2))'*dc; % discretization of chirp rate
chop = max(round(chirpReso/dc), 1);
crate = crate(1:chop:N);
crate = crate(crate<=chirpMax & crate>=-chirpMax);
% if N=10, crate = [-4 -3 -2 -1 0 1 2 3 4]/100
% if N=11, crate = [-5 -4 -3 -2 -1 0 1 2 3 4]/121 ====== not sure what it is.

xi = f; % Output frequency
lamb = crate; % Output chirp rate

% Length
t_len = length(t); % Input time length
fLen = length(f); % number of bins of frequency
cLen = length(crate); % number of bins of chirp rate
tau_Len = length(tau) ; % Output time length
xi_len = length(xi); % Output frequency length
lamb_len = length(lamb); % Output chirp rate length
%====================================================================
%% check input signals
if (xcol~=1)
    error('X must have only one column');
%elseif highFreq > 0.5
%    error('TopFreq must be a value in [0, 0.5]');
elseif (tDS < 1) || (rem(tDS,1))
    error('tDS must be an integer value >= 1');
end

[hrow,hcol] = size(h); 
if (hcol~=1)||(rem(hrow,2)==0) % ||= or
    error('H must be a smoothing window with odd length');
end
ht =(-Q:Q)' ; % ht=[-Lh -Lh+1 ... Lh-1 Lh]
%====================================================================
%% run STFT and reassignment rule
tfc = zeros(cLen, fLen, tau_Len); 	% chirplet transform
tfrsq = zeros(lamb_len, xi_len, tau_Len); % synchrosqueezed chirplet transform
tfrsq_proposed = zeros(lamb_len, xi_len, tau_Len); % Proposed algorithm
tfrtic = f ; % positive frequency
tcrtic = crate;                   % Chirp axis

Ex = mean(abs(x(min(t/dt+1):max(t/dt))).^2);
Threshold = 1.0e-8*Ex;  % threshold 臨界點

% Running reassignment:
fprintf(['Chirp-rate total: ',num2str(cLen), '; now:     ']) ;
for cidx = 1:cLen
    fprintf('\b\b\b\b') ;	tmp = sprintf('%4d',cidx) ; fprintf(tmp) ;
    chirp = crate(cidx);
    for tau_idx = 1:tau_Len
        % tidx is the current time in output time tick:
        tidx = tau_idx*tDS+1;
        indices  = mod((tidx-Q:tidx+Q), t_len)+1; % Truncate from n-Q to n+Q
        % tau is the relevant index associated with ti
        %tau = -min([round(N/2)-1,Lh,tidx-1]):min([round(N/2)-1,Lh,xrow-tidx]);
        % indices is the absolute index in the "evaluation window"
        %indices= rem(N+tau,N)+1; % remainder

        tf0 = zeros(N, 1) ; tf1 = zeros(N, 1) ; tf2 = zeros(N, 1) ;
        tf0x1 = zeros(N, 1) ; tf1x1 = zeros(N, 1) ; tf0x2 = zeros(N,1);
        tf3 = zeros(N, 1); tf2x1 = zeros(N, 1); tf1x2 = zeros(N, 1);
        tf0x3 = zeros(N, 1); tf0x4 = zeros(N,1); tf1x3 = zeros(N, 1);
        tf2x2 = zeros(N, 1);
        % tfaxb = chirplet transform with x^(a)*[(d/dx)^(b)g]
        tf0(1:N) = x(indices).*h.*exp(-pi*1i*chirp.*...
            (ht).^2*dt^2); % for CT with window g
        tf1(1:N) = x(indices).*Dh.*exp(-pi*1i*chirp.*...
            (ht).^2*dt^2); % for CT with window g'
        tf2(1:N) = x(indices).*DDh.*exp(-pi*1i*chirp.*...
            (ht).^2*dt^2); % for CT with window g''
        tf0x1(1:N) = x(indices).*h.*(ht*dt).*...
            exp(-pi*1i*chirp.*(ht).^2*dt^2); % for CT with window xg
        tf1x1(1:N) = x(indices).*Dh.*(ht*dt).*...
            exp(-pi*1i*chirp.*(ht).^2*dt^2); % for CT with window xg'
        tf0x2(1:N) = x(indices).*h.*(ht*dt).^2.*...
            exp(-pi*1i*chirp.*(ht).^2*dt^2); % for CT with window x^2g
        tf3(1:N) = x(indices).*DDDh.*exp(-pi*1i*chirp.*...
            (ht).^2*dt^2); % for CT with window g'''
        tf2x1(1:N) = x(indices).*DDh.*(ht*dt).*...
            exp(-pi*1i*chirp.*(ht).^2*dt^2); % for CT with window xg''
        tf1x2(1:N) = x(indices).*Dh.*(ht*dt).^2.*...
            exp(-pi*1i*chirp.*(ht).^2*dt^2); % for CT with window x^2g'
        tf0x3(1:N) = x(indices).*h.*(ht*dt).^3.*...
            exp(-pi*1i*chirp.*(ht).^2*dt^2); % for CT with window x^3g
        tf0x4(1:N) = x(indices).*h.*(ht*dt).^4.*...
            exp(-pi*1i*chirp.*(ht).^2*dt^2); % for CT with window x^4g
        tf1x3(1:N) = x(indices).*Dh.*(ht*dt).^3.*...
            exp(-pi*1i*chirp.*(ht).^2*dt^2); % for CT with window x^3g'
        tf2x2(1:N) = x(indices).*DDh.*(ht*dt).^2.*...
            exp(-pi*1i*chirp.*(ht).^2*dt^2); % for CT with window x^2g''
        % g''', xg'', x^2g', x^3g

        % Fast Fourier transform with N points: 
        tf0 = fft(tf0, N) ; %tf0 = tf0(1:N/2) ;
        tf1 = fft(tf1, N) ; %tf1 = tf1(1:N/2) ;
        tf2 = fft(tf2, N) ; %tf2 = tf2(1:N/2) ;
        tf0x1 = fft(tf0x1, N) ; %tf0x1 = tf0x1(1:N/2) ;
        tf1x1 = fft(tf1x1, N) ; %tf1x1 = tf1x1(1:N/2) ;
        tf0x2 = fft(tf0x2, N) ; %tf0x2 = tf0x2(1:N/2) ;
        tf3 = fft(tf3, N) ; %tf3 = tf3(1:N/2) ;
        tf2x1 = fft(tf2x1, N) ; %tf2x1 = tf2x1(1:N/2) ;
        tf1x2 = fft(tf1x2, N) ; %tf1x2 = tf1x2(1:N/2) ;
        tf0x3 = fft(tf0x3, N) ; %tf0x3 = tf0x3(1:N/2) ;
        tf0x4 = fft(tf0x4, N) ; %tf0x4 = tf0x4(1:N/2) ;
        tf1x3 = fft(tf1x3, N) ; %tf1x3 = tf1x3(1:N/2) ;
        tf2x2 = fft(tf2x2, N) ; %tf2x2 = tf2x2(1:N/2) ;

        % Shift to single spectrum and select nonnegative frequencies :
        tf0 = dt*exp(1i*2*pi*Q*(0:Q)'/N).*tf0(1:(Q+1));
        tf1 = dt*exp(1i*2*pi*Q*(0:Q)'/N).*tf1(1:(Q+1));
        tf2 = dt*exp(1i*2*pi*Q*(0:Q)'/N).*tf2(1:(Q+1));
        tf3 = dt*exp(1i*2*pi*Q*(0:Q)'/N).*tf3(1:(Q+1));
        tf0x1 = dt*exp(1i*2*pi*Q*(0:Q)'/N).*tf0x1(1:(Q+1));
        tf1x1 = dt*exp(1i*2*pi*Q*(0:Q)'/N).*tf1x1(1:(Q+1));
        tf2x1 = dt*exp(1i*2*pi*Q*(0:Q)'/N).*tf2x1(1:(Q+1));
        tf0x2 = dt*exp(1i*2*pi*Q*(0:Q)'/N).*tf0x2(1:(Q+1));
        tf1x2 = dt*exp(1i*2*pi*Q*(0:Q)'/N).*tf1x2(1:(Q+1));
        tf2x2 = dt*exp(1i*2*pi*Q*(0:Q)'/N).*tf2x2(1:(Q+1));
        tf0x3 = dt*exp(1i*2*pi*Q*(0:Q)'/N).*tf0x3(1:(Q+1));
        tf1x3 = dt*exp(1i*2*pi*Q*(0:Q)'/N).*tf1x3(1:(Q+1));
        tf0x4 = dt*exp(1i*2*pi*Q*(0:Q)'/N).*tf0x4(1:(Q+1));


        % Frequency components:
        %jcol_vec = ((1:1:N/2)')./N;
        % q0 and q1:
        q0 = (tf0x1.*tf1 - tf1x1.*tf0 + 2.0*pi*1i*chirp*(tf0x2.*tf0 - tf0x1.*tf0x1))./(tf0.*tf0);
        q1 = (tf0x2.*tf1 - tf1x2.*tf0 + 2.0*pi*1i*chirp*(tf0x3.*tf0 - tf0x2.*tf0x1))./(tf0.*tf0);

        Dtf0 = -tf1 + 2.0*pi*1i*f.*tf0 + 2.0*pi*1i*chirp*tf0x1; % (d/dt)T^g_f
        DDtf0 = tf2 - 4*pi*1i*f.*tf1 - 2.0*pi*1i*chirp*tf1x1 + (2.0*pi*1i*f).^2.*tf0 +...
                        2.0*(2.0*pi*1i*chirp)*(2.0*pi*1i*f).*tf0x1 ...
                        - 2.0*pi*1i*chirp*(tf0+tf1x1) + (2.0*pi*1i*chirp)^2*tf0x2; % (d/dt)^2 T^g_f
        DDDtf0 = -tf3 + 6.0*pi*1i*f.*tf2 + 6.0*pi*1i*chirp*tf2x1 + tf1.*(6.0*pi*1i*chirp - 12.0*(pi*1i*f).^2)...
                + tf1x1.*(-24.0*(pi*1i*f)*(pi*1i*chirp)) + tf0.*(8.0*(pi*1i*f).^3 - 12.0*(pi*1i*chirp)*(pi*1i*f))...
                + tf1x2*(-12.0*(pi*1i*chirp)^2) + tf0x1.*(24.0*(pi*1i*f).*(pi*1i*f)*(pi*1i*chirp) - 12.0*(pi*1i*chirp)^2)...
                + tf0x2.*(24.0*(pi*1i*chirp)^2*(pi*1i*f)) + tf0x3*(8.0*(pi*1i*chirp)^3); % (d/dt)^3 T^g_f
        D_DTT = (DDtf0.*tf0 - Dtf0.*Dtf0)./(tf0.^2); % d/dt(d/dt(T^g_f)/T^g_f)
        DD_DTT = ((DDDtf0.*tf0 - DDtf0.*Dtf0).*tf0.*tf0 - 2.0*tf0.*Dtf0.*(DDtf0.*tf0 - Dtf0.*Dtf0))./(tf0.^4);% (d/dt)^2(d/dt(T^g_f)/T^g_f)
        
        Dtf0x1 =  -tf0 - tf1x1 + 2*pi*1i*f.*tf0x1 + 2*pi*1i*chirp*tf0x2;
        Dtf0x2 = -2*tf0x1 - tf1x2 + 2*pi*1i*f.*tf0x2 + 2*pi*1i*chirp*tf0x3;
        Dtf0x3 = -3*tf0x2 - tf1x3 + 2*pi*1i*f.*tf0x3 + 2*pi*1i*chirp*tf0x4;
        
        Dtf1x0 =          -tf2 + 2*pi*1i*f.*tf1 + 2*pi*1i*chirp*tf1x1;
        Dtf1x1 = -tf1 - tf2x1 + 2*pi*1i*f.*tf1x1 + 2*pi*1i*chirp*tf1x2;
        Dtf1x2 = -2*tf1x1 - tf2x2 + 2*pi*1i*f.*tf1x2 + 2*pi*1i*chirp*tf1x3;

        q0_up = tf0x1.*tf1 - tf1x1.*tf0 + 2.0*pi*1i*chirp*(tf0x2.*tf0 - tf0x1.*tf0x1);
        Dq0_up_dq0up = Dtf0x1.*tf1 + tf0x1.*Dtf1x0 - Dtf1x1.*tf0 - tf1x1.*Dtf0 +...
            2*pi*1i*chirp*(Dtf0x2.*tf0 + tf0x2.*Dtf0 - 2*tf0x1.*Dtf0x1);
        Dq0_up_dq0down = 2*tf0.*Dtf0;
        Dq0 = (Dq0_up_dq0up.*tf0.^2 - Dq0_up_dq0down.*q0_up)./(tf0.^4);
        
        q1_up = tf0x2.*tf1 - tf1x2.*tf0 + 2.0*pi*1i*chirp*(tf0x3.*tf0 - tf0x2.*tf0x1);
        Dq1_up_dq1up = Dtf0x2.*tf1 + tf0x2.*Dtf1x0 - Dtf1x2.*tf0 - tf1x2.*Dtf0 +...
            2*pi*1i*chirp*(Dtf0x3.*tf0 + tf0x3.*Dtf0 - Dtf0x2.*tf0x1 - tf0x2.*Dtf0x1);
        Dq1_up_dq1down = 2*tf0.*Dtf0;
        Dq1 = (Dq1_up_dq1up.*tf0.^2 - Dq1_up_dq1down.*q1_up)./(tf0.^4);

%         % theoratical derivatives for q_j s.
%         Dq0 = (tf1x1.*tf0.*tf0.*(-tf1 + 2.0*pi*1i*jcol_vec.*tf0 + 2.0*pi*1i*chirp*tf0x1)...
%              + (-2.0*pi*1i*chirp)*tf0x2.*tf0.*tf0.*(-tf1 + 2.0*pi*1i*jcol_vec.*tf0 + 2.0*pi*1i*chirp*tf0x1)...
%              + tf0.*tf0.*(tf1 - 2.0*pi*1i*chirp*tf0x1).*(-tf0 - tf1x1 + 2.0*pi*1i*jcol_vec.*tf0x1 + 2.0*pi*1i*chirp*tf0x2)...
%              + tf0.*tf0.*tf0.*(-2.0*pi*1i*chirp*2.0*tf0x1 - 4.0*pi*1i*chirp*tf1x2 - 4*pi*pi*chirp*chirp*tf0x3 ...
%              + tf1 - 4*pi*pi*chirp*jcol_vec.*tf0x2 - 2.0*pi*1i*jcol_vec.*tf1x1 + tf2x1)...
%              + tf0x1.*tf0.*(-tf0.*tf2 - 2.0*pi*1i*jcol_vec.*tf0.*tf1 + 4.0*pi*1i*chirp*tf1x1.*tf0 + 2.0*pi*1i*chirp*tf0.*tf0 ...
%              - 4.0*pi*pi*chirp*jcol_vec.*tf0.*tf0x1 + 4.0*pi*pi*chirp*chirp*tf0.*tf0x2 + 2.0*tf1.*tf1 - 8.0*pi*1i*chirp*tf1.*tf0x1 ...
%              - 8.0*pi*pi*chirp*chirp*tf0x1.*tf0x1))./(tf0.*tf0.*tf0.*tf0);
%         Dq1 = (tf1x2.*tf0.*tf0.*(-tf1 + 2.0*pi*1i*jcol_vec.*tf0 + 2.0*pi*1i*chirp*tf0x1)...
%              + (-2.0*pi*1i*chirp)*tf0x3.*tf0.*tf0.*(-tf1 + 2.0*pi*1i*jcol_vec.*tf0 + 2.0*pi*1i*chirp*tf0x1)...
%              + tf0.*tf0.*(tf1 - 2.0*pi*1i*chirp*tf0x1).*(-2.0*tf0x1 - tf1x2 + 2.0*pi*1i*jcol_vec.*tf0x2 + 2.0*pi*1i*chirp*tf0x3)...
%              + tf0.*tf0.*tf0.*(-2.0*pi*1i*chirp*3.0*tf0x2 - 4.0*pi*1i*chirp*tf1x3 - 4*pi*pi*chirp*chirp*tf0x4 ...
%              + 2.0*tf1x1 - 4*pi*pi*chirp*jcol_vecsqSTCT.*tf0x3 - 2.0*pi*1i*jcol_vec.*tf1x2 + tf2x2)...
%              + tf0x2.*tf0.*(-tf0.*tf2 - 2.0*pi*1i*jcol_vec.*tf0.*tf1 + 4.0*pi*1i*chirp*tf1x1.*tf0 + 2.0*pi*1i*chirp*tf0.*tf0 ...
%              - 4.0*pi*pi*chirp*jcol_vec.*tf0.*tf0x1 + 4.0*pi*pi*chirp*chirp*tf0.*tf0x2 + 2.0*tf1.*tf1 - 8.0*pi*1i*chirp*tf1.*tf0x1 ...
%              - 8.0*pi*pi*chirp*chirp*tf0x1.*tf0x1))./(tf0.*tf0.*tf0.*tf0);

        den_q0 = (2.0*q0.*q0 + Dq1.*q0 - q1.*Dq0);
        %den_q0_new = Check_zero(den_q0,thres_q0);
        den_q0_new = den_q0;

        % CT:
        tfc(cidx,:,tau_idx) = tf0;

        % SCT :
        theta0 = 0;
        %lambda0 = (D_DTT - pi*1i*theta0.*q1)./Check_zero(q0,thres_q0)./(2.0*pi*1i);
        lambda0 = (D_DTT - pi*1i*theta0.*q1)./q0./(2.0*pi*1i);
        lambda = real(lambda0) ;
        %omega = real((Dtf0 - 2.0*pi*1i*lambda0.*tf0x1 - pi*1i*theta0.*tf0x2)./(2.0*pi*1i*Check_zero(tf0,thres_q0))); 
        omega = real((Dtf0 - 2.0*pi*1i*lambda0.*tf0x1 - pi*1i*theta0.*tf0x2)./(2.0*pi*1i*tf0));  
        omega_idx = round((omega-xi(1))/(xi(end)-xi(1))*(xi_len-1))+1; 
        lambda_idx = round((lambda-lamb(1))/(lamb(end)-lamb(1))*(lamb_len-1))+1;
        for jcol_ = 1: fLen
            jhat_ = omega_idx(jcol_);
            khat_ = lambda_idx(jcol_);     
            if (abs(tf0(jcol_)) > Threshold)
                if (jhat_ <= xi_len) && (jhat_ >= 1) && (khat_ >= 1) && (khat_ <= lamb_len) 
                    tfrsq(khat_,jhat_,tau_idx) =...
                        tfrsq(khat_,jhat_,tau_idx) + tf0(jcol_) ;
                end
            end
        end 

        % Proposed algorithm:
        theta0_proposed = real(((DD_DTT.*q0 - D_DTT.*Dq0)./den_q0_new./(1i*pi)));%/(N^(beta)); %10,sqrt(N/2)
        %lambda0_proposed = (D_DTT - pi*1i*theta0_proposed.*q1)./Check_zero(q0,thres_q0)./(2.0*pi*1i);
        lambda0_proposed = (D_DTT - pi*1i*theta0_proposed.*q1)./q0./(2.0*pi*1i);
        lambda_proposed = real(lambda0_proposed) ;
        %omega_proposed = real((Dtf0 - 2.0*pi*1i*lambda0_proposed.*tf0x1 - pi*1i*theta0_proposed.*tf0x2)./(2.0*pi*1i*Check_zero(tf0,thres_q0)));
        omega_proposed = real((Dtf0 - 2.0*pi*1i*lambda0_proposed.*tf0x1 - pi*1i*theta0_proposed.*tf0x2)./(2.0*pi*1i*tf0));
        omega_pro_idx = round((omega_proposed-xi(1))/(xi(end)-xi(1))*(xi_len-1))+1; 
        lambda_pro_idx = round((lambda_proposed-lamb(1))/(lamb(end)-lamb(1))*(lamb_len-1))+1;
        for jcol = 1: fLen
            jhat = omega_pro_idx(jcol);
            khat = lambda_pro_idx(jcol);     
            if (abs(tf0(jcol)) > Threshold)
                if (jhat <= xi_len) && (jhat >= 1) && (khat >= 1) && (khat <= lamb_len) 
                    tfrsq_proposed(khat,jhat,tau_idx) =...
                        tfrsq_proposed(khat,jhat,tau_idx) + tf0(jcol) ;
                end
            end
        end       
    end
end
fprintf('\n') ;
end