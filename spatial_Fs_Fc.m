function [Fs, Fc, Fc_beta, Fs_beta, Feq_beta, phi_beta] = spatial_Fs_Fc(alpha_s, beta, N)

% Use safe variable names
syms alpha_s beta N 

%% ---------- S1 ----------
S1 = ( alpha_s/2 + beta ) * cos( beta - pi/N ) + ( alpha_s/2 - beta ) * cos( beta + pi/N ) - ( 2*pi/N - alpha_s/2 + beta ) * cos( alpha_s - beta - pi/N ) + ( 3*alpha_s/2 - beta - 2*pi/N ) * cos( alpha_s - beta - 3*pi/N );
%% ---------- S2 ----------
S2 = ( 3*alpha_s/2 - beta - 4*pi/N ) * cos( 2*alpha_s - beta - 3*pi/N ) + ( 5*alpha_s/2 - beta - 4*pi/N ) * cos( 2*alpha_s - beta - 5*pi/N ) - ( 6*pi/N - 5*alpha_s/2 + beta ) * cos( 3*alpha_s - beta - 5*pi/N ) + ( 7*alpha_s/2 - beta - 6*pi/N ) * cos( 3*alpha_s - beta - 7*pi/N );
%% ---------- S3 ----------
S3 = ( 7*alpha_s/2 - beta - 8*pi/N ) * cos( 4*alpha_s - beta - 7*pi/N ) + ( 9*alpha_s/2 - beta - 8*pi/N ) * cos( 4*alpha_s - beta - 9*pi/N ) - ( 10*pi/N - 9*alpha_s/2 + beta ) * cos( 5*alpha_s - beta - 9*pi/N ) + ( 11*alpha_s/2 - beta - 10*pi/N ) * cos( 5*alpha_s - beta - 11*pi/N );
%% ---------- Fs and Fc ----------
Fs = simplify( S1 - (S2 + S3)/2 );
Fc = simplify( sqrt(sym(3))/2 * ( S3 - S2 ) );

%% Here we convert it to a system function
S1 = symfun(S1,[alpha_s, beta, N]);
S1 = matlabFunction(S1);
%% Here we convert it to a system function
S2 = symfun(S2,[alpha_s, beta, N]);
S2 = matlabFunction(S2);
%% Here we convert it to a system function
S3 = symfun(S3,[alpha_s, beta, N]);
S3 = matlabFunction(S3);
%% Here we convert it to a system function
Fs = symfun(Fs,[alpha_s, beta, N]);
%% Here we convert it to a system function
Fc = symfun(Fc,[alpha_s, beta, N]);
% the derivatives of Fc and Fs with respect to beta
Fc_beta = symfun(diff(Fc, beta),[alpha_s, beta, N]);
Fs_beta = symfun(diff(Fs, beta),[alpha_s, beta, N]);
%% ---------- Equivalent amplitude and phase ----------
Feq_beta = simplify( sqrt(Fs_beta^2 + Fc_beta^2) );
phi_beta = simplify( atan2(Fc_beta, Fs_beta) );
% Here we convert them into Matlab functions
Fc = matlabFunction(Fc);
Fs = matlabFunction(Fs);
Fc_beta = matlabFunction(Fc_beta);
Fs_beta = matlabFunction(Fs_beta);
Feq_beta = matlabFunction(Feq_beta);
phi_beta = matlabFunction(phi_beta);


end



