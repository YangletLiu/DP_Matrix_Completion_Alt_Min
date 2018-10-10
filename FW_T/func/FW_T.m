function [ output ] = FW_T(par)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% this function implements Frank-Wolfe-thresholding method
% 0.5*norm(Omega.*(L+S-M), 'fro')^2 + lambda_1*||L||_* + lambda_2*||S||_1
%
% please refer to the following paper:
%     "Scalable robust matrix recovery: Frank-Wolfe meets proximal methods"
% Input
%   par.M: data matrix (observations)
%   par.lambda_1: lambda_1
%   par.lambda_2: lambda_2
%   par.Omega: sampling operator
%   par.showvideo
%   par.framesize
%
% Output:
%   output.L: the low-rank component recovered
%   output.S: the sparse component recovered
%   output.hist: objective values in history
%   output.iter: number of iterations used
%   output.time: time elapsed
%  Cun Mu, Apr. '16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% setup
M = par.M; % data matrix
[m,n] = size(M);


% stopping criterion by default (if no specification)
epsilon = 10^-3; showvideo = 0;
if isfield(par, 'epsilon') epsilon = par.epsilon; end
if isfield(par, 'showvideo') showvideo = par.showvideo; end


lambda_1 = par.lambda_1;
lambda_2 = par.lambda_2;
iter = par.iter;
Omega = par.Omega; % obeational pattern
frameSize = par.framesize;

flag_T = false;
if m > n
    
    M = M';
    Omega = Omega';
    [m,n] = size(M);
    flag_T = true;
  
end

%% FW-T method

%% initialization
L = zeros(m,n); S = zeros(m,n);
t_1 = 0; t_2 = 0;

U_1 = 0.5*norm(M,'fro')^2/lambda_1; % initial rough guess
U_2 = 0.5*norm(M,'fro')^2/lambda_2; % initial rough guess

history = 0; % to store the values of each iteration
time = cputime;
fprintf('%6s       %3s         %10s     %5s \n', ...
    'iter.', 'obj. ', 'rel. err.', 'count');

temp = Omega.*(L + S - M);  % gradient
history(1) = norm(temp,'fro')^2/2+lambda_1*t_1+lambda_2*t_2;
fprintf('|  %2d |   %10.5d   |   %7.5d       |%2.1d | \n', ...
    0, history(1), inf, 0);
count = 0;

%% loop
for k = 1: iter
    
    
    %------------------linearization subproblem-----------------------%
    [U,ev,V] = lansvd(temp,1,'l'); % top singular pair    
    
    D_L = -U*V';
    
    if lambda_1 >= ev
        V_L = zeros(m,n); V_t_1 = 0;
    else
        V_L = U_1*D_L; V_t_1 = U_1;
    end
    
    % [mag ind] = max(reshape(abs(temp), numel(abs(temp)), 1));
    [mag ind] = max(abs(temp(:)));
    j = floor((ind-1)/m)+1; i = mod(ind-1,m)+1;
    sign_ = sign(temp(i,j));
    D_S = zeros(m,n);
    D_S(i,j) = -sign_;
    
    if lambda_2 >= mag
        V_S = 0; V_t_2 = 0;
    else
        V_S = U_2*D_S; V_t_2 = U_2;
    end
    
    
    %--------------------- use QP (exact search) ---------------------%
    H = zeros(2,2);
    temp_1 = Omega.*(V_L-L); temp_2 = Omega.*(V_S-S); % temp = L + S - M;
    H(1,1) = norm(temp_1, 'fro')^2;
    H(2,2) = norm(temp_2, 'fro')^2;
    H(1,2) = temp_1(:)'*temp_2(:);
    H(2,1) = H(1,2);
    f = zeros(1,2);
    f(2) = temp(:)'*temp_2(:);
    f = f + [lambda_1*(V_t_1-t_1), lambda_2*(V_t_2-t_2)];
    
    % using QP solvers
    lb = zeros(2,1);
    ub = [1;1];
    f(1) = temp(:)'*temp_1(:);
    options.Display = 'off';
    options.TolFun = 1e-12;
    x = quadprog(H,f,[],[],[],[],lb,ub,[],options);    
    
    %----------------------- update L and S --------------------------%
    alpha = x(1);
    beta = x(2);
    L = (1-alpha)*L + alpha*V_L; t_1 = t_1 + alpha*(V_t_1-t_1);
    S = (1-beta)*S + beta*V_S; t_2 = t_2 + beta*(V_t_2-t_2);
    
    %------------------------ thresholding----------------------------%
    temp_3 = S-Omega.*(L+S-M);
    S = max(temp_3 - lambda_2, 0);
    S = S + min(temp_3 + lambda_2, 0);
    t_2 = norm(S(:),1);
    
    %----------- update U_1 and U_2 to a better esitmate -------------%
    temp = Omega.*(L + S - M);  % gradient
    history(k+1) = norm(temp,'fro')^2/2+lambda_1*t_1+lambda_2*t_2;
    U_1 = min(history(k+1)/lambda_1, U_1);
    U_2 = min(history(k+1)/lambda_2, U_2);
    
    rel_err =  abs(history(k+1)-history(k))/history(k);
    
    if rel_err <epsilon
        count = count+1;
    else
        count=0;
    end
    if count == 5
        fprintf('|  %2d |   %10.5d   |   %10.5d   |%2.1d | \n', ...
            k, history(k+1), rel_err, count);
        break;
    end
    
    fprintf('|  %2d |   %10.5d   |   %10.5d   |%2.1d | \n', ...
        k, history(k+1), rel_err, count);
    
end
time_elapsed = cputime - time;

if flag_T  
   
    % transpose back
    M = M';
    L = L';
    S = S';
    Omega = Omega';
    
end


%% output
output.L = L;
output.S = S;
output.hist = history;
output.iter = k;
output.time = time_elapsed;

% summary of output
fprintf('**************************************************************\n')
fprintf('Summary of Output: \n')
fprintf('%20s  %10s  \n', 'time', 'iter.')
fprintf('%5s  %15.4d  %6d  \n', 'FW-T', output.time, output.iter)
fprintf('**************************************************************\n')


if showvideo
    
    % original video
    toplay = M;
    
    no_frames = size(toplay,2);
    width = frameSize(1);
    length = frameSize(2);
    video_ori = zeros(width,length,no_frames); % width*length*frames
    for i=1:no_frames
        video_ori(:,:,i) = reshape(toplay(:,i),width,length);
    end
        
    % background term video
    toplay = output.L;
    no_frames = size(toplay,2);
    width = frameSize(1);
    length = frameSize(2);
    video_L = zeros(width,length,no_frames); % width*length*frames
    for i=1:no_frames
        video_L(:,:,i) = reshape(toplay(:,i),width,length);
    end
    
    % sparse term video
    toplay = output.S;
    no_frames = size(toplay, 2);
    width = frameSize(1);
    length = frameSize(2);
    video_S = zeros(width,length,no_frames); % width*length*frames
    for i=1:no_frames
        video_S(:,:,i) = reshape(toplay(:,i),width,length);
    end
    
    figure();
    
    video_combine = [video_ori, video_L, abs(video_S)];
    insert = zeros(width,10,no_frames)+ max(max(max(video_combine)));
    video_combine = [video_ori, insert, video_L, insert, abs(video_S)];
    
   % display video now 
   for i=1:no_frames
        imagesc(video_combine(:,:,i)); 
        colormap('gray'); axis off;  % title('FW-T');
        pause(.01)
   end
    
end
end