% code by Brian Hunt, Zhixin Lu, Jaideep Pathak & Artur Perevalov
% prediction based on Jaeger reservoir model

%% Organizing the data for the prediction

% generating some data for the example
T = 20000;
time = [1:T]'/100;  % simple sin signal
% creating a 2D signal that contains an offset, usual sin signal, nonlinear sin and some noize
%data = [10 + sin(5*time).*sin(sqrt(13)*time) + tanh(5*sin(17*time)) + 0.1*rand(T,1), (4*sin(10*time))] ;
data = dd;
size(data)

% some preproduction for data
data = data - mean(data); % substracting the mean value might be useful
data = data/std(data(:));

% splitting the data in two parts - one for training one for predicting
start = 1;  % starting point of training
training_length = 19000;  % length of training
after_training = size(data,1)-1-training_length;  % length of prediction

Utr = data(start:start+training_length,:);                                      % training data
Ute = data(start+training_length+1:start+training_length+after_training,:);     % test data

%% Setting up parameters 

% parameters of the RC
rng(0);             % seed random number generator
alpha =1;           % slowing parameter
trans = 1000;       % transient time 
reg = 1e-6;         % regularization parameter
rho =   1;          % spectral radius
p = 0.0195;         % connection probability
N = 2^11;           % amount of neurons
scale_rho = 0.05;  % input strength


train = size(Utr,1);              % training time
dim = size(Utr,2);                % number of coordinates
scale = scale_rho*ones(1,dim);    % input strength for each coordinate
coord = 1;                        % coordinate to graph

U = Utr;                     

A = sparse(double(rand(N) < p));        % apply connection probability
A(A~=0) = 2*rand(1,length(find(A)))-1;  % randomize nonzero entries
A = rho*A/abs(eigs(A,1));               % scale spectral radius

Win = (2*rand(N,dim)-1)*diag(scale);    % input weights

X = zeros(train+1, N);                  % preallocate
X(1,:) = 2*rand(1,N)-1;
                            
                            
%% initialize reservoir
for n = 1:train                                 % listening phase
    X(n+1,:) = (1-alpha) * X(n,:) + alpha*tanh(X(n,:)*A' + U(n,:)*Win');   % reservoir with input
end

init = X(train+1,:);        % initial random state for prediction
X = X((trans+1):train,:);   % throw away transient part
U = U((trans+1):train,:);

%% Here is the magic
% solving system of equations
Wout = (X'*X + reg*eye(N))\(X'*U);  % fit X*Wout to U with regularization


%% Plotting training phase 
plots = 1; %1 if you need to visualise plots

if plots == 1  % plotting training part
    figure(1)
    tvals = (trans+1):train;
    subplot(2,1,1)              % graph training signal versus reservoir output
    plot(tvals, U(:,coord), 'b', tvals, X*Wout(:,coord), 'r'); axis tight
    subplot(2,1,2)              % graph normalized RMS error
    plot(tvals, sqrt(sum((X*Wout-U).^2,2)./(sum((X*Wout).^2,2)+sum(U.^2,2))))

    axis tight
end

%% Prediction phase

test = size(Ute,1);             % testing time
U = Ute;

XX = zeros(test, N);
XX(1,:) = init;

% prediction
for n = 1:(test-1)      % predict using autonomous reservoir with feedback
    XX(n+1,:) = (1-alpha) * XX(n,:) + alpha*tanh(XX(n,:)*A' + (XX(n,:)*Wout)*Win');
end

%% Plotting prediction phase
if plots == 1
    figure(2)
    tvals = 1:test;
    subplot(2,1,1)              % graph test signal versus reservoir prediction
    plot(tvals, U(:,coord), 'b', tvals, XX*Wout(:,coord), 'r'); axis tight
    subplot(2,1,2)              % graph normalized RMS error
    plot(tvals, sqrt(sum((XX*Wout-U).^2,2)./(sum((XX*Wout).^2,2)+sum(U.^2,2))))

    axis tight
end

%% Saving the prediction phase data
predictions = XX*Wout; 
experiment = U;

%% Clearing the variables
clearvars A after_training coord dim init
clearvars n N out p reg rho scale shift start test train training_length
clearvars trans  Ute Utr Win  i scale_rho Ustd data_norm_coef k fid ans timestep
clearvars tvals a b windowsize plots error_type data
clearvars Wout X U XX mean_d s stdd
clearvars alpha er_matr er_vec res_array time T