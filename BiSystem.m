function obj=BiSystem(Input)

%% load parameters
K=Input.K;
N1=Input.N1;
N2=Input.N2;
N3=Input.N3;
nuw1=Input.nuw1;
nuw2=Input.nuw2;
epsilon=Input.epsilon;


%% Generate system 
X  = randn(N1, K);
H1 = randn(N2, N1);
H2 = randn(N3, N2);

Z1=H1*X;
X2=Z1+sqrt(nuw1)*randn(N2,K);

%% valid entries of observation 
omega = false(N2,K);
ind = randperm(N2*K);
omega(ind(1:ceil(epsilon*N2*K))) = true;
X2(~omega) = 0;


Z2=H2*X2;
Y=Z2+sqrt(nuw2)*randn(N3, K);





%% save parameters
obj.X=X;
obj.Y=Y;
obj.H1=H1;
obj.H2=H2;
obj.Z1=Z1;
obj.omega=omega;


end