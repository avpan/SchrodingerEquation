clear all;clc;clf;
n = 1; %n = 1,2,3 eigenstates
N = 256; %N initialize
L = 2;
h = L/(N);
xo = -1;
x = xo;

%find wave function, psi for N where psi = psi_n
for i = 1:N+1
    if(i == 1 || i == N+1)
        psi(i) = 0;
    else
        psi(i) = sin(n*pi*x/2);
        X(i-1) = x;
    end
    x = x + h;
end
x = x-h;
psi = transpose(psi); 
X = transpose(X);
%initialize potential matrix
V = zeros(N);
V(1,1) = Inf;
V(N+1,N+1) = Inf;

%using 6th order difference we create the d^2/dx^2
for i = 1:N+1
    if(i-3 > 0 && i+3 < N+1)
        d2dx2(i,i-3) = 1/90;
        d2dx2(i,i-2) = -3/20;
        d2dx2(i,i-1) = 3/2;
        d2dx2(i,i) = -49/18;
        d2dx2(i,i+1) = 3/2;
        d2dx2(i,i+2) = -3/20;
        d2dx2(i,i+3) = 1/90;
    elseif(i - 2 > 0 && i+3 < N+1)
        d2dx2(i,i-2) = -3/20;
        d2dx2(i,i-1) = 3/2;
        d2dx2(i,i) = -49/18;
        d2dx2(i,i+1) = 3/2;
        d2dx2(i,i+2) = -3/20;
        d2dx2(i,i+3) = 1/90;
    elseif(i-1 > 0 && i+3 < N+1)
        d2dx2(i,i-1) = 3/2;
        d2dx2(i,i) = -49/18;
        d2dx2(i,i+1) = 3/2;
        d2dx2(i,i+2) = -3/20;
        d2dx2(i,i+3) = 1/90;
    elseif(i+1 > N+1)
        d2dx2(i,i-3) = 1/90;
        d2dx2(i,i-2) = -3/20;
        d2dx2(i,i-1) = 3/2;
        d2dx2(i,i) = -49/18;
    elseif(i+2 > N+1)
        d2dx2(i,i-3) = 1/90;
        d2dx2(i,i-2) = -3/20;
        d2dx2(i,i-1) = 3/2;
        d2dx2(i,i) = -49/18;
        d2dx2(i,i+1) = 3/2;
    elseif(i+3 > N+1)
        d2dx2(i,i-3) = 1/90;
        d2dx2(i,i-2) = -3/20;
        d2dx2(i,i-1) = 3/2;
        d2dx2(i,i) = -49/18;
        d2dx2(i,i+1) = 3/2;
        d2dx2(i,i+2) = -3/20;
    elseif(i+3 == N+1)
        d2dx2(i,i-3) = 1/90;
        d2dx2(i,i-2) = -3/20;
        d2dx2(i,i-1) = 3/2;
        d2dx2(i,i) = -49/18;
        d2dx2(i,i+1) = 3/2;
        d2dx2(i,i+2) = -3/20;
        d2dx2(i,i+3) = 1/90;
    elseif(i==1)
        d2dx2(i,i) = -49/18;
        d2dx2(i,i+1) = 3/2;
        d2dx2(i,i+2) = -3/20;
        d2dx2(i,i+3) = 1/90;
    end
end
d2dx2 = d2dx2/h^2;
%to find eigenvalues for infinite sqwell
%I must delete first row/colum and last row/column
subpsi = removerows(psi,N+1);
subpsi = removerows(subpsi,1);
psi;
subpsi;

subV = removerows(V,N+1);
subV = removerows(subV.',N+1)'; %delete columns and convert back to its original 
subV = removerows(subV.',1)';
subV = removerows(subV,1);

%add perturbation
P = 5;
for i = 1:(N-1)
    subVprime(i,i) = subV(i,i) + X(i)/2;
end
subVprime;

subd2dx2 = removerows(d2dx2,N+1);
subd2dx2 = removerows(subd2dx2.',N+1)';
subd2dx2 = removerows(subd2dx2.',1)';
subd2dx2 = removerows(subd2dx2,1);
d2dx2;
-subd2dx2;

%construct the Hamiltonian 
%satisfys H|v> = lambda|v>
H = -subd2dx2/2 + subV;
[v,d] = eig(H);
lambda = eig(H);
n = 1:size(lambda); %n = 1,2,3...
Energy = transpose(n.^2*pi^2/8);

%normalize eigenvectors
for n = 1:4
    v(:,n) = v(:,n)/norm(v(:,n).*v(:,n));
end

% for n = 60:63
%     v(:,n) = v(:,n)/norm((v(:,n).*v(:,n)));
% end

%construct perturb Hamiltonian
Hprime = -subd2dx2/2 + subVprime;
[vprime,dprime] = eig(Hprime);
lambda2 = eig(Hprime);

for n = 1:4
    vprime(:,n) = vprime(:,n)/norm(vprime(:,n).*vprime(:,n));
end

%Show that perturb eigenvalues En(1) = <wave(0)n|H'|wave(0)n>
% M1 = v(:,1)';
% M2 = v(:,1);
% E = M1*Hprime*M2
% dprime(1,1)

Estat = 0;
size_E = size(lambda);
for i = 1:size_E
    M1 = v(:,i)';
    M2 = v(:,i);
    E_p(i) = M1*Hprime*M2;   
    percent = abs(E_p(i)-lambda2(i))/lambda2(i) * 100;
    %E_p(i)/lambda2(i)
    if(percent < 1)
        %E_p(i) = M1*Hprime*M2;   
        Estat = Estat+1;
    end    
end
'The Eigenvalue less than 1% with the accuracy of perturbation theory is'
Estat

k = 0;
mode = 0;
for i = 1:size(lambda)
    bool = abs(lambda(i)-Energy(i))/Energy(i)*100;
    Error(i) = abs(lambda(i)-Energy(i))/Energy(i)*100;
    if(bool < 1)
        mode(i) = i;
        k = k+1;
    end
end

mode;
'the amount of eigenvalues less than 1% is'
k

%plot eigenvectors
hold on;
plot(X,v(:,1),'--r*','LineWidth',1.5);
plot(X,v(:,2),'--b*','LineWidth',1.5);
plot(X,v(:,3),':g*','LineWidth',1.5);
plot(X,v(:,4),'-.black*','LineWidth',1.5);
legend('\Psi_1(x)','\Psi_2(x)','\Psi_3(x)','\Psi_4(x)');
xlabel('x [a]'),ylabel('\Psi(x)');

%plot eigenvectors with Exact
n = 1;
ExactO = cos(n*pi*X/2);
ExactE = sin(n*pi*X/2);
hold on;
%plot(X,ExactE,'-black','LineWidth',1.5);
plot(X,ExactO,'-black','LineWidth',1.5);
plot(X,v(:,n),'--black*','LineWidth',1.5);
legend('Exact','\Psi_1(x)');
xlabel('x [a]'),ylabel('\Psi(x)');

% %spurious modes
n = 63;
plot(X,v(:,n),'-black*','LineWidth',1.5);
xlabel('x [a]'),ylabel('\Psi_{63}(x)');

%plot Eigenvalues
n = 1:size(lambda);
plot(n,Error,'black*');
xlabel('n [n=1,2,3...]'),ylabel('Error [%]');

%plot perturbed and unperturbed
hold on;
for i = 1:4
plot(X,v(:,i),'-black','LineWidth',1.5);
plot(X,vprime(:,i),'--black','LineWidth',1.5);
end
xlabel('x [a]'),ylabel('\Psi(x)');

%plot the perturbed and unperturbed energy
loglog(Energy,E_p,'blacko','LineWidth',1.5);
xlabel('E [hbar^{2}m^{-1}a^{2}]'), ylabel('E [hbar^{2}m^{-1}a^{2}]');
