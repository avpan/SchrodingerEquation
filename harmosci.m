clear all;clc;clf;
n = 1; %n = 1,2,3 eigenstates
N = 64; %N initialize
L = 20;
h = L/(N);
xo = -L/2;
x = xo;

%find wave function, psi for N where psi = psi_n
for i = 1:N+1
    if(i == 1 || i == N+1)
        X(i) = x;
    else
        X(i) = x;
    end
    x = x + h;
end
x = x-h;
X = transpose(X);

%initialize potential matrix
for i = 1:N+1
    V(i,i) = X(i)^2;
end

% %using 6th order difference we create the d^2/dx^2
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

% %I must delete first row/colum and last row/column
%X = removerows(X,N+1);
%X = removerows(X,1);
subV = removerows(V,N+1);
subV = removerows(subV.',N+1)'; %delete columns and convert back to its original 
subV = removerows(subV.',1)';
subV = removerows(subV,1);
V;
subV;

%add perturbation
for i = 1:size(V)
    Vprime(i,i) = V(i,i) + X(i)^4;
end
Vprime;

subd2dx2 = removerows(d2dx2,N+1);
subd2dx2 = removerows(subd2dx2.',N+1)';
subd2dx2 = removerows(subd2dx2.',1)';
subd2dx2 = removerows(subd2dx2,1);
d2dx2;
-subd2dx2;

%construct the Hamiltonian 
%satisfys H|v> = lambda|v>
H = -d2dx2/2 + V/2;
[v,d] = eig(H);
lambda = eig(H);
n = 0:(size(lambda)-1);
Energy = n+1/2;

%normalize the eigenvector
for n = 1:4
    v(:,n) = v(:,n)/norm((v(:,n).*v(:,n)));
end

%construct perturb Hamiltonian
Hprime = -d2dx2/2 + Vprime;
[vprime,dprime] = eig(Hprime);
lambda2 = eig(Hprime);

%normalize the eigenvector
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
    percent(i) = abs(E_p(i)-lambda2(i))/lambda2(i) * 100;
    if(percent(i) < 10)
        Estat = Estat+1;
    end    
end
'The Eigenvalue less than 1% with the accuracy of perturbation theory is'
Estat

k = 0;
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

%plot first four eigenvectors
% hold on;
% plot(X,v(:,1),'blacko','LineWidth',1.5);
% plot(X,v(:,2),'blackx','LineWidth',1.5);
% plot(X,v(:,3),'black+','LineWidth',1.5);
% plot(X,v(:,4),'black*','LineWidth',1.5);
% legend('\Psi_0(x)','\Psi_1(x)','\Psi_2(x)','\Psi_3(x)')
% xlabel('x [a]'),ylabel('\Psi(x)');

%Hermite Polynomial exact solution
% wave_0 = exp(-X.^2/2);
% wave_1 = sqrt(2)*X.*exp(-X.^2/2);

%plot Hermite poly. 
% hold on;
% % plot(X,wave_0,'-black','LineWidth',1.5); 
% % plot(X,v(:,1)/pi^(.25),'--blacko','LineWidth',1.5);
% 
% plot(X,wave_1,'-black','LineWidth',1.5); 
% plot(X,v(:,2)/pi^(.25),'--blacko','LineWidth',1.5);
% xlabel('x [a]'),ylabel('\Psi(x) [hbar(m*\omega_{0})^{-1/4}]');

% %plot Eigenvalues
% n = 0:(size(lambda)-1);
% plot(n,Error,'black*');
% xlabel('n [n=0,1,2...]'),ylabel('Error [%]');

% %plot perturbed and unperturbed
% hold on;
% for i = 1:3
% plot(X,v(:,i),'-black','LineWidth',1.5);
% plot(X,vprime(:,i),'--black','LineWidth',1.5);
% end
% xlabel('x [a]'),ylabel('\Psi(x)');

%plot the perturbed and unperturbed energy
% loglog(Energy,E_p,'blacko','LineWidth',1.5);
% xlabel('E [hbar\omega_{0}]'), ylabel('E [hbar\omega_{0}]');