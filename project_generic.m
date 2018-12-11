% pendulum example
% masses m1,m2 okay
% springs' internal energy e1,e2
% positions r1,r2 okay
% momenta q1,q2
% springs' entropy sa,sb
% springs' natural length lamb_a0,lamb_b0
% masses' kinetic energy kin1,kin2
% springs' length lamb_a,lamb_b
% springs' absolute temperature theta_a,theta_b
% springs' free energy psi_a,psi_b
% conductivity constant k

global m1;
global m2;
global ca;
global cb;
global lamb_a0;
global lamb_b0;
global k;
global mod;
global entropy;
global dedlambda;
global count1;
count1=0;
m1=1;
m2=2;
ca=0.1;
cb=1;
lamb_a0=2;
lamb_b0=1;
k=300;
% time interval T
T = 0:0.3:500;
r1 = zeros(2,1);
r2 = zeros(2,1);
p1 = zeros(2,1);
p2 = zeros(2,1);
r1(:,1) = [12;0];
p1(:,1) = [0;2]; %quantità di moto iniziale massa 1
r2(:,1) = [2.2;0];
p2(:,1) = [1;-2]; %quantità di moto iniziale massa 2
theta_a = zeros(1,length(T));
theta_b = zeros(1,length(T));
theta_a(1) = 290;
theta_b(1) = 442;

dedlambda = @(modq,c,lamb_0,s) (1/(modq))*(c*log(modq/lamb_0)+60-60*...
    exp((s-0.2*log(modq/lamb_0))/5));

mod = @(x,y) sqrt((x^2)+(y^2));
entropy = @(s,modq,lamb_0) s-0.2*log(modq/lamb_0);
thetafun = @(s,modq,lamb_0) 300.*exp((s-0.2*log(modq/lamb_0))./5);
sa = 5*log(theta_a(1)/300)+0.2*log(mod(r1(1,1),r1(2,1))/lamb_a0);
sb = 5*log(theta_b(1)/300)+0.2*log(mod(r2(1,1)-r1(1,1),r2(2,1)-r1(2,1))/lamb_b0);

e = @(modq,c,lamb_0,s) (c/2).*log(modq/lamb_0).*log(modq/lamb_0)...
	+60*log(modq/lamb_0)+1500.*(exp((s-0.2.*log(modq/lamb_0))./5)-1);

%psi_ = @(modq,c,lamb_0,theta) (c/2).*log(modq/lamb_0).*log(modq/lamb_0) ...
    %-0.2*(theta-300)*log(modq/lamb_0)+5*(theta-300-theta*log(theta/300));
    
%variabili di stato
z=[r1;r2;p1;p2;sa;sb];

%[t,y]=ode45(@different,[0 25],z(:,1));

%z=monolitichTC(@differententropTC,z(:,1),T);

z=mid_p_r_entr_var(@differententrop,z(:,1),T);
    
lamb_a = zeros(1,length(T));
lamb_b = zeros(1,length(T));
Etot = zeros(1,length(T));

%definisco variabili da plottare
for i = 1:length(T)
	lamb_a(i) = mod(z(1,i),z(2,i));
	lamb_b(i) = mod(z(3,i)-z(1,i),z(4,i)-z(2,i));
	theta_a(i) = thetafun(z(9,i),mod(z(1,i),z(2,i)),lamb_a0);
    theta_b(i) = thetafun(z(10,i),mod(z(3,i)-z(1,i),z(4,i)-z(2,i)),lamb_b0);
	Etot(i) = ((mod(z(5,i),z(6,i))^(2))./(2*m1)) + ((mod(z(7,i),z(8,i))^(2))./(2*m2)) ...
		+ e(lamb_a(i),ca,lamb_a0,z(9,i)) + e(lamb_b(i),cb,lamb_b0,z(10,i)); %...
        %+ psi_(lamb_a(i),ca,lamb_a0,theta_a(i)) + psi_(lamb_b(i),cb,lamb_b0,theta_b(i));
end

% rappresentazione grafica

subplot(2,2,1);   %energy
hold on;
grid on;
axis([0 500 10 5000]);
plot(T,Etot);
xlabel('time[t]');
ylabel('Energy[e]');

subplot(2,2,2);   %entropy
hold on;
grid on;
axis([0 500 0 10]);
plot(T,z(9,:)+z(10,:));
xlabel('time[t]');
ylabel('Entropy[s]');

subplot(2,2,3);   %temperatures
hold on;
grid on;
axis([0 500 100 500]);
plot(T,theta_a);
plot(T,theta_b);
xlabel('time[t]');
ylabel('Temperatures[theta]'); % ?

subplot(2,2,4);   %log strains
hold on;
grid on;
axis([0 500 0 10]);
plot(T,log(lamb_a));
plot(T,log(lamb_b));  % - 0.55.*ones(1,length(T))
xlabel('time[t]');
ylabel('Log strains');

%plot eventuale delle posizioni

%{
figure
hold on;
grid on;
for i=1:15 %length(T)-1
    title('Panoramica posizioni');
    xplot=[0;z(1,i);z(3,i)];
    yplot=[0;z(2,i);z(4,i)];
    plot(z(1,i),z(2,i),'.','Color','b','MarkerSize',30); %massa 1
    plot(z(3,i),z(4,i),'.','Color','g','MarkerSize',30); %massa 2
    plot(xplot,yplot,'k'); %molle  
end
%}