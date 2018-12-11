function zout = differententropTC(~,zn,znplus)

global m1;
global m2;
global ca;
global cb;
global lamb_a0;
global lamb_b0;
global k;
global count1;

%[q1x;q1y;q2x;q2y;p1x;p1y;p2x;p2y;sa;sb] = z;

zout=ones(10,1);
middlezn = 0.5*(zn+znplus);

zout(1) = middlezn(5)./m1;
zout(2) = middlezn(6)./m1;
zout(3) = middlezn(7)./m2;
zout(4) = middlezn(8)./m2;

ea_cap_di_pi = @(pi1,pi5) (ca/2).*(log(sqrt(pi1)/lamb_a0)).*(log(sqrt(pi1)/lamb_a0)) ...
        +60.*log(sqrt(pi1)/lamb_a0)+1500.*(exp((pi5-0.2.*log(sqrt(pi1)/lamb_a0)./5))-1);

eb_cap_di_pi = @(pi2,pi6) (cb/2).*(log(sqrt(pi2)/lamb_b0)).*(log(sqrt(pi2)/lamb_b0)) ...
        +60.*log(sqrt(pi2)/lamb_b0)+1500.*(exp((pi6-0.2.*log(sqrt(pi2)/lamb_b0)./5))-1);

D1e = @(pi12,c,lamb_0,pi56) 0.5*(1/(pi12)).*(c*log(sqrt(pi12)/lamb_0)+60-60 ...
    *exp((pi56-0.2*log(sqrt(pi12)/lamb_0))./5));

D2e = @(pi56,pi12,lamb_0) 300.*exp((pi56-0.2*log(sqrt(pi12)/lamb_0))/5); % a seconda che io metta o meno 5 il grafico risulta costante o meno

PIn = simmetryvariable(zn);
PInplus = simmetryvariable(znplus);

if abs(PInplus(1)-PIn(1)) ~= 0 % a seconda della tolleranza impostata cambia il grafico
	Dpi1ea = 0.5*(ea_cap_di_pi(PInplus(1),PInplus(5))-ea_cap_di_pi(PIn(1),PInplus(5)))./(PInplus(1)-PIn(1)) ...
		 + 0.5*(ea_cap_di_pi(PInplus(1),PIn(5))-ea_cap_di_pi(PIn(1),PIn(5)))/(PInplus(1)-PIn(1));

	%Dpi1eb = 0.5*(eb_cap_di_pi(PInplus(1),PInplus(6))-eb_cap_di_pi(PIn(1),PInplus(6)))/(PInplus(1)-PIn(6)) ...
		 %+ 0.5*(eb_cap_di_pi(PInplus(1),PIn(6))-eb_cap_di_pi(PIn(1),PIn(6)))/(PInplus(1)-PIn(1));		 

else
	Dpi1ea = 0.5.*(D1e((PIn(1)+PInplus(1))*0.5,ca,lamb_a0,PInplus(5)) ...
		+ D1e((PIn(1)+PInplus(1))*0.5,ca,lamb_a0,PIn(5))); count1=count1+1; %controllo
end

if abs(PInplus(2)-PIn(2)) ~= 0
	Dpi2eb = 0.5*(eb_cap_di_pi(PInplus(2),PInplus(6))-eb_cap_di_pi(PIn(2),PInplus(6)))/(PInplus(2)-PIn(2)) ...
		 + 0.5*(eb_cap_di_pi(PInplus(2),PIn(6))-eb_cap_di_pi(PIn(2),PIn(6)))/(PInplus(2)-PIn(2));
else
	Dpi2eb = 0.5*(D1e((PIn(2)+PInplus(2))*0.5,cb,lamb_b0,PInplus(6)) ...
		+D1e((PIn(2)+PInplus(2))*0.5,cb,lamb_b0,PIn(6)));
end

zout([5 6]) = -middlezn([1 2]).*2.*Dpi1ea -(middlezn([1 2])-middlezn([3 4])).*2.*Dpi2eb;

zout([7 8]) = -(middlezn([3 4])-middlezn([1 2])).*2*Dpi2eb;

if abs(PInplus(5)-PIn(5)) ~= 0
	theta_a_star = 0.5*(ea_cap_di_pi(PInplus(1),PInplus(5))-ea_cap_di_pi(PInplus(1),PIn(5)))/(PInplus(5)-PIn(5)) ...
		 + 0.5*(ea_cap_di_pi(PIn(1),PInplus(5))-ea_cap_di_pi(PIn(1),PIn(5)))/(PInplus(5)-PIn(5));
		 
else
	theta_a_star = 0.5*(D2e(0.5*(PInplus(5)+PIn(5)),PInplus(1),lamb_a0) ...
		+D2e(PInplus(5)+PIn(5),PIn(1),lamb_a0));
end

if abs(PInplus(6)-PIn(6)) ~= 0
	theta_b_star = 0.5*(eb_cap_di_pi(PInplus(2),PInplus(6))-eb_cap_di_pi(PInplus(2),PIn(6)))/(PInplus(6)-PIn(6)) ...
		 + 0.5*(eb_cap_di_pi(PIn(2),PInplus(6))-eb_cap_di_pi(PIn(2),PIn(6)))/(PInplus(6)-PIn(6));
else
	theta_b_star = 0.5*(D2e(0.5*(PInplus(6)+PIn(6)),PInplus(2),lamb_b0) ...
		+D2e(PInplus(6)+PIn(6),PIn(2),lamb_b0));
end

zout(9)=k*((theta_b_star/theta_a_star)-1);

zout(10)=k*((theta_a_star/theta_b_star)-1);


end
