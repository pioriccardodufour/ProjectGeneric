function PI = simmetryvariable(z)
	
	%global m1 m2;
	%e = @(lamb,s) 0.05.*log((lamb)/lamb_0).*log((lamb)/lamb_0)...
        %+60*log((lamb)/lamb_0)+1500*(exp((s-0.2.*log((lamb)/lamb_0)/5))-1);

	%q1x = z(1);
	%q1y = z(2);
	%q2x = z(3);
	%q2y = z(4);
	%p1x = z(5);
	%p1y = z(6);
	%p2x = z(7);
	%p2y = z(8);
	%sa = z(9);
	%sb = z(10);

	%[q1x;q1y;q2x;q2y;p1x;p1y;p2x;p2y;sa;sb] = z;

	pi1 = dot(z([1 2]),z([1 2]));
	pi2 = dot(z([3 4])-z([1 2]),z([3 4])-z([1 2]));
	pi3 = dot(z([5 6]),z([5 6]));
	pi4 = dot(z([7 8]),z([7 8]));
	pi5 = z(9); %sa
	pi6 = z(10);
	pi7 = dot(z([1 2]),z([5 6])); 
	pi8 = dot(z([3 4]),z([7 8])); 
	pi9 = dot(z([1 2]),z([3 4]));

	PI = [pi1;pi2;pi3;pi4;pi5;pi6;pi7;pi8;pi9];

	%gradPI(1,:) = [2*q1x 2*q1y zeros(1,8)];
	%gradPI(2,:) = [2*(q1x - q2x) 2*(q1y - q2y) 2*(q2x - q1x) 2*(q2y - q1y) zeros(1,6)];
	%gradPI(3,:) = [0 0 0 0 2*p1x 2*p1y 0 0 0 0];
	%gradPI(4,:) = [0 0 0 0 0 0 2*p2x 2*p2y 0 0];
	%gradPI(5,:) = [0 0 0 0 0 0 0 0 1 0];
	%gradPI(6,:) = [0 0 0 0 0 0 0 0 0 1];
	%gradPI(7,:) = [p1x p1y 0 0 0 q1x q1y 0 0 0 0];
	%gradPI(8,:) = [0 0 p1x p1y 0 0 q2x q2y 0 0];
	%gradPI(9,:) = [q2x q2y q1x q1y 0 0 0 0 0 0];

	%gradPItrasp = gradPI';

	%Etilde = pi3/(2*m1) + pi4/(2*m2) + feval(e,pi1,pi5) + feval(e,pi2,pi6);

	%Stilde = pi5 + pi6;

	

end	