function [znew] = differententrop(~,zold)

global m1;
global m2;
global ca;
global cb;
global lamb_a0;
global lamb_b0;
global k;
global dedlambda;
global mod;
global entropy;

znew=zeros(10,1);
znew(1)=zold(5)./m1;
znew(2)=zold(6)./m1;
znew(3)=zold(7)./m2;
znew(4)=zold(8)./m2;

znew([5 6]) = -dedlambda(mod(zold(1),zold(2)),ca,lamb_a0,zold(9)).*zold([1 2])./(mod(zold(1),zold(2))) ...
	-dedlambda(mod(zold(3)-zold(1),zold(4)-zold(2)),cb,lamb_b0,zold(10)).*(zold([1 2])-zold([3 4])) ...
	./mod(zold(3)-zold(1),zold(4)-zold(2));

znew([7 8])= -dedlambda(mod(zold(3)-zold(1),zold(4)-zold(2)),cb,lamb_b0,zold(10)) ...
	.*(zold([3 4])-zold([1 2]))./mod(zold(3)-zold(1),zold(4)-zold(2));

znew(9)=k*(exp((entropy(zold(10),mod(zold(3)-zold(1),zold(4)-zold(2)),lamb_b0)-...
    entropy(zold(9),mod(zold(1),zold(2)),lamb_a0))/5)-1);

znew(10)=k*(exp((entropy(zold(9),mod(zold(1),zold(2)),lamb_a0)-...
    entropy(zold(10),mod(zold(3)-zold(1),zold(4)-zold(2)),lamb_b0))/5)-1);


end

