function y = monolitichTC(odefun, y0, t)

h = diff(t);
y(:,1) = y0;

for n = 1:length(t)-1
	y(:,n+1) = fsolve(@(z) -z+y(:,n)+h(n)*odefun(t(n),y(:,n),z), y(:,n), ...
		optimset('Display', 'off', 'MaxIter', 1000));
end

end
