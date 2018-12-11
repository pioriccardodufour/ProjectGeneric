function y = mid_p_r_entr_var(odefun, y0, t)

h = diff(t);
y(:,1) = y0;

    for n = 1:length(t)-1
        y(:,n+1) = fsolve(@(z) -z+y(:,n)+h(n)*odefun(t(n),y(:,n).*0.5 + z.*0.5), y(:,n), ...
            optimset('Display', 'off', 'MaxIter', 1000));
    end

end
