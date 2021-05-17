function xi = CasadiDyn(qddFunc,xi,ui, dt) 

q = xi(1:ceil(length(xi)/2));
qd = xi(ceil(length(xi)/2)+1:end);
qdd = full(qddFunc(q,qd,ui));
xi = xi + [qd;qdd]*dt;


end