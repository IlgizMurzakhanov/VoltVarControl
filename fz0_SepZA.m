function [z0,a0,sol0] = fz0_SepZA(G,qmax,GEN,X,epsilon,IncremFlag)

clear('yalmip');
z = sdpvar(4*G,1);
a = sdpvar(G,1);
cost = sum(z); %% just some objective which still leaves the obtained solution feasible

% z(1:G) -> v_bar
con = (0.95*ones(G,1) <= z(1:G));
con = con + (z(1:G) <= 1.05*ones(G,1));

% z(G+1:2*G) -> delta
con = con + (zeros(G,1) <= z(G+1:2*G));
con = con + (z(G+1:2*G) <= 0.03*ones(G,1));

% z(2*G+1:3*G) -> sigma
con = con + (z(G+1:2*G) + 0.02*ones(G,1) <= z(2*G+1:3*G));
con = con + (z(2*G+1:3*G) <= 0.18*ones(G,1));

con = con + (z(2*G+1:3*G)-z(G+1:2*G) <= (diag(qmax))*z(3*G+1:4*G)); % z(3*G+1:4*G) -> с
con = con + (z(3*G+1:4*G) >= sum(X(GEN,GEN))'/(1-epsilon)); 
if IncremFlag == 0
    con = con + ((X(GEN,GEN))*a <= (1-epsilon)*ones(G,1)); 
    con = con + (diag(a)*z(3*G+1:4*G) >= ones(G,1)); 
end
settings = sdpsettings('verbose',0,'solver','gurobi');
sol0 = optimize(con,cost,settings);
if ~sol0.problem % checking if solver had any issues
    disp('solved z0')
    z0 = double(z);
    a0 = double(a);
end

end
