function [Qs,Vs,sol1,Cost,StartQP,EndQP,qU_hat] = fInPr(G,S,tV,c,X,GEN,delta,qmax,Agen,bv,Vr,i,Qs_prev,sigma,StartQP,EndQP)
% fInPr(G,T,VTildeExcSla,c,XBusExcSla,IndGenExcSlaInt,delta,qU,Agen);

    Cost = [];

    % initialize matrices with equilibrium injections and voltages (solutions of inner problem)
    Qs = zeros(G,S);

    % Defining qmax as qhat (new formulation)
    qU_hat = (ones(G,1)./c).*(sigma - delta);
    qmax = qU_hat;
    
    for scenario = 1:S
        vtilde = tV(:,scenario);
        clear('yalmip');

        %% Variables
        q = sdpvar(G,1);

        %% Cost function
        cost = q'*(diag(c)+X(GEN,GEN))*q/2 + norm(diag(delta/2)*q,1) + q'*(vtilde(GEN)-bv); 

        %% Constraints
        con = (-qmax <= q <= qmax);

        %% Settings and run
    %     settings = sdpsettings('verbose',0,'solver','sdpt3');

%         if i == 1
            settings = sdpsettings('verbose',0,'solver','gurobi');

            %     disp('Started solving S QPs')
            StartQP(scenario*i) = toc;

            sol1 = optimize(con,cost,settings);

            %     disp('Completed solving S QPs')
            EndQP(scenario*i) = toc;

%         elseif i > 1
%             settings = sdpsettings('verbose',0,'solver','gurobi','usex0',1);
%             assign(q,Qs_prev(:,scenario));
%             sol1 = optimize(con,cost,settings);
%         end
                  
        if ~sol1.problem % checking if solver had any issues
            Qs(:,scenario) = double(q);
            Cost = [Cost; double(cost)];
        end
    end

    Vs = X*Agen*Qs + tV;
    
end

