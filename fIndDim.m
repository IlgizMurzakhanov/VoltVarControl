function [NBusIncSla,NGenIncSla,NBusExcSla,NGenExcSla,IndGenExcSla,IndGenExcSlaOrd] = fIndDim(pf_caseA)
    % Defining dimensions
    NBusIncSla = size(pf_caseA.bus,1);
    NGenIncSla = size(pf_caseA.gen,1);

    NBusExcSla = NBusIncSla-1;
    NGenExcSla = NGenIncSla-1;

    IndGenExcSla = pf_caseA.gen(2:end,1); % generator indices
    IndGenExcSlaOrd = sort(IndGenExcSla);
end

