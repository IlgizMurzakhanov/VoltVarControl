function [PgExcSla_orig,QgExcSla_orig,Qg_orig,PcExcSla_orig,QcExcSla_orig,...
    VmPFExcSla_orig,pfCaseA] = fOrigPgQgPcQcVmPF(NBusExcSla,NBusIncSla,caseA,IndGenExcSlaInt,G,mpopt)
    % Creating vectors
    PgExcSla_orig = zeros(NBusExcSla,1);
    QgExcSla_orig = zeros(NBusExcSla,1);

    % Generation
    PgExcSla_orig(IndGenExcSlaInt) = caseA.gen(2:G+1,2);
    QgExcSla_orig(IndGenExcSlaInt) = caseA.gen(2:G+1,3);
    Qg_orig = QgExcSla_orig(IndGenExcSlaInt);

    % Loads
    PcExcSla_orig = caseA.bus(2:NBusIncSla,3);
    QcExcSla_orig = caseA.bus(2:NBusIncSla,4);

    % Power flow voltages in buses except slack
    pfCaseA = runpf(caseA,mpopt);
    VmPFExcSla_orig = pfCaseA.bus(2:end,8); 
end

