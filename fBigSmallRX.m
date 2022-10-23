function [RBusExcSla,XBusExcSla,IndGenExcSlaInt,R,X,XColRed] = fBigSmallRX(caseA,IndGenExcSlaOrd)

% Via approach 1: via YBUS built-in function
% Extract resistance and reactance matrices
[YBUS, YF, YT] = makeYbus(caseA);
YBusExcSla = YBUS(2:end,2:end);
ZBusExcSla = inv(YBusExcSla);

RBusExcSla = real(ZBusExcSla);
XBusExcSla = imag(ZBusExcSla);

% Via Vassilis Kekatos's lecture, pages 11, 4
% https://www.faculty.ece.vt.edu/kekatos/pdsa/Lecture11.pdf
% Branch-bus incidence matrix
Ainc = makeIncidence(caseA);

% Reduced branch-bus incidence matrix
A = Ainc(:,2:end); % if get mistake in this line, then case has switched-off lines

% Helping matrix F
F = inv(A);

% Extract lines' r and x
r = caseA.branch(:,3);
x = caseA.branch(:,4);

% Create diagonal matrices
Dr = diag(r);
Dx = diag(x);

% Create R and X matrices
RBusExcSla_viaF = F*Dr*transpose(F);
XBusExcSla_viaF = F*Dx*transpose(F);

% APPROACHES 1 AND 2 ARE EQUIVALENT
RBusExcSla == RBusExcSla_viaF;
XBusExcSla == XBusExcSla_viaF;

% Internal index (numbering without slack bus)
IndGenExcSlaInt = IndGenExcSlaOrd - 1; 

% Reduced R, X matrices
R = RBusExcSla(IndGenExcSlaInt,IndGenExcSlaInt);
X = XBusExcSla(IndGenExcSlaInt,IndGenExcSlaInt);

R_viaF = RBusExcSla_viaF(IndGenExcSlaInt,IndGenExcSlaInt);
X_viaF = XBusExcSla_viaF(IndGenExcSlaInt,IndGenExcSlaInt);

% Partly reduced (in columns) X
XColRed = XBusExcSla(:,IndGenExcSlaInt);

end

