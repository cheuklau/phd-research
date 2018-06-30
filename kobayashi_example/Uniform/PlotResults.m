clear;
clc;

% Reference solution GLC (263168 directions)
Reference = 3.3822878E-05;

% LS solutions
LS = [8.8151441308E-05, ...
0.0000000000E+00, ...
0.0000000000E+00, ...
2.0581266835E-04, ...
1.1287227730E-04, ...
5.6653601168E-05];

% Calculate LS relative error
LSerr = abs(LS - Reference) / Reference;

% Calculate LS mesh size
LSmesh = [24, 80, 168, 288, 440, 624];
for i = 1 : size(LSmesh, 2)
    LSmesh(i) = 1 / sqrt(LSmesh(i));
end

% GC solutions
GC = [0.0000000000E+00, ...
3.5356839205E-05, ...
3.2568070259E-05, ...
3.3695009664E-05, ...
3.3803432064E-05, ...
3.3883981204E-05, ...
3.3839185490E-05];

% GC relative error
GCerr = abs(GC - Reference) / Reference;

% GC mesh size
GCmesh = [80, ...
1088, ...
4224, ...
5928, ...
16640, ...
37248, ...
66048];

for i = 1 : size(GCmesh, 2)
    GCmesh(i) = 1 / sqrt(GCmesh(i));    
end

% QR solutions
QR = [1.6126298032E-04, ...
0.0000000000E+00, ...
9.7859397268E-06, ...
1.0481314523E-05, ...
1.1005904967E-04, ...
2.6657165639E-05];

% QR relative error
QRerr = abs(QR - Reference) / Reference;

% QR mesh size
QRmesh = [24, ...
80, ...
168, ...
288, ...
440, ...
624];
for i = 1 : size(QRmesh, 2)
    QRmesh(i) = 1 / sqrt(QRmesh(i));
end

% LDFE-ST solutions
LDFEST = [0.0000000000E+00, ...
5.9413788125E-06, ...
5.2162578229E-05, ...
1.8184284573E-05, ...
3.5895422534E-05, ...
3.5098913465E-05, ...
3.4034194573E-05, ...
3.3786160734E-05];

% LDFE-ST relative error
LDFESTerr = abs(LDFEST - Reference) / Reference;

% LDFE-ST mesh size
LDFESTmesh = [32, 128, 512, 2048, 8192, 32768, 131072, 524288];
for i = 1 : size(LDFESTmesh, 2)
   LDFESTmesh(i) = 1 / sqrt(LDFESTmesh(i)); 
end

%{
% LDFE-SQ-sa solutions
LDFESQSA = [0.0000000000E+00, ...
1.2568327554E-04, ...
6.2427745876E-05, ...
5.0862222832E-05, ...
3.5182641978E-05, ...
1.6686526493E-05, ...
4.2474608971E-05, ...
3.6295860644E-05, ...
3.2262680821E-05, ...
3.7326280409E-05, ...
3.1789441544E-05, ...
3.3760430056E-05, ...
3.3882649307E-05, ...
3.3639217921E-05];
%}

% LDFE-SQ-sa solutions
LDFESQSA = [0.0000000000E+00, ... y
6.2427745876E-05, ...y
5.0862222832E-05, ...y
3.2262680821E-05, ...y
3.3760430056E-05, ...y
3.3639217921E-05];

% LDFE-SQ-sa relative error
LDFESQSAerr = abs(LDFESQSA - Reference) / Reference;

% LDFE-SQ mesh size
LDFESQmesh = [96, ...y
864, ...y
1536, ...y
7776, ...y
42336, ...y
161376];
for i = 1 : size(LDFESQmesh, 2)
   LDFESQmesh(i) = 1 / sqrt(LDFESQmesh(i));     
end
%{
% QDFE-SQ-SA solutions
QDFESQSA = [1.5312032073E-04, ...
8.4165989231E-05, ...
4.0811020753E-05, ...
3.9690853051E-05, ...
3.8219698633E-05, ...
3.4684836085E-05, ...
3.5967733255E-05, ...
3.1765746194E-05, ...
3.5695805155E-05, ...
3.4923080232E-05, ...
3.2335767585E-05, ...
3.2705191733E-05, ...
3.3870390920E-05, ...
3.2838293518E-05];
%}
% QDFE-SQ-SA solutions
QDFESQSA = [1.5312032073E-04, ...y
8.4165989231E-05, ...y
4.0811020753E-05, ...y
3.4684836085E-05, ...y
3.2335767585E-05, ...y
3.3870390920E-05, ...y
3.2838293518E-05];

% QDFE-SQ-SA relative error
QDFESQSAerr = abs(QDFESQSA - Reference) / Reference;

% Calculate QDFE-SQ mesh size
QDFESQmesh = [216, ...y
864, ...y
1944, ...y
7776, ...y
26136, ...y
95256, ...y
146016];
for i = 1 : size(QDFESQmesh, 2)
    QDFESQmesh(i) = 1 / sqrt(QDFESQmesh(i));
end

% Fourth order mesh size
fourthMesh = [32, 216, 1944, 5400, 10584, 17496, 26136, 55296, 95256, 146016];
for i = 1 : size(fourthMesh, 2)
    fourthMesh(i) = 1 / sqrt(fourthMesh(i));
end

% Plot results for (35, 95, 35) position using QDFE as reference solution
slope = 2;
intercept = log((GCerr(1)+LSerr(1)+QRerr(1)+LDFESTerr(1)+LDFESQSAerr(1)+QDFESQSAerr(1))/6);
figure
LSplot = loglog(LSmesh, LSerr, 'm-s');
set(LSplot, 'LineWidth', 1.25, 'MarkerSize', 10, 'color', [0 0.5 0]);
hold on
GCplot = loglog(GCmesh, GCerr, 'm-+');
set(GCplot, 'LineWidth', 1.25, 'MarkerSize', 10, 'color', [0 0.5 0]);
hold on
QRplot = loglog(QRmesh, QRerr, 'm-*');
set(QRplot, 'LineWidth', 1.25, 'MarkerSize', 10, 'color', [0 0.5 0]);
hold on
LDFESTplot = loglog(LDFESTmesh, LDFESTerr, 'b-s');
set(LDFESTplot, 'LineWidth', 1.25, 'MarkerSize', 10);
hold on
LDFESQSAplot = loglog(LDFESQmesh, LDFESQSAerr, 'k-s');
set(LDFESQSAplot, 'LineWidth', 1.25, 'MarkerSize', 10);
hold on
QDFESQSAplot = loglog(QDFESQmesh, QDFESQSAerr, 'k-+');
set(QDFESQSAplot, 'LineWidth', 1.25, 'MarkerSize', 10);
hold on
fourthOrder = loglog(fourthMesh, exp(slope .* log(fourthMesh / LDFESTmesh(1)) + intercept), 'r');
set(fourthOrder, 'LineWidth', 1.25);
grid on
set(gca,'FontSize', 14)
xlabel('Mesh Length', 'FontSize', 18)
ylabel('|Relative Error|', 'FontSize', 18)
legend('LS', 'GC', 'QR', 'LDFE-ST', 'LDFE-SQ', ...
     'QDFE-SQ', '2nd-order', 'Location', 'Southeast')