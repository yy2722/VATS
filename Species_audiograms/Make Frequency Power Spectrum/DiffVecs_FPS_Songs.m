%Manually load vectors to run Diff Analyses
handles.save = 'C:\Users\Woolley Lab-MB\Desktop\MATLAB Testing\Figs\Analyses\CrossCorrelations\FPS\diffVec\'; %where cross correlation results will be saved to

%Manually assign files to be run in Diff Analyses
handles.C = 'CB_Calls Only';
%handles.C = 'CB_Songs Only';
handles.D = 'GW_Vocalizations';
handles.G = 'LF_Calls Only';
%handles.G = 'LF_Songs Only';
handles.J = 'LW_Calls Only';
%handles.J = 'LW_Songs Only';
handles.M = 'OF_Calls Only';
%handles.M = 'OF_Songs Only';
handles.P = 'RF_Calls Only';
%handles.P = 'RF_Songs Only';
handles.S = 'SF_Calls Only';
%handles.S = 'SF_Songs Only';
handles.V = 'ZF_Calls Only';
%handles.V = 'ZF_Songs Only';

%Manually loads and assigns variables to vetors in Excel for cross-correlation
C = xlsread('C:\Users\Woolley Lab-MB\Desktop\Multispecies Vocalizations\FPS\Multispecies FPS Vectors.xlsx','All Species','B2:XQ2'); %CB_Calls Only
%C = xlsread('C:\Users\Woolley Lab-MB\Desktop\Multispecies Vocalizations\FPS\Multispecies FPS Vectors.xlsx','All Species','B3:XQ3'); %CB_Songs Only
D = xlsread('C:\Users\Woolley Lab-MB\Desktop\Multispecies Vocalizations\FPS\Multispecies FPS Vectors.xlsx','All Species','B4:XQ4'); %GW_Vocalizations
G = xlsread('C:\Users\Woolley Lab-MB\Desktop\Multispecies Vocalizations\FPS\Multispecies FPS Vectors.xlsx','All Species','B5:XQ5'); %LF_Calls Only
%G = xlsread('C:\Users\Woolley Lab-MB\Desktop\Multispecies Vocalizations\FPS\Multispecies FPS Vectors.xlsx','All Species','B6:XQ6'); %LF_Songs Only
J = xlsread('C:\Users\Woolley Lab-MB\Desktop\Multispecies Vocalizations\FPS\Multispecies FPS Vectors.xlsx','All Species','B7:XQ7'); %LW_Calls Only
%J = xlsread('C:\Users\Woolley Lab-MB\Desktop\Multispecies Vocalizations\FPS\Multispecies FPS Vectors.xlsx','All Species','B8:XQ8'); %LW_Song Only
M = xlsread('C:\Users\Woolley Lab-MB\Desktop\Multispecies Vocalizations\FPS\Multispecies FPS Vectors.xlsx','All Species','B9:XQ9'); %OF_Calls Only
%M = xlsread('C:\Users\Woolley Lab-MB\Desktop\Multispecies Vocalizations\FPS\Multispecies FPS Vectors.xlsx','All Species','B10:XQ10'); %OF_Songs Only
P = xlsread('C:\Users\Woolley Lab-MB\Desktop\Multispecies Vocalizations\FPS\Multispecies FPS Vectors.xlsx','All Species','B11:XQ11'); %RF_Calls Only
%P = xlsread('C:\Users\Woolley Lab-MB\Desktop\Multispecies Vocalizations\FPS\Multispecies FPS Vectors.xlsx','All Species','B12:XQ12'); %RF_Songs Only
S = xlsread('C:\Users\Woolley Lab-MB\Desktop\Multispecies Vocalizations\FPS\Multispecies FPS Vectors.xlsx','All Species','B13:XQ13'); %SF_Calls Only
%S = xlsread('C:\Users\Woolley Lab-MB\Desktop\Multispecies Vocalizations\FPS\Multispecies FPS Vectors.xlsx','All Species','B14:XQ14'); %SF_Songs Only
V = xlsread('C:\Users\Woolley Lab-MB\Desktop\Multispecies Vocalizations\FPS\Multispecies FPS Vectors.xlsx','All Species','B15:XQ15'); %ZF_Calls Only
%V = xlsread('C:\Users\Woolley Lab-MB\Desktop\Multispecies Vocalizations\FPS\Multispecies FPS Vectors.xlsx','All Species','B16:XQ16'); %ZF_Songs Only


%Semi-auto make diff vectors for songs
dCD = C - D;
dCG = C - G;
dCJ = C - J;
dCM = C - M;
dCP = C - P;
dCS = C - S;
dCV = C - V;

dDG = D - G;
dDJ = D - J;
dDM = D - M;
dDP = D - P;
dDS = D - S;
dDV = D - V;

dGJ = G - J;
dGM = G - M;
dGP = G - P;
dGS = G - S;
dGV = G - V;

dJM = J - M;
dJP = J - P;
dJS = J - S;
dJV = J - V;

dMP = M - P;
dMS = M - S;
dMV = M - V;

dPS = P - S;
dPV = P - V;

dSV = S - V;

%Square values
sdCD = dCD.^2;
sdCG = dCG.^2;
sdCJ = dCJ.^2;
sdCM = dCM.^2;
sdCP = dCP.^2;
sdCS = dCS.^2;
sdCV = dCV.^2;

sdDG = dDG.^2;
sdDJ = dDJ.^2;
sdDM = dDM.^2;
sdDP = dDP.^2;
sdDS = dDS.^2;
sdDV = dDV.^2;

sdGJ = dGJ.^2;
sdGM = dGM.^2;
sdGP = dGP.^2;
sdGS = dGS.^2;
sdGV = dGV.^2;

sdJM = dJM.^2;
sdJP = dJP.^2;
sdJS = dJS.^2;
sdJV = dJV.^2;

sdMP = dMP.^2;
sdMS = dMS.^2;
sdMV = dMV.^2;

sdPS = dPS.^2;
sdPV = dPV.^2;

sdSV = dSV.^2;

%Take sum of squared values
SSsdCD = sum(sdCD);
SSsdCG = sum(sdCG);
SSsdCJ = sum(sdCJ);
SSsdCM = sum(sdCM);
SSsdCP = sum(sdCP);
SSsdCS = sum(sdCS);
SSsdCV = sum(sdCV);

SSsdDG = sum(sdDG);
SSsdDJ = sum(sdDJ);
SSsdDM = sum(sdDM);
SSsdDP = sum(sdDP);
SSsdDS = sum(sdDS);
SSsdDV = sum(sdDV);

SSsdGJ = sum(sdGJ);
SSsdGM = sum(sdGM);
SSsdGP = sum(sdGP);
SSsdGS = sum(sdGS);
SSsdGV = sum(sdGV);

SSsdJM = sum(sdJM);
SSsdJP = sum(sdJP);
SSsdJS = sum(sdJS);
SSsdJV = sum(sdJV);

SSsdMP = sum(sdMP);
SSsdMS = sum(sdMS);
SSsdMV = sum(sdMV);

SSsdPS = sum(sdPS);
SSsdPV = sum(sdPV);

SSsdSV = sum(sdSV);


% %Take sqrt of sum of squared values
% sqrtSS = sqrt(SS);
sqrtSSsdCD = sqrt(SSsdCD);
sqrtSSsdCG = sqrt(SSsdCG);
sqrtSSsdCJ = sqrt(SSsdCJ);
sqrtSSsdCM = sqrt(SSsdCM);
sqrtSSsdCP = sqrt(SSsdCP);
sqrtSSsdCS = sqrt(SSsdCS);
sqrtSSsdCV = sqrt(SSsdCV);

sqrtSSsdDG = sqrt(SSsdDG);
sqrtSSsdDJ = sqrt(SSsdDJ);
sqrtSSsdDM = sqrt(SSsdDM);
sqrtSSsdDP = sqrt(SSsdDP);
sqrtSSsdDS = sqrt(SSsdDS);
sqrtSSsdDV = sqrt(SSsdDV);

sqrtSSsdGJ = sqrt(SSsdGJ);
sqrtSSsdGM = sqrt(SSsdGM);
sqrtSSsdGP = sqrt(SSsdGP);
sqrtSSsdGS = sqrt(SSsdGS);
sqrtSSsdGV = sqrt(SSsdGV);

sqrtSSsdJM = sqrt(SSsdJM);
sqrtSSsdJP = sqrt(SSsdJP);
sqrtSSsdJS = sqrt(SSsdJS);
sqrtSSsdJV = sqrt(SSsdJV);

sqrtSSsdMP = sqrt(SSsdMP);
sqrtSSsdMS = sqrt(SSsdMS);
sqrtSSsdMV = sqrt(SSsdMV);

sqrtSSsdPS = sqrt(SSsdPS);
sqrtSSsdPV = sqrt(SSsdPV);

sqrtSSsdSV = sqrt(SSsdSV);

 
%Save
save([handles.save '_DiffVecResults.mat'],'SSsd*','sqrtSSsd*'); %saves Sum of Squares and sqrt of SS
