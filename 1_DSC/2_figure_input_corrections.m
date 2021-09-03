
%% Type S sensor

inputs_corr{1,11}='Z-CorrectionS/CorrS-6_1-50Kmin-PtY2O3-Up-';
inputs_corr{2,11}='S';
inputs_corr{3,11}='Correction-S-6';
inputs_corr{4,11}=5;

inputs_corr{1,10}='Z-CorrectionS/CorrS-5_1-50Kmin-PtY2O3-Up-';
inputs_corr{2,10}='S';
inputs_corr{3,10}='Correction-S-5';
inputs_corr{4,10}=5;

inputs_corr{1,9}='Z-CorrectionS/CorrS-4_1-50Kmin-PtY2O3-Up-';
inputs_corr{2,9}='S';
inputs_corr{3,9}='Correction-S-4';
inputs_corr{4,9}=5;

inputs_corr{1,8}='Z-CorrectionS/CorrS-3_1-50Kmin-PtY2O3-Up-';
inputs_corr{2,8}='S';
inputs_corr{3,8}='Correction-S-3';
inputs_corr{4,8}=5;

inputs_corr{1,7}='Z-CorrectionS/CorrS-2_1-50Kmin-PtY2O3-Up-';
inputs_corr{2,7}='S';
inputs_corr{3,7}='Correction-S-2';
inputs_corr{4,7}=5;

%% Type P sensor

inputs_corr{1,6}='Z-CorrectionP/Corr3-8_1-50Kmin-PtY2O3-Up-';
inputs_corr{2,6}='P';
inputs_corr{3,6}='Correction-P-3-8';
inputs_corr{4,6}=5;

inputs_corr{1,5}='Z-CorrectionP/Corr3-3_1-50Kmin-PtY2O3-Up-';
inputs_corr{2,5}='P';
inputs_corr{3,5}='Correction-P-3-3';
inputs_corr{4,5}=5;

inputs_corr{1,4}='Z-CorrectionP/Corr3-2_1-50Kmin-PtY2O3-Up-';
inputs_corr{2,4}='P';
inputs_corr{3,4}='Correction-P-3-2';
inputs_corr{4,4}=5;

inputs_corr{1,3}='Z-CorrectionP/Corr3-1_1-50Kmin-PtY2O3-Up-';
inputs_corr{2,3}='P';
inputs_corr{3,3}='Correction-P-3-1';
inputs_corr{4,3}=5;

inputs_corr{1,2}='Z-CorrectionP/2-Correction-2-PtRh+Y2O3-Up-';
inputs_corr{2,2}='P';
inputs_corr{3,2}='Correction-P-2';
inputs_corr{4,2}=5;

%% Row labels

inputs_corr{1,1}='path';
inputs_corr{2,1}='sensor';
inputs_corr{3,1}='graphtitle';
inputs_corr{4,1}='runs';

save('figure_inputs_corr.mat','inputs_corr');
