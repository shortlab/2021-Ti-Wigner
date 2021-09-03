clear

load('figure_inputs_corr.mat')

%% Set temperature limits and save figures?

% Select temperatures to fit cubic baseline to [Celsius]

initial=300;
peak1initial=350;

peak1final=450;
peak2initial=500;

peak2final=575;
final=600;

% Save figures?

tosave='N';

%% Import all the Corrections

% Columns are Temp (C), Time (min), DSC (uV!), Sensitivity (uV/mW)

    %% Type S sensor

    initialplotCH(1,'linear','linear')

    % Import all corrections and calculate sensitivity
    
    % Sapphire
    P0=mean([288.29999,290.00000,289.10001]);
    P1=mean([961.77399,957.28418,955.18536]);
    P2=mean([0.84398,0.83493,0.89982]);
    P3=mean([-0.40450,-0.43120,-0.38942]);
    P4=mean([-0.44141,-0.42635,-0.38837]);
    P5=mean([0.65314,0.65436,0.58712]);
    
    for i=7:11
        for j=1:5
            cmd=strcat('wc -l',{' <'},inputs_corr{1,i},num2str(j),'.csv');
            [~,s] = system(cmd{1});
            w=str2num(s);
            range=strcat('41:',num2str(w-1));
            correction_s_all(j,i-6)={readmatrix(strcat(inputs_corr{1,i},num2str(j),'.csv'),'FileType','text','Range',range)};
            
            for k=1:size(correction_s_all{j,i-6},1)
                z=(correction_s_all{j,i-6}(k,1)-P0)/P1;
                correction_s_all{j,i-6}(k,4)=(P2+P3*z+P4*z^2+P5*z^3)*exp(-z^2);
            end
        end
    end
    
    for i=1:5
        for j=1:5
            [rowinitial] = find(correction_s_all{j,i}(:,1)>initial,1);
            [rowfinal] = find(correction_s_all{j,i}(:,1)>final,1);
            
            set(gca,'ColorOrderIndex',j)
            scatter(correction_s_all{j,i}(rowinitial:rowfinal,1),correction_s_all{j,i}(rowinitial:rowfinal,3)./correction_s_all{j,i}(rowinitial:rowfinal,4),'.','DisplayName',strcat('S-',num2str(i),'-',num2str(j)))
        end
    end
    
    % Plot all corrections raw
    
    title('Type S Corrections')
    ylabel('$DSC~(mW)$','FontSize',20,'Interpreter','latex')
    if tosave == 'Y'
        savefig('Corrections-S,P/Type-S-Corrections-All')
        saveas(gcf,'Corrections-S,P/Type-S-Corrections-All','tif')
    end
    
    
    %% Fit cubic baselines to everything and subtract
    
    correction_s_mW={};
    correction_s_length=[];
    initialplotCH(2,'linear','linear')
    
    for i=1:5               % Correction number
        for j=2:5           % Run number
            
            % Find indices of data windows
            [rowinitial] = find(correction_s_all{j,i}(:,1)>initial,1);
            [rowpeak1initial] = find(correction_s_all{j,i}(:,1)>peak1initial,1);
            [rowpeak1final] = find(correction_s_all{j,i}(:,1)>peak1final,1);
            [rowpeak2initial] = find(correction_s_all{j,i}(:,1)>peak2initial,1);
            [rowpeak2final] = find(correction_s_all{j,i}(:,1)>peak2final,1);
            [rowfinal] = find(correction_s_all{j,i}(:,1)>final,1);
            subset=[rowinitial:rowfinal,rowpeak1initial,rowpeak1final:rowpeak2initial,rowpeak2final:rowfinal];
            correction_s_length(j-1,i)=size(rowinitial:rowfinal,2);
            
            constructorfit=fit(correction_s_all{j,i}(subset,1),correction_s_all{j,i}(subset,3)./correction_s_all{j,i}(subset,4),'poly3');
            correction_s_temp(j-1,i)={correction_s_all{j,i}(rowinitial:rowfinal,1)};
            correction_s_time(j-1,i)={correction_s_all{j,i}(rowinitial:rowfinal,2)};
            baseline=constructorfit(correction_s_temp{j-1,i});
            correction_s_mW(j-1,i)={correction_s_all{j,i}(rowinitial:rowfinal,3)./correction_s_all{j,i}(rowinitial:rowfinal,4)-baseline};
            set(gca,'ColorOrderIndex',j)
            scatter(correction_s_temp{j-1,i},correction_s_mW{j-1,i},'.','DisplayName',strcat('S-',num2str(i),'-',num2str(j)))
        end
    end
    
    % Plot all correction subtractions
    
    title('Type S Subtractions')
    ylabel('$DSC~(mW)$','FontSize',20,'Interpreter','latex')
    if tosave == 'Y'
        savefig('Corrections-S,P/Type-S-Corrections-Sub')
        saveas(gcf,'Corrections-S,P/Type-S-Corrections-Sub','tif')
    end
    
    
    %% Average the subtractions
    
    cat_correction_s_temp=[];
    cat_correction_s_time=[];
    cat_correction_s_mW=[];
    interpolated_s_temp={};
    interpolated_s_time={};
    interpolated_s_mW={};
    
    [minlength,minlength_index]=min(correction_s_length,[],'all','linear');
    minlength_x=ceil(minlength_index/size(correction_s_length,1));
    minlength_y=mod(minlength_index,size(correction_s_length,1));
    
    for i=1:5
        for j=1:4
            interpolated_s_temp{j,i}=correction_s_temp{minlength_y,minlength_x};
            interpolated_s_time{j,i}=correction_s_time{minlength_y,minlength_x};
            interpolated_s_mW{j,i}=interp1(correction_s_temp{j,i},correction_s_mW{j,i},correction_s_temp{minlength_index});
            cat_correction_s_temp=cat(2,cat_correction_s_temp,interpolated_s_temp{j,i});
            cat_correction_s_time=cat(2,cat_correction_s_time,interpolated_s_time{j,i});
            cat_correction_s_mW=cat(2,cat_correction_s_mW,interpolated_s_mW{j,i});
        end
    end
    
    correction_s_mean_temp=cat(2,mean(cat_correction_s_temp,2),std(cat_correction_s_temp,0,2)/sqrt(20));
    correction_s_mean_time=cat(2,mean(cat_correction_s_time,2),std(cat_correction_s_time,0,2)/sqrt(20));
    correction_s_mean_mW=cat(2,mean(cat_correction_s_mW,2),std(cat_correction_s_mW,0,2)/sqrt(20));

    % Plot average of subtractions
    
    initialplotCH(3,'linear','linear')
    scatter(correction_s_mean_temp(:,1),correction_s_mean_mW(:,1),'.','DisplayName','S-Mean')
    errorinc=50;
    errorbar(gca,correction_s_mean_temp(1:errorinc:end,1),correction_s_mean_mW(1:errorinc:end,1),correction_s_mean_mW(1:errorinc:end,2),'Color',[0, 0, 0],'DisplayName','± std-err','LineStyle','none')
    title('Type S Corrections')
    ylabel('$DSC~(mW)$','FontSize',20,'Interpreter','latex')

    if tosave == 'Y'
        savefig('Corrections-S,P/Type-S-Corrections-Avg')
        saveas(gcf,'Corrections-S,P/Type-S-Corrections-Avg','tif')
    end

    % Export variables
    
    save('correction_s.mat','correction_s_mean_temp','correction_s_mean_mW');
    
    close all
    
    
%% Integrate signal over - two ROIs
%     
%     peak1initial=380;
%     peak1final=470;
%     peak2initial=500;
%     peak2final=590;
% 
%     [rowpeak1initial] = find(correction_s_mean_temp(:,1)>peak1initial,1);
%     [rowpeak1final] = find(correction_s_mean_temp(:,1)>peak1final,1);
%     [rowpeak2initial] = find(correction_s_mean_temp(:,1)>peak2initial,1);
%     [rowpeak2final] = find(correction_s_mean_temp(:,1)>peak2final,1);
% 
%     S_peak1=trapz(60*correction_s_mean_time(rowpeak1initial:rowpeak1final),correction_s_mean_mW(rowpeak1initial:rowpeak1final,1)); % mJ
%     S_peak2=trapz(60*correction_s_mean_time(rowpeak2initial:rowpeak2final),correction_s_mean_mW(rowpeak2initial:rowpeak2final,1)); % mJ
% 
%     S_upper=correction_s_mean_mW(:,1)+correction_s_mean_mW(:,2);
%     S_lower=correction_s_mean_mW(:,1)-correction_s_mean_mW(:,2);
%     S_upper_peak1=trapz(60*correction_s_mean_time(rowpeak1initial:rowpeak1final),S_upper(rowpeak1initial:rowpeak1final)); % mJ
%     S_upper_peak2=trapz(60*correction_s_mean_time(rowpeak2initial:rowpeak2final),S_upper(rowpeak2initial:rowpeak2final)); % mJ
%     S_lower_peak1=trapz(60*correction_s_mean_time(rowpeak1initial:rowpeak1final),S_lower(rowpeak1initial:rowpeak1final)); % mJ
%     S_lower_peak2=trapz(60*correction_s_mean_time(rowpeak2initial:rowpeak2final),S_lower(rowpeak2initial:rowpeak2final)); % mJ
% 
%     S_peak1_uncertainty=S_upper_peak1-S_peak1; % mJ
%     S_peak2_uncertainty=S_upper_peak2-S_peak2; % mJ
    
    
    %% Integrate signal over - for all temperatures
    
    peak1initial=380;
    peak2final=590;

    [rowpeak1initial] = find(correction_s_mean_temp(:,1)>peak1initial,1);
    [rowpeak2final] = find(correction_s_mean_temp(:,1)>peak2final,1);

    S_peak=trapz(60*correction_s_mean_time(rowpeak1initial:rowpeak2final),correction_s_mean_mW(rowpeak1initial:rowpeak2final,1)); % mJ

    S_upper=correction_s_mean_mW(:,1)+correction_s_mean_mW(:,2);
    S_lower=correction_s_mean_mW(:,1)-correction_s_mean_mW(:,2);
    S_upper_peak=trapz(60*correction_s_mean_time(rowpeak1initial:rowpeak2final),S_upper(rowpeak1initial:rowpeak2final)); % mJ
    S_lower_peak=trapz(60*correction_s_mean_time(rowpeak1initial:rowpeak2final),S_lower(rowpeak1initial:rowpeak2final)); % mJ

    S_peak_uncertainty=S_upper_peak-S_peak; % mJ
    

    
%% Type P sensor

initial=300;
peak1initial=350;

peak1final=450;
peak2initial=500;

peak2final=575;
final=600;

    initialplotCH(4,'linear','linear')

    % Import all corrections and calculate sensitivity
  
    % Sapphire
    P0=mean([223.39999,223.39999,222.50000]);
    P1=mean([689.52338,681.08417,679.69708]);
    P2=mean([4.13972,4.25021,4.55513]);
    P3=mean([-0.02857,0.21801,0.33153]);
    P4=mean([-5.51378,-5.89944,-6.27144]);
    P5=mean([4.58082,4.66168,4.87595]);
    
    for i=2:6
        for j=1:5
            cmd=strcat('wc -l',{' <'},inputs_corr{1,i},num2str(j),'.csv');
            [~,s] = system(cmd{1});
            w=str2num(s);
            range=strcat('41:',num2str(w-1));
            correction_p_all(j,i-1)={readmatrix(strcat(inputs_corr{1,i},num2str(j),'.csv'),'FileType','text','Range',range)};
            
            for k=1:size(correction_p_all{j,i-1},1)
                z=(correction_p_all{j,i-1}(k,1)-P0)/P1;
                correction_p_all{j,i-1}(k,4)=(P2+P3*z+P4*z^2+P5*z^3)*exp(-z^2);
            end
        end
    end
    
    for i=1:5
        for j=1:5
            [rowinitial] = find(correction_p_all{j,i}(:,1)>initial,1);
            [rowfinal] = find(correction_p_all{j,i}(:,1)>final,1);
            
            set(gca,'ColorOrderIndex',j)
            scatter(correction_p_all{j,i}(rowinitial:rowfinal,1),correction_p_all{j,i}(rowinitial:rowfinal,3)./correction_p_all{j,i}(rowinitial:rowfinal,4),'.','DisplayName',strcat('P-',num2str(i),'-',num2str(j)))
        end
    end
    
    % Plot all corrections raw
    
%     title('Type P Corrections')
    ylabel('Power (mW)','FontSize',20,'Interpreter','tex')
    legend off
    
    if tosave == 'Y'
        savefig('Corrections-S,P/Type-P-Corrections-All')
        saveas(gcf,'Corrections-S,P/Type-P-Corrections-All','tif')
    end
    
    
    %% Fit cubic baselines to everything and subtract
    
    correction_p_mW={};
    correction_p_length=[];
    initialplotCH(5,'linear','linear')
    
    for i=1:5               % Correction number
        for j=2:5           % Run number
            
            % Find indices of data windows
            [rowinitial] = find(correction_p_all{j,i}(:,1)>initial,1);
            [rowpeak1initial] = find(correction_p_all{j,i}(:,1)>peak1initial,1);
            [rowpeak1final] = find(correction_p_all{j,i}(:,1)>peak1final,1);
            [rowpeak2initial] = find(correction_p_all{j,i}(:,1)>peak2initial,1);
            [rowpeak2final] = find(correction_p_all{j,i}(:,1)>peak2final,1);
            [rowfinal] = find(correction_p_all{j,i}(:,1)>final,1);
            subset=[rowinitial:rowfinal,rowpeak1initial,rowpeak1final:rowpeak2initial,rowpeak2final:rowfinal];
            correction_p_length(j-1,i)=size(rowinitial:rowfinal,2);
            
            constructorfit=fit(correction_p_all{j,i}(subset,1),correction_p_all{j,i}(subset,3)./correction_p_all{j,i}(subset,4),'poly3');
            correction_p_temp(j-1,i)={correction_p_all{j,i}(rowinitial:rowfinal,1)};
            correction_p_time(j-1,i)={correction_p_all{j,i}(rowinitial:rowfinal,2)};
            baseline=constructorfit(correction_p_temp{j-1,i});
            correction_p_mW(j-1,i)={correction_p_all{j,i}(rowinitial:rowfinal,3)./correction_p_all{j,i}(rowinitial:rowfinal,4)-baseline};
            set(gca,'ColorOrderIndex',j)
            scatter(correction_p_temp{j-1,i},correction_p_mW{j-1,i},'.','DisplayName',strcat('P-',num2str(i),'-',num2str(j)))
        end
    end
    
    % Plot all correction subtractions
    
%     title('Type P Subtractions')
    ylabel('Power (mW)','FontSize',20,'Interpreter','tex')
    legend off
    
    if tosave == 'Y'
        savefig('Corrections-S,P/Type-P-Corrections-Sub')
        saveas(gcf,'Corrections-S,P/Type-P-Corrections-Sub','tif')
    end
    
    
    %% Average the subtractions
    
    cat_correction_p_temp=[];
    cat_correction_p_time=[];
    cat_correction_p_mW=[];
    interpolated_p_temp={};
    interpolated_p_time={};
    interpolated_p_mW={};
    
    [minlength,minlength_index]=min(correction_p_length,[],'all','linear');
    minlength_x=ceil(minlength_index/size(correction_p_length,1));
    minlength_y=mod(minlength_index,size(correction_p_length,1));
    
    for i=1:5
        for j=1:4
            interpolated_p_temp{j,i}=correction_p_temp{minlength_y,minlength_x};
            interpolated_p_time{j,i}=correction_p_time{minlength_y,minlength_x};
            interpolated_p_mW{j,i}=interp1(correction_p_temp{j,i},correction_p_mW{j,i},correction_p_temp{minlength_index});
            cat_correction_p_temp=cat(2,cat_correction_p_temp,interpolated_p_temp{j,i});
            cat_correction_p_time=cat(2,cat_correction_p_time,interpolated_p_time{j,i});
            cat_correction_p_mW=cat(2,cat_correction_p_mW,interpolated_p_mW{j,i});
        end
    end
    
    correction_p_mean_temp=cat(2,mean(cat_correction_p_temp,2),std(cat_correction_p_temp,0,2)/sqrt(20));
    correction_p_mean_time=cat(2,mean(cat_correction_p_time,2),std(cat_correction_p_time,0,2)/sqrt(20));
    correction_p_mean_mW=cat(2,mean(cat_correction_p_mW,2),std(cat_correction_p_mW,0,2)/sqrt(20));

    % Plot average of subtractions
    
    initialplotCH(6,'linear','linear')
    scatter(correction_p_mean_temp(:,1),correction_p_mean_mW(:,1),'.','DisplayName','P-Mean')
    errorinc=50;
    errorbar(gca,correction_p_mean_temp(1:errorinc:end,1),correction_p_mean_mW(1:errorinc:end,1),correction_p_mean_mW(1:errorinc:end,2),'Color',[0, 0, 0],'DisplayName','± std-err','LineStyle','none')
%     title('Type P Corrections')
    ylabel('Power (mW)','FontSize',20,'Interpreter','tex')
    legend off

    if tosave == 'Y'
        savefig('Corrections-S,P/Type-P-Corrections-Avg')
        saveas(gcf,'Corrections-S,P/Type-P-Corrections-Avg','tif')
    end

    % Export variables
    
    save('correction_p.mat','correction_p_mean_temp','correction_p_mean_mW');

    close all
    
    
%% Integrate signal over - two ROIs

%     peak1initial=380;
%     peak1final=470;
%     peak2initial=500;
%     peak2final=590;
% 
%     [rowpeak1initial] = find(correction_p_mean_temp(:,1)>peak1initial,1);
%     [rowpeak1final] = find(correction_p_mean_temp(:,1)>peak1final,1);
%     [rowpeak2initial] = find(correction_p_mean_temp(:,1)>peak2initial,1);
%     [rowpeak2final] = find(correction_p_mean_temp(:,1)>peak2final,1);
% 
%     P_peak1=trapz(60*correction_p_mean_time(rowpeak1initial:rowpeak1final),correction_p_mean_mW(rowpeak1initial:rowpeak1final,1)); % mJ
%     P_peak2=trapz(60*correction_p_mean_time(rowpeak2initial:rowpeak2final),correction_p_mean_mW(rowpeak2initial:rowpeak2final,1)); % mJ
% 
%     P_upper=correction_p_mean_mW(:,1)+correction_p_mean_mW(:,2);
%     P_lower=correction_p_mean_mW(:,1)-correction_p_mean_mW(:,2);
%     P_upper_peak1=trapz(60*correction_p_mean_time(rowpeak1initial:rowpeak1final),P_upper(rowpeak1initial:rowpeak1final)); % mJ
%     P_upper_peak2=trapz(60*correction_p_mean_time(rowpeak2initial:rowpeak2final),P_upper(rowpeak2initial:rowpeak2final)); % mJ
%     P_lower_peak1=trapz(60*correction_p_mean_time(rowpeak1initial:rowpeak1final),P_lower(rowpeak1initial:rowpeak1final)); % mJ
%     P_lower_peak2=trapz(60*correction_p_mean_time(rowpeak2initial:rowpeak2final),P_lower(rowpeak2initial:rowpeak2final)); % mJ
% 
%     P_peak1_uncertainty=P_upper_peak1-P_peak1; % mJ
%     P_peak2_uncertainty=P_upper_peak2-P_peak2; % mJ
    
%% Average sensor uncertainty - for two ROIs
   
% load('inputs.mat')
% 
% cat_sample_mass=[];
% 
% for i=2:size(inputs,2)
%     cat_sample_mass=cat(2,cat_sample_mass,inputs{7,i});
% end
% 
% mean_sample_mass=mean(cat_sample_mass,2); % mg
% 
% I_peak1_uncertainty=(3*S_peak1_uncertainty+6*P_peak1_uncertainty)/(9*mean_sample_mass); % J/g
% I_peak2_uncertainty=(3*S_peak2_uncertainty+6*P_peak2_uncertainty)/(9*mean_sample_mass); % J/g
% U_peak1_uncertainty=(9*P_peak1_uncertainty)/(9*mean_sample_mass); % J/g
% U_peak2_uncertainty=(9*P_peak2_uncertainty)/(9*mean_sample_mass); % J/g
% 
% save('mean+std_err-corrections-ROIs.mat','I_peak1_uncertainty','I_peak2_uncertainty','U_peak1_uncertainty','U_peak2_uncertainty');
% 
% 

    %% Integrate signal over - for all temperatures

    peak1initial=380;
    peak2final=590;

    [rowpeak1initial] = find(correction_p_mean_temp(:,1)>peak1initial,1);
    [rowpeak2final] = find(correction_p_mean_temp(:,1)>peak2final,1);

    P_peak=trapz(60*correction_p_mean_time(rowpeak1initial:rowpeak2final),correction_p_mean_mW(rowpeak1initial:rowpeak2final,1)); % mJ

    P_upper=correction_p_mean_mW(:,1)+correction_p_mean_mW(:,2);
    P_lower=correction_p_mean_mW(:,1)-correction_p_mean_mW(:,2);
    P_upper_peak=trapz(60*correction_p_mean_time(rowpeak1initial:rowpeak2final),P_upper(rowpeak1initial:rowpeak2final)); % mJ
    P_lower_peak=trapz(60*correction_p_mean_time(rowpeak1initial:rowpeak2final),P_lower(rowpeak1initial:rowpeak2final)); % mJ

    P_peak_uncertainty=P_upper_peak-P_peak; % mJ


%% Average sensor uncertainty - for all temperatures
    
load('inputs.mat')

cat_sample_mass=[];

for i=2:size(inputs,2)
    cat_sample_mass=cat(2,cat_sample_mass,inputs{7,i});
end

mean_sample_mass=mean(cat_sample_mass,2); % mg

I_peak_uncertainty=(3*S_peak_uncertainty+6*P_peak_uncertainty)/(9*mean_sample_mass); % J/g
U_peak_uncertainty=(9*P_peak_uncertainty)/(9*mean_sample_mass); % J/g

save('mean+std_err-corrections-alltemps.mat','I_peak_uncertainty','U_peak_uncertainty');

 