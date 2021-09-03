clear

load('figure_inputs.mat')
load('figure_inputs_corr.mat')
load('correction_s.mat')
load('correction_p.mat')

% Save figures?
              
inputs(5,:)={'N'};          


%% Perform these calculations for the chosen samples

% Choose samples to loop over (see inputs to determine sample numbers)

for v=2:size(inputs,2)
    
    % Select temperatures to fit cubic baseline to [Celsius]

    initial=300;
    peak1initial=350;

    peak1final=450;
    peak2initial=500;

    peak2final=575;
    final=600;
    
    % Put file paths within 'inputs' into 'data'
    
    data=cell(inputs{3,v},2);      
    for i = 1:size(data,1)          
        data(i,1)={strcat(inputs{1,v},'/',inputs{1,v},inputs{2,v},num2str(i))};  
    end
        
    %% Import DSC data for all runs
    
    % Columns are Temp (C), Time (min), DSC (uV/mg), Sensitivity (uV/mW)

    for i=1:size(data,1)
            cmd=strcat('wc -l',{' <'},data(i,1),'.csv');
            [~,s] = system(cmd{1});
            w=str2num(s);
            range=strcat('41:',num2str(w-1));
            data(i,2)={readmatrix(char(data(i,1)),'FileType','text','Range',range)};
    end
    
    %% Undo temperature correction for Type S runs
    
    if inputs{4,v} == 'S'               
        
        B0=-5077.1;
        B1=608.7;
        B2=-446.0;

        for i=1:size(data,1)
            for j=1:size(data{i,2},1)
                T=data{i,2}(j,1);
                data{i,2}(j,1)=T-(1E-3*B0+1E-5*B1*T+1E-8*B2*T.^2);
            end
        end
        
    else    
    end
    
    %% Calculate sensitivity and store in data{i,2}(:,4)
    
    if inputs{4,v} == 'S'           % Sapphire standard

        P0=mean([288.29999,290.00000,289.10001]);
        P1=mean([961.77399,957.28418,955.18536]);
        P2=mean([0.84398,0.83493,0.89982]);
        P3=mean([-0.40450,-0.43120,-0.38942]);
        P4=mean([-0.44141,-0.42635,-0.38837]);
        P5=mean([0.65314,0.65436,0.58712]);

        for i=1:size(data,1)
            for j=1:size(data{i,2},1)
                z=(data{i,2}(j,1)-P0)/P1;
                data{i,2}(j,4)=(P2+P3*z+P4*z^2+P5*z^3)*exp(-z^2);
            end
        end
        
    elseif inputs{4,v} == 'P'       % Sapphire standard

        P0=mean([223.39999,223.39999,222.50000]);
        P1=mean([689.52338,681.08417,679.69708]);
        P2=mean([4.13972,4.25021,4.55513]);
        P3=mean([-0.02857,0.21801,0.33153]);
        P4=mean([-5.51378,-5.89944,-6.27144]);
        P5=mean([4.58082,4.66168,4.87595]);

        for i=1:size(data,1)
            for j=1:size(data{i,2},1)
                z=(data{i,2}(j,1)-P0)/P1;
                data{i,2}(j,4)=(P2+P3*z+P4*z^2+P5*z^3)*exp(-z^2);
            end
        end
        
    else
    end
    
    %% Restore corrected data (Type P) to 'raw'
    
    if inputs{4,v} == 'P'               
        
        % Import Correction files
        
        for i=1:5
                cmd=strcat('wc -l',{' <'},inputs{6,v},num2str(i),'.csv');
                [~,s] = system(cmd{1});
                w=str2num(s);
                range=strcat('41:',num2str(w-1));
                correction(i,1)={readmatrix(strcat(inputs{6,v},num2str(i),'.csv'),'FileType','text','Range',range)};
        end

        % Add Correction back onto data

        for i=1:size(data,1)
            data{i,2}(1:size(data{i,2},1),3)=data{i,2}(1:size(data{i,2},1),3)+(correction{i,1}(1:size(data{i,2},1),3)/inputs{7,v});
        end        

    else
    end
    
    %% Plot all raw DSC

    initialplotCH(1,'linear','linear')

    for i=1:size(data,1)
        
        % Find indices of data windows
        [rowinitial] = find(data{i,2}(:,1)>initial,1);
        [rowfinal] = find(data{i,2}(:,1)>final,1);
        subset=rowinitial:rowfinal;
        
        set(gca,'ColorOrderIndex',i)
        scatter(data{i,2}(subset,1),data{i,2}(subset,3)./data{i,2}(subset,4),'.','DisplayName',strcat(num2str(i),'-rawDSC'))
        title(inputs{1,v})
        ylabel('\it DSC (mW/mg)','FontSize',20)
        axis([initial final -4 -1])
        
    end

    if inputs{5,v} == 'Y'
        savefig(strcat(inputs{1,v},'/',inputs{1,v},inputs{2,v},'1subset-rawDSC'))
        saveas(gcf,strcat(inputs{1,v},'/',inputs{1,v},inputs{2,v},'1subset-rawDSC'),'tif')
        saveas(gcf,strcat('./1-rawDSC','/',inputs{1,v},inputs{2,v},'1subset-rawDSC'),'tif')
    else
    end
        
    %% %%   Fit cubic polynomial to raw data and subtract
    
    initialplotCH(2,'linear','linear')

    for i=1:size(data,1)
        
        % Find indices of data windows
        [rowinitial] = find(data{i,2}(:,1)>initial,1);
        [rowpeak1initial] = find(data{i,2}(:,1)>peak1initial,1);
        [rowpeak1final] = find(data{i,2}(:,1)>peak1final,1);
        [rowpeak2initial] = find(data{i,2}(:,1)>peak2initial,1);
        [rowpeak2final] = find(data{i,2}(:,1)>peak2final,1);
        [rowfinal] = find(data{i,2}(:,1)>final,1);
        subset=[rowinitial:rowpeak1initial,rowpeak1final:rowpeak2initial,rowpeak2final:rowfinal];

        % Fit cubic polynomial to data
        [baseline]=fit(data{i,2}(subset,1),data{i,2}(subset,3)./data{i,2}(subset,4),'poly3');
        
        % Subtract cubic polynomial from data
        data(i,3)={cat(2,data{i,2}(rowinitial:rowfinal,1),data{i,2}(rowinitial:rowfinal,2),data{i,2}(rowinitial:rowfinal,3)./data{i,2}(rowinitial:rowfinal,4)-baseline(data{i,2}(rowinitial:rowfinal,1)),data{i,2}(rowinitial:rowfinal,4))};
        data_length(i,1)=size(data{i,3},1);
        
        set(gca,'ColorOrderIndex',i)
        plotCH(data{i,3}(:,1),data{i,3}(:,3),'DisplayName',strcat(num2str(i),'-rawsub'))
        title(inputs{1,v})
        ylabel('\it DSC mW/mg','FontSize',20)
        axis([initial final -15E-3 10E-3])

    end

    if inputs{5,v} == 'Y'
        savefig(strcat(inputs{1,v},'/',inputs{1,v},inputs{2,v},'2subset-rawsub-All'))
        saveas(gcf,strcat(inputs{1,v},'/',inputs{1,v},inputs{2,v},'2subset-rawsub-All'),'tif')
        saveas(gcf,strcat('./2-rawsub-All','/',inputs{1,v},inputs{2,v},'2subset-rawsub-All'),'tif')
    else
    end
    
    %% Subtract correction from all runs

    for i=1:size(data,1)
        
        [min_data_length,min_data_length_index]=min(data_length,[],'all','linear');
                
        data{i,4}=cat(2,data{min_data_length_index,3}(:,1),data{min_data_length_index,3}(:,2),interp1(data{i,3}(:,1),data{i,3}(:,3),data{min_data_length_index,3}(:,1)),data{min_data_length_index,3}(:,4));
        
        if inputs{4,v} == 'S'
            interp_corr_s_mean_temp=data{min_data_length_index,3}(:,1);
            interp_corr_s_mean_mW=interp1(correction_s_mean_temp(:,1),correction_s_mean_mW(:,1),data{min_data_length_index,3}(:,1));
            data{i,4}(:,3)=data{i,4}(:,3)-interp_corr_s_mean_mW./inputs{7,v};
        elseif inputs{4,v} == 'P'
            interp_corr_p_mean_temp=data{min_data_length_index,3}(:,1);
            interp_corr_p_mean_mW=interp1(correction_p_mean_temp(:,1),correction_p_mean_mW(:,1),data{min_data_length_index,3}(:,1));
            data{i,4}(:,3)=data{i,4}(:,3)-interp_corr_p_mean_mW./inputs{7,v};
        else
        end
        
    end
    
    %% Plot all Corrected data

    initialplotCH(3,'linear','linear')
    errorincrement=25;
    
    if inputs{4,v} == 'S'
        std_err_corr=interp1(correction_s_mean_temp(:,1),correction_s_mean_mW(:,2)./inputs{7,v},data{1,4}(:,1));
    elseif inputs{4,v} == 'P'
        std_err_corr=interp1(correction_p_mean_temp(:,1),correction_p_mean_mW(:,2)./inputs{7,v},data{1,4}(:,1));
    else
    end

    for i=1:size(data,1)
        set(gca,'ColorOrderIndex',i)
        scatter(data{i,4}(:,1),data{i,4}(:,3),'.','DisplayName',strcat(num2str(i),'-corrsub'))
        set(gca,'ColorOrderIndex',i)
        errorbar(gca,data{i,4}(1:errorincrement:end,1),data{i,4}(1:errorincrement:end,3),std_err_corr(1:errorincrement:end),'LineStyle','none','DisplayName','std\_err\_corr')
        title(inputs{1,v})
        
    end
    axis([initial final -8E-3 12E-3])
    ylabel('\it DSC (mW/mg)','FontSize',20)

    if inputs{5,v} == 'Y'
        savefig(strcat(inputs{1,v},'/',inputs{1,v},inputs{2,v},'3subset-corrsub-All'))
        saveas(gcf,strcat(inputs{1,v},'/',inputs{1,v},inputs{2,v},'3subset-corrsub-All'),'tif')
        saveas(gcf,strcat('./3-corrsub-All','/',inputs{1,v},inputs{2,v},'3subset-corrsub-All'),'tif')
    else
    end
    
     %% Calculate mean of subs for runs 2+

     subs=[];
     
    for i=2:size(data,1)
        subs=cat(2,subs,data{i,4}(:,3));    
    end

    averagedsubs=mean(subs,2);
    std_err_subs=std(subs,0,2)/sqrt(size(subs,2));
    
    std_err_all=sqrt(std_err_corr.^2+std_err_subs.^2);

    %% Plot subs for run 1 and mean of runs 2-5
    
    initialplotCH(4,'linear','linear')
    set(gca,'ColorOrderIndex',1)
    plotCH(data{1,4}(:,1),data{1,4}(:,3),'DisplayName','1-corrsub')
    set(gca,'ColorOrderIndex',1)
    errorbar(gca,data{1,4}(1:errorincrement:end,1),data{1,4}(1:errorincrement:end,3),std_err_corr(1:errorincrement:end),'LineStyle','none','DisplayName','std\_err\_corr')    
    set(gca,'ColorOrderIndex',2)
    plotCH(data{1,4}(:,1),averagedsubs,'DisplayName','(2+)-corrsub')
    set(gca,'ColorOrderIndex',2)
    errorbar(gca,data{1,4}(1:errorincrement:end,1),averagedsubs(1:errorincrement:end),std_err_all(1:errorincrement:end),'Color',[0.9, 0, 0],'DisplayName','std\_err\_all','LineStyle','none')
    title(inputs{1,v})
    ylabel('\it DSC (mW/mg)','FontSize',20)
    axis([initial final -8E-3 12E-3])

    if inputs{5,v} == 'Y'
        savefig(strcat(inputs{1,v},'/',inputs{1,v},inputs{2,v},'4subset-corrsub-Avg'))
        saveas(gcf,strcat(inputs{1,v},'/',inputs{1,v},inputs{2,v},'4subset-corrsub-Avg'),'tif')
        saveas(gcf,strcat('./4-corrsub-Avg','/',inputs{1,v},inputs{2,v},'4subset-corrsub-Avg'),'tif')
    else
    end

    %% Plot difference between run 1 and mean of runs 2-5

    diff_sub=data{1,4}(:,3)-averagedsubs;
    
    initialplotCH(5,'linear','linear')
    plotCH(data{1,4}(:,1),diff_sub,'DisplayName','(1-2+)-diff')
    errorbar(gca,data{1,4}(1:errorincrement:end,1),diff_sub(1:errorincrement:end),std_err_all(1:errorincrement:end),'Color',[0, 0, 0],'DisplayName','std\_err\_all','LineStyle','none')
    title(inputs{1,v})
    ylabel('\it DSC (mW/mg)','FontSize',20)
    axis([initial final -4E-3 10E-3])

    if inputs{5,v} == 'Y'
        savefig(strcat(inputs{1,v},'/',inputs{1,v},inputs{2,v},'5subset-newsub-Diff'))
        saveas(gcf,strcat(inputs{1,v},'/',inputs{1,v},inputs{2,v},'5subset-newsub-Diff'),'tif')
        saveas(gcf,strcat('./5-newsub-Diff','/',inputs{1,v},inputs{2,v},'5subset-newsub-Diff'),'tif')
    else
    end
    
    close all
    
    inputs{8,1}='diff_sub';
    inputs{8,v}=cat(2,data{1,4}(:,1),data{1,4}(:,2),diff_sub);
    inputs{9,1}='errorbar';
    inputs{9,v}=cat(2,data{1,4}(:,1),data{1,4}(:,2),std_err_all);
    
    %% Integrate signal over - two ROIs
%     
%     peak1initial=380;
%     peak1final=470;
%     peak2initial=500;
%     peak2final=590;
% 
%     [rowpeak1initial] = find(data{1,4}(:,1)>peak1initial,1);
%     [rowpeak1final] = find(data{1,4}(:,1)>peak1final,1);
%     [rowpeak2initial] = find(data{1,4}(:,1)>peak2initial,1);
%     [rowpeak2final] = find(data{1,4}(:,1)>peak2final,1);
%     
%     peak1=trapz(60*data{1,4}(rowpeak1initial:rowpeak1final,2),diff_sub(rowpeak1initial:rowpeak1final)); % J/g
%     peak2=trapz(60*data{1,4}(rowpeak2initial:rowpeak2final,2),diff_sub(rowpeak2initial:rowpeak2final)); % J/g
% 
% %     lim=axis;
% %     coordpeak1x=[peak1initial,peak1initial,peak1final,peak1final];
% %     coordpeak1y=[lim(3),lim(4),lim(4),lim(3)];
% %     coordpeak2x=[peak2initial,peak2initial,peak2final,peak2final];
% %     coordpeak2y=[lim(3),lim(4),lim(4),lim(3)];
% %     shadepeak1=patch(coordpeak1x,coordpeak1y,[0 0 0],'DisplayName','Peak 1','FaceAlpha',0.05,'EdgeColor','none');
% %     shadepeak2=patch(coordpeak2x,coordpeak2y,[0 0 0],'DisplayName','Peak 2','FaceAlpha',0.05,'EdgeColor','none');
% 
%     inputs{10,1}='peak1,peak2';
%     inputs{10,v}=cat(1,peak1,peak2);
%     save('inputs.mat','inputs');

    %% Integrate signal over - for all temperatures
    
    peak1initial=380;
    peak2final=590;

    [rowpeak1initial] = find(data{1,4}(:,1)>peak1initial,1);
    [rowpeak2final] = find(data{1,4}(:,1)>peak2final,1);
    
    peak=trapz(60*data{1,4}(rowpeak1initial:rowpeak2final,2),diff_sub(rowpeak1initial:rowpeak2final)); % J/g

%     lim=axis;
%     coordpeak1x=[peak1initial,peak1initial,peak1final,peak1final];
%     coordpeak1y=[lim(3),lim(4),lim(4),lim(3)];
%     shadepeak1=patch(coordpeak1x,coordpeak1y,[0 0 0],'DisplayName','Peak 1','FaceAlpha',0.05,'EdgeColor','none');

    inputs{10,1}='peak';
    inputs{10,v}=peak;
    save('inputs.mat','inputs');
    
end % of v loop = each sample


%% Mean of integrals (peak1, peak2) from each sample
 
% peak1_I=[];
% peak2_I=[];
% peak1_U=[];
% peak2_U=[];
% 
% for i=2:10 % I
%         peak1_I=cat(2,peak1_I,inputs{10,i}(1,1));
%         peak2_I=cat(2,peak2_I,inputs{10,i}(2,1));
% end
% 
% for i=11:19 % U
%         peak1_U=cat(2,peak1_U,inputs{10,i}(1,1));
%         peak2_U=cat(2,peak2_U,inputs{10,i}(2,1));
% end
% 
% meanpeak1_I=mean(peak1_I,2);
% meanpeak2_I=mean(peak2_I,2);
% meanpeak1_U=mean(peak1_U,2);
% meanpeak2_U=mean(peak2_U,2);
% 
% std_err_meanpeak1_I=std(peak1_I,0,2)/sqrt(size(peak1_I,2));
% std_err_meanpeak2_I=std(peak2_I,0,2)/sqrt(size(peak2_I,2));
% std_err_meanpeak1_U=std(peak1_U,0,2)/sqrt(size(peak1_U,2));
% std_err_meanpeak2_U=std(peak2_U,0,2)/sqrt(size(peak2_U,2));
% 
% save('mean+std_err-samples-ROIs.mat','meanpeak1_I','meanpeak1_U','meanpeak2_I','meanpeak2_U','std_err_meanpeak1_I','std_err_meanpeak1_U','std_err_meanpeak2_I','std_err_meanpeak2_U');


%% Mean of integral (peak) from each sample
 
peak_I=[];
peak_U=[];

for i=2:10 % I
        peak_I=cat(2,peak_I,inputs{10,i}(1,1));
end

for i=11:19 % U
        peak_U=cat(2,peak_U,inputs{10,i}(1,1));
end

meanpeak_I=mean(peak_I,2);
meanpeak_U=mean(peak_U,2);

std_err_meanpeak_I=std(peak_I,0,2)/sqrt(size(peak_I,2));
std_err_meanpeak_U=std(peak_U,0,2)/sqrt(size(peak_U,2));

save('mean+std_err-samples-alltemps.mat','meanpeak_I','meanpeak_U','std_err_meanpeak_I','std_err_meanpeak_U');


%% legend for cell array 'inputs'

% inputs{1,1}='sample';
% inputs{2,1}='name';
% inputs{3,1}='runs';
% inputs{4,1}='thermocouple';
% inputs{5,1}='tosave';
% inputs{6,1}='correction';
% inputs{7,1}='mass';
% inputs{8,1}='diff_sub';
% inputs{9,1}='std_err_both';
% inputs{10,1}='peak1,peak2';

%% legend for cell array 'data'

% Columns = Temp (C), Time (min), DSC (*/mg), Sensitivity (uV/mW)

% COLUMN 1 = File Names
% COLUMN 2 = Raw data, DSC (uV/mg)
% COLUMN 3 = Raw data minus cubic fit = 'sub', DSC (mW/mg)
% COLUMN 4 = Corrected data minus the cubic fit, DSC (mW/mg)
