clear

load('figure_inputs_corr.mat')

%% Set temperature limits and save figures?

% Select temperatures [Celsius]

initial=300;
final=990;

% Save figures?

tosave='N';

%% Import all the Corrections

% Columns are Temp (C), Time (min), DSC (uV/mg), Sensitivity (uV/mW)

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
        for j=5:5
            cmd=strcat('wc -l',{' <'},inputs_corr{1,i},num2str(j),'.csv');
            [~,s] = system(cmd{1});
            w=str2num(s);
            range=strcat('41:',num2str(w-1));
            correction_s_all(j-4,i-6)={readmatrix(strcat(inputs_corr{1,i},num2str(j),'.csv'),'FileType','text','Range',range)};
            
            for k=1:size(correction_s_all{j-4,i-6},1)
                z=(correction_s_all{j-4,i-6}(k,1)-P0)/P1;
                correction_s_all{j-4,i-6}(k,4)=(P2+P3*z+P4*z^2+P5*z^3)*exp(-z^2);
            end
        end
    end
    
    for i=1:5
        for j=1
            [rowinitial] = find(correction_s_all{j,i}(:,1)>initial,1);
            [rowfinal] = find(correction_s_all{j,i}(:,1)>final,1);
            
            scatter(correction_s_all{j,i}(rowinitial:rowfinal,1),correction_s_all{j,i}(rowinitial:rowfinal,3)./correction_s_all{j,i}(rowinitial:rowfinal,4),'.','DisplayName',strcat('S-',num2str(i),'-',num2str(j+4)))
        end
    end
    
    % Plot all corrections raw
    
    title('Type S DSC')
    ylabel('$DSC~(mW)$','FontSize',20,'Interpreter','latex')
    if tosave == 'Y'
        savefig('6-AB/Correction-AB/Type-S-AB-Corrections-DSC')
        saveas(gcf,'6-AB/Correction-AB/Type-S-AB-Corrections-DSC','tif')
    end
    
    
    %%
    
    correction_s_mW={};
    correction_s_length=[];
    
    for i=1:5               % Correction number
        for j=1             % Run number
            
            % Find indices of data windows
            [rowinitial] = find(correction_s_all{j,i}(:,1)>initial,1);
            subset=[rowinitial:rowfinal];
            correction_s_length(j,i)=size(rowinitial:rowfinal,2);
            
            correction_s_temp(j,i)={correction_s_all{j,i}(rowinitial:rowfinal,1)};
            correction_s_mW(j,i)={correction_s_all{j,i}(rowinitial:rowfinal,3)./correction_s_all{j,i}(rowinitial:rowfinal,4)};%-baseline};
        end
    end
       
    %% Average the subtractions
    
    cat_correction_s_temp=[];
    cat_correction_s_mW=[];
    interpolated_s_temp={};
    interpolated_s_mW={};
    
    [minlength,minlength_index]=min(correction_s_length,[],'all','linear');
    minlength_x=ceil(minlength_index/size(correction_s_length,1));
    
    for i=1:5
        for j=1:1
            interpolated_s_temp{j,i}=correction_s_temp{1,minlength_x};
            interpolated_s_mW{j,i}=interp1(correction_s_temp{j,i},correction_s_mW{j,i},correction_s_temp{minlength_index});
            cat_correction_s_temp=cat(2,cat_correction_s_temp,interpolated_s_temp{j,i});
            cat_correction_s_mW=cat(2,cat_correction_s_mW,interpolated_s_mW{j,i});
        end
    end
    
    correction_s_mean_temp=cat(2,mean(cat_correction_s_temp,2),std(cat_correction_s_temp,0,2)/sqrt(20));
    correction_s_mean_mW=cat(2,mean(cat_correction_s_mW,2),std(cat_correction_s_mW,0,2)/sqrt(20));

    % Plot average of subtractions
    
    initialplotCH(3,'linear','linear')
    scatter(correction_s_mean_temp(:,1),correction_s_mean_mW(:,1),'.','DisplayName','S-Mean')
    errorinc=50;
    errorbar(gca,correction_s_mean_temp(1:errorinc:end,1),correction_s_mean_mW(1:errorinc:end,1),correction_s_mean_mW(1:errorinc:end,2),'Color',[0, 0, 0],'DisplayName','± std-err','LineStyle','none')
    title('Type S Averaged')
    ylabel('$DSC~(mW)$','FontSize',20,'Interpreter','latex')

    if tosave == 'Y'
        savefig('6-AB/Correction-AB/Type-S-AB-Corrections-Avg')
        saveas(gcf,'6-AB/Correction-AB/Type-S-AB-Corrections-Avg','tif')
    end

    % Export variables
    
    save('correction_ab_s.mat','correction_s_mean_temp','correction_s_mean_mW');

%% Type P sensor

    initialplotCH(4,'linear','linear')

    % Import all corrections and calculate sensitivity
  
%     Sapphire
    P0=mean([223.39999,223.39999,222.50000]);
    P1=mean([689.52338,681.08417,679.69708]);
    P2=mean([4.13972,4.25021,4.55513]);
    P3=mean([-0.02857,0.21801,0.33153]);
    P4=mean([-5.51378,-5.89944,-6.27144]);
    P5=mean([4.58082,4.66168,4.87595]);
    
    for i=2:6
        for j=5:5
            cmd=strcat('wc -l',{' <'},inputs_corr{1,i},num2str(j),'.csv');
            [~,s] = system(cmd{1});
            w=str2num(s);
            range=strcat('41:',num2str(w-1));
            correction_p_all(j-4,i-1)={readmatrix(strcat(inputs_corr{1,i},num2str(j),'.csv'),'FileType','text','Range',range)};
            
            for k=1:size(correction_p_all{j-4,i-1},1)
                z=(correction_p_all{j-4,i-1}(k,1)-P0)/P1;
                correction_p_all{j-4,i-1}(k,4)=(P2+P3*z+P4*z^2+P5*z^3)*exp(-z^2);
            end
        end
    end
    
    for i=1:5
        for j=1:1
            [rowinitial] = find(correction_p_all{j,i}(:,1)>initial,1);
            [rowfinal] = find(correction_p_all{j,i}(:,1)>final,1);
            
            scatter(correction_p_all{j,i}(rowinitial:rowfinal,1),correction_p_all{j,i}(rowinitial:rowfinal,3)./correction_p_all{j,i}(rowinitial:rowfinal,4),'.','DisplayName',strcat('P-',num2str(i),'-',num2str(j+4)))
        end
    end
    
    % Plot all corrections raw
    
    title('Type P DSC')
    ylabel('$DSC~(mW)$','FontSize',20,'Interpreter','latex')
    if tosave == 'Y'
        savefig('6-AB/Correction-AB/Type-P-AB-Corrections-DSC')
        saveas(gcf,'6-AB/Correction-AB/Type-P-AB-Corrections-DSC','tif')
    end
    
    
    %% Fit cubic baselines to everything and subtract
    
    correction_p_mW={};
    correction_p_length=[];
    
    for i=1:5               % Correction number
        for j=1:1           % Run number
            
            % Find indices of data windows
            [rowinitial] = find(correction_p_all{j,i}(:,1)>initial,1);
            [rowfinal] = find(correction_p_all{j,i}(:,1)>final,1);
            subset=[rowinitial:rowfinal];
            correction_p_length(j,i)=size(rowinitial:rowfinal,2);
            
            correction_p_temp(j,i)={correction_p_all{j,i}(rowinitial:rowfinal,1)};
            correction_p_mW(j,i)={correction_p_all{j,i}(rowinitial:rowfinal,3)./correction_p_all{j,i}(rowinitial:rowfinal,4)};%-baseline};

         end
    end    
    
    %% Average the subtractions
    
    cat_correction_p_temp=[];
    cat_correction_p_mW=[];
    interpolated_p_temp={};
    interpolated_p_mW={};
    
    [minlength,minlength_index]=min(correction_p_length,[],'all','linear');
    minlength_x=ceil(minlength_index/size(correction_p_length,1));
        
    for i=1:5
        for j=1:1
            interpolated_p_temp{j,i}=correction_p_temp{1,minlength_x};
            interpolated_p_mW{j,i}=interp1(correction_p_temp{j,i},correction_p_mW{j,i},correction_p_temp{minlength_index});
            cat_correction_p_temp=cat(2,cat_correction_p_temp,interpolated_p_temp{j,i});
            cat_correction_p_mW=cat(2,cat_correction_p_mW,interpolated_p_mW{j,i});
        end
    end
    
    correction_p_mean_temp=cat(2,mean(cat_correction_p_temp,2),std(cat_correction_p_temp,0,2)/sqrt(20));
    correction_p_mean_mW=cat(2,mean(cat_correction_p_mW,2),std(cat_correction_p_mW,0,2)/sqrt(20));

    % Plot average of subtractions
    
    initialplotCH(6,'linear','linear')
    scatter(correction_p_mean_temp(:,1),correction_p_mean_mW(:,1),'.','DisplayName','P-Mean')
    errorinc=50;
    errorbar(gca,correction_p_mean_temp(1:errorinc:end,1),correction_p_mean_mW(1:errorinc:end,1),correction_p_mean_mW(1:errorinc:end,2),'Color',[0, 0, 0],'DisplayName','± std-err','LineStyle','none')
    title('Type P Corrections')
    ylabel('$DSC~(mW)$','FontSize',20,'Interpreter','latex')

    if tosave == 'Y'
        savefig('6-AB/Correction-AB/Type-P-AB-Corrections-Avg')
        saveas(gcf,'6-AB/Correction-AB/Type-P-AB-Corrections-Avg','tif')
    end

    % Export variables
    
    save('correction_ab_p.mat','correction_p_mean_temp','correction_p_mean_mW');

%     close all
    