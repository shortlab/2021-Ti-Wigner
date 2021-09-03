clear

load('figure_inputs.mat')
load('figure_inputs_corr.mat')
load('correction_ab_s.mat')
load('correction_ab_p.mat')

%% Set temperature limits and save figures?

% Select temperatures to fit baseline to [Celsius]

initial=600;
peak1initial=800;
peak1final=970;
final=989;

% Save figures?
              
tosave='N'; 

avg_aB=[];

%% Perform these calculations for the chosen samples

% Choose samples to loop over (see inputs to determine sample numbers)

for v=2:size(inputs,2)
    
    % Put file paths within 'inputs' into 'data'
    
    if inputs{4,v} == 'S'
        offset=4;   % The first run to start evaluating = 2nd 1000C heat.
        runs=2;     % Number of runs to 1000C to evaluate.
    elseif inputs{4,v} == 'P'
        offset=6;
        runs=3;
    else
    end
    
    if strcmp(inputs{1,v},'UA2-3')
        offset=5;
        runs=1;
    else
    end
    
    data=cell(runs,2);      
    for i = offset:offset+runs-1    
        data(i-offset+1,1)={strcat(inputs{1,v},'/',inputs{1,v},inputs{2,v},num2str(i))};  
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
        
        for i=5:5
                cmd=strcat('wc -l',{' <'},inputs{6,v},num2str(i),'.csv');
                [~,s] = system(cmd{1});
                w=str2num(s);
                range=strcat('41:',num2str(w-1));
                correction(i-4,1)={readmatrix(strcat(inputs{6,v},num2str(i),'.csv'),'FileType','text','Range',range)};
        end

        % Add Correction back onto data

        for i=1:size(data,1)
            data{i,2}(1:size(data{i,2},1),3)=data{i,2}(1:size(data{i,2},1),3)+(correction{1,1}(1:size(data{i,2},1),3)/inputs{7,v});
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
        
%         set(gca,'ColorOrderIndex',i)
%         plotCH(data{i,2}(subset,1),data{i,2}(subset,3)./data{i,2}(subset,4),'DisplayName',strcat(num2str(i),'-raw'))
%         title(inputs{1,v})
%         ylabel('\it DSC (mW/mg)','FontSize',20)
%         axis([initial final -2 0])
        
    end

    %% %%   Save subset of data within temperature window
    
    for i=1:size(data,1)
        
        % Find indices of data windows
        [rowinitial] = find(data{i,2}(:,1)>initial,1);
        [rowfinal] = find(data{i,2}(:,1)>final,1);
        subset=[rowinitial:rowfinal];

        data(i,3)={cat(2,data{i,2}(rowinitial:rowfinal,1),data{i,2}(rowinitial:rowfinal,2),data{i,2}(rowinitial:rowfinal,3)./data{i,2}(rowinitial:rowfinal,4),data{i,2}(rowinitial:rowfinal,4))};
        data_length(i,1)=size(data{i,3},1);
        
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

    for i=1:size(data,1)
        
%         plotCH(data{i,4}(:,1),data{i,4}(:,3),'DisplayName',strcat(num2str(i),'-corr'))
        
    end    
    
    %% Fit linear baseline to subset and subtract from data
    
    for i=1:size(data,1)
    
        [rowinitial] = find(data{i,4}(:,1)>initial,1)+1;
        [rowpeak1initial] = find(data{i,4}(:,1)>peak1initial,1);
        [rowpeak1final] = find(data{i,4}(:,1)>peak1final,1);
        [rowfinal] = find(data{i,4}(:,1)>final,1)-1;
        subset=[rowinitial:rowpeak1initial,rowpeak1final:rowfinal];
        
        [baseline]=fit(data{i,4}(subset,1),data{i,4}(subset,3),'poly1');

        data(i,5)={cat(2,data{i,4}(rowinitial:rowfinal,1),data{i,4}(rowinitial:rowfinal,2),data{i,4}(rowinitial:rowfinal,3)-baseline(data{i,4}(rowinitial:rowfinal,1)),data{i,4}(rowinitial:rowfinal,4))};
        set(gca,'ColorOrderIndex',i+1)
        plotCH(data{i,5}(:,1),data{i,5}(:,3),'DisplayName',strcat(num2str(i),'-sub'))
    
    end
    
    %% Integrate aB
    
    integral_initial=820;
    integral_final=980;
    axis([initial final -2.5 0.5])
    aB_peak=[];
    
    for i=1:size(data,1)
    
        [rowintinitial] = find(data{i,5}(:,1)>integral_initial,1);
        [rowintfinal] = find(data{i,5}(:,1)>integral_final,1);
        
        newcolors=[0, 0, 0;0.9, 0, 0;0, 0, 0.9;0.4, 0.8, 0.2;0.6, 0.2, 0.8;];
        aB_peak(i,1)=trapz(60*data{i,5}(rowintinitial:rowintfinal,2),data{i,5}(rowintinitial:rowintfinal,3)); % J/g
        shadeintegral=patch(data{i,5}(rowintinitial:rowintfinal,1),data{i,5}(rowintinitial:rowintfinal,3),newcolors(i+1,:),'DisplayName',strcat(num2str(i+1),'-\alpha/\beta'),'FaceAlpha',0.2);

        txt=strcat('\DeltaH_{\alpha/\beta}^{',num2str(offset-1+i),'} =',{' '},num2str(round(aB_peak(i,1),1)),' J/g');
        lim=axis;
        X=lim(1)+0.1*(lim(2)-lim(1));
        Y=lim(3)+(0.8-0.1*i)*(lim(4)-lim(3));

        text(X,Y,txt,'FontSize',15,'Interpreter','tex','Color',newcolors(i+1,:))
        
    end
    
    mean_aB_peak=cat(2,mean(aB_peak),std(aB_peak));
    txt=strcat('\DeltaH_{\alpha/\beta}^{avg.}',' =',{' '},num2str(round(mean_aB_peak(1),1)),'\pm',num2str(round(mean_aB_peak(2),1)),{' '},' J/g');
    Y=lim(3)+(0.8-0.1*(i+2))*(lim(4)-lim(3));
    text(X,Y,txt,'FontSize',15,'Interpreter','tex','Color',newcolors(1,:),'FontWeight','normal')
    
    title(inputs{1,v},'FontSize',20,'FontWeight','normal')
    legend off
    
    
    if tosave == 'Y'
%         savefig(strcat(inputs{1,v},'/',inputs{1,v},inputs{2,v},'AB'))
%         saveas(gcf,strcat(inputs{1,v},'/',inputs{1,v},inputs{2,v},'AB'),'tif')
        saveas(gcf,strcat('./6-AB','/',inputs{1,v},inputs{2,v},'AB'),'tif')
    else
    end
    
    name_aB{v-1,1}=inputs{1,v};
    avg_aB=cat(1,avg_aB,mean_aB_peak);
    
%     close all
    
end

%% Plot all averages

I=9;
U=9; % excluding UA2-3 as only heated to 1000C once.

I_avg_aB=round(mean(avg_aB(1:I,1)),1);
U_avg_aB=round(mean(avg_aB(I+1:I+U,1)),1);
I_avg_aB_u=round(sqrt((std(avg_aB(1:I,1)))^2+(mean(avg_aB(1:I,2)))^2),1);
U_avg_aB_u=round(sqrt((std(avg_aB(I+1:I+U,1)))^2+(mean(avg_aB(I+1:I+U,2)))^2),1);

initialplotCH(2,'linear','linear')
plot(avg_aB(:,1),'x','DisplayName','mean\_aB','LineWidth',1.5)
errorbar(avg_aB(:,1),avg_aB(:,2),'Color',[0, 0, 0],'DisplayName','± sigma','LineStyle','none')

plot([0:size(avg_aB,1)+2]',-87.11*ones(size(avg_aB,1)+3,1),'--','Color',0.6*ones(1,3),'LineWidth',1,'DisplayName','Ref.')
patch([0,0,20,20],-[91.47,82.75,82.75,91.47],'k','FaceAlpha',0.05,'DisplayName','± error','EdgeColor','none')
text(19.5,-87.11,'\DeltaH_{\alpha/\beta}^{Ref.}','FontSize',15,'Interpreter','tex','Color',0*ones(1,3));

xticks(1:size(avg_aB))
xtickangle(90)
set(gca,'xticklabel',name_aB,'XMinorTick','off')
ylabel('Specific enthalpy \alpha/\beta (J/g)','FontSize',20)
xlabel('')
lim=axis;
axis([lim(1) lim(2)-1 -105 -65])
lim=axis;
legend off
% title('   Irradiated                    Unirradiated','FontWeight','normal','FontSize',15,'HorizontalAlignment','center')
% plot(((I+U+1)/2)*ones(1,size(lim(3):lim(4),2)),lim(3):lim(4),'.','Color',0.6*ones(1,3),'LineWidth',1,'DisplayName','Ref.')


% table1=strcat('\begin{tabular}{lll}','&\bf J/g & $\pm$ ',' \\ I&',num2str(I_avg_aB),'&',num2str(I_avg_aB_u),'\\ U&',num2str(U_avg_aB),'&',num2str(U_avg_aB_u),'\\ \end{tabular}');
% annotation('textbox', [0.73 0.73 0.18 0.2], 'string', table1,'Interpreter','latex','FontSize',13,'BackgroundColor','w','FitBoxToText','on','LineWidth',1,'Margin',3)

 if tosave == 'Y'
    savefig('./6-AB/avg-AB')
    saveas(gcf,'./6-AB/avg-AB','tif')
 else
 end

 %% Same as above but percent uncertainty!
 
 for i=1:size(avg_aB,1)
     avg_aB_percent(i,1)=100*((avg_aB(i,1)--87.11)/87.11);
 end
 
initialplotCH(3,'linear','linear')
plot(avg_aB_percent,'x','LineWidth',1.5)
xticks(1:size(avg_aB))
xtickangle(90)
set(gca,'xticklabel',name_aB,'XMinorTick','off')
ylabel('Uncertainty \DeltaH_{\alpha/\beta}^{avg.} (%)','FontSize',20)
xlabel('')
legend off

 
% inputs{1,1}='sample';
% inputs{2,1}='name';
% inputs{3,1}='runs';
% inputs{4,1}='thermocouple';
% inputs{5,1}='tosave';
% inputs{6,1}='correction';
% inputs{7,1}='mass';

%% legend for cell array 'data'

% Columns = Temp (C), Time (min), DSC (*/mg), Sensitivity (uV/mW)

% COLUMN 1 = File Names
% COLUMN 2 = Raw data, all temperature range, DSC (uV/mg)
% COLUMN 3 = Raw data, subset temperature, DSC (mW/mg)
% COLUMN 4 = Corrected data, subset temperature, DSC (mW/mg)
