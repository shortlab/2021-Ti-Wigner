clear

load('inputs.mat');

initial=300; % C
final=600;   % C
tosave='N';             % Save figures?

for v=2:size(inputs,2)
    
    templength(v,1)=size(inputs{8,v}(:,1),1);
    
end

[mintemplength,mintemplength_index]=min(templength(templength > 0),[],'all','linear');

U=0;
I=0;
U_diff_sub={};
I_diff_sub={};

for v=2:size(inputs,2)
    
    if strncmp(inputs{1,v},'U',1)
        U=U+1;
        U_diff_sub=cat(2,U_diff_sub,{U;cat(2,inputs{8,1+mintemplength_index}(:,1),inputs{8,1+mintemplength_index}(:,2),interp1(inputs{8,v}(:,1),inputs{8,v}(:,3),inputs{8,1+mintemplength_index}(:,1)));cat(2,inputs{9,1+mintemplength_index}(:,1),inputs{9,1+mintemplength_index}(:,2),interp1(inputs{9,v}(:,1),inputs{9,v}(:,3),inputs{9,1+mintemplength_index}(:,1)))});
    else
        I=I+1;
        I_diff_sub=cat(2,I_diff_sub,{I;cat(2,inputs{8,1+mintemplength_index}(:,1),inputs{8,1+mintemplength_index}(:,2),interp1(inputs{8,v}(:,1),inputs{8,v}(:,3),inputs{8,1+mintemplength_index}(:,1)));cat(2,inputs{9,1+mintemplength_index}(:,1),inputs{9,1+mintemplength_index}(:,2),interp1(inputs{9,v}(:,1),inputs{9,v}(:,3),inputs{9,1+mintemplength_index}(:,1)))});
    end
    
end

cat_U_temp=[];

for i=1:mintemplength
    for j=1:U
        cat_U_temp(i,j)=U_diff_sub{2,j}(i,1);
        cat_U_time(i,j)=U_diff_sub{2,j}(i,2);
        cat_U_diff_sub(i,j)=U_diff_sub{2,j}(i,3);
        cat_U_ste(i,j)=U_diff_sub{3,j}(i,3);
    end
    for j=1:I
        cat_I_temp(i,j)=I_diff_sub{2,j}(i,1);
        cat_I_time(i,j)=I_diff_sub{2,j}(i,2);
        cat_I_diff_sub(i,j)=I_diff_sub{2,j}(i,3);
        cat_I_ste(i,j)=I_diff_sub{3,j}(i,3);
    end
end

mean_U_temp=mean(cat_U_temp,2);
mean_U_diff_sub=mean(cat_U_diff_sub,2);
mean_U_time=mean(cat_U_time,2);

mean_I_temp=mean(cat_I_temp,2);
mean_I_diff_sub=mean(cat_I_diff_sub,2);
mean_I_time=mean(cat_I_time,2);

% Uncertainty (imported) from corrections x2 + averaging runs 2:5 in figure_crank 
mean_U_ste=mean(cat_U_ste,2);
mean_I_ste=mean(cat_I_ste,2);

% Uncertainty due to the averaging of multiple U/I samples
std_U=std(cat_U_diff_sub,0,2);
std_I=std(cat_I_diff_sub,0,2);
ste_U=std_U./sqrt(U);
ste_I=std_I./sqrt(I);

ste_U_tot=sqrt(mean_U_ste.^2+ste_U.^2);
ste_I_tot=sqrt(mean_I_ste.^2+ste_I.^2);

%% Plot mean of I and U

initialplotCH_portrait(1,'linear','linear')
axisposition=[0.137545944931073 0.15+(140/560) 0.766317365269462 0.557339416620955]; % left bottom width height
errorinc=25;

UnIrr=plot(mean_U_temp,1000*mean_U_diff_sub,'DisplayName','Unirradiated','LineWidth',2,'Color',0.6*ones(1,3));
errorbar(gca,mean_U_temp(1:errorinc:end),1000*mean_U_diff_sub(1:errorinc:end),1000*ste_U_tot(1:errorinc:end),'Color',0.6*ones(1,3),'LineStyle','none')

Irr=plot(mean_I_temp,1000*mean_I_diff_sub,'DisplayName','Irradiated','LineWidth',2,'Color',[0, 0, 0]);
errorbar(gca,mean_I_temp(1:errorinc:end),1000*mean_I_diff_sub(1:errorinc:end),1000*ste_I_tot(1:errorinc:end),'Color',[0, 0, 0],'LineStyle','none')

axis([initial final -1E0 7E0])
ylabel('Specific Power (\muW/mg)','FontSize',20)
xlabel('Temperature (\circC)','FontSize',20)

txt=strcat('\uparrow',{' '},'exo.');
lim=axis;
X=lim(1)+0.075*(lim(2)-lim(1));
Y=6.2E0;
Yexo=Y+0.1E0;
text(X,Yexo,txt,'FontSize',15,'Interpreter','tex','Color',0.3*ones(1,3))

if tosave == 'Y'
    savefig('mean_plot')
    saveas(gcf,'mean_plot','tif')
    savefig('./Figures/mean_plot')
    saveas(gcf,'./Figures/mean_plot','tif')
else
end

%% Plot error contributions

% initialplotCH(2,'linear','linear')
% plot(mean_I_temp,mean_I_ste,':','DisplayName','I-(corr+2:5)','LineWidth',1,'Color',[0, 0, 0])
% hold on
% plot(mean_I_temp,ste_I,'--','DisplayName','I-(X)','LineWidth',1,'Color',[0, 0, 0])
% plot(mean_I_temp,ste_I_tot,'DisplayName','I-tot','LineWidth',1,'Color',[0, 0, 0])
% plot(mean_U_temp,mean_U_ste,':','DisplayName','U-(corr+2:5)','LineWidth',1,'Color',[0.9, 0, 0])
% plot(mean_U_temp,ste_U,'--','DisplayName','U-(X)','LineWidth',1,'Color',[0.9, 0, 0])
% plot(mean_U_temp,ste_U_tot,'DisplayName','U-tot','LineWidth',1,'Color',[0.9, 0, 0])
% % axis([initial final 0E-3 1.2E-3])
% title('Error Contributions')
% 
% if tosave == 'Y'
%     savefig('mean_errors')
%     saveas(gcf,'mean_errors','tif')
%      savefig('./Figures/mean_errors')
%     saveas(gcf,'./Figures/mean_errors','tif')
% else
% end
% 
% close

%% Integrate between ROI limits

peak1initial=380;
peak1final=470;
peak2initial=500;
peak2final=590;

[rowpeak1initial] = find(mean_I_temp>peak1initial,1);
[rowpeak1final] = find(mean_I_temp>peak1final,1);
[rowpeak2initial] = find(mean_I_temp>peak2initial,1);
[rowpeak2final] = find(mean_I_temp>peak2final,1);

I_peak1=trapz(60*mean_I_time(rowpeak1initial:rowpeak1final),mean_I_diff_sub(rowpeak1initial:rowpeak1final)); % J/g
I_peak2=trapz(60*mean_I_time(rowpeak2initial:rowpeak2final),mean_I_diff_sub(rowpeak2initial:rowpeak2final)); % J/g

I_upper=mean_I_diff_sub+ste_I_tot;
I_lower=mean_I_diff_sub-ste_I_tot;
I_upper_peak1=trapz(60*mean_I_time(rowpeak1initial:rowpeak1final),I_upper(rowpeak1initial:rowpeak1final)); % J/g
I_upper_peak2=trapz(60*mean_I_time(rowpeak2initial:rowpeak2final),I_upper(rowpeak2initial:rowpeak2final)); % J/g
I_lower_peak1=trapz(60*mean_I_time(rowpeak1initial:rowpeak1final),I_lower(rowpeak1initial:rowpeak1final)); % J/g
I_lower_peak2=trapz(60*mean_I_time(rowpeak2initial:rowpeak2final),I_lower(rowpeak2initial:rowpeak2final)); % J/g

[rowpeak1initial] = find(mean_U_temp>peak1initial,1);
[rowpeak1final] = find(mean_U_temp>peak1final,1);
[rowpeak2initial] = find(mean_U_temp>peak2initial,1);
[rowpeak2final] = find(mean_U_temp>peak2final,1);

U_peak1=trapz(60*mean_U_time(rowpeak1initial:rowpeak1final),mean_U_diff_sub(rowpeak1initial:rowpeak1final)); % J/g
U_peak2=trapz(60*mean_U_time(rowpeak2initial:rowpeak2final),mean_U_diff_sub(rowpeak2initial:rowpeak2final)); % J/g

U_upper=mean_U_diff_sub+ste_U_tot;
U_lower=mean_U_diff_sub-ste_U_tot;
U_upper_peak1=trapz(60*mean_U_time(rowpeak1initial:rowpeak1final),U_upper(rowpeak1initial:rowpeak1final)); % J/g
U_upper_peak2=trapz(60*mean_U_time(rowpeak2initial:rowpeak2final),U_upper(rowpeak2initial:rowpeak2final)); % J/g
U_lower_peak1=trapz(60*mean_U_time(rowpeak1initial:rowpeak1final),U_lower(rowpeak1initial:rowpeak1final)); % J/g
U_lower_peak2=trapz(60*mean_U_time(rowpeak2initial:rowpeak2final),U_lower(rowpeak2initial:rowpeak2final)); % J/g

tstat_1=(I_peak1-U_peak1)/(I_peak1-I_lower_peak1);
pvalue_1=1-normcdf(tstat_1,0,1);

tstat_2=(I_peak2-U_peak2)/(I_peak2-I_lower_peak2);
pvalue_2=1-normcdf(tstat_2,0,1);

I1=sprintf('%.2f',I_peak1);
U1=sprintf('%.2f',abs(U_peak1));
I2=sprintf('%.2f',I_peak2);
U2=sprintf('%.2f',U_peak2);

I1_u=num2str(round(I_peak1-I_lower_peak1,2)); % J/g
U1_u=num2str(round(U_peak1-U_lower_peak1,2)); % J/g
I2_u=num2str(round(I_peak2-I_lower_peak2,2)); % J/g
U2_u=num2str(round(U_peak2-U_lower_peak2,2)); % J/g
I1_p=num2str(round(pvalue_1,3));
I2_p=num2str(round(pvalue_2,3));

% table1=strcat('\begin{tabular}{lllll}', '\textsf{( J/g )}& \textsf{1} & $\pm$ & \textsf{2} & $\pm$ \\[5pt]', ...
%                                         '\textsf{Irradiated} &\textsf{',I1,'}&\textsf{',I1_u,'}&\textsf{',I2,'}&\textsf{',I2_u,'}\\[5pt]',...
%                                         '\textsf{Unirradiated} &\textsf{',U1,'}&\textsf{',U1_u,'}&\textsf{',U2,'}&\textsf{',U2_u,'}\\[5pt]',...
%                                         '\textsf{p-value} & \textsf{',I1_p,'}&&\textsf{',I2_p,'}&','\end{tabular}');
% annotation('textbox', [axisposition(1) axisposition(2)-0.35 axisposition(3) 0.2], 'string', table1,'Interpreter','latex','FontSize',16,'BackgroundColor','w','FitBoxToText','off','LineWidth',1,'Margin',5)

annotation('rectangle',[axisposition(1) axisposition(2)-0.33 axisposition(3) 0.16],'FaceColor','w','LineWidth',1)

table0='\newline Irradiated\newline Unirradiated\newline \itp-value\newline';
annotation('textbox', [axisposition(1) axisposition(2)-0.33 axisposition(3)*(80/300) 0.16],'string', table0,'Interpreter','tex','FontSize',15,'BackgroundColor','w','FitBoxToText','off','LineWidth',1,'Margin',5,'LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle')

table1=strcat('  \int ROI 1\newline',I1,'\pm',I1_u,'\newline',U1,'\pm',U1_u,'\newline','{     }\it',I1_p);
annotation('textbox', [axisposition(1)+axisposition(3)*(80/300) axisposition(2)-0.33 axisposition(3)*(90/300) 0.16],'string', table1,'Interpreter','tex','FontSize',15,'BackgroundColor',0.95*ones(1,3),'FitBoxToText','off','LineWidth',1,'Margin',5,'LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle')

table2=strcat('  \int ROI 2\newline',I2,'\pm',I2_u,'\newline',U2,'\pm',U2_u,'\newline','{     }\it',I2_p);
annotation('textbox', [axisposition(1)+axisposition(3)*(200/300) axisposition(2)-0.33 axisposition(3)*(90/300) 0.16],'string', table2,'Interpreter','tex','FontSize',15,'BackgroundColor',0.95*ones(1,3),'FitBoxToText','off','LineWidth',1,'Margin',5,'LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle')

annotation('rectangle',[axisposition(1) axisposition(2)-0.33 axisposition(3) 0.16],'FaceColor','w','LineWidth',1,'FaceAlpha',0)

annotation('textbox',[axisposition(1) axisposition(2)-0.16 axisposition(3) 0.05],'string','Stored energy (J/g)','FontSize',20,'LineStyle','none')


coordpeak1x=[peak1initial,peak1initial,peak1final,peak1final];
coordpeak1y=[lim(3),lim(4),lim(4),lim(3)];
coordpeak2x=[peak2initial,peak2initial,peak2final,peak2final];
coordpeak2y=[lim(3),lim(4),lim(4),lim(3)];

shadepeak1=patch(coordpeak1x,coordpeak1y,[0 0 0],'DisplayName','Peak 1','FaceAlpha',0.05,'EdgeColor','none');
shadepeak2=patch(coordpeak2x,coordpeak2y,[0 0 0],'DisplayName','Peak 2','FaceAlpha',0.05,'EdgeColor','none');
legend([Irr UnIrr])


arrow1x=[axisposition(1)+((peak1initial-initial)/(final-initial))*axisposition(3), axisposition(1)+((peak1final-initial)/(final-initial))*axisposition(3)];
arrow1y=[axisposition(2)+((Y-lim(3))/(lim(4)-lim(3)))*axisposition(4), axisposition(2)+((Y-lim(3))/(lim(4)-lim(3)))*axisposition(4)];

arrow1=annotation('arrow',arrow1x, arrow1y,'HeadStyle','vback3','Color',0.3*ones(1,3));
arrow1b=annotation('arrow',[arrow1x(2),arrow1x(1)], arrow1y,'HeadStyle','vback3','Color',0.3*ones(1,3));

arrow2x=[axisposition(1)+((peak2initial-initial)/(final-initial))*axisposition(3), axisposition(1)+((peak2final-initial)/(final-initial))*axisposition(3)];
arrow2y=[axisposition(2)+((Y-lim(3))/(lim(4)-lim(3)))*axisposition(4), axisposition(2)+((Y-lim(3))/(lim(4)-lim(3)))*axisposition(4)];

arrow2=annotation('arrow',arrow2x, arrow2y,'HeadStyle','vback3','Color',0.3*ones(1,3));
arrow2b=annotation('arrow',[arrow2x(2),arrow2x(1)], arrow2y,'HeadStyle','vback3','Color',0.3*ones(1,3));

peak1numX=peak1initial+0.35*(peak1final-peak1initial);
peak1numY=Y+0.25E0;
peak1num=text(peak1numX,peak1numY,'ROI 1','FontSize',15,'Interpreter','tex','Color',0.3*ones(1,3));

peak2numX=peak2initial+0.35*(peak2final-peak2initial);
peak2numY=Y+0.25E0;
peak2num=text(peak2numX,peak2numY,'ROI 2','FontSize',15,'Interpreter','tex','Color',0.3*ones(1,3));

if tosave == 'Y'
    savefig('mean_plot')
    saveas(gcf,'mean_plot','tif')
    savefig('./Figures/mean_plot')
    saveas(gcf,'./Figures/mean_plot','tif')
else
end

set(legend,...
    'Position',[0.201042309604859 0.749987154878241 0.256756756756757 0.0888979370249728],...
    'LineWidth',1,...
    'FontSize',15,...
    'EdgeColor',[0.7 0.7 0.7]);


%% Integrate over full temperature range

[rowpeak1initial] = find(mean_I_temp>peak1initial,1);
[rowpeak2final] = find(mean_I_temp>peak2final,1);

I_all=trapz(60*mean_I_time(rowpeak1initial:rowpeak2final),mean_I_diff_sub(rowpeak1initial:rowpeak2final)); % J/g
I_upper_all=trapz(60*mean_I_time(rowpeak1initial:rowpeak2final),I_upper(rowpeak1initial:rowpeak2final)); % J/g
I_lower_all=trapz(60*mean_I_time(rowpeak1initial:rowpeak2final),I_lower(rowpeak1initial:rowpeak2final)); % J/g

[rowpeak1initial] = find(mean_U_temp>peak1initial,1);
[rowpeak2final] = find(mean_U_temp>peak2final,1);

U_all=trapz(60*mean_U_time(rowpeak1initial:rowpeak2final),mean_U_diff_sub(rowpeak1initial:rowpeak2final)); % J/g
U_upper_all=trapz(60*mean_U_time(rowpeak1initial:rowpeak2final),U_upper(rowpeak1initial:rowpeak2final)); % J/g
U_lower_all=trapz(60*mean_U_time(rowpeak1initial:rowpeak2final),U_lower(rowpeak1initial:rowpeak2final)); % J/g


%% Add integral standard errors (sample + correction) in quadrature - for two ROIs

load('mean+std_err-samples-ROIs.mat')
load('mean+std_err-corrections-ROIs.mat')

std_err_peak1_I=sqrt((std_err_meanpeak1_I)^2+(I_peak1_uncertainty)^2);
std_err_peak2_I=sqrt((std_err_meanpeak2_I)^2+(I_peak2_uncertainty)^2);
std_err_peak1_U=sqrt((std_err_meanpeak1_U)^2+(U_peak1_uncertainty)^2);
std_err_peak2_U=sqrt((std_err_meanpeak2_U)^2+(U_peak2_uncertainty)^2);


%% Format these standard errors and plot

tstat_1=(I_peak1-U_peak1)/(std_err_peak1_I);
pvalue_1=1-normcdf(tstat_1,0,1);
I1_p=num2str(round(pvalue_1,4));

tstat_2=(I_peak2-U_peak2)/(std_err_peak2_I);
pvalue_2=1-normcdf(tstat_2,0,1);
I2_p=num2str(round(pvalue_2,4));

I1_u=num2str(round(std_err_peak1_I,2)); % J/g
U1_u=num2str(round(std_err_peak1_U,2)); % J/g
I2_u=num2str(round(std_err_peak2_I,2)); % J/g
U2_u=num2str(round(std_err_peak2_U,2)); % J/g

table0='\newline Irradiated\newline Unirradiated\newline \itp-value\newline';
annotation('textbox', [axisposition(1) axisposition(2)-0.33 axisposition(3)*(80/300) 0.16],'string', table0,'Interpreter','tex','FontSize',15,'BackgroundColor','w','FitBoxToText','off','LineWidth',1,'Margin',5,'LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle')

table1=strcat('  \int ROI 1\newline',I1,'\pm',I1_u,'\newline',U1,'\pm',U1_u,'\newline','{   }\it',I1_p);
annotation('textbox', [axisposition(1)+axisposition(3)*(80/300) axisposition(2)-0.33 axisposition(3)*(90/300) 0.16],'string', table1,'Interpreter','tex','FontSize',15,'BackgroundColor',0.95*ones(1,3),'FitBoxToText','off','LineWidth',1,'Margin',5,'LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle')

table2=strcat('  \int ROI 2\newline',I2,'\pm',I2_u,'\newline',U2,'\pm',U2_u,'\newline','{   }\it',I2_p);
annotation('textbox', [axisposition(1)+axisposition(3)*(200/300) axisposition(2)-0.33 axisposition(3)*(90/300) 0.16],'string', table2,'Interpreter','tex','FontSize',15,'BackgroundColor',0.95*ones(1,3),'FitBoxToText','off','LineWidth',1,'Margin',5,'LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle')

annotation('rectangle',[axisposition(1) axisposition(2)-0.33 axisposition(3) 0.16],'FaceColor','w','LineWidth',1,'FaceAlpha',0)

annotation('textbox',[axisposition(1) axisposition(2)-0.16 axisposition(3) 0.05],'string','Stored energy (J/g)','FontSize',20,'LineStyle','none')

% close all


%% Add integral standard errors (sample + correction) in quadrature - for all temperatures

load('mean+std_err-samples-alltemps.mat')
load('mean+std_err-corrections-alltemps.mat')

std_err_peak_I=sqrt((std_err_meanpeak_I)^2+(I_peak_uncertainty)^2);
std_err_peak_U=sqrt((std_err_meanpeak_U)^2+(U_peak_uncertainty)^2);
