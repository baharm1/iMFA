function cout = plot_metabolite(name,i,MID_sim,MID,SD,time_points,t)
cout = 0;

% Set the file path
global folder
filename = [folder,'\plots\',char(name),'.png'];

%% Plot the data
figure;
% calculate the configuration of the time
if(length(i) == 2)
    tiledlayout(1,2);
elseif(length(i) <= 4)
    tiledlayout(2,2);
elseif(length(i) <= 6)
    tiledlayout(2,3);
elseif(length(i) <= 9)
    tiledlayout(3,3);
end

% Create plots
for n = i'
    nexttile;
    plot(t, MID_sim(:,find(i == n)));
    hold on;
    scatter(time_points, MID(find(i == n),:));
    hold on;
    errorbar(time_points,MID(find(i == n),:),SD(find(i == n),:), 'LineStyle','none');
    title([char(name),' M+',num2str(n)]);
    xlabel('time (h)');
    ylabel('MID');
end

saveas(gcf, filename);
close
