%% Vintage Peak monte carlo simulation
% Adam Hawkins, last major edit: Aug 2022
% This code block relies on Matlab code provided within Balco et al. (2008),
% Borchers et al. (2016), and Lifton et al.(2014)

close all; 
% profile on
% profile viewer
clear all;
% site = [1 2 3 4]

snow_z = (220./2); %Snow depth in cm; snow depth avg depth from nearby pillow.

% z = [0.9.*ice_z 0.9*ice_z(length(ice_z))+((0:1:10).*2.65)];%ice cover
z = [0.38.*snow_z 0.38*snow_z(length(snow_z))+((1:1:100).*2.65)]; %mass depth in cm
for site = 4
    if site == 4
        out = depth(z,50.26192000,-124.30493000,1479); %VPC4
        obsz = [0.0250    0.1250    0.2250    0.3250    0.4250].*2.65.*100;
        obs = [188600 179400  151100  134700  146400];
        err = [203100   187000  158200  141800  153700]-[158600   164300  136400  120100  131200];
        err_2sig = [209800  191200  162500  145600  157700]-[115000 142500  115000  99200   109500];
        n14_min = 7754; %5890; with snow
        N10 = 1.669E+05;
        N10_err = 5.925E+03;
    elseif  site == 3
        out = depth(z,50.25635000, -124.30430000,1606);
        obsz = [0.0250    0.0750    0.1250    0.2250    0.3250].*2.65.*100;%VPC3  
        obs = [81000 74100 70900 59900 85300];
        err = [93200 81200 78000 66900 92300]-[55100 59000 55700 44700 70400];
        err_2sig = [98500   84400   81300   70200   95800]-[18300   37000   33500   22900   48800];
        n14_min = 2325; %1870; %2400 with snow
        N10 = 6.104E+04;
        N10_err = 3.090E+03;
    elseif site == 2
        out = depth(z,50.25346000, -124.30145000,1726);%VPC2
        obsz = [0.0250    0.0750    0.1250    0.1750    0.2750].*2.65.*100;
        obs = [84000 79000 80200 78000 63700];
        err = [95100 86100 87300 85000 70800]-[59900 64100 65400 63000 48800];
        err_2sig = [100200  89400   90500   88200   73900]-[25400   42600   43500   41300   27100];
        n14_min = 2199; %1770; %2270 with snow
        N10 = 5.694E+04;
        N10_err = 2.241E+03;
    elseif site == 1
        out = depth(z,50.25861000, -124.28943000,1726); %VPC1
        obsz = [0.0250    0.0750    0.1750    0.2750    0.3750].*2.65.*100;%in meters
        obs = [140100 130700 110600 120900 122900];
        err = [151700 137800 117700 128200 134800]-[115800 116000 95700 106000 98100];
        err_2sig = [157100  141700  121300  132000  140300]-[80100  94200   73600   84300   62000];
        n14_min = 4078; %3210; % 4230 with snow
        N10 = 1.179E+05;
        N10_err = 3.351E+03;
    end
    
        l14 = log(2)/5730;
        l10 = log(2)/1.38e6;
        P14 = out.P14_tot;
        P10 = out.P10_tot;

        num_iter=10000;

        e = [0 0.005 0.01 0.02]; %subglacial erosion (cm/yr)
        e_air = 0.00; %cm/yr subaerial erosion
        tinh = 3000; %inheritence from previous exposure

        for i = 1:num_iter
            te = unifrnd(n14_min, 15000);  
            tb = unifrnd(0, 15000-n14_min); %Set te/tb for the iteration
            for j = 1:length(obsz)
                obs_z = 0.38*max(snow_z) + obsz(j); % mass depth of the observed 14C measurement under snowcover
                idx = interp1(z, 1:numel(z), obs_z, 'nearest'); %returns index value of z nearest to the depth of the observation. 
                idx_10 = interp1(z,1:numel(z), (0.38*max(snow_z) + (2.5*2.65)), 'nearest'); 

                P14_sam = P14(idx); %Find P14 for each depth(z) in obs_z
                P10_sam = P10(idx_10);
                for k = 1:length(e) %with same te/tb, run iteration three times with dif e values

                    out.N14(i,j,k) = P14_sam./l14.*(1-exp(-te.*l14)).*exp(-l14.*tb).*exp(-2.65.*(e(k).*tb./160)).*exp(-2.65.*(e_air.*te./160));
%                     out.N10(i,1,k) = P10_sam./l10.*(1-exp(-te.*l10)).*exp(-l10.*tb).*exp(-2.65.*(e(k).*tb./160)).*exp(-2.65.*(e_air.*te./160));
                    out.N10(i,1,k) = (P10_sam./l10.*(1-exp(-tinh.*l10)))+ P10_sam./l10.*(1-exp(-te.*l10)).*exp(-l10.*tb).*exp(-2.65.*(e(k).*tb./160)).*exp(-2.65.*(e_air.*te./160));
%                     out.N14(i,j,k) = P14_sam./l14.*(1-exp(-te.*l14)).*exp(-l14.*tb).*exp(-2.65.*(e(k).*tb./160));
%                     out.N10(i,1,k) = P10_sam./l10.*(1-exp(-te.*l10)).*exp(-l10.*tb).*exp(-2.65.*(e(k).*tb./160));
                    out.te(i,j,k) = te;
                    out.tb(i,j,k) = tb;
                end
            end
        end
     
        % Filter results so that only modeled profiles within 2 sigma range of observations along entire core are captured
        obs_mat = repmat(obs,length(out.N14(:,1,:)),1,length(e));
        err_mat = repmat(err_2sig,length(out.N14(:,1,:)),1,length(e));
        out.N14(out.N14<=(obs_mat-err_mat) | out.N14 >=(obs_mat+err_mat)) = NaN;
        out.N10(out.N10<=(N10-(2*N10_err)) | out.N10 >=(N10+(2*N10_err))) = NaN; 

        out.Bechi2 = ((N10 - out.N10).^2./out.N10);
        out.Bechisum = sum(out.Bechi2,2);
        out.chi2 = ((obs - out.N14).^2./out.N14); 
        out.chisum = sum(out.chi2,2);

        tb = out.tb(:,1,:);
        te = out.te(:,1,:);
        chisum = out.chisum(:,1,:);
        bechisum = out.Bechisum(:,1,:);
        chi_total = chisum + bechisum;
        out.results = cat(2, te, tb, chisum, bechisum, chi_total);
        out.surface = cat(2, te, tb, out.N14(:,1,:), out.N10(:,1,:));

        obs_z = 0.3*max(snow_z) + obsz;
end

%% Plotting the observed vs modelled nuclide profiles
% fig=figure(1);
% t = tiledlayout(2,length(e)/2);
% title(t, strcat('Core', {' '}, num2str(site)), 'FontWeight', 'bold');
% xlabel(t, 'Concentration (at g^{-1})');
% ylabel(t, 'Mass Depth (g cm^{-2})');
% for i = 1:length(e)
%     nexttile
%     N14 = out.N14(:,:,i);
%     N14 = N14((all((~isnan(N14)),2)),:); %filtering out results with measurements outside of 2 sig bounds along any length of core. 
%     N10mod = out.N10(:,:,i);
%     N10mod = N10mod((all((~isnan(N10mod)),2)),:);
% 
%     boundedline(obs,obs_z, err_2sig, 'orientation', 'horiz', 'alpha');
% %     xlim([0 max(obs)]);
%     ylim([(min(obs_z)-5) max(obs_z)]);
%     title(e(i) + " cm yr^{-1} erosion");
%     hold on;
%     errorbar(N10, (0.3*max(snow_z)+(0.0250*2.65*100)), N10_err,'horiz', 'Marker', 'square','Color','r');
%     hold on;
% 
%     for k = 1:length(N14(:,1,1))
%         h = semilogx(N14(k,:),obs_z,'LineWidth',1.5, 'Color','k');
%     end
%     hold on;
%     for l = 1:length(N10mod(:,1,1))
%         m = plot(N10mod(l,:),(0.3*max(snow_z)+(0.0250*2.65*100)),'Marker','o', 'Color','k');
%     end
%     hold on;
%     grid off;
%     ax = gca;
%     ax.YDir = 'reverse';
% %     ax.FontAngle = 'italic'; ax.FontSize = 12; ax.FontName = 'Arial'; ax.FontWeight = 'bold';
%     txt = ['n sims = ' sprintf('%d',size(N14,1))];
%     xlim(xlim+[-10,10]);ylim(ylim+[-10,10]);%use xlim and ylim to fix location
%     text(max(xlim), max(ylim)*0.90,txt, 'Horiz','right', 'Vert','top')
% end
% filename = ['Monte_Core' num2str(site) '_' num2str(tinh) 'yrInheritance' '_SmallE.jpg'];
% fpath = '/Users/adamhawkins/Documents/Vintage_Peak/Figures/';
% saveas(t, fullfile(fpath, filename));
% 
% close
%% 

fig = figure('Position',[10 10 1000 200]);
% t = tiledlayout(2,length(e)/2);
t = tiledlayout('horizontal', 'TileSpacing','Compact')
% title(t, strcat('Core', {' '}, num2str(site)), 'FontWeight', 'bold');

for i= 1:length(e)
    nexttile
    data = out.results(:,:,i);
    data = data((all((~isnan(data)),2)),:);
    x = data(:,1);
    y = data(:,2);
    z = data(:,5); % column 5 is summed 14C and 10Be chi
    [X,Y] = meshgrid(min(x):max(x),min(y):max(y));
    Z=griddata(x,y,z,X,Y);
    contourf(X,Y,Z,20); %Make contour plot
    % pointsize=10;
    % scatter(x,y,pointsize,z, 'filled');
    xlim([n14_min 15000]);
    ylim([0 15000-n14_min]);
    %colorbar;
    xlabel('exposure duration (kyr)');
    if i == 1
        label = ['Core ',num2str(site)];
        ylabel({['\bf' label];['\rm' 'burial duration (kyr)']});
    else
    end
    if i == 4
        contourcbar
    else
    end
    % zlabel('chi-squared');
    txt = ['n sims = ' sprintf('%d',size(x,1))];
    xlim(xlim+[-10,10]);ylim(ylim+[-10,10]);%use xlim and ylim to fix location
    text(min(xlim)*1.05, max(ylim)*0.70,txt, 'Horiz','left', 'Vert','bottom')
%     text((15000-n14_min)/2, 5000, txt, 'LineWidth',1);
    if site == 1
        title(e(i) + " cm yr^{-1} erosion");
    else
    end
    hold on;
    minValue = min(data(:,5));
    [row, column] = find(data == minValue);
    min_chi = data(row,:);
    if not(isempty(data))
        plot(min_chi(1,1), min_chi(1,2),'+', 'Color', 'red', 'MarkerSize', 15 , 'LineWidth',4);
        exp_txt = ['Exposure duration: ' num2str(round(min_chi(1,1))) ' kyr'];
        bur_txt = ['Burial duration: ' num2str(round(min_chi(1,2))) ' kyr'];
        % text(min(xlim)*1.05, max(ylim)*0.85,exp_txt, 'Color', 'white', 'FontWeight', 'bold', 'Horiz','left', 'Vert','bottom');
        text(min(xlim)*1.05, max(ylim)*0.90,exp_txt, 'Horiz','left', 'Vert','bottom');
        text(min(xlim)*1.05, max(ylim)*0.80,bur_txt, 'Horiz','left', 'Vert','bottom');
    else
    end
end
% ylabel('burial duration (kyr)');
filename = ['Monte_Core' num2str(site) '_' num2str(tinh) 'yrInheritance' '_SmallEContourPlot.jpg'];
fpath = '/Users/adamhawkins/Documents/Vintage/Matlab_output/';
saveas(t, fullfile(fpath, filename));

%Let's look at how much having our depth profiles help to narrow acceptable
%exposure histories:
sum(sum(~isnan(out.surface),2)==4)
sum(sum(~isnan(out.results),2)==5)