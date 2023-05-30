% Compute curvature from a series of (x,y) positions.
%
%
% The computation uses orthogonal regression due to uncertainties in both
% the x and y directions (but note that it's assumed the variance of the
% uncertainties in x and y are equal).
%
% Written by Jacob Notbohm, University of Wisconsin-Madison, 2013-2019
% 

clear;
close all;
clc;

% --- USER INPUTS ---
domainname = 'domain.tif'; % Name of multipage tif file with domains
downsample_factor = 1;
makeplots = 1; % Set to 1 to make plots; set to zero to skip plotting
N = 51; % Number of points in window
% N = 51; % default: Number of points in window
% -------------------
% Curvature limits for plotting
% minkappa = -0.03;%default
% maxkappa = 0.03;
minkappa = -0.04;
maxkappa = 0.04;

% Data info
info = imfinfo(domainname);
num_images = length(info);

% Preallocate cell of curvature and angle values
xyktheta = cell(num_images,1);

for k=1:num_images
        domain = imread(domainname,k);
    domain = logical( domain/max(domain(:)) ); % Convert to logical array
    BNDRY = bwboundaries(domain); % This should be 1 cell with boundary coordinates
    
    if ~isempty(BNDRY)
        BNDRY = BNDRY{1};
        x_bndry = BNDRY(:,2); % x coordinates are the columns
        y_bndry = BNDRY(:,1); % y coordinates are the rows
        
        % Downsampling can reduce noise
%         x = downsample(x_bndry,downsample_factor);
%         y = downsample(y_bndry,downsample_factor);
%         
        % Need to specify location of monolayer edge that you want to
        % compute

        seq=1100:2000; 
        % data you are interested in analyzing, varies by
%         image
        x = downsample(x_bndry(seq),downsample_factor);
        y = downsample(y_bndry(seq),downsample_factor);
        
        % data you are interested in analyzing, varies by
%         image
%         x = downsample(x_bndry(:),downsample_factor);
%         y = downsample(y_bndry(:),downsample_factor);
%

        if length(x) > N
            
            % % Smooth data to reduce noise
%             span=11;
%             x=smooth(x,span); y=smooth(y,span);
            
            % --- Compute the angle of a line of best fit to a window of N
            % points using the orthogonal regression ---
            
            
            % Wrap data so that curvatures can be computed for endpoints of vectors x
            % and y
            xwrap = [ x(length(x)-(N-1)/2:length(x)) ; x ; x(1:(N+1)/2) ];
            ywrap = [ y(length(y)-(N-1)/2:length(y)) ; y ; y(1:(N+1)/2) ];
            
            % Preallocate variables
            s_tot = zeros(length(x),1)*nan; % Path length
            kappa = s_tot; % Curvature
            theta0_save = s_tot; % Angles
            
            % For each window
            for win=1:length(x)
                % -1th window indices
                idx = win + (0:N-1);
                % x and y data
                winx = xwrap(idx);
                winy = ywrap(idx);
                % Centroids
                xc_minus = mean(winx);
                yc_minus = mean(winy);
                % Complex notation
                winz = winx + 1i*winy;
                zc = xc_minus + 1i*yc_minus;
                % Compute angle
                W = sum( (winz-zc).^2 );
                % Use atan2 instead of angle to achieve 4-quadrant arctangent
                theta_minus = atan2(imag(W),real(W))/2;
                
                % 0th window indices
                idx = win+1 + (0:N-1);
                % x and y data
                winx = xwrap(idx);
                winy = ywrap(idx);
                % Centroids
                xc_0 = mean(winx);
                yc_0 = mean(winy);
                % Complex notation
                winz = winx + 1i*winy;
                zc = xc_0 + 1i*yc_0;
                % Compute angle
                W = sum( (winz-zc).^2 );
                theta_0 = atan2(imag(W),real(W))/2;
                
                % The angle is only in (-pi,pi] due to the division by 2.
                % Can try to account for this by determining whether
                % displacement in x is to right or left. If it's to left,
                % then angle theta_0 should be theta_0+pi.
                dx = winx(length(winx)) - winx(1);
                if dx<0
                    theta_0_fullcircle = theta_0 + pi;
                else
                    theta_0_fullcircle = theta_0;
                end
                
                
                % +1st window indices
                idx = win+2 + (0:N-1);
                % x and y data
                winx = xwrap(idx);
                winy = ywrap(idx);
                % Centroids
                xc_plus = mean(winx);
                yc_plus = mean(winy);
                % Complex notation
                winz = winx + 1i*winy;
                zc = xc_plus + 1i*yc_plus;
                % Compute angle
                W = sum( (winz-zc).^2 );
                theta_plus = atan2(imag(W),real(W))/2;
                
                % There is a branch cut between quadrants 2 and 3 (at -pi -> pi). After
                % dividing by 2, this means the branch cut is between -pi/2 -> pi/2.
                % The half quadrants potentially affected by the branch cut are given
                % by -pi/2 -> -pi/4 and pi/4 -> pi/2. If angles are in these half
                % quadrants, add pi to the angles in quadrant 4 to move the branch cut.
                thetas = [theta_minus theta_0 theta_plus]; % Vector of angles
                clear idx_q3 idx_q4;
                idx_q2 = max( thetas>pi/4 & thetas<=pi/2 );
                idx_q3 = max( thetas>=-pi/2 & thetas<=-pi/4 );
                if idx_q2==1 && idx_q3==1
                    % Find angles in 4th quadrant
                    clear idx;
                    idx = thetas>-pi/2 & thetas<0;
                    % Add pi to move branch cut
                    thetas(idx) = thetas(idx)+pi;
                end
                
                % Compute arclengths
                s_minus = -sqrt( (xc_minus-xc_0)^2 + (yc_minus-yc_0)^2 );
                s_0 = 0;
                s_plus = sqrt( (xc_plus-xc_0)^2 + (yc_plus-yc_0)^2 );
                s = [s_minus s_0 s_plus]; % Vector of arclengths
                
                % Compute d theta / d s with linear regression to get curvature
                coeffs = polyfit(s,thetas,1);
                % Save total arclength
                if win==1
                    s_tot(win,1) = 0;
                else
                    s_tot(win,1) = s_tot(win-1,1) + abs(s_minus);
                end
                % Store curvature
                kappa(win,1) = coeffs(1);
                
                % Store angles
                
                theta0_save(win,:) = theta_0_fullcircle; % This is computed with atan2 (4 quadrant angle) and before correcting for branch cut
                
            end
            
            % Output mean and standard dev of curvatures
            disp(['Mean curvature is ',num2str(mean(kappa))]);
            disp(['Stand dev curvature is ',num2str(std(kappa))]);
             disp(['Magnitude (norm) curvature is ',num2str(norm(kappa))]);
             disp(['RMS curvature is ',num2str(rms(kappa))]);
            % --- PLOTS ---
            if makeplots == 1
                
                % % % --- Angle plot (For debugging only) ---
                
% %                 hf0 = figure;
% %                 pos = get(hf0,'position');
% %                 set(hf0,'color','w','position',[0.25*pos(1) 0.10*pos(2) 1.5*pos(3) 1.5*pos(4)]);
% %                 hold on
% %                 xplot = 1:length(theta0_save);
% %                 % plot(xplot,thetas_save(:,1)*180/pi,'ro');
% %                 plot(xplot,theta0_save(:,2)*180/pi,'ko');
% %                 % plot(xplot,thetas_save(:,3)*180/pi,'bo');
% %                 xlabel('Data point number','fontsize',11);
% %                 ylabel('Angle (deg)','fontsize',11);
                
                
                % --- Plot curvature vs path length ---
                
                hf1 = figure;
                pos = get(hf1,'position');
                set(hf1,'color','w','position',[0.2*pos(1) 0.1*pos(2) 1.5*pos(3) 0.8*pos(4)]);
                
                subplot(1,2,1)
                plot(s_tot,kappa,'-')
                xlabel('Path distance, s','fontsize',11);
                ylabel('Curvature, \kappa','fontsize',11);
                set(gca,'fontsize',11);
                axis square
                box off
                
                subplot(1,2,2)
                
                fs = 1/10; %(1.5*downsample_factor); % Approximate sampling frequency (1/pix)
                [kappa_psd, freq_vec] = pwelch(kappa);
                % plot(freq_vec*fs,kappa_psd,'b-')
                xx = (1:length(kappa_psd))/(length(kappa_psd))*fs;
                plot(xx,kappa_psd,'b-');
                
                xlabel('Wavenumber, (1/pix)','fontsize',11);
                ylabel('PSD','fontsize',11);
                set(gca,'fontsize',11)
                xlim([0 .1])
                axis square
                box off
                
                % Save plot
                set(hf1,'PaperPositionMode','auto');
                print(hf1,'-dpng','-r450','curvatures_s');
                
                
                % --- Color plot of curvatures at (x,y) positions ---
                
                hf2 = figure;
                set(hf2,'color','w','position',[1.1*pos(1) 0.1*pos(2) 1.5*pos(3) 1.5*pos(4)]);
                
                % Set colormap
                colors = [0 0 189 ; ...
                    0 0 255 ; ...
                    0 66 255 ; ...
                    0 132 255 ; ...
                    0 189 255 ; ...
                    0 255 255 ; ...
                    66 255 189 ; ...
                    132 255 132 ; ...
                    189 255 66 ; ...
                    255 255 0 ; ...
                    255 189 0 ; ...
                    255 132 0 ; ...
                    255 66 0 ; ...
                    255 0 0 ; ...
                    189 0 0 ; ...
                    132 0 0 ]/255 ;
                % Set spacing
                spacing=(maxkappa-minkappa)/length(colors);
                
                % Plot with color specified
                for kk=1:length(colors)
                    idx = (kappa>=(kk-1)*spacing+minkappa) & (kappa<=kk*spacing+minkappa) ;
                    vec_color=colors(kk,:);
                    plot(x(idx),y(idx),'o','Color',vec_color,'MarkerFaceColor',vec_color,'MarkerSize',2.25);
                    hold on
                end
                % Plot data above and below color limits with an x
                idx = kappa<minkappa;
                plot(x(idx),y(idx),'x','Color',[0 0 189]/255,'MarkerSize',6);
                idx = kappa>maxkappa;
                plot(x(idx),y(idx),'x','Color',[132 0 0]/255,'MarkerSize',6);
                
                % Plot cell island outline
                plot(x,y,'w-');
                
                xlabel('X Position','fontsize',11);
                ylabel('Y Position','fontsize',11);
                set(gca,'fontsize',11);
                axis xy
                axis equal
                
                set(gca,'color','k');
                posa = get(gca,'position');
                set(gca,'position',[posa(1) posa(2) 0.9*posa(3) 0.9*posa(4)]);
                
                % Plot white box at position s=0 (Note positive direction is CCW)
                hold on
                plot(x(1),y(1),'sw');
                
                % % % For checking path length - plot numbers at various positions along path
                % % threshold_vals = [0 500 1000 1500 2000 2500];
                % % for k=1:length(threshold_vals);
                % %     idxs(k) = find(s_tot>threshold_vals(k),1,'first');
                % %     text(x(idxs(k)),y(idxs(k)),num2str(threshold_vals(k)/1000),'color',[0.8 0.8 0.8]);
                % % end
                
                % Colorbar
%                  caxis([-0.01 0.01]);
                caxis([minkappa maxkappa]);
                hc = colorbar;
                %colormap(redblue);
                posc = get(hc,'Position');
                % set(hc,'Position',[posc(1) posc(2) 0.5*posc(3) posc(4)]);
                set(hc,'Position',[1.1*posa(3) posa(2) 0.6*posc(3) posc(4)]);
                % Create textbox
                annotation(hf2,'textbox',[0.97*posa(3) 0.7721 0.3 0.039],...
                    'String',{'Curvature, \kappa'},'FitBoxToText','off','LineStyle','none','fontsize',11);
                
                % Save plot
                set(hf2,'PaperPositionMode','auto');
                set(hf2,'InvertHardCopy','off');
                 print(hf2,'-dpng','-r450','curvatures_xy');
                
                % --- Histogram ---
                
                hf3 = figure;
                pos = get(hf3,'position');
                set(hf3,'color','w','position',[0.2*pos(1) 1.2*pos(2) 0.8*pos(3) 0.8*pos(4)]);
                
                bincenters = linspace(minkappa,maxkappa,51);
                [num, ~] = hist(kappa,bincenters);
                
                bar(bincenters,num/length(kappa)*100);
                xlabel('Curvature, \kappa','fontsize',11);
                ylabel('Frequency (%)','fontsize',11);
                set(gca,'fontsize',11);
                box off
                xlim([min(bincenters) max(bincenters)])
                % Put mean and standard deviation on plot
                annotation(hf3,'textbox',[0.1483 0.7767 0.5 0.15],'fontsize',11,'linestyle','none',...
                    'String',{['\kappa_{mean} = ',num2str(mean(kappa))],['\kappa_{stdv} = ',num2str(std(kappa))],['\kappa_{rms} = ',num2str(rms(kappa))]});
                % Save plot
                set(hf3,'PaperPositionMode','auto');
                set(hf3,'InvertHardCopy','off');
                 print(hf3,'-dpng','-r450','histogram_xy');
              
               % close all;
            end
            
            % Save x, y, kappa, and theta to k-th cell
            xyktheta{k} = [x, y, kappa, theta0_save];
            
        end
    end
end

% Save curvatures
save('curvatures_angles.mat','xyktheta');


