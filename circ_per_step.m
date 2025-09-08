clc; clear; close all;

timeSteps = [53.3, 53.9, 54.4, 54.9, 55.4, 55.9, 56.4]; % time steps for the half cycle
xq = -1:0.003:2;   % grid resolution for structured grid in x direction (0.003m)
yq = -2:0.003:2;   % grid resolution for structured grid in y direction (0.003m)
[Xq, Yq] = meshgrid(xq, yq); % mesh grid function to create 2D grid of x y coords

neg_vorticity_threshold = -10;  % vorticity threshold to filter out weaker vortices
pad = 0.20; % space around airfoil for plotting 
timeSteps = timeSteps(1:end-1); % ignore last time step since it is periodic
numSteps = numel(timeSteps); % time steps left
centroidX_history = NaN(numSteps,1); % store centroid location in x
centroidY_history = NaN(numSteps,1); % store centroid location in y
circulation_history = NaN(numSteps,1); % circulation storage

for t = 1:numSteps
    vortFilename    = sprintf('1.5m/data-%.1f.csv', timeSteps(t)); % open vorticity field data files
    airfoilFilename = sprintf('1.5m/airfoil-position/airfoil-%.1f.csv', timeSteps(t)); % open airfoil position files

    if ~isfile(vortFilename) || ~isfile(airfoilFilename) % skip if file missing
        continue;
    end

    raw = readmatrix(vortFilename); %% raw vorticity data
    if any(isnan(raw(1,:))) || any(isinf(raw(1,:))) % for first rows (some corrupted so this allows code to skip without throwing error)
        raw = raw(2:end,:);
    end
    x_raw = raw(:,1); y_raw = raw(:,2); vort_z = raw(:,12); % x,y and vorticity z values
    ok = ~isnan(x_raw) & ~isnan(y_raw) & ~isnan(vort_z); % removing bad points
    x_raw = x_raw(ok); y_raw = y_raw(ok); vort_z = vort_z(ok); % corrected data files with NaN removed

    Vq = griddata(x_raw, y_raw, vort_z, Xq, Yq, 'linear'); % interpolate vorticity field onto the structured grid

    A = readmatrix(airfoilFilename); % aerofoil coords
    if any(isnan(A(1,:))) || any(isinf(A(1,:))) % keeping valid coords (same as vorticity fields above)
        A = readmatrix(airfoilFilename,'Range','A2');
    end
    xa = A(:,1); ya = A(:,2); % x and y coords 
    valid = ~isnan(xa) & ~isnan(ya); % same as above
    xa = xa(valid); ya = ya(valid);
    XYuniq = unique([xa, ya], 'rows', 'stable'); % ensures only unique points stay to remove duplicate points
    xa = XYuniq(:,1); ya = XYuniq(:,2); % aerofoil with duplicates and NaN's removed

    % CHATGPT FIX FOR AEROFOIL PLOTTING
    if numel(xa) < 4
        k = convhull(xa, ya);
        xb = xa(k); yb = ya(k);
    else
        k = boundary(xa, ya, 0.2); % shrink factor, tweak if needed
        xb = xa(k); yb = ya(k);
    end
    if ~(xb(1)==xb(end) && yb(1)==yb(end))
        xb(end+1)=xb(1); yb(end+1)=yb(1);
    end

    % identifies bounds of aerofoil and replaces all points inside with NaN
    airfoil_mask = inpolygon(Xq, Yq, xb, yb);
    Vq(airfoil_mask) = NaN;
    Vq(Vq > neg_vorticity_threshold) = NaN;
% identifying vortex regiosn
    vortexMask = ~isnan(Vq); % valid vortex points
    if ~any(vortexMask(:)), continue; end % skips if no vortex can be detected
    CC = bwconncomp(vortexMask); % finds connected regions of vortex
    [~, largestID] = max(cellfun(@numel, CC.PixelIdxList)); % selects largest vortex
    largestMask = false(size(Vq)); 
    largestMask(CC.PixelIdxList{largestID}) = true;
    Vq_lev = NaN(size(Vq)); %keep ONLY largest vortex
    Vq_lev(largestMask) = Vq(largestMask);

    % centroid and circulation calculation
    dx = xq(2)-xq(1); dy = yq(2)-yq(1); % grid spacing for area 
    valid = ~isnan(Vq_lev);
    xV = Xq(valid); yV = Yq(valid); % valid x and y vortex points
    w  = abs(Vq_lev(valid)); % weight points by vorticity magnitude
    centroidX = sum(xV.*w)/sum(w); % weighted x pos
    centroidY = sum(yV.*w)/sum(w); % weighted y pos
    circulationVal = sum(Vq_lev(valid))*dx*dy; % circulation of vortex

    % saving results
    centroidX_history(t) = centroidX;
    centroidY_history(t) = centroidY;
    circulation_history(t) = circulationVal;

    % plotting
    figure;

    contourf(Xq, Yq, Vq_lev, 40, 'LineColor','none'); % fill contour plot of vorticity
    colormap(turbo); 
    colorbar;
    
    % adjust colour range for clarity
    vmax = max(abs(Vq_lev(:)), [], 'omitnan');
    if isempty(vmax) || vmax==0
        vmax = 1;
    end
    caxis([-45 -10]); 
    
    hold on; 
    plot(xb, yb, 'k-', 'LineWidth', 2);       % airfoil outline
    plot(centroidX, centroidY, 'ro', 'MarkerSize', 8, 'LineWidth', 2); % centroid marker
    hold off;
    
    axis equal;
    xlabel('x (m)'); ylabel('y (m)');
    xlim([min(xb)-pad, max(xb)+pad]);
    ylim([min(yb)-pad, max(yb)+pad]);
    title(sprintf('t = %.1f', timeSteps(t)));  
end


figure; plot(timeSteps, centroidX_history, '-or'); hold on;
plot(timeSteps, centroidY_history, '-ob'); hold off;
xlabel('Time (s)'); ylabel('Centroid (m)'); legend('X','Y'); grid on;

figure; plot(timeSteps, abs(circulation_history), '-k'); 
xlabel('Time (s)'); ylabel('Circulation (m^2/s)'); grid on;

figure; plot(centroidX_history, centroidY_history,'-o'); 
xlabel('X centroid (m)'); ylabel('Y centroid (m)'); axis equal; grid on;

results = table(timeSteps(:), centroidX_history, centroidY_history, circulation_history, ...
    'VariableNames', {'Time_s','CentroidX_m','CentroidY_m','Circulation_m2s'});
writetable(results, '1.5m-centroid-track.csv');
%%


data05 = readtable('0.5m-centroid-track.csv');
data10 = readtable('1m-centroid-track.csv');
data15 = readtable('1.5m-centroid-track.csv');

T = 2*pi;   % period

% --- non-dimensional time
tau05 = mod(data05.Time_s / T, 1);
tau10 = mod(data10.Time_s / T, 1);
tau15 = mod(data15.Time_s / T, 1);


figure;
plot(tau05, abs(data05.Circulation_m2s), '-o', 'DisplayName', '$z/c = 0.5$'); hold on;
plot(tau10, abs(data10.Circulation_m2s), '-s', 'DisplayName', '$z/c = 1$');
plot(tau15, abs(data15.Circulation_m2s), '-^', 'DisplayName', '$z/c = 1.5$'); hold off;

xlabel('$t/T$');
ylabel('$\Gamma$ ($m^2/s$)');
legend('Location','best');
grid on;

% show later half of the cycle (6/12 -> 12/12)
xlim([0.5 1]);
xticks(0.5:0.1:1.0);
