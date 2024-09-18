clear;
close all;

c0 = 343; % Speed of sound (m/s)
rho0 = 1.2; % Air density (kg/m^3)
freq = 10:10:500; % Define frequencies of interest
beta = 0.001;
W=2*pi*freq;
left = 1;
right = 0;
q_y = 1.5;
s = sqrt(2);

%% Input position of sources and receivers(listeners' ear)

% low frequency range sources positions
q_number_low = input('number of sources in low freqency :\n ');
q_low = [];
q_spacing = input('spacing of sources :\n');
% q_y = input('y coordinate of sources :\n');
q_each = [q_number_low,2];
if mod(q_number_low, 2) == 0
    for i = 1:q_number_low/2
        qn = [-(i-1/2)*q_spacing,q_y];
        q_each(i,:) = qn;
    end
    for i = 1:q_number_low/2
        qm = [(i-1/2)*q_spacing,q_y];
        q_each(i+q_number_low/2,:) = qm;
    end
else
    for i = 1:(q_number_low-1)/2
        qn = [-i*q_spacing,q_y];
        q_each(i,:) = qn;
    end
    for i = 1:(q_number_low-1)/2
        qn = [i*q_spacing,q_y];
        q_each(i+(q_number_low-1)/2+1,:) = qn;
    end
    q_each((q_number_low+1)/2,:)=[0,q_y];
end
q_low =  cat(1,q_low,q_each);
q_low = sortrows(q_low, 1);

% mid frequency sources positions
q_number_mid = input('number of sources in mid freqency :\n ');
q_mid = [];
q_spacing = input('spacing of sources :\n');
% q_y = input('y coordinate of sources :\n');
q_each = [q_number_mid,2];
if mod(q_number_mid, 2) == 0
    for i = 1:q_number_mid/2
        qn = [-(i-1/2)*q_spacing,q_y];
        q_each(i,:) = qn;
    end
    for i = 1:q_number_mid/2
        qm = [(i-1/2)*q_spacing,q_y];
        q_each(i+q_number_mid/2,:) = qm;
    end
else
    for i = 1:(q_number_mid-1)/2
        qn = [-i*q_spacing,q_y];
        q_each(i,:) = qn;
    end
    for i = 1:(q_number_mid-1)/2
        qn = [i*q_spacing,q_y];
        q_each(i+(q_number_mid-1)/2+1,:) = qn;
    end
    q_each((q_number_mid+1)/2,:)=[0,q_y];
end
q_mid =  cat(1,q_mid,q_each);
q_mid = sortrows(q_mid, 1);

% high frequency range sources positions
q_number_high = input('number of sources in high freqency :\n ');
q_high = [];
q_spacing = input('spacing of sources :\n');
% q_y = input('y coordinate of sources :\n');
q_each = [q_number_high,2];
if mod(q_number_high, 2) == 0
    for i = 1:q_number_high/2
        qn = [-(i-1/2)*q_spacing,q_y];
        q_each(i,:) = qn;
    end
    for i = 1:q_number_high/2
        qm = [(i-1/2)*q_spacing,q_y];
        q_each(i+q_number_high/2,:) = qm;
    end
else
    for i = 1:(q_number_high-1)/2
        qn = [-i*q_spacing,q_y];
        q_each(i,:) = qn;
    end
    for i = 1:(q_number_high-1)/2
        qn = [i*q_spacing,q_y];
        q_each(i+(q_number_high-1)/2+1,:) = qn;
    end
    q_each((q_number_high+1)/2,:)=[0,q_y];
end
q_high =  cat(1,q_high,q_each);
q_high = sortrows(q_high, 1);

% listeners' positions
p_number = input('number of listeners :\n');
p_left =  zeros(p_number, 2);
p_right =  zeros(p_number, 2);
for i = 1:p_number
    fprintf('No.%d listener position :\n', i);
    x = input('x-coordinate :\n');
    p_left(i, :) = [x-0.1 0];
    p_right(i, :) = [x+0.1 0];
end
p = [p_left; p_right];

% p = [-0.45 0; -0.25 0; 0.25 0; 0.45 0];
q_total = [q_low; q_mid; q_high];


%% plot position of sources and receivers

figure;
hold on;
plot(q_low(:,1), q_low(:,2), 'ro', 'MarkerFaceColor','r');
plot(q_mid(:,1), q_mid(:,2), 'r*', 'MarkerFaceColor','r');
plot(q_high(:,1), q_high(:,2), 'rx', 'MarkerFaceColor','r');
plot(p(:,1),p(:,2), 'b<', 'MarkerFaceColor', 'b');
xlabel('X Coordinate/m');
ylabel('Y Coordinate/m');
title('Position of Sound Sources and Listeners');
legend('low_freq sources','mid_freq sources','high_freq sources', 'Listeners Ear');
axis equal;
grid on;
hold off;


%% Plot condition number vs freqency

if mod(q_number_low, 2) == 1
    [~,C_low]=function_odd(c0,rho0,p,q_low,freq,left,right,beta,s);
else
    [~,C_low]=function_even(c0,rho0,p,q_low,freq,left,right,beta);
end
figure('Position',[0,0,400,300]);
hold on;
semilogy(freq,C_low, 'LineWidth', 1.5);
title('Condition Number vs Mid Frequency');
xlabel('Frequency (Hz)');
ylabel('Condition Number');
set(gca, 'YScale', 'log');
grid on;
ax = gca;
ax.Box = 'on'; 
ax.LineWidth = 1; 
ax.TickDir = 'in';
hold off;

if mod(q_number_mid, 2) == 1
    [~,C_mid]=function_odd(c0,rho0,p,q_mid,freq,left,right,beta,s);
else
    [~,C_mid]=function_even(c0,rho0,p,q_mid,freq,left,right,beta);
end
figure;
hold on;
semilogy(freq,C_mid, 'LineWidth', 1.5);
title('Condition Number vs Mid Frequency');
xlabel('Frequency (Hz)');
ylabel('Condition Number');
set(gca, 'YScale', 'log');
grid on;
hold off;

if mod(q_number_high, 2) == 1
    [~,C_high]=function_odd(c0,rho0,p,q_high,freq,left,right,beta,s);
else
    [~,C_high]=function_even(c0,rho0,p,q_high,freq,left,right,beta);
end
figure;
hold on;
semilogy(freq,C_high, 'LineWidth', 1.5);
title('Condition Number vs High Frequency'); 
xlabel('Frequency (Hz)');
ylabel('Condition Number');
set(gca, 'YScale', 'log');
grid on;
hold off;

figure('Position',[0,0,400,300]);
hold on;
% semilogy(freq(freq <= 2500), C_low(freq <= 2500), 'LineWidth', 1.5, 'Color', 'b');
% semilogy(freq(freq >= 2400 & freq <= 7200), C_mid(freq >= 2400 & freq <= 7200), 'LineWidth', 1.5, 'Color', 'r');
% semilogy(freq(freq >= 7000), C_high(freq >= 7000), 'LineWidth', 1.5, 'Color', 'g');
semilogy(freq,C_low, 'LineWidth', 1.5);
semilogy(freq,C_mid, 'LineWidth', 1.5);
% semilogy(freq,C_high, 'LineWidth', 1.5);
% title('Condition Number vs Frequency');
xlabel('Frequency (Hz)');
ylabel('Condition Number');
set(gca, 'YScale', 'log');
grid on;
legend('400-2500 Hz', '2400-7200 Hz', 'Over 7000 Hz');
ax = gca;
ax.Box = 'on';
ax.LineWidth = 1;
ax.TickDir = 'in';
hold off;


%% Construct a room plane

n = 500;  % sample point
x = linspace(-2, 2, n)';
y = linspace(-1, 2, n)';
[X, Y] = meshgrid(x, y); 
figure;
plot(x(194),y(167),'x');
hold on
plot(-0.45, 0,'x');
hold off


%% Plot sound field in whole room

value_rep = zeros(4,length(W));

for ii = 1:length(W)% for all frequency
    % in low, mid and high frequency range, the code only change parameters:
    % 'q_number' and 'q'
    % and color range for low is ()
    f_tmp = freq(ii); % This is the frequency we are investigating at the moment
    if f_tmp<2100 % low frequency f<2600Hz
        if mod(q_number_low, 2) == 1
            [v,~]=function_odd(c0,rho0,p,q_low,f_tmp,left,right,beta,s);
        else
            [v,~]=function_even(c0,rho0,p,q_low,f_tmp,left,right,beta);
        end
        nump = size(q_low, 1);
        C = zeros(n, n, nump);
        r = zeros(n, n, nump);% Dimension 3 is for each source (Here is 4 sources)
        for jj = 1:nump
            r(:, :, jj) = sqrt((X - q_low(jj, 1)).^2 + (Y - q_low(jj, 2)).^2);
            r(r == 0) = 1e-6;
        end
        p_room = zeros(n,n);
        for kk = 1:n % for x-coordinate
            for jj =1:n % for y-coordinate
                C(kk, jj, :) = rho0*exp(-1i*(W(ii)/c0)*r(kk,jj,:))./(4*pi*r(kk,jj,:));
                u = size(nump);
                for pp =1:nump % Superposition the sound pressure of each sound source to a point
                    u(pp) = C(kk, jj, pp)*v(pp);
                end
                p_room(kk, jj) = sum(u);
            end
        end
        value_rep(1,ii) = p_room(167,194);
        value_rep(2,ii) = p_room(167,219);
        value_rep(3,ii) = p_room(167,281);
        value_rep(4,ii) = p_room(167,306);
        figure(6);
        surf(X, Y, 20*log10(abs(p_room)), 'EdgeColor', 'none');
        hold on
        maxZ = max(20*log10(abs(p_room)), [], 'all') + 10;
        h1 = plot3(q_low(:,1), q_low(:,2),repmat(maxZ, size(q_low, 1),1),'*m');
        h2 = plot3(p(:,1),p(:,2),repmat(maxZ, size(p, 1),1), 'squarer','LineWidth', 2,'MarkerSize', 10);
        view(2);
        colorbar;
        clim([-20 20]); 
        colormap('jet'); 
        % title(sprintf('Frequency: %d Hz', freq(ii)));
        xlabel('X Coordinate');
        ylabel('Y Coordinate');
        % legend([h1, h2], {'Sources', 'Receivers'}, 'Location', 'northeastoutside');
        hold off;
        pause(0.00001);
    elseif f_tmp>=2100 && f_tmp<7000  % mid frequency 2600Hz<f<7300Hz
        if mod(q_number_mid, 2) == 1
            [v,~]=function_odd(c0,rho0,p,q_mid,f_tmp,left,right,beta,s);
        else
            [v,~]=function_even(c0,rho0,p,q_mid,f_tmp,left,right,0);
        end
        nump = size(q_mid, 1);
        C = zeros(n, n, nump);
        r = zeros(n, n, nump);
        for jj = 1:nump
            r(:, :, jj) = sqrt((X - q_mid(jj, 1)).^2 + (Y - q_mid(jj, 2)).^2);
            r(r == 0) = 1e-6;
        end
        p_room = zeros(n,n);
        for kk = 1:n
            for jj =1:n
                C(kk, jj, :) = rho0*exp(-1i*(W(ii)/c0)*r(kk,jj,:))./(4*pi*r(kk,jj,:));
                u = size(nump);
                for pp =1:nump
                    u(pp) = C(kk, jj, pp)*v(pp);
                end
                p_room(kk, jj) = sum(u);
            end
        end
        value_rep(1,ii) = p_room(167,194);
        value_rep(2,ii) = p_room(167,219);
        value_rep(3,ii) = p_room(167,281);
        value_rep(4,ii) = p_room(167,306);        
        figure(6);
        surf(X, Y, 20*log10(abs(p_room)), 'EdgeColor', 'none');
        hold on
        maxZ = max(20*log10(abs(p_room)), [], 'all') + 10;
        h1 = plot3(q_mid(:,1), q_mid(:,2),repmat(maxZ, size(q_mid, 1),1),'*m');
        h2 = plot3(p(:,1),p(:,2),repmat(maxZ, size(p, 1),1), 'squarer','LineWidth', 2,'MarkerSize', 10);
        view(2);
        colorbar;
        clim([-20 20]); 
        colormap('jet'); 
        % title(sprintf('Frequency: %d Hz', freq(ii)));
        xlabel('X Coordinate');
        ylabel('Y Coordinate');
        % legend([h1, h2], {'Sources', 'Receivers'}, 'Location', 'northeastoutside');
        hold off;
        pause(0.00001);
    else % high frequency f>7300Hz
        if mod(q_number_high, 2) == 1
            [v,~]=function_odd(c0,rho0,p,q_high,f_tmp,left,right,beta,s);
        else
            [v,~]=function_even(c0,rho0,p,q_high,f_tmp,left,right,0);
        end
        nump = size(q_high, 1);
        C = zeros(n, n, nump);
        r = zeros(n, n, nump);
        for jj = 1:nump
            r(:, :, jj) = sqrt((X - q_high(jj, 1)).^2 + (Y - q_high(jj, 2)).^2);
            r(r == 0) = 1e-6;
        end
        p_room = zeros(n,n);
        for kk = 1:n
            for jj =1:n
                C(kk, jj, :) = rho0*exp(-1i*(W(ii)/c0)*r(kk,jj,:))./(4*pi*r(kk,jj,:));
                u = size(nump);
                for pp =1:nump
                    u(pp) = C(kk, jj, pp)*v(pp);
                end
                p_room(kk, jj) = sum(u);
            end
        end
        value_rep(1,ii) = p_room(167,194);
        value_rep(2,ii) = p_room(167,219);
        value_rep(3,ii) = p_room(167,281);
        value_rep(4,ii) = p_room(167,306);        
        figure(6);
        surf(X, Y, 20*log10(abs(p_room)), 'EdgeColor', 'none');
        hold on
        maxZ = max(20*log10(abs(p_room)), [], 'all') + 10;
        h1 = plot3(q_high(:,1), q_high(:,2),repmat(maxZ, size(q_high, 1),1),'*m');
        h2 = plot3(p(:,1),p(:,2),repmat(maxZ, size(p, 1),1), 'squarer','LineWidth', 2,'MarkerSize', 10);
        view(2);
        colorbar;
        clim([-30 30]);
        colormap('jet');
        % title(sprintf('Frequency: %d Hz', freq(ii)));
        xlabel('X Coordinate');
        ylabel('Y Coordinate');
        % legend([h1, h2], {'Sources', 'Receivers'}, 'Location', 'northeastoutside');
        hold off;
        pause(0.00001);
    end
end


%% error

figure('Position',[0,0,400,300]);
hold on;
for ii = 1:size(value_rep, 1)
    plot(freq,abs(value_rep(ii,:)).^2, 'LineWidth', 1.5);
end

ylabel('Error');
legend('left_1','right_1', 'left_2','right_2');
grid on;
ax = gca;
ax.Box = 'on'; 
ax.LineWidth = 1; 
ax.TickDir = 'in';
hold off; 
figure('Position',[0,0,400,300]);
hold on;
for ii = 1:2:size(value_rep,1)
    plot(freq,10*log10(abs(value_rep(ii,:)).^2./abs(value_rep(ii+1,:)).^2), 'LineWidth', 1.5);
end
% title('cross-talk cancellation performance');
xlabel('Frequency (Hz)');
ylabel('cross-talk cancellation');
set(gca, 'YScale', 'log');
legend('listener 1','listener 2');
grid on;
ax = gca;
ax.Box = 'on'; 
ax.LineWidth = 1; 
ax.TickDir = 'in';
hold off;
