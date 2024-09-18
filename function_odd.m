function [v_opt,C] = function_odd(c0,rho0,p_temp,q_temp,freq,value_left,value_right,beta,s)
    W=2*pi*freq;
    p = sortrows(p_temp, 1);
    q = sortrows(q_temp, 1);

    %% calculate distance between ...
    numv = size(p, 1); 
    nump = size(q, 1); 
    r = zeros(numv, nump); 

    for ii = 1:numv
        for jj = 1:nump
            r(ii, jj) = sqrt((p(ii, 1) - q(jj, 1))^2 + (p(ii, 2) - q(jj, 2))^2);
        end
    end

    B_temp = rho0./(4*pi.*r);
    [~, cols] = size(B_temp);
    midcols = ceil(cols / 2);
    B_temp(:, midcols) = B_temp(:, midcols) * s;

    %% Assign a value to the microphone

    d = [value_left; value_right];%左右耳理想值desired value in left and right ear
    D = zeros(numv, 2);
    for ii = 1:numv
        D(ii,1) = mod(ii, 2);
        D(ii,2) = mod(ii-1, 2);
    end

    %% Reproduced signal
    v_opt = zeros(nump,length(W));
    C = zeros(1,length(W));
    I = eye(numv);
    for ii = 1:length(W)
        B = B_temp.*exp(-1i*(W(ii)/c0).*r);
        v_opt(:,ii) = B'* pinv(B* B'+ beta* I)* D* d;
        C(ii) = cond(B); 
    end
    
end