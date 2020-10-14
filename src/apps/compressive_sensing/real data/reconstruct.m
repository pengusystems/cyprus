function spectra = reconstruct(measurements, measurement_matrix, reconstruction_method)
% The length of the signal to estimate.
element_length = size(measurement_matrix, 2);

% Number of samples we got from compress.
n_cs_samples = size(measurement_matrix, 1);

% Create an inverse of the sparsing dictionary recalling that:
% 1. W = dftmtx(element_length); 
% 2. c = W*x_t
W_inv = (dftmtx(element_length)'/element_length); % x_t = W_inv*c


% Reconstruct.
switch reconstruction_method
    %%%%%%%%%%%%%%%%%%%%
    % Matching Pursuit %
    %%%%%%%%%%%%%%%%%%%%
    case 'MP'
        % Compressive sensing scheme: y = A*W_inv*c -> y = B*c
        % Normalize the columns of B (an over complete dictionary).
        y = measurements;
        A = measurement_matrix;
        B = A*W_inv;        
        B_cell = num2cell(B,1);
        B_cols_norm = cellfun(@norm, B_cell);

        % Initialize the coefficients containers.
        k = ceil(n_cs_samples/log(element_length));
        coeffs = zeros(k,1);
        c = zeros(element_length,1);
        elements_used = zeros(k,1);    
        
        % Iterate MP.
        ii = 1;    
        while (ii <= k)
            % Project the measurements on the dictionary B
            % and select the coefficient with the largest abs imaginery value.
            prj = B'*y ./ B_cols_norm';
            [~, I] = sort(abs(imag(prj)), 'descend');
            dup = find(elements_used == I(1));
            
            % If element hasn't been used before store it
            % and subtract its contribution (projection normalized).
            % Otherwise, add the new contribution to the one already
            % stored.
            if isempty(dup)
                coeffs(ii) = prj(I(1)) / B_cols_norm(I(1));
                c(I(1)) = coeffs(ii);                
                y = y - coeffs(ii) * B(:,I(1));
                elements_used(ii) = I(1);
                ii = ii + 1;
            else
                temp = prj(I(1)) / B_cols_norm(I(1));
                y = y - temp * B(:,I(1));
                coeffs(dup) = coeffs(dup) + temp;
            end
        end
        spectra = c;   
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Matching Pursit - symmetric %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'MPSymmetryEven'
        if (mod(element_length, 2) ~= 0)
            error('MPSymmetryEven requires a measurement matrix with an even number of columns');
        end
        % Compressive sensing scheme: y = A*W_inv*c -> y = B*c
        % Normalize the columns of B (an over complete dictionary).        
        y = measurements;
        A = measurement_matrix;
        B = A*W_inv;        
        B_cell = num2cell(B,1);
        B_cols_norm = cellfun(@norm, B_cell);

        % Initialize the coefficients containers.
        k = 2 * ceil(n_cs_samples/log(element_length));
        coeffs = zeros(k,1);
        c = zeros(element_length,1);
        elements_used = zeros(k,1);    
        
        % Iterate MP.
        ii = 1;    
        while (ii <= k)
            % Project the measurements on the dictionary B
            % and select the coefficient with the largest abs imaginery value.
            prj = B'*y ./ B_cols_norm';
            [~, I] = sort(abs(imag(prj)), 'descend');
            coeff1_index = I(1);
            
            % If it is the DC or center component.
            if ((coeff1_index == 1) || (coeff1_index == element_length / 2))
                dup = find(elements_used == coeff1_index);
                if isempty(dup)
                    coeffs(ii) = prj(coeff1_index) / B_cols_norm(coeff1_index);
                    c(coeff1_index) = coeffs(ii);                
                    y = y - coeffs(ii) * B(:,coeff1_index);
                    elements_used(ii) = coeff1_index;
                    ii = ii + 1;
                else
                    temp = prj(coeff1_index) / B_cols_norm(I(1));
                    y = y - temp * B(:,coeff1_index);
                    coeffs(dup) = coeffs(dup) + temp;
                end                
                
            % Else take advantage of the symmetry.
            else
                coeff2_index = element_length - coeff1_index + 2;                            
                dup1 = find(elements_used == coeff1_index);
                dup2 = find(elements_used == coeff2_index);            
                if (isempty(dup1) && isempty(dup2))
                    coeffs(ii)   = prj(coeff1_index) / B_cols_norm(coeff1_index);
                    coeffs(ii+1) = prj(coeff2_index) / B_cols_norm(coeff2_index);                
                    c(coeff1_index) = coeffs(ii); 
                    c(coeff2_index) = coeffs(ii+1);                                
                    y = y - coeffs(ii) * B(:,coeff1_index) - coeffs(ii+1) * B(:,coeff2_index);
                    elements_used(ii)   = coeff1_index;
                    elements_used(ii+1) = coeff2_index;                
                    ii = ii + 2;
                else
                    temp1 = prj(coeff1_index) / B_cols_norm(coeff1_index);
                    temp2 = prj(coeff2_index) / B_cols_norm(coeff2_index);                
                    coeffs(dup1) = coeffs(dup1) + temp1;
                    coeffs(dup2) = coeffs(dup2) + temp2;
                    y = y - temp1 * B(:,coeff1_index) - temp2 * B(:,coeff2_index);                
                end
            end
        end
        spectra = c;          
end