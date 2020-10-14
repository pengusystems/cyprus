function spectra = reconstruct(measurements, measurement_matrix, reconstruction_method)
element_length = size(measurement_matrix, 2);
n_cs_samples = size(measurement_matrix, 1);
W = dftmtx(element_length); % c = W*x_t
W_inv = (W'/element_length); % x_t = W_inv * c

switch reconstruction_method
    %%%%%%%%%%
    % L1 Min %
    %%%%%%%%%%    
    case 'L1'
        y = measurements;
        A = measurement_matrix;
        B = A*W_inv;
        epsilon1 = 0.1;        
        epsilon2 = 0.5;
        
        cvx_begin
            variable c(element_length) complex;         
            minimize(norm(c,1)) ;
            subject to 
                norm(y-B*c,2) <= epsilon1;
                %abs(c(1:floor(element_length/2)) + c(ceil(element_length/2)+1:end)) <= epsilon2;
        cvx_end
        spectra = c;

    %%%%%%%%%%%%%%%%%%%%
    % Matching Pursuit %
    %%%%%%%%%%%%%%%%%%%%
    case 'MP'
        y = measurements;
        A = measurement_matrix;
        B = A*W_inv;
        B_cell = num2cell(B,1);
        B_rows_norm = cellfun(@norm, B_cell);

        k = ceil(n_cs_samples/log(element_length));
        coeffs = zeros(k,1);
        c = zeros(element_length,1);
        elements_used = zeros(k,1);    
        ii = 1;    
        while (ii <= k)
            prj = B'*y ./ B_rows_norm';
            [~, I] = sort(abs(prj), 'descend');
            dup = find(elements_used == I(1));
            if isempty(dup)
                coeffs(ii) = prj(I(1)) / B_rows_norm(I(1));
                c(I(1)) = coeffs(ii);
                y = y - coeffs(ii) * B(:,I(1));
                elements_used(ii) = I(1);
                ii = ii + 1;
            else
                temp = prj(I(1)) / B_rows_norm(I(1));
                y = y - temp * B(:,I(1));
                coeffs(dup) = coeffs(dup) + temp;
            end
        end
        spectra = c; 
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Iterative Hard Thresholding For Structured Sparsity %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'ITH'    
        exit; % This is not working. Need to add the structuring.
        y = measurements;
        A = measurement_matrix;
        B = A*W_inv;
        B_cell = num2cell(B,1);
        B_rows_norm = cellfun(@norm, B_cell);

        k = ceil(n_cs_samples/log(element_length));
        c = zeros(element_length,1);
        ii = 1;    
        while (ii <= k)
            % Compute signal estimate.
            prj = B'*y ./ B_rows_norm';
            
            % FIXME.            
            % Prune estimate according to structure.
            %c = M(c_before_structuring, k);
            [~, I] = sort(abs(prj), 'descend');
            prj(I(3:end)) = 0;
            c = c + prj;
            
            % Update measurement residuals.
            temp = prj(I(1)) / B_rows_norm(I(1));
            y = y - temp * B(:,I(1));    
            temp = prj(I(2)) / B_rows_norm(I(2));
            y = y - temp * B(:,I(2));  
            ii = ii + 2;
        end
        spectra = c;         
end