function [x, detA] = Cholesky(main_diag, upper_diag1, upper_diag2, lower_diag1, lower_diag2, b)
    % -Cholesky method for solving Ax = b
    % -(main_diag, upper_diag1, upper_diag2, lower_diag1, lower_diag2) vectors are the main
    % diagnoals of the symmetric positive pentadiagonal definite matrix (A)
    % -If the size of the matrix is n, then main_diag vector has n elements,
    % upper_diag1 has n-1 elements, upper_diag2 has n-2, and lower_diags
    % are relativly equal to upper_diags
    % -b: right-hand side vector
    % -x: vector output
    % -detA : is the determinent of A
    % -A is represented by 5 vectors and L is represented by 3 vectors
    n = length(main_diag);  % Size of the matrix
    
    % Cholesky factorization
    L1 = zeros(n, 1);
    L2 = zeros(n-1, 1);
    L3 = zeros(n-2, 1);
    detA = 1; % initialize determinant of L to 1
    
    % Calculate L1(1)
    L1(1) = sqrt(main_diag(1));
    detA = detA * L1(1); % update determinant of L
    
    % Calculate L2(1) and L1(2)
    L2(1) = upper_diag1(1) / L1(1);
    L1(2) = sqrt(main_diag(2) - L2(1)^2);
    detA = detA * L1(2); % update determinant of L
    
    % Calculate L3(1), L3(2),L1(3) and L2(2)
    L3(1) = upper_diag2(1) / L1(1);
    L3(2) = (upper_diag2(2) - L1(1) * 0) / L1(2);
    L2(2) = (upper_diag1(2) - L2(1)*L3(1)) /L1(2);
    L1(3) = sqrt(main_diag(3) - L3(1)^2 - L2(2)^2);
    detA = detA * L1(3); % update determinant of L
   
    

    % Calculate remaining elements of L
    
    for j = 4:n
        L3(j-2) = (upper_diag2(j-2)) / L1(j-2);
        L2(j-1) = (lower_diag1(j-1) - L3(j-2) * L2(j-2)) / L1(j-1);
        L1(j) = sqrt(main_diag(j) - L3(j-2)^2 - L2(j-1)^2);
        
        
        
        detA = detA * L1(j); % update determinant of L
        
    end
    %adding the last item of L2
    
    detA = detA * detA;
    
 

    % Forward substitution to solve Ly = b
    y = zeros(n, 1);
    y(1) = b(1) / L1(1);
    
    for i = 2:n
        if i == 2
            y(i) = (b(i) - L2(i-1) * y(i-1)) / L1(i);
        else
            y(i) = (b(i) - L2(i-1) * y(i-1) - L3(i-2) * y(i-2)) / L1(i);
        end
    end
    
    % Backward substitution to solve L'x = y
    x = zeros(n, 1);
    x(n) = y(n) / L1(n);
    
    for i = n-1:-1:1
        if i == n-1
            x(i) = (y(i) - L2(i) * x(i+1)) / L1(i);
        elseif i == n-2
            x(i) = (y(i) - L2(i) * x(i+1) - L3(i) * x(i+2)) / L1(i);
        else
            x(i) = (y(i) - L2(i) * x(i+1) - L3(i) * x(i+2)) / L1(i);
        end
    end
end
