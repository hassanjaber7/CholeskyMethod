function testCholesky(n)
isPositive=false;

while (isPositive == false)
     Q = randi([-10, 10],n);
    P = Q'*Q;
    disp(Q');
    disp(P);
   A = zeros(n);
   A = A + diag(diag(P)) + diag(diag(P,1), 1) + diag(diag(P,2), 2) + diag(diag(P,-1), -1) + diag(diag(P,-2), -2);
    main_diag = diag(P);
    upper_diag1 = diag(P,1);
    upper_diag2 = diag(P,2);
    lower_diag1 = diag(P,-1);
    lower_diag2 = diag(P,-2);
    
    
    
    [R,flag] = chol(A);
    
    % Display the result
    if flag ==0
        isPositive=true;
        disp('The matrix is positive definite.');
        disp(A)
    
    end


end
 b = randi([1, 1000],n,1);
 % Start the timer
tic;

 % Call the Cholesky function to solve Ax = b
 
   [x, detA] = Cholesky(main_diag, upper_diag1, upper_diag2, lower_diag1, lower_diag2, b);
% Stop the timer
elapsedTime = toc;
 disp(chol(A))
 % Compute the condition number of A
   cond_A = norm(A)*norm(inv(A));

   expected_x = linsolve(A,b); 
   
    % Compute the relative error
    error1 = norm(x - expected_x) / norm(expected_x);
    %Compute the forward stability error
    error2 = norm(x - expected_x) / (norm(expected_x)*cond_A);

    %Compute the backward stability error
    error3 = norm(b - A*expected_x)/(norm(A)*norm(expected_x));
    
    fprintf('Condition Number (cond(A)): %.6f\n', cond_A);
    fprintf('Relative Error: %.20f\n', error1);
    fprintf('Forward Stability Error: %.20f\n', error2);
    fprintf('Backward stability Error: %.20f\n', error3);
    % Display the elapsed time
    disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);
