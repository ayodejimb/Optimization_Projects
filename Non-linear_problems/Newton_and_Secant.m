%%              OPTIMAL SOLUTION FOR A NON-LINEAR PROGRAMMING PROBLEM              %%
%%                   // WITH  NEWTON'S  METHOD  //                       %%

% This program uses NEWTON'S method to find the minimizer to a non-linear 
% programming problem

%% (a) Newton's Method
clear all; %% Clearing workspace
clc


% Using x as a variable
syms x

%%%%% Objective function
f = @(x) 2*x^4 - 5*x^3 + 100*x^2 + 30*x - 75;

% The initial condition
xo = 2;

% The first and Second derivative of the function
f_prime = diff(f,x); 
f_doub_prime = diff(f_prime,x);

% Using the inline function to allow computation of fprime(xo) and
% f_doub_prime(xo)
f_prime_x = inline(f_prime);
f_doub_prime_x = inline(f_doub_prime);

% The Newton first Iteration
x = xo - f_prime_x(xo)/f_doub_prime_x(xo);

% Computing the check condition
f_prime_xx = f_prime_x(xo);

% The stopping criterion 
while (f_prime_xx >= 1*10^-4) % Do this repeatedly until the stopping criterion is true
    
    % make xo = x
    xo = x;

    % Do the next iteration
    x = xo - f_prime_x(xo)/f_doub_prime_x(xo);
    
    % The check condition
    f_prime_xx = f_prime_x(xo);

end 

% Minimized cost
f_min = f(x);

%%%% Displaying the results
fprintf ('The minimized cost is %d \n', f_min)
fprintf ('The minimized variable is %d \n', x)


% Plotting the function
ezplot(f,[-2 2]);
grid
title('Plot of the function 2x^4 - 5x^3 + 100x^2 + 30x - 75')
hold on
plot(x, f_min,'*')
 

%% (b) Secant's Method
clear all; %% Clearing workspace
clc
format short;

% Using x as a variable
syms x

%%%%% Objective function
f = @(x) 2*x^4 - 5*x^3 + 100*x^2 + 30*x - 75;

% The initial conditions
x0 = 2;
x1 = 2.1;

% The first derivative of the function
f_prime = diff(f,x); 

% Using the inline function to allow computation of fprime(xo)
f_prime_x = inline(f_prime);

% To be used in the iteration of the secant
f_prime_x0 = f_prime_x(x0);

% The Check condition
f_prime_x1 = f_prime_x(x1);

% The stopping criterion 
while (f_prime_x1 >= 1*10^-4) % Do this repeatedly until the stopping criterion is true
    
    % The Secant Iteration
    x2 = x1 - (( (x1-x0)/(f_prime_x1-f_prime_x0) ) * f_prime_x1);

     x0 = x1;

     x1 = x2;
    
    % The Check condition
    f_prime_x1 = f_prime_x(x1);

end 

% Minimized cost
f_min = f(x1);

%%%% Displaying the results
fprintf ('The minimized cost is %d \n', f_min)
fprintf ('The minimized variable is %d \n', x1)

% Plotting the Graph
ezplot(f,[-2 2]);
grid
title('Plot of the function 2x^4 - 5x^3 + 100x^2 + 30x - 75')
ylabel("f(x)")
xlabel("x")
hold on
plot(x1, f_min,'*')





