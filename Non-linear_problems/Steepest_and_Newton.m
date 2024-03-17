%% The Rosenbrock function using steepest descent method but Newton's Method for the Line Search %%

clear;
close; clc;

% The function (this allows using a variable (alpha) in the function)
syms f(x1,x2)
f(x1,x2) = (1 - x1)^2 + 100*(- x1^2 + x2)^2;

% The gradient of the function
g = gradient(f);
gf = inline(g);

% The Hessian of the function
h = hessian(f);
hf = inline(h);

syms alphaa
syms argmin_k

% To plot the result
fsurf(f,[-100 100 -100 100]);
hold on;

% To store the x1, x2 and f(x1,x2) computed for each steepest descent iteration
% It is useful to know the number of iterations performed before the stopping criterion is reached.  
func_values = [];
x1_values = [];
x2_values = [];

% Initial alpha given
alpha = 0.1;

% Initial xo given
xo = [-1;2];

% Steepest Descent computation
x = xo - alpha*(gf(xo(1),xo(2)));

% To store the computation result for plotting
title('The Plot of Rosenbrock iteraions using Steepest Descent','FontSize',12)
xlabel('x1','FontSize',12)
ylabel('x2','FontSize',12)
zlabel('f(x1,x2)','FontSize',12)

plot_f = f(x(1),x(2));
plot3(x(1),x(2),double(plot_f),'o-','linewidth',3)
hold on;

func_values = [func_values, double(plot_f)];
x1_values = [x1_values , x(1)];
x2_values = [x2_values, x(2)];

% Steepest Descent Check condition
check_condition = (norm(x -xo))/norm(xo);

while (check_condition > 10^-3)
    % Assign the new x to the previous one
    xo = x;

    % Substituting xo into the function
    % One-dimensional search
    find_alpha = xo - alphaa*(gf(xo(1),xo(2)));
    argmin_k= f(find_alpha(1), find_alpha(2));

    % Convert from Symbolic to Function
    f_newton = matlabFunction(argmin_k);

    % The first and Second derivative of the function
    f_prime = diff(f_newton, alphaa); 
    f_doub_prime = diff(f_prime,alphaa);
    
    % Using the inline function to do the first and second derivative for Newton
    f_prime_alpha = inline(f_prime);
    f_doub_prime_alpha = inline(f_doub_prime);
    
    % Subsituting alpha into the first and second derivative
    f_prime_alphaa = f_prime_alpha(alpha);
    f_doub_prime_alphaa = f_doub_prime_alpha(alpha);
    
    % Computing the new alpha
    alpha_k_plus_one = alpha - f_prime_alphaa/f_doub_prime_alphaa;

    % Newton iteration check condition 
    check_condition_2 = (abs(alpha_k_plus_one - alpha))/alpha;

    % The Newton one_dimensional search for optimal alpha
    while check_condition_2 > 1*10^-4  % Do this repeatedly until the stopping criterion is reached
        
        alpha = alpha_k_plus_one;
        
        f_prime_alphaa = f_prime_alpha(alpha);
        f_doub_prime_alphaa = f_doub_prime_alpha(alpha);
        
        alpha_k_plus_one = alpha - f_prime_alphaa/f_doub_prime_alphaa;

        % Newton check condition 
        check_condition_2 = (abs(alpha_k_plus_one - alpha))/alpha;

    end

        % Do the Steepest Descent with the new alpha
        x = xo - alpha_k_plus_one*(gf(xo(1),xo(2)));

        % Steepest Descent Check condition
        check_condition = (norm(x -xo))/norm(xo);

        % Assign the new alpha to the previous one
        alpha = alpha_k_plus_one;


        % For plotting
        plot_f = f(x(1),x(2));
        plot3(x(1),x(2),double(plot_f),'o-','linewidth',3)
        hold on
        
        func_values = [func_values, double(plot_f)];
        x1_values = [x1_values , x(1)];
        x2_values = [x2_values, x(2)];

end

% Count the number of iterations and print to window
iteration_count = size(x1_values,2);
fprintf('%d iterations completed until the stopping criteria is reached \n', iteration_count);