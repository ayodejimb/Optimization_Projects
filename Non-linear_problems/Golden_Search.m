
%%              OPTIMAL SOLUTION FOR A NON-LINEAR PROGRAMMING PROBLEM              %%
%%                   // WITH GOLDEN SECTION   METHOD  //                %%

% This programm use GOLDEN_SECTION method to find the minimizer to a non-linear 
% programming problem of interval x to y

%% PART 2.1 
clear all; %% Clearing workspace
clc

%%%%% Objective function
f = @(x) x.^4 + 4*x.^3 + 9*x.^2 + 6*x + 6; 

%%%% Interval [x y] and the stopping criteria (eps)
x = -2; %%% Lower limit
y = 2;  %%% Upper limit
eps = 10^-2;

%%%%%%% Golden ratio, constant
phi = 0.618;

%%%%%%%%%% Defining interval
x1 = y - (y - x)*phi;  
x2 = x + (y - x)*phi;
i = 0;

%%%%%%%%%%%% Loop till tolerance isn't satisfied.
while abs(x1 - x2) > eps

    %%% Check function at x1 and x2
    fx1 = f(x1);
    fx2 = f(x2);

    %%%%% Update interval
   if fx1 < fx2
       y = x2;
       x2 = x1;
       x1 = y - (y - x)*phi;

   else
       x = x1;
       x1 = x2;
       x2 = x + (y - x)*phi;
   end
end


%%%%%%%%% Plotting the graph for visualization
figure;
hold on;
a = -0.65:0.02: -0.25;
plot(a, f(a));
plot(x2, f(x2), "x", "color", "red");
plot(x1, f(x1), "x", "color", "red");
title("Graph of F(x) against x")
xlabel("x")
ylabel("F(x)")
legend("F(x)", "Interval (x2 - x1)")
hold off;

%%%%%% Optimal interval (x1 and x2)
range = [x1, x2];

%%%%%%% Computing the average of the two limits
xmin = (x1 + x2)/2;

%%%%%%% Computing the optimal cost function
fmin = f(xmin);


%%%%% Displaying the result
fprintf ('\n')
disp (['The range is [' num2str(range) ']'])  
disp (['The Minimizer is [' num2str(xmin) ']'])
disp (['The optimal cost is [' num2str(fmin) ']'])
       