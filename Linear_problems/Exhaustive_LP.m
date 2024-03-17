
%%           //   OPTIMAL SOLUTION FOR A LINEAR PROGRAMMING PROBLEM   //   %%
%%                  // WITH EXHAUSTIVE SEARCH METHOD  //                   %%

% This function solves linear programming problem by first finding all 
% basic feasible solutions using fundamental theorems and then comparing them to find the optimal one 
% for the particular problem.
% Note: The function solves optimization problem defined in standard form

%% PROGRAM STARTS PART 1.1

% Function declaration
function [f,x,B] = Exhaustive_LP(A,b,c)



% The function takes as argument A, b and c.
% Number of variables = n (no of columns in the A matrix)
% Number of constraints = m (no of rows in the b matrix)


% number of variables and number of constraints
  n = size(A,2);
  m = size (b,1);
  

% If condition to check if n > m 
   if n>m

 % number of possibilities for basic feasible solution
  n_possib = nchoosek(n,m);

 % pairs of the basic feasible solutions
    n_pairs = nchoosek(1:n,m);

 % (4) Computation of basic feasible solutions loop
   all_bfs = [];
    
    for i=1:n_possib
        
        Bx = A(:,n_pairs(i,:));  % Extraction from the basis
        
        if det(Bx) ~= 0

            xb = (Bx)\b;         % Inverse of B multiplied with b
    
            if all(xb >=0 ) && all(xb ~= inf) && all(xb~= -inf) % Conditions for feasible solutions
    
                bfs = zeros(n,1);
                bfs(n_pairs(i,:),:) = xb;  % bfs matrix
                    
                all_bfs = [all_bfs, bfs];  % all basic feasible solutions
           
            else
                 disp(['The solution xb = [' num2str(xb(:).') ']'' will constitute a basic infeasible solution'])
                 fprintf('\n');
            end

        else
            disp ('')
           
        end

    end

 else 
     
    disp("Number of variables must be greater than the number of constraints");


   end

   % (5) Finding the optimal cost, optimal BFS and the corresponding basis
   cost_values = c'*all_bfs;
    
   % Optimal cost
   [f,i] = max(cost_values); %%%%%%%%% Calculating the Maximum

   % Optimal basic feasible solution
   x = all_bfs(:,i);

   % corresponding basis B
   B = find(x');

   disp (['The optimal basic feasible solution is [' num2str(x(:).') ']'])
   fprintf ('The optimal cost is %d \n', f)
   disp (['The corresponding basis is [' num2str(B(:).') ']'])

end

