
%%     //         OPTIMAL SOLUTION FOR A LINEAR PROGRAMMING PROBLEM        //   %%
%%                   // WITH SIMPLEX METHOD  //                       %%

% The one-phase simplex method is based on the assumption that the we have
% an initial basic feasible solution from the standard form

% A = A matrix containing the coefficient of the variables (m x n matrix)
% b = The right hand side of the constraint (Column vector)
% v = indices of the columns of the basis corresponding to an initial BFS (Row vector)
% c = coefficients of the objective function  (Column vector)
%% PROGRAM STARTS  PART 1.2

% Function declaration
function [f,x,B] = LP_Simplex(A,b,c,v)


% Forming the simplex tableau
costs = -c'; %%%%% Cost = -c to solve maximization problem, Cost = c to solve minimization problem
tableau = [costs 0; A, b];


% Determining the number of non-basic initial variables
 nbivrbles = size(A,2) - size(v,2); 

% Determine the basic variables values in the objective function using
% their indices in the tableau matrix
for i=1:size(v,2)

    basics_var(i) = costs(v(i));

end

% Compute the reduced cost rq
for j = 1:size(tableau,2)       %No of columns in the simplex tableau

    add_val = 0;

  for k = 1:size(basics_var,2)
       
      val = tableau(k+1,j) * basics_var(1,k);  

      add_val = add_val + val;   %addition loop (Zq)

  end 

    rq(j) = tableau(1,j) - add_val;  %rq = Cq - Zq

end

% Adding the reduced cost rq result to the tableau table
    tableau = [tableau; rq];

% Putting into a table
    v_current = v;
    simplex_tab = array2table(tableau);
    
    for i = 1:size(simplex_tab,2)

    simplex_tab.Properties.VariableNames(i) = strcat ('x', {num2str(i)});

    end

    simplex_tab.Properties.VariableNames(i) = {'b'} %replacing the last column title with var name 'b' 


% Condition check for Optimal BFS (no negative value in the rq? stop! Else, proceed)
%     The simplex Algorithm
    if any(rq(:,1:end-1) < 0) 
        
     while any(rq(:,1:end-1) < 0) % Do this as long as the reduced cost has a negative value

       disp ('=== The current BFS is not optimal === ')
       disp ("   === The Next Simplex Tableau === ")

       [q_value, q_index] = min(rq(:,1:end-1)); %index and value (highest -ve) of the pivot column

       % Selecting the pivot row
       b_extract = tableau(2:end-1, end); % extracting the b values only
       val_extract = tableau(2:end-1, q_index); %extracting the values of the index column

       for n = 1:size(b_extract,1)

            % Check for division by a negative value/negative of the
            % numerator
            if b_extract(n,:) && val_extract(n,:) > 0 

                pivot_rows_div(n) = b_extract(n,:)/val_extract(n,:);  

            else 

                pivot_rows_div(n) = inf;
                
            end
           
       end

       % Check if the division is not unbounded 
       if all(pivot_rows_div) ~= inf

         [p_value, p_index] = min(pivot_rows_div); %index and value of the pivot row

       % Index of the leaving and entering variable into the basis
        fprintf (['The column index of the leaving and entering variable into the basis are %d and %d ' ...
               'respectively'],v(p_index),q_index)

       % Calculating the current corresponding Basis B
        outgoing_basics = v_current(p_index);
         
         % Replacing the outgoing basis with the incoming one
         v_current(v_current==outgoing_basics) = q_index;
        
       % Selecting the pivot/key element 
           key_element = tableau(p_index+1, q_index);
        
       % Computing the unit pivot
           unit_pivot = tableau(p_index+1,:)./key_element;
        
       % Filling the simplex tableau
           for m = 1:size(tableau,1)-1  %total number of rows in the tableau
            
               if (m ~= p_index)    % is m equal to the index of the pivot row?

                   tableau(m+1,:) = tableau(m+1,:) - (unit_pivot.*tableau(m+1,q_index));

               else 
                   tableau(m+1,:) = unit_pivot;  %replace the row with the unit pivot row
               end

           end

           % Putting into a table
            simplex_tab = array2table(tableau);
            
            for i = 1:size(simplex_tab,2)
        
               simplex_tab.Properties.VariableNames(i) = strcat ('x', {num2str(i)});
        
            end
        
            simplex_tab.Properties.VariableNames(i) = {'b'} %replacing the last column title with var name 'b'

           % Check the reduced cost
           rq = tableau(end,:);

            % Selecting the optimal basic feasible solution x
                for p = 1:size(tableau,1)-2
                    x(p) = tableau(p+1,end);
                end
        
               % Selecting the optimal cost f
               f = tableau(end,end);

        else 
                 disp('The problem is unbounded')          
       end

    end
            % Print this onto the screen when the while condition is false
            % Which implies that we have arrived at an optimal solution
            disp (['The optimal basic feasible solution is [' num2str(x(:).') ']'])
            fprintf ('The optimal cost is %d \n', f)
            disp (['The corresponding basis is [' num2str(v_current(:).') ']'])
    
    else

        % Selecting the optimal basic feasible solution x
        for p = 1:size(tableau,1)-2
            x(p) = tableau(p+1,end);
        end

       % Selecting the optimal cost f
       f = tableau(end,end);

       % The corresponding is basis is choosen using v index since no
       % simplex iteration is performed

       disp (['The optimal basic feasible solution is [' num2str(x(:).') ']'])
       fprintf ('The optimal cost is %d \n', f)
       disp (['The corresponding basis is [' num2str(v(:).') ']'])       

    end


end

