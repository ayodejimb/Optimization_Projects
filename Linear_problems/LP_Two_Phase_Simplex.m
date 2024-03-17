
%%              OPTIMAL SOLUTION FOR A LINEAR PROGRAMMING PROBLEM      %%
%%                // TWO PHASE SIMPLEX METHOD  //                  %%

% This method is use to solve linear programming problem which doesn't have 
% obvious basic feasible soluion when expressed in Standard Form

% A = A matrix containing the coefficient of the variables (m x n matrix)
% b = The right hand side of the constraint (Column vector)
% c = Coefficients of the objective function (Column vector)
%% PROGRAM STARTS PART 1.3

% Function declaration
function [f,x,B] = LP_Two_Phase_Simplex(A,b,c)



% ====================  Phase I =========================================
disp (" ====================== Phase I ============================= ")

A_column = size(A,2);
n = size(A,2); % no of columns in the A matrix
m = size(A,1); % no of rows in the A matrix

id_matrix = eye(m);
id_col = size(id_matrix,2);  % No of columns in the identity matrix

% Check for identity component in each of the columns of A
    bfs_index = [];
    for j = 1:id_col

           for i=1:n
        
            if  isequal((A(:,i)), id_matrix(:,j))
                bfs_index = [bfs_index,i];    %index of BFS present in the initial A matrix
    
            else
    
                 A(:,end+1) = id_matrix(:,j);
    
            end

          end

    end

 % Remove all duplicates to obtain the A matrix with a BFS
 A = unique(A.','rows','stable').';

 % index of the BFS in the constructed A matrix with BFS
 for k = n+1:size(A,2)
     bfs_index = [bfs_index,k];
 end
% The initial cost
c_initial = [zeros(1,n), ones(1,k-n)];

% Now, the Simplex Tableau
tableau = [c_initial 0; A, b];

% Calculating the basic variables values in the initial_cost using
% their indices in the tableau matrix


% Rearrange the index of bfs based on the appearance of the diagonal
    bfs_mat = [];
for l = 1:size(bfs_index,2)
    
    bfs = tableau(2:m+1,bfs_index(1,l));

   bfs_mat = [bfs_mat bfs];    % initial bfs matrix in the Tableau

end

% Checking loop to rearrange the order in the basis based on the identity
% column of the initial BFS 

 for l = 1:size(bfs_index,2)

    if l==1 && bfs_mat(1,l) == 1   % first column with leading '1'
       
       bfs_index(1,1) = bfs_index(l);  % first in basis
      
    elseif l~=1 && bfs_mat(1,l) == 1  % column has 1 as the first element

       temp = bfs_index(1,l-1);  % take the index before
       bfs_index(1,1) = bfs_index(l);
       bfs_index(l) = temp;      % save temp to l current position

    elseif l~=1 && bfs_mat(l,l) == 1   % (2,2,3,3, etc) % column has corresponding diagonal

        bfs_index(1,l) = bfs_index(l);

    elseif l~=1 && bfs_mat(l-1,l-1) == 1   % column takes the diagonal of the preceeding one
        
        temp = bfs_index(1,l-1);
        bfs_index(1,l-1) = bfs_index(l);
        bfs_index(l) = temp; 
        
    else                                 % none of the condition listed above

        temp = bfs_index(1,end);         % store the last index in temp
        bfs_index(1,end) = bfs_index(l); % move the current index to the end
        bfs_index(l) = temp;             % bring temp to the current index
        
                       
    end
 end


for i=1:size(bfs_index,2)

    basics_var(i) = c_initial(bfs_index(i));

end

% (4) Compute the reduced cost rq

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
    B = bfs_index;
    simplex_tab = array2table(tableau);
    
    for i = 1:size(simplex_tab,2)

    simplex_tab.Properties.VariableNames(i) = strcat ('x', {num2str(i)});

    end

    disp (' === The Canonical Tableau ==== ')

    simplex_tab.Properties.VariableNames(i) = {'b'} %replacing the last column title with var name 'b'


   % (5) Condition check for Optimal BFS (no negative value in the rq? stop! Else, proceed)

    if any(rq(:,1:end-1) < 0) 
        
     while any(rq(:,1:end-1) < 0) % Do this as long as the reduced cost has a negative value

       disp ('=== The current BFS is not optimal === ')
       disp ("   === The Next Canonical Tableau === ")

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
               'respectively'],bfs_index(p_index),q_index)

       % Calculating the current corresponding Basis B
        outgoing_basics = B(p_index);
         
         % Replacing the outgoing basis with the incoming one
         B(B==outgoing_basics) = q_index;
        
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
            % disp (['The optimal basic feasible solution is [' num2str(x(:).') ']'])
            fprintf ('The optimal cost is %d \n', f)
            disp (['The original problem has a BFS corresponding to [' num2str(B(:).') ']'])
            fprintf ('\n')
    
    else

        % Selecting the optimal basic feasible solution x
        for p = 1:size(tableau,1)-2
            x(p) = tableau(p+1,end);
        end

       % Selecting the optimal cost f
       f = tableau(end,end);

       % The corresponding basis is choosen using bfs_index since no
       % simplex iteration is performed

       % disp (['The optimal basic feasible solution is [' num2str(x(:).') ']'])
       fprintf ('The optimal cost is %d \n', f)
       disp (['The original problem has a BFS corresponding to [' num2str(bfs_index(:).') ']']) 
       fprintf ('\n')

    end


% ====================  Phase II =========================================
    fprintf('\n')
    disp (" ====================== Phase II ============================= ")

% (1) Removing the Artificial Variables Column and Cost
 
  
  tableau(1,:) = [];
  tableau(end,:) = [];

  sol = tableau(:,end);

  tableau = tableau(:,1:A_column);
  tableau = [tableau,sol];

 % Adding the original cost to the tableau
  % Computing the number of zeros to be added to the cost function to form
  % the tableau
  siz_tab = size(tableau); 
  siz_c = size(c);
  siz_tab = siz_tab(:,2);
  siz_c = siz_c(:,1);

  subtrac = siz_tab - siz_c;
  subtrac = zeros(1,subtrac);

  c = [c',subtrac];
  tableau = [c; tableau];

  
  % Getting the values of the current basis in the original cost
    for i=1:size(B,2)

        basics_var(i) = c(B(i));
    
    end


  % Calculating the reduced cost
  rq = [1,size(tableau,2)];       % Modifying the size of the reduced cost to the original problem
  for j = 1:size(tableau,2)       % No of columns in the simplex tableau

    add_val = 0;

  for k = 1:size(basics_var,2)
       
      val = tableau(k+1,j) * basics_var(1,k);  

      add_val = add_val + val;   % addition loop (Zq)

  end 

    rq(j) = tableau(1,j) - add_val;  %rq = Cq - Zq

  end

  % Adding the reduced cost rq result to the tableau table
    tableau = [tableau; rq];

  % Putting into a table
    B;
    simplex_tab = array2table(tableau);
    
    for i = 1:size(simplex_tab,2)

    simplex_tab.Properties.VariableNames(i) = strcat ('x', {num2str(i)});

    end
    
    fprintf('\n')
    disp (' === The Canonical Tableau ==== ')

    simplex_tab.Properties.VariableNames(i) = {'b'} %replacing the last column title with var name 'b'



  % (2) Repeat the Simplex Algorithm to find the optimal BFS

    if any(rq(:,1:end-1) < 0) 
        
     while any(rq(:,1:end-1) < 0) % Do this as long as the reduced cost has a negative value
        
       fprintf ('\n')  
       disp ('=== The current BFS is not optimal === ')
       disp ("   === The Next Canonical Tableau === ")

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
               'respectively'],bfs_index(p_index),q_index)

       % Calculating the current corresponding Basis B
        outgoing_basics = B(p_index);
         
         % Replacing the outgoing basis with the incoming one
         B(B==outgoing_basics) = q_index;
        
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
            % B = current_basis;

            % Print this onto the screen when the while condition is false
            % Which implies that we have arrived at an optimal solution
            disp (['The optimal basic feasible solution is [' num2str(x(:).') ']'])
            fprintf ('The optimal cost is %d \n', f)
            disp (['The corresponding basis are [' num2str(B(:).') ']'])
    
    else

        % Selecting the optimal basic feasible solution x
        for p = 1:size(tableau,1)-2
            x(p) = tableau(p+1,end);
        end

       % Selecting the optimal cost f
       f = tableau(end,end);

       % The corresponding basis is choosen using bfs_index since no
       % simplex iteration is performed
        
       fprintf ('\n')
       disp (['The optimal basic feasible solution is [' num2str(x(:).') ']'])
       fprintf ('The optimal cost is %d \n', f)
       disp (['The corresponding basis are [' num2str(B(:).') ']'])       

    end

end

