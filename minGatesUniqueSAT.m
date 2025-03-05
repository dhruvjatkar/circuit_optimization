function [min_gates, min_W, min_A, W, A] = minGatesUniqueSAT(n_, F, rmax, order, solver)%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%
%close all
%clear all



%order = 'MSBfirst';

%n = 4;

m = 0:2^n_-1;
U = dec2bin(m', n_) == '1'; %rows are different input vectors

if strcmp(order, 'MSBfirst')
    U = U;
elseif strcmp(order, 'LSBfirst')
    U = fliplr(U);
else
    error('unknown option for order')
end

U = [U, zeros(2^n_, 1)];

n = n_ + 1;

%F = [0, 0, 0, 0, 0, 1, 0, 0];

%rng(192834);

%F = randi(2, 1, 2^n)-1;

for r = 1:rmax
    fprintf('Trying %d gates...', r)
    sat = satProblem();


    sat.addDecisionVariable('x', [zeros(n, 2^n_); ones(r, 2^n_)]);

    s_mask = zeros(n+r,n+r, n+r);

    for i = n+1:n+r
        for k = 1:i-1
            for j = 1:k-1
                s_mask(i,j,k) = 1;
            end
        end
    end


    sat.addDecisionVariable('s', s_mask);
    %add clauses for each NOR gate
    for i = n+1:n+r
        for k = 1:i-1
            for j = 1:k-1
                for t = 1:2^n_
                    %
                    if (j > n) && (k > n)
                        %both inputs are coming from NOR gates
                        sat.addClause([-sat.var('s', i, j, k),  sat.var('x', j, t),  sat.var('x', k, t),  sat.var('x', i, t)]);
                        sat.addClause([-sat.var('s', i, j, k),  sat.var('x', j, t), -sat.var('x', k, t), -sat.var('x', i, t)]);
                        sat.addClause([-sat.var('s', i, j, k), -sat.var('x', j, t),  sat.var('x', k, t), -sat.var('x', i, t)]);
                        sat.addClause([-sat.var('s', i, j, k), -sat.var('x', j, t), -sat.var('x', k, t), -sat.var('x', i, t)]);
                    elseif (j > n) && (k <= n)
                        %j from gate, k from input
                        if U(t, k)
                            %x_kt is true
                            %sat.addClause([-sat.var('s', i, j, k),  sat.var('x', j, t),  sat.var('x', k, t),  sat.var('x', i, t)]);
                            sat.addClause([-sat.var('s', i, j, k),  sat.var('x', j, t),                      -sat.var('x', i, t)]);
                            %sat.addClause([-sat.var('s', i, j, k), -sat.var('x', j, t),  sat.var('x', k, t), -sat.var('x', i, t)]);
                            sat.addClause([-sat.var('s', i, j, k), -sat.var('x', j, t),                      -sat.var('x', i, t)]);
                        else
                            %x_kt is false
                            sat.addClause([-sat.var('s', i, j, k),  sat.var('x', j, t),                       sat.var('x', i, t)]);
                            %sat.addClause([-sat.var('s', i, j, k),  sat.var('x', j, t), -sat.var('x', k, t), -sat.var('x', i, t)]);
                            sat.addClause([-sat.var('s', i, j, k), -sat.var('x', j, t),                      -sat.var('x', i, t)]);
                            %sat.addClause([-sat.var('s', i, j, k), -sat.var('x', j, t), -sat.var('x', k, t), -sat.var('x', i, t)]);
                        end
                    elseif (j <= n) && (k > n)
                        %j from input, k from gate
                        if U(t, j)
                            %x_jt is true
                            %sat.addClause([-sat.var('s', i, j, k),  sat.var('x', j, t),  sat.var('x', k, t),  sat.var('x', i, t)]);
                            %sat.addClause([-sat.var('s', i, j, k),  sat.var('x', j, t), -sat.var('x', k, t), -sat.var('x', i, t)]);
                            sat.addClause([-sat.var('s', i, j, k),                        sat.var('x', k, t), -sat.var('x', i, t)]);
                            sat.addClause([-sat.var('s', i, j, k),                       -sat.var('x', k, t), -sat.var('x', i, t)]);
                        else
                            %x_jt is false
                            sat.addClause([-sat.var('s', i, j, k),                        sat.var('x', k, t),  sat.var('x', i, t)]);
                            sat.addClause([-sat.var('s', i, j, k),                       -sat.var('x', k, t), -sat.var('x', i, t)]);
                            %sat.addClause([-sat.var('s', i, j, k), -sat.var('x', j, t),  sat.var('x', k, t), -sat.var('x', i, t)]);
                            %sat.addClause([-sat.var('s', i, j, k), -sat.var('x', j, t), -sat.var('x', k, t), -sat.var('x', i, t)]);
                        end
                    else
                        %j from input, k from input
                        if U(t, j) && U(t, k)
                            %x_jt is true, x_kt is true
                            %sat.addClause([-sat.var('s', i, j, k),  sat.var('x', j, t),  sat.var('x', k, t),  sat.var('x', i, t)]);
                            %sat.addClause([-sat.var('s', i, j, k),  sat.var('x', j, t), -sat.var('x', k, t), -sat.var('x', i, t)]);
                            %sat.addClause([-sat.var('s', i, j, k), -sat.var('x', j, t),  sat.var('x', k, t), -sat.var('x', i, t)]);
                            sat.addClause([-sat.var('s', i, j, k),                                            -sat.var('x', i, t)]);
                        elseif U(t, j) && ~U(t, k)
                            %x_jt is true, x_kt is false
                            %sat.addClause([-sat.var('s', i, j, k),  sat.var('x', j, t),  sat.var('x', k, t),  sat.var('x', i, t)]);
                            %sat.addClause([-sat.var('s', i, j, k),  sat.var('x', j, t), -sat.var('x', k, t), -sat.var('x', i, t)]);
                            sat.addClause([-sat.var('s', i, j, k),                                            -sat.var('x', i, t)]);
                            %sat.addClause([-sat.var('s', i, j, k), -sat.var('x', j, t), -sat.var('x', k, t), -sat.var('x', i, t)]);
                        elseif ~U(t, j) && U(t, k)
                            %x_jt is false, x_kt is true
                            %sat.addClause([-sat.var('s', i, j, k),  sat.var('x', j, t),  sat.var('x', k, t),  sat.var('x', i, t)]);
                            sat.addClause([-sat.var('s', i, j, k),                                            -sat.var('x', i, t)]);
                            %sat.addClause([-sat.var('s', i, j, k), -sat.var('x', j, t),  sat.var('x', k, t), -sat.var('x', i, t)]);
                            %sat.addClause([-sat.var('s', i, j, k), -sat.var('x', j, t), -sat.var('x', k, t), -sat.var('x', i, t)]);
                        else
                            %x_jt is false, x_kt is false
                            sat.addClause([-sat.var('s', i, j, k),                                             sat.var('x', i, t)]);
                            %sat.addClause([-sat.var('s', i, j, k),  sat.var('x', j, t), -sat.var('x', k, t), -sat.var('x', i, t)]);
                            %sat.addClause([-sat.var('s', i, j, k), -sat.var('x', j, t),  sat.var('x', k, t), -sat.var('x', i, t)]);
                            %sat.addClause([-sat.var('s', i, j, k), -sat.var('x', j, t), -sat.var('x', k, t), -sat.var('x', i, t)]);
                        end
                    end
                end
            end
        end
    end

    %clauses to make each gate have two inputs
    for i = n+1:n+r
        C = [];
        for k = 1:i-1
            for j = 1:k-1
                C(end+1) = sat.var('s', i, j, k);
            end
        end
        sat.addClause(C);
    end

    %clauses to make us select exactly no more than one s_ijk for each i
    for i = n+1:n+r
        for k = 1:i-1
            for j = 1:k-1
                %
                for k_ = 1:i-1
                    for j_ = 1:k_-1
                        if (j_ ~= j) || (k_ ~= k)
                            sat.addClause([-sat.var('s', i, j, k), -sat.var('s', i, j_, k_)]);
                        end
                    end
                end
            end
        end
    end


    %clauses to make output match F
    for t = 1:2^n_
        if F(t)
            %output is true
            sat.addClause(sat.var('x', n+r, t));
        else
            %output is false
            sat.addClause(-sat.var('x', n+r, t));
        end
    end

    %clauses to make representation canonical
    for i = n+1:n+r-1
        %
        for k = 1:i-1
            for k_ = 1:k-1
                for j = 1:k-1
                    for j_ = 1:k_-1
                        if (sat.var('s', i, j, k) > 0) && (sat.var('s', i+1, j_, k_) > 0)
                            sat.addClause([-sat.var('s', i, j, k), -sat.var('s', i+1, j_, k_)]);
                        end
                    end
                end
            end
            k_ = k;
            for j = 1:k-1
                % if i == 5 && k == 4 && j == 2
                %      pause
                % end
                for j_ = 1:j
                    if (sat.var('s', i, j, k) > 0) && (sat.var('s', i+1, j_, k_) > 0)
                        sat.addClause([-sat.var('s', i, j, k), -sat.var('s', i+1, j_, k_)]);
                    end
                end
            end
        end

    end

    s = sat.isSat(solver);

    if s.sat == 1
        %found circuit with r gates
        min_gates = r;
        fprintf('done. Feasible\n')
        break;
    elseif r == rmax
        %uh oh
        min_gates = -1;
        fprintf('done. Not feasible\n')
        fprintf('went up to rmax and did not find any circuits\n')
        return
    else
        fprintf('done. Not feasible\n')
    end
end


%clauses to make representation canonical
% for i = n+1:n+r-1
%     %
%     for k = 1:i-1
%         for k_ = 1:k-1
%             for j = 1:k-1
%                 for j_ = 1:k_-1
%                     if (sat.var('s', i, j, k) > 0) && (sat.var('s', i+1, j_, k_) > 0)
%                         sat.addClause([-sat.var('s', i, j, k), -sat.var('s', i+1, j_, k_)]);
%                     end
%                 end
%             end
%         end
%         k_ = k;
%         for j = 1:k-1
%             % if i == 5 && k == 4 && j == 2
%             %      pause
%             % end
%             for j_ = 1:j
%                 if (sat.var('s', i, j, k) > 0) && (sat.var('s', i+1, j_, k_) > 0)
%                     sat.addClause([-sat.var('s', i, j, k), -sat.var('s', i+1, j_, k_)]);
%                 end
%             end
%         end
%     end
% 
% end

%s = sat.isSat(solver);

[min_W, min_A] = sat2dag(s, n, r);

W = {min_W};
A = {min_A};
fprintf('Finding additional circuits...')
Ncircuits = 1;
while s.sat == 1
    sat.addClause(-s.assignment);
    s = sat.isSat(solver);
    if s.sat()
        Ncircuits = Ncircuits + 1;
        fprintf('%d, ', Ncircuits)
        [W{end+1}, A{end+1}] = sat2dag(s, n, r);
    else
        fprintf('\n')
        break
    end
end

end

function [W, A] = sat2dag(s, n, r)
S = s.values('s');
S = S{1};

WA = zeros(n+r, r);
%A = zeros(r, r);
for i = n+1:n+r
    [~, J, K] = ind2sub(size(S(i,:,:)), find(S(i,:,:) == 1));

    if ~isscalar(J)
        error('J not scalar')
    end
    if ~isscalar(K)
        error('K not scalar')
    end
    WA(J, i) = 1;
    WA(K, i) = 1;
end


W = WA(1:n, n+1:end);
A = WA(n+1:end, n+1:end);
end