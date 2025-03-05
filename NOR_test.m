%
close all
clear all

r = 12;

order = 'MSBfirst';

n = 4;

m = 0:2^n-1;
U = dec2bin(m', n) == '1'; %rows are different input vectors

if strcmp(order, 'MSBfirst')
    U = U;
elseif strcmp(order, 'LSBfirst')
    U = fliplr(U);
else
    error('unknown option for order')
end

%F = [0, 0, 0, 0, 0, 1, 0, 0];

rng(192834);

F = randi(2, 1, 2^n)-1;

sat = satProblem();


sat.addDecisionVariable('x', [zeros(n, 2^n); ones(r, 2^n)]);

s_mask = zeros(r, r, r);

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
            for t = 1:2^n
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
                    elseif U(t, j) && U(t, k)
                        %x_jt is true, x_kt is false
                        %sat.addClause([-sat.var('s', i, j, k),  sat.var('x', j, t),  sat.var('x', k, t),  sat.var('x', i, t)]);
                        %sat.addClause([-sat.var('s', i, j, k),  sat.var('x', j, t), -sat.var('x', k, t), -sat.var('x', i, t)]);
                        sat.addClause([-sat.var('s', i, j, k),                                            -sat.var('x', i, t)]);
                        %sat.addClause([-sat.var('s', i, j, k), -sat.var('x', j, t), -sat.var('x', k, t), -sat.var('x', i, t)]);
                    elseif U(t, j) && U(t, k)
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


%clauses to make output match F
for t = 1:2^n
    if F(t)
        %output is true
        sat.addClause(sat.var('x', r, t));
    else
        %output is false
        sat.addClause(-sat.var('x', r, t));
    end
end

%clauses to make representation canonical
for i = n+1:n+r-1
    %
    for k = 1:i-1
        for k_ = 1:k-1
            for j = 1:k-1
                for j_ = 1:k_-1
                    sat.addClause([-sat.var('s', i, j, k), -sat.var('s', i+1, j_, k_)]);
                end
            end
        end
        k_ = k;
        for j = 1:k-1
            for j_ = 1:j
                sat.addClause([-sat.var('s', i, j, k), -sat.var('s', i+1, j_, k_)]);
            end
        end
    end

end


s = sat.isSat()





