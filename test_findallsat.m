%
close all
clear all

%r = 12;

order = 'MSBfirst';

n = 4;

%F = [0, 0, 0, 0, 0, 1, 0, 0];

%rng(192834);
%rng(98734);
%rng(9324);
F = randi(2, 1, 2^n)-1;

[min_gates, min_W, min_A, W, A] = minGatesUniqueSAT(n, F, 15, order, 'glucose');


%%

Ncircuits = length(W);

for i = 1:Ncircuits
    if any(sum(W{i}, 1) + sum(A{i}, 1) > 2)
        error('more than two inputs to a gate')
    end
    if ~isCanon(W{i}, A{i})
        error('found non canonical representation')
    end
end



function c = isCanon(W, A)

WA = [W;A];

r = size(A, 1);
n = size(W, 1);

v = 2.^[0:(n+r-1)];

for i = 1:size(W, 2)-1
    if v*WA(:,i) >= v*WA(:,i+1)
        c = 0;
        return
    end
end
c = 1;


end