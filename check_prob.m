A = zeros(10000,80);
r = rand(size(A));
p = exp(log(0.7)/N);
for ii = 1:size(A,1)
    for jj = 1:80
        if r(ii,jj) <  p
            A(ii,jj) = 1;
        end
    end
end

sum((sum(A,2)/80) == 1)/size(A,1) - 0.7