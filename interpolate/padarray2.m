function z = padarray2(z,N)
    [m,n] = size(z);
    z = [zeros(ceil(N)) zeros(ceil(N),n) zeros(ceil(N),floor(N));
        zeros(m,ceil(N)) z zeros(m,floor(N));
        zeros(floor(N),ceil(N)) zeros(floor(N),n) zeros(floor(N))];
end