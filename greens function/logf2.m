function val = logf2(x)
    val = log(x);
    val(x < 0) = val(x < 0) - 2i*pi; 
end