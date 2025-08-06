function numberedgrid = get_grid(X,Y)

    numbers = 1:length(X(:));
    numberedgrid = reshape(numbers, size(X));

end
