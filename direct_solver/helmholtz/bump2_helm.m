function coefs = bump2_helm(X,Y,amp,width)

    coefs = cell(1);
    coefs{1} = amp*exp(-(X.^2 + Y.^2)/(2*width^2));

end