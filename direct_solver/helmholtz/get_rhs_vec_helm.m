function [rhs_vec, rhs] = get_rhs_vec_helm(coefs,zk,phiinc)

V = coefs{1};

rhs = - zk^2*V.*phiinc;
rhs_vec = rhs(:);

end