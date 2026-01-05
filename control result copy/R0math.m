syms gamma delta1 delta2 delta3 p1 p2 p3 beta1 beta2 beta3 S sigma V N lambda zeta birth alpha
dvdx = [gamma+zeta 0 0 0; 
    -gamma (delta1+p1+zeta) 0 0;
    0 -p1 (delta2+p2+zeta) 0;
    0 0 -p2 (delta3+p3+zeta)];
inv_dvdx = inv(dvdx)
dfdx = [0 beta1*(S+sigma*V) beta2*(S+sigma*V) beta3*(S+sigma*V); 
    0 0 0 0;
    0 0 0 0;
    0 0 0 0];
m = dfdx*inv(dvdx);
n = det(m)

M1 = (beta1*gamma*(S + V*sigma))/((gamma + zeta)*(delta1 + p1 + zeta)) + (beta2*gamma*p1*(S + V*sigma))/((gamma + zeta)*(delta1 + p1 + zeta)*(delta2 + p2 + zeta)) + (beta3*gamma*p1*p2*(S + V*sigma))/((gamma + zeta)*(delta1 + p1 + zeta)*(delta2 + p2 + zeta)*(delta3 + p3 + zeta));
M2 = (beta1*(S + V*sigma))/(delta1 + p1 + zeta) + (beta2*p1*(S + V*sigma))/((delta1 + p1 + zeta)*(delta2 + p2 + zeta)) + (beta3*p1*p2*(S + V*sigma))/((delta1 + p1 + zeta)*(delta2 + p2 + zeta)*(delta3 + p3 + zeta)); 
M3 = (beta2*(S + V*sigma))/(delta2 + p2 + zeta) + (beta3*p2*(S + V*sigma))/((delta2 + p2 + zeta)*(delta3 + p3 + zeta));
M4 = (beta3*(S + V*sigma))/(delta3 + p3 + zeta)


S = birth/(alpha+zeta)
V = alpha*birth/(zeta*(alpha+zeta))

N1 = (beta1*gamma*(S + V*sigma))/((gamma + zeta)*(delta1 + p1 + zeta)) + (beta2*gamma*p1*(S + V*sigma))/((gamma + zeta)*(delta1 + p1 + zeta)*(delta2 + p2 + zeta)) + (beta3*gamma*p1*p2*(S + V*sigma))/((gamma + zeta)*(delta1 + p1 + zeta)*(delta2 + p2 + zeta)*(delta3 + p3 + zeta));
N2 = (beta1*(S + V*sigma))/(delta1 + p1 + zeta) + (beta2*p1*(S + V*sigma))/((delta1 + p1 + zeta)*(delta2 + p2 + zeta)) + (beta3*p1*p2*(S + V*sigma))/((delta1 + p1 + zeta)*(delta2 + p2 + zeta)*(delta3 + p3 + zeta)); 
N3 = (beta2*(S + V*sigma))/(delta2 + p2 + zeta) + (beta3*p2*(S + V*sigma))/((delta2 + p2 + zeta)*(delta3 + p3 + zeta));
N4 = (beta3*(S + V*sigma))/(delta3 + p3 + zeta)



N = [N1 N2 N3 N4; 
     0 0 0 0;
     0 0 0 0;
     0 0 0 0];
eig(N)

M = [M1 M2 M3 M4; 
     0 0 0 0;
     0 0 0 0;
     0 0 0 0];

J = [-alpha-zeta 0 0 0 0 0 0 0; 
     alpha -zeta 0 0 0 0 0 0;
     0 0 -gamma-zeta 0 0 0 0 0;
     0 0 gamma -delta1-p1-zeta 0 0 0 0;
     0 0 0 p1 -delta2-p2-zeta 0 0 0;
     0 0 0 0 p2 -delta3-p3-zeta 0 0;
     0 0 0 delta1 delta2 delta3 -zeta 0;
     0 0 0 0 0 p3 0 0];

j = eig(J)

JJ = [-alpha-zeta-lambda 0 0 0 0 0 0 0; 
     alpha -zeta-lambda 0 0 0 0 0 0;
     0 0 -gamma-zeta-lambda 0 0 0 0 0;
     0 0 gamma -delta1-p1-zeta-lambda 0 0 0 0;
     0 0 0 p1 -delta2-p2-zeta-lambda 0 0 0;
     0 0 0 0 p2 -delta3-p3-zeta-lambda 0 0;
     0 0 0 delta1 delta2 delta3 -zeta-lambda 0;
     0 0 0 0 0 p3 0 0-lambda];

det(JJ)