function res = derivs1(t,states,params,D,u)
%Unpack states
x = states(1);
theta = states(2);
y = states(3);
xdot = states(4);
thetadot = states(5);
ydot = states(6);

%Unpack params
m1 = params(1);
m2 = params(2);
m3 = params(3);
h2 = params(4);
R = params(5);
g = params(6);

%Unpack D
d = D(1);
ddot = D(2);
dddot = D(3);

%Unpack u
Fc = u(1);
Fm = u(2);

%Calculate l and ldot
l = (1./((m1 + m2).*cos(theta))).*(m1.*(y - d - h2.*cos(theta)) + h2.*(m1 + m2).*cos(theta));
ldot_term1 = (1./(m1 + m2)).*(sin(theta)./((cos(theta)).^2)).*(m1.*(y - d - h2.*cos(theta)) +  h2.*(m1 + m2).*cos(theta));
ldot_term2 = (1./((m1 + m2).*cos(theta))).*(m1.*(ydot - ddot + h2.*thetadot.*sin(theta)) - h2.*(m1 + m2).*thetadot.*sin(theta));
ldot = ldot_term1 + ldot_term2;

xddot = (1./m3).*Fc - (sin(theta)./m3).*Fm;
thetaddot = ((dddot + g).*sin(theta) - 2.*ldot.*thetadot)./l - (1./(m3.*l)).*Fc + (sin(theta)./(m3.*l)).*Fm;

%long stuff for yddot
term1 = dddot;
term2 = -(h2 - ((m1 + m2)./m1).*(l + h2)).*(((dddot + g).*sin(theta) - 2.*ldot.*thetadot)./l).*sin(theta);
term3 = (h2.*(1 + (m1 + m2)./m1)).*thetadot.^2.*cos(theta);
term4 = -((m1 + m2)./m1).*(dddot + g).*(cos(theta)).^2;
term5 = -(2.*thetadot.*(m1 + m2)./m1).*ldot.*sin(theta);
coeffFc = (sin(theta)./(m3.*l)).*(h2 - ((m1 + m2)./m1).*(l + h2)) - ((m1 + m2)./(m1.*m3)).*sin(theta).*cos(theta);
coeffFm = ((m3 + (m1 + m2).*(sin(theta)).^2.*cos(theta))./(m1.*m3)) - ((sin(theta)).^2./(m3.*l)).*(h2 - ((m1 + m2)./m1).*(l + h2));

yddot = term1 +  term2 + term3 + term4 + term5 + coeffFc.*Fc + coeffFm.*Fm;


res = [xdot;thetadot;ydot;xddot;thetaddot;yddot];


end