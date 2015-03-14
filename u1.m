function res = u1(states,statesd,statesd_accel,params,sensors,controlparams)

    %h1 is assumed to be measured from sensors

    %Unpack states
    x = states(1);
    theta = states(2);
    y = states(3);
    xdot = states(4);
    thetadot = states(5);
    ydot = states(6);
    
    %Form state vector
    q = states(2:3);
    qdot = states(5:6);
    
    %Unpack desired states
    qd = statesd(2:3);
    qd_dot = statesd(5:6);
    qd_ddot = statesd_accel(2:3);
    
    %Unpack params
    m1 = params(1);
    m2 = params(2);
    m3 = params(3);
    h2 = params(4);
    R = params(5);
    g = params(6);
    
    %Unpack sensors
    h1 = sensors(1);
    ldot = sensors(2);
    
    %Unpack controlparams
    eta = controlparams(1);
    lambda = controlparams(2);
    phi = controlparams(3);
    dddot_bound = controlparams(4);
    ddot_bound = controlparams(5);
    
    %error variables
    qtilda = q - qd;
    qtildadot = qdot - qd_dot;
    s = qtildadot + lambda.*qtilda;
    
    %Calculate l using h1
    l = (m2.*h2 + m1.*(h1 + h2))./(m1 + m2);
    
    %Calculate fhat
    fhat(1,1) = (g.*sin(theta) - 2.*ldot.*thetadot)./l;
        term1 = 0;
        term2 = 2.*thetadot.^2.*(sin(theta)).^2.*y./(cos(theta)).^2 + 2.*thetadot.*ydot.*sin(theta)./cos(theta);
        term3 = ((g.*sin(theta) - 2.*ldot.*thetadot)./l).*(h1 + h2).*sin(theta);
        term4 = (((m1 + m2)./m1).*(l.*thetadot.^2 - g.*cos(theta)) - (h1 + h2).*thetadot.^2.*cos(theta)).*cos(theta);
    fhat(2,1) = term1 + term2 + term3 + term4;
    
    %Calculate B
    B(1,1) = -1./(m3.*l);
    B(1,2) = sin(theta)./(m3.*l);
    B(2,1) = -((m1 + m2).*sin(theta).*cos(theta)./(m1.*m3) + (h1 + h2).*sin(theta)./(m3.*l));
    B(2,2) = (m3 + (m1 + m2).*(sin(theta)).^2.*cos(theta))./(m1.*m3) + (h1 + h2).*(sin(theta)).^2./(m3.*l);
    
    uhat = B\(qd_ddot - fhat - lambda.*qtildadot);
    
    %Calculate F
    F(1,1) = ((dddot_bound + g).*sin(theta) - 2.*ldot.*thetadot)./l;
        termA = dddot_bound;
        termB = 2.*thetadot.^2.*(sin(theta)).^2.*(y - dddot_bound)./(cos(theta)).^2 + 2.*thetadot.*(ydot - ddot_bound).*sin(theta)./cos(theta);
        termC = ((dddot_bound + g.*sin(theta) - 2.*ldot.*thetadot)./l).*(h1 + h2).*sin(theta);
        termD = (((m1 + m2)./m1).*(l.*thetadot.^2 - (dddot_bound + g).*cos(theta)) - (h1 + h2).*thetadot.^2.*cos(theta)).*cos(theta);
    F(2,1) = termA + termB + termC + termD;
    
    K(1,1) = F(1,1) + eta;
    K(1,2) = 0;
    K(2,1) = 0;
    K(2,2) = F(2,1) + eta;
        
    u = uhat + K*sat(s./phi);
    
%     u = uhat;

    res = u;
end