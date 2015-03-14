
%2.152 project: inverted pendulum with active suspension
%here, I switch to d(x) instead of d(t)

clear all;
close all;

%state is [x; theta; y; xdot; thetadot; ydot];

%Simulation parameters
dt = 1e-4;
tend = 10;
t = 0:dt:tend;

%Physical parameters
m1 = 1;
m2 = 1;
m3 = 1;
h2 = 1;
R = 0.5;
g = 9.8;
    params = [m1; m2; m3; h2; R; g];
    
%initialize the disturbances 
A = 0.1;
omega = 25;
d(1:length(t)) = 0;
ddot(1:length(t)) = 0;
dddot(1:length(t)) = 0;
% d(1:length(t)) = 0;
% ddot(1:length(t)) = 0;
% dddot(1:length(t)) = 0;
    D = [d;ddot;dddot];

%set the desired trajectories
v = 2.5;
    xd = v.*t;
    thetad(1:length(t)) = pi/16;
    yd(1:length(t)) = 2;
    xd_dot(1:length(t)) = v;
    thetad_dot(1:length(t)) = 0;
    yd_dot(1:length(t)) = 0;
        statesd = [xd; thetad; yd; xd_dot; thetad_dot; yd_dot];
    xd_ddot(1:length(t)) = 0;
    thetad_ddot(1:length(t)) = 0;
    yd_ddot(1:length(t)) = 0;
        statesd_accel = [xd_ddot; thetad_ddot; yd_ddot];

%Controller parameters
eta = 0.1;
lambda = 1;
phi = 0.5;
xddot_upperbound = 17;
dddot_bound = A.*omega.*xddot_upperbound + A.*omega.^2.*v.^2;
ddot_bound = A.*omega.*v;
controlparams = [eta; lambda; phi;dddot_bound;ddot_bound];
%Velocity control parameters
P = 0.01;
Pi = 0.05;
maxtheta = pi/8;

%initial conditions
x0 = 0;
theta0 = 0;
y0 = 3;
xdot0 = 0;
thetadot0 = 0;
ydot0 = 0;
    ICs = [x0; theta0; y0; xdot0; thetadot0; ydot0];

%initialize states
states = zeros(6,length(t));
    x(1:length(t)) = 0;
    theta(1:length(t)) = 0;
    y(1:length(t)) = 0;
    xdot(1:length(t)) = 0;
    thetadot(1:length(t)) = 0;
    ydot(1:length(t)) = 0;
states(:,1) = ICs;
    x(1:length(t)) = x0;
    theta(1:length(t)) = theta0;
    y(1:length(t)) = y0;
    xdot(1:length(t)) = xdot0;
    thetadot(1:length(t)) = thetadot0;
    ydot(1:length(t)) = ydot0;


%initialize h1, which is an intermediate redundant variable that is useful
h1(1:length(t)) = 0;

%initialize control input vectors
Fc(1:length(t)) = 0;
Fm(1:length(t)) = 0;
u = zeros(2,length(t));

%intialize other stuff that needs initializing
verrorI = 0;
xddot = 0;
h1(1:length(t)) = 0;

for index = 1:length(t) - 1
    %calculate intermediate variables l, h1, ldot
    l = (1./((m1 + m2).*cos(theta(index)))).*(m1.*(y(index) - d(index) - h2.*cos(theta(index))) + h2.*(m1 + m2).*cos(theta(index)));
    h1(index) = (m1 + m2).*(l - h2)./m1;
    ldot_term1 = (1./(m1 + m2)).*(sin(theta(index))./((cos(theta(index))).^2)).*(m1.*(y(index) - d(index) - h2.*cos(theta(index))) +  h2.*(m1 + m2).*cos(theta(index)));
    ldot_term2 = (1./((m1 + m2).*cos(theta(index)))).*(m1.*(ydot(index) - ddot(index) + h2.*thetadot(index).*sin(theta(index))) - h2.*(m1 + m2).*thetadot(index).*sin(theta(index)));
    ldot = ldot_term1 + ldot_term2;
    sensors = [h1(index); ldot];
    
    theaccels = derivs2(t(index),states(:,index),params,D(:,index),u(:,index+1));
    xddot = theaccels(4);
    %calculate the disturbances
    d(index+1) = A.*sin(omega.*x(index));
    ddot(index+1) = A.*omega.*xdot(index).*cos(omega.*x(index));
    dddot(index+1) = A.*omega.*xddot.*cos(omega.*x(index)) - A.*omega.^2.*(xdot(index)).^2.*sin(omega.*x(index));
    D(:,index+1) = [d(index+1); ddot(index+1); dddot(index+1)];
    
    %CONTROLLER
    %     Fc(index) = 0;
    %     Fm(index) = 0;
    %     u(:,index) = [Fc(index); Fm(index)];   
        %outer loop for xdot control (by controlling theta)
%         verrorI = verrorI + dt.*(xdot(index) - xd_dot(index)); %rectangular approximation for integral
        if index == 1
            verrorI = verrorI;
        else
            verrorI = verrorI + (dt./2).*(xdot(index) - xd_dot(index) + xdot(index-1) - xd_dot(index-1));
        end
        statesd(2,index) = -P.*(xdot(index) - xd_dot(index)) - Pi.*verrorI;
        %cap the desired angle
        if abs(statesd(2,index)) > maxtheta
            statesd(2,index) = sign(statesd(2,index)).*maxtheta;
        end
        %inner loop for angle and height control
        u(:,index+1) = u1(states(:,index),statesd(:,index),statesd_accel(:,index),params,sensors,controlparams);    
    
        
    %RK4
    k1 = dt.*derivs2(t(index),states(:,index),params,D(:,index),u(:,index+1));
    k2 = dt.*derivs2(t(index)+dt./2,states(:,index)+0.5.*k1,params,D(:,index),u(:,index+1));
    k3 = dt.*derivs2(t(index)+dt./2,states(:,index)+0.5.*k2,params,D(:,index),u(:,index+1));
    k4 = dt.*derivs2(t(index)+dt,states(:,index)+k3,params,D(:,index),u(:,index+1));
    states(:,index+1) = states(:,index) + (1/6).*(k1 + 2.*k2 + 2.*k3 + k4);
    
    x(index+1) = states(1,index+1);
    theta(index+1) = states(2,index+1);
    y(index+1) = states(3,index+1);
    xdot(index+1) = states(4,index+1);
    thetadot(index+1) = states(5,index+1);
    ydot(index+1) = states(6,index+1);
end

figure;
hold on;
plot(t,x,'k.-','LineWidth',2);
plot(t,xdot,'b.-','LineWidth',2);
legend('x','xdot');

figure;
hold on;
plot(t,theta,'k.-','LineWidth',2);
plot(t,thetadot,'b.-','LineWidth',2);
legend('theta','thetadot');

figure;
hold on;
plot(t,y,'k.-','LineWidth',2);
plot(t,ydot,'b.-','LineWidth',2);
legend('y','ydot');

figure;
hold on;
plot(t,u(1,:),'k.-','LineWidth',2);
plot(t,u(2,:),'b.-','LineWidth',2);
legend('Fc','Fm');


%animate
figure;
for tindex=1:1/(100.*dt):length(t)-1
    rot = [cos(-theta(tindex)) -sin(-theta(tindex)); sin(-theta(tindex)) cos(-theta(tindex))];
    xmin = x(tindex) - 4;
    xmax = x(tindex) + 4;
    ymin = -3; 
    ymax = 5;
    clf;
    plot(x(1:tindex),d(1:tindex)-R/2,'LineWidth',2); %draw the road
    hold on;
    axis([xmin xmax ymin ymax]);
    rectangle('Curvature',[1,1],'Position',[x(tindex)-R/2,d(tindex)-R/2,R,R],'FaceColor','k'); %draw the wheel
    plot([x(tindex) x(tindex)+h2.*sin(theta(tindex))],[d(tindex) d(tindex)+h2.*cos(theta(tindex))],'k.-','LineWidth',5); %draw the rigid link
        platecoords1 = rot*[-0.5; 0] + [x(tindex)+h2.*sin(theta(tindex)); d(tindex)+h2.*cos(theta(tindex))];
        platecoords2 = rot*[0.5; 0] + [x(tindex)+h2.*sin(theta(tindex)); d(tindex)+h2.*cos(theta(tindex))];
        platecoords = [platecoords1 platecoords2];
%     rectangle('Curvature',[1,1],'Position',[x(tindex)+h2.*sin(theta(tindex))-R/2,d(tindex)+h2.*cos(theta(tindex))-R/2,R,R],'FaceColor','k'); %draw the solid support
    plot(platecoords(1,:),platecoords(2,:),'k.-','LineWidth',6);
    plot([x(tindex)+h2.*sin(theta(tindex)) x(tindex)+(h1(tindex)+h2).*sin(theta(tindex))],[d(tindex)+h2.*cos(theta(tindex)) y(tindex)],'r.-','LineWidth',3); %draw the motor
    plot(x(1:tindex)+(h1(1:tindex)+h2).*sin(theta(1:tindex)),y(1:tindex),'g.:'); %draw the path of the load
    rectangle('Curvature',[1,1],'Position',[x(tindex)+(h1(tindex)+h2).*sin(theta(tindex))-R/2,y(tindex)-R/2,R,R],'FaceColor','c','EdgeColor','c'); %draw the load
    plot(x(tindex)+(h1(tindex)+h2).*sin(theta(tindex)),y(tindex),'rx','LineWidth',3); %draw the load center
    title(t(tindex));
    M(tindex) = getframe;
end
