%projanimate.m

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