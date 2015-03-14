clear all;
close all;

m1 = 5;
m2 = 3;
h1 = 45;
h2 = 43;

l = (m2*h2 + m1*(h1 + h2)) / (m1 + m2);

disp((m1 + m2)*l^2);

disp((m1^2 / (m1 + m2))*h1^2 + (m1 + m2)*h2^2 + 2*m1*h1*h2);

disp((m2*h2^2 + m1*(h1^2 + h2^2) + 2*m1*h1*h2));