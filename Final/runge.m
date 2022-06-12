%Runge Kutta 4th order
h=input("h : 0.001 or 0.0005 or 0.00025 :: ");
if(h==0.001)
    max_iter=5000;
elseif(h==0.0005)
    max_iter=10000;
else
    max_iter=20000;
end
x=zeros(max_iter+1,1);y=zeros(max_iter+1,1);
x(1)=3;y(1)=5;
for j=1:max_iter;
    k1x=4*x(j)-48*x(j)*y(j);k1y=-3*y(j)+39*x(j)*y(j);
    k2x=4*(x(j)+k1x*h/2)-48*(x(j)+k1x*h/2)*y(j);k2y=-3*(y(j)+k1y*h/2)+39*x(j)*(y(j)+k1y*h/2);
    k3x=4*(x(j)+k2x*h/2)-48*(x(j)+k2x*h/2)*y(j);k3y=-3*(y(j)+k2y*h/2)+39*x(j)*(y(j)+k2y*h/2);
    k4x=4*(x(j)+k3x*h)-48*(x(j)+k3x*h)*y(j);k4y=-3*(y(j)+k3y*h)+39*x(j)*(y(j)+k3y*h);
    x(j+1)=x(j)+h*(k1x+k2x+k3x+k4x)/6;
    y(j+1)=y(j)+h*(k1y+k2y+k3y+k4y)/6;
end
t=0:h:5;
subplot(1,2,1);
xlabel('t');
ylabel('x,y');
p1=plot(t,x);
p1.LineWidth=2;
hold on
p2=plot(t,y);
p2.LineWidth=2;
legend('prey','predator')
hold off
subplot(1,2,2);
xlabel('x');
ylabel('y');
p3=plot(x,y);
legend('xlabel:prey, ylabel:pradator')
p3.LineWidth=2;
