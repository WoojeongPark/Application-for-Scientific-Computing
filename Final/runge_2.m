%Runge Kutta 2nd order
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
    k1x=h*(4*x(j)-48*x(j)*y(j));k1y=h*(-3*y(j)+39*y(j)*x(j));
    k2x=h*(4*(x(j)+k1x/2)-48*(x(j)+k1x/2)*y(j));k2y=h*(-3*(y(j)+k1y/2)+39*(y(j)+k1y/2)*x(j));
    x(j+1)=x(j)+k2x;
    y(j+1)=y(j)+k2y;
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

