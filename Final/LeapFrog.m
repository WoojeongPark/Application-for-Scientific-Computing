%Implicit Euler
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
x(2)=x(1)+(4*x(1)-48*x(1)*y(1))*h;
y(2)=y(1)+(-3*y(1)+39*x(1)*y(1))*h;
for j=2:max_iter
    x(j+1)=x(j-1)+(4*x(j-1)-48*x(j-1)*y(j-1))*h*2;
    y(j+1)=y(j-1)+(-3*y(j-1)+39*x(j-1)*y(j-1))*h*2;
end
t=0:h:5;
subplot(2,1,1);
xlabel('t');
ylabel('x,y');
p1=plot(t,x);
p1.LineWidth=2;
hold on
p2=plot(t,y);
p2.LineWidth=2;
legend('prey','predator')
hold off
subplot(2,1,2);
xlabel('x');
ylabel('y');
p3=plot(x,y);
legend('xlabel:prey, ylabel:pradator')
p3.LineWidth=2;