type lotka
t0=0;
tfinal=5;
y0=[3;5];   
[t,y] = ode45(@lotka,[t0 tfinal],y0);
plot(t,y)
title('Predator/Prey Populations Over Time')
xlabel('t')
ylabel('Population')
legend('Prey','Predators','Location','North')
