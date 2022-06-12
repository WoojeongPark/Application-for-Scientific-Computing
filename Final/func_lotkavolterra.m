function yp = lotka(t,y)
yp = diag([4-48*y(2), -3+39*y(1)])*y;
end