function theta = calcNormal(p1,p2)

x1 = p1(2);
y1 = p1(1);
z1 = 0;

x2 = p2(2);
y2 = p2(1);
z2 = 0;

x3 = 0;
y3 = 0;
z3 = 1;

v = [x2-x1, y2-y1];
vec = [v(2), -v(1)];

theta = atan(abs(vec(2))/abs(vec(1)));

if sign(vec(1)) == 1 && sign(vec(2)) == 1
    theta = 2*pi-theta;
elseif sign(vec(1)) == -1 && sign(vec(2)) == 1
    theta = pi + theta;
elseif sign(vec(1)) == -1 && sign(vec(2)) == -1
    theta = pi - theta;
elseif sign(vec(1)) == 1 && sign(vec(2)) == -1
    theta = theta;
elseif sign(vec(1)) == -1 && vec(2) == 0
    theta = pi;
elseif vec(1) == 0 && sign(vec(2)) == -1
    theta = pi/2;
elseif vec(1) == 0 && sign(vec(2)) == 1
    theta = (3*pi)/2;
end




end