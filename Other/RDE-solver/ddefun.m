function dydt = ddefun(t,y,Z)
ydelay = Z(:,1);
dydt =  - y + ydelay(1) + t;
end