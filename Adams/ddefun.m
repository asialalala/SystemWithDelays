% implementation for ddends

%  SECOND NDDE
% function yp = ddefun(t,y,ydel,ypdel) 
%     yp = y + ydel - 2*ypdel;
% end

%  FIRST NDDE
% function yp = ddefun(t,y,ydel,ypdel) 
%     yp = y + ydel - 1/4*ypdel;
% end


function yp = ddefun(t,y,ydel) 
    A = 1;
    yp = A*ydel;
end