function [ ps ] = root_find(pg,tol,g,WL,WR)
%ROOT_FIND Finds the root of function f_ps using Newton-Raphson

rerr = 1;
while rerr > tol
    pnew = pg - f_ps(pg,g,WL,WR)/fp_ps(pg,g,WL,WR);
    rerr = abs(f_ps(pnew,g,WL,WR) - f_ps(pg,g,WL,WR));
    pg = pnew;
end
ps = pnew;
end