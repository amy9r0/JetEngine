function [root,e]=FalsePosition(p0,p1,err,f,it)
i=      2;
q0=     f(p0);
q1=     f(p1);
root=   NaN;
e=      NaN;
if q0*q1>0
    warning('The interval given does not contain a root.')
    return
end
while i<=it
    p=  p1-q1*(p1-p0)/(q1-q0);
    if abs(p-p1)<err
        e=      abs(p-p1);
        root=   p;
        return
    end
    i=  i+1;
    q=  f(p);
    if q*q1<0
        p0= p1;
        q0= q1;
    end
    p1= p;
    q1= q;
end
if i==it
    warning('The number of iterations reached maximum before a root was determined.')
end