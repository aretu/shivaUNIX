function T=Tmeas_2(X) %X=mcV
keyboard
D=X*104857593.75;
%thermo linear
temp=(10000*D)./(2^23-D);
c=(1/(1.2873851E-3 + 2.3575235E-4*log(temp) + 9.49786E-8*(log(temp))^3) - 273.85);
% Tc T-->V

c1=3.9450128025E1;
c2=2.3622373598E-2;
c3=-3.2858906784E-4;
c4=-4.9904828777E-6;
c5=-6.7509059173E-8;
c6=-5.7410327428E-10;
c7=-3.1088872894E-12;
c8=-1.0451609365E-14;
c9=-1.9889266878E-17;
c10=-1.6322697486E-20;

d1=-1.7600413686E1;
d2=3.8921204975E1;
d3=1.8558770032E-2;
d4=-9.9457592874E-5;
d5=3.1840945719E-7;
d6=-5.6072844889E-10;
d7=5.6075059059E-13;
d8=-3.2020720003E-16;
d9=9.7151147152E-20;
d10=-1.2104721275E-23;
d11=1.185976E2;
d12=-1.183432E-4;
d13=-126.9686;

% per 0 - 1372;
V=d1 + d2*c + d3*c.^2 + d4*c.^3 + d5*c.^4 + d6*c.^5 + ...
    d7*c.^6 + d8*c.^7 + d9*c.^8 + d10*c.^9 + ...
    d11*exp((d12)*(c-d13)*(c-d13));

% per -200 - 0
I0=find(c<0); 

if ~isempty(I0) %�C
V(I0)=c1*c(I0)+ c2*c(I0).^2 + c3*c(I0).^3 + c4*c(I0).^4 + c5*c(I0).^5 + c6*c(I0).^6 + ...
    c7*c(I0).^7+ c8*c(I0).^8 + c9*c(I0).^9 + c10*c(I0)^.10 ;
end
    
% Tc V-->T

I1=V>=0;
I2=V>=20644.0;
I=I1+I2;

e1=2.5173462E-2;  %0 -200  - 0
e2=-1.1662878E-6;
e3=-1.0833638E-9;
e4=-8.9773540E-13;
e5=-3.7342377E-16;
e6=-8.6632643E-20;
e7=-1.450598E-23;
e8=-5.1920577E-28;  

f1=-1.318058E+2;  %2 500 - 1372 
f2=4.830222E-2;
f3=-1.646031E-6;
f4=5.464731E-11;
f5=-9.650715E-16;
f6=8.802193E-21;
f7=-3.110810E-26;

g1=2.508355E-2;  %1 0 -500
g2=7.860106E-8;
g3=-2.503131E-10;
g4=8.315270E-14;
g5=-1.228034E-17;
g6=9.804036E-22;
g7=-4.413030E-26;
g8=1.057734E-30;
g9=-1.052755E-35;

I0=find(I==0);
if ~isempty(I0)
    
T(I0)= v(I0)*e1+ v(I0).^2*e2 + v(I0).^3*e3+ v(I0).^4*e4 + v(I0).^5*e5 + ...
    v(I0).^6*e6 + v(I0).^7*e7 + v(I0).^8*e8 + v(I0).^9*e9;
end

I0=find(I==1);
if ~isempty(I0)
   
T(I0)=v(I0)*g1+ v(I0).^2*g2 + v(I0).^3*g3+ v(I0).^4*g4 + v(I0).^5*g5 + ...
    v(I0).^6*g6 + v(I0).^7*g7 + v(I0).^8*g8 + v(I0).^9*g9 +  v(I0).^10*g10;
end

I0=find(I==2);
if ~isempty(I0)

T(I0)=f1+ v(I0)*f2 + v(I0).^2*f3 + v(I0).^3*f4+ v(I0).^4*f5 + v(I0).^5*f6 + ...
v(I0).^6*f7;

end
    

