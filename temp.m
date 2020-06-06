tI=1;
%tI=find(time>=0,1,'first');
%t=time(tI:dn:end)/1000;
dt=stamp(tI:dn:end)/1000;
%dt=diff(t); dt(end+1)=dt(end); dt=abs(dt);
s=shear(tI:dn:end).*1E6;
v=vel(tI:dn:end);

Gamma = (Rho*c*sqrt(Kappa*pi)); 

% c=950; Rho=2900; Kappa=0.48E-6; %Kappa=3.44E-007; %gabbro 
% Gamma = Rho*c*sqrt(Kappa*pi); 

% %cp=(50*1167/100 + 40*1164/100 + 10*913/100); ro=(50*2765/100 + 40*3279/100 + 10*5150/100); key= 0.72E-006; %key=1.27E-006; %carbonates
% %cp = 880; ro = 2700; key = 1.48E-6; %new basalt 5/11/2012
% %cp=(cal/kg C)--> 1J=0.239*cal; %key= (m^2/s) %rho=kg/m^3

M=length(s);
Temp(1:2)=0;


for m=3:M
    aa=0;
for n=2:m-1
    bb = (v(n) .* s(n)) * dt(n) / sqrt(dt(n)*(m-n)) /Gamma ;
%    bb = (v(n) .* s(n)) * dt(n)  ;

    aa=aa+bb;
end
    Temp(m) = aa;
    %H(m) = H(1)*exp(-T(m)/200); 
    %alpha(m)=min(sn./H + 1*s(m),1); 
    %Tau(m) = alpha(m)*y0*exp(-T(m)/Tc);
    

end

cw=4186;
Rhow=999.97;
for i=1:length(Temp)
    Tf(i)=25 + c*Rho*Temp(i)/(cw*Rhow);
end
    

figure(101);

subplot(1,2,1); 
plot(t,s.*v,'r')
xlabel('Time (s)')
ylabel('Shear Stress (Pa) * Slip Rate (m/s)')

subplot(1,2,2);
plot(t,Temp,'--k')
xlabel('Time #')
ylabel('T (^oC)')


