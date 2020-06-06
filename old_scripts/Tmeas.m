function T=Tmeas(X) %X=mcV


%slave di TempEstimation
a0=-1.318058E2;
a1=4.830222E-2;
a2=-1.646031E-6;
a3=5.464731E-11;
a4=-9.650715E-16;
a5=8.802193E-21;
a6=-3.110810E-26;

T1=a0+a1.*X+a2*X.^2+a3*X.^3+a4*X.^4+a5*X.^5+a6*X.^6;

b0=-1.31E2;
b1=4.830E1;
b2=-1.646E0;
b3=5.46E-2;
b4=-9.65E-4;
b5=8.89E-6;
b6=-3.11E-8;
T2=b0+b1.*X+b2.*(X.^2)+b3.*(X.^3)+b4*(X.^4)+b5*(X.^5)+b6*(X.^6)+0*(X.^7)+0*X.^8;

T=size(T1);
%T=T1;
I=find(T1<1000); if ~isempty(I);  T(I)=T1(I); end
I=find(T2>=1000); if ~isempty(I); T(I)=T2(I); end 
%I=find(T<0); %offset
J=find(min(T));
T=T-T(J);
X=X-min(X);
T=1400/55*X/1000;
T=T+20;