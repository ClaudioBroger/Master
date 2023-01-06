function auxs=fiteoexp(vpl)

%w = [0.0012 0.0005 0.0007 0.0018 0.0014 0.0018 0.0023 0.0013 0.0013 0.0014 0.0002 0.0002 0.0002 0.0002];
%tp = [0  478.2000  556.6000  608.0000  648.2000  680.0000  706.6000  725.0000  734.9400  756.6000 771.2000  795.0000  825.0000  863.8000];
%exp= [0.0073 0.0069 0.0072 0.0065 0.0058 0.0058 0.0046 0.0018 0.0016 0.0013 0.0382e-3 0.0757e-3 0.0441e-3 0.2482e-3];

tp = [740.2000  747.2000  747.7000  748.2000  748.7000  749.2000  749.7000  750.2000  750.7000 ...
    751.2000  751.7000  752.2000  752.7000  753.2000  754.2000  755.2000  757.2000  760.2000 ...
    765.2000  780.2000  795.2000  810.2000  840.2000];

w = (1e-3)*[ 0.5950    0.9304    0.5023    0.3522    0.4558    0.7142    0.3529    0.8285    0.7643 ...
      0.1523    0.6873    0.1702    0.4453    0.3324    0.2457    0.4229    0.5552    0.4805 ...
      0.2170    0.4954    0.3861    0.4209    0.4057];
  
expo = [  0.0073    0.0063    0.0061    0.0060    0.0065    0.0073    0.0063    0.0062    0.0059 ...
    0.0062    0.0058    0.0054    0.0055    0.0053    0.0052    0.0046    0.0046    0.0039 ...
    0.0030    0.0015    0.0004    0.0002    0.0002];



vp=10.^vpl;


eval=zeros(1,length(tp));

sum=0;

for i=1:length(tp)

    
    eval(i) = vp(1)*exp(-vp(2)*tp(i));
    sum=sum+0.5*((eval(i)-expo(i))/w(i))^2;

end

auxs=sum;


return