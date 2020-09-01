%All the length unit here is nm
clear all
wvl = linspace(450,950,1000) % wavelength 500nm - 800nm 



%% Set up layer properties layer order top-down air-bn-sample-bn-graphite-SiO2-Si
d1=30;  %top bn thickness 
d3=30;  %bot bn thickness 
f1=0.3;
f2=0.2
g1=0.015;
g2=0.02;
eng = 1243./wvl
e2 = 5+f1./(-eng.^2+1.68^2-1i*g1*eng)+f2./(-eng.^2+2^2-1i*g2*eng); %sample layer dielectric function, assume a peak at 740
n0 = 1; % air refractive indewvl
n1 = 1.8; %top BN refractive indewvl
n2 = e2.^(1/2); %sample layer refractive indewvl
n3 = 1.8; %bottom bn refractive indewvl
e1=2.5; 
e3=1.66;
e4=e1+1i*e3*wvl/500;  %graphite dielectric function (rough)
n4=e4.^(1/2); %graphite refractive indewvl
d4 = 2; % graphite thickness
n5 = 1.46; %SiO2 refractive indewvl
n6 = sqrt(15.590+0.21635*1i); % Si refractive indewvl
d5 = 90;  %SiO2 thickness
d2 = 1.3; %sample thickness



%% Simulate reflection
newDefaultColors = jet(20);
figure;
title(sprintf('RC for different BN thickness, graphite %d nm',d4))
set(gca,'Ydir','reverse')
xlabel('Wavelength (nm)')
ylabel('-RC')
set(gca, 'ColorOrder', newDefaultColors, 'NextPlot', 'replacechildren');

%%
for d1 = 30:10:60
    for d3 = 30:10:60
        r01 = (n0-n1)/(n0+n1);
        r10 = (n1-n0)/(n0+n1);
        t01 = 2*n0/(n0+n1);
        t10 = 2*n1/(n0+n1);
        r12 = (n1-n2)./(n1+n2);
        r21 = (n2-n1)./(n1+n2);
        r23 = (n2-n3)./(n2+n3);
        r32 = (n3-n2)./(n2+n3);
        t12 = 2*n1./(n1+n2);
        t21 = 2*n2./(n1+n2);
        t23 = 2*n2./(n2+n3);
        t32 = 2*n3./(n2+n3);
        r13 = (n1-n3)/(n1+n3);
        r31 = (n3-n1)/(n1+n3);
        t13 = 2*n1/(n1+n3);
        t31 = 2*n3/(n1+n3);
        r34 = (n3-n4)./(n3+n4);
        r43 = (n4-n3)./(n3+n4);
        t34 = 2*n3./(n3+n4);
        t43 = 2*n4./(n3+n4);
        t45 = 2*n4./(n4+n5);
        t54 = 2*n5./(n4+n5);
        r45 = (n4-n5)./(n4+n5);
        r54 = (n5-n4)./(n4+n5);
        r56 = (n5-n6)./(n5+n6);
        t56=2*n5./(n5+n6);
        t65=2*n6./(n5+n6);
        rrrr = r45+t45.*r56.*t54.*exp(2*2*pi.*n5*d5./wvl*1i)./(1-r54.*r56.*exp(2*2*pi.*n5*d5./wvl*1i));
        rrr = r34+t34.*rrrr.*t43.*exp(2*2*pi*n4*d4./wvl*1i)./(1-r43.*rrrr.*exp(2*2*pi*n4*d4./wvl*1i));
        rr = r23+t23.*rrr.*t32.*exp(2*2*pi*n3*d3./wvl*1i)./(1-r32.*rrr.*exp(2*2*pi*n3*d3./wvl*1i));
        r = r12+t12.*rr.*t21.*exp(2*2*pi*n2*d2./wvl*1i)./(1-r21.*rr.*exp(2*2*pi*n2*d2./wvl*1i));
        r2= r01+t01.*r.*t10.*exp(2*2*pi*n1*d1./wvl*1i)./(1-r10*r.*exp(2*2*pi*n1*d1./wvl*1i));
        %ra= t01*(1+r).*exp(2*pi*n1*d1./wvl*1i)./(1-r10*r.*exp(2*2*pi*n1*d1./wvl*1i));
        r = r13+t13*rrr.*t31.*exp(2*2*pi*n3*d3./wvl*1i)./(1-r31*rrr.*exp(2*2*pi*n3*d3./wvl*1i));
        r1= r01+t01.*r*t10.*exp(2*2*pi*n1*d1./wvl*1i)./(1-r10*r.*exp(2*2*pi*n1*d1./wvl*1i));
        R=abs(r2).^2;
        y=R./abs(r1).^2-1;
        hold all
        plot(wvl,y,'DisplayName', sprintf('top BN = %d nm, bottom BN = %d nm',d1,d3),'LineWidth', 3);
    end
end
