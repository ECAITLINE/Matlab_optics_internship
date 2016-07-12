# Matlab_optics_internship
# The newest version
clc; clear all;
close all;

lambda = 633e-9; % wavelength of illlumination
n_au = 0.197 + 1i*3.0908; % refractive index of gold at wavelength of 633 nm
epsilon_au = n_au^2; % relative permittivity of gold
n_air = 1; % refractive index of air
epsilon_air = 1; % relative permittivity of air
lambda_spp = real(lambda*sqrt((epsilon_au + epsilon_air)/(epsilon_au*epsilon_air))); % wavelength of spp
kspp = 2*pi/lambda_spp; % wavevector of spp
J0=0;J1=1;J2=2;J3=3;J4=4;
x = linspace(-1e-6,1e-6,500); y = linspace(-1e-6,1e-6,500);
[X,Y] = meshgrid(x,y); R = sqrt(X.^2 + Y.^2); theta = atan2(Y,X);

Espp0 = 1i^J0*besselj(J0,kspp*R).*exp(1i*J0*theta);
Espp1 = 1i^J1*besselj(J1,kspp*R).*exp(1i*J1*theta);
Espp2 = 1i^J2*besselj(J2,kspp*R).*exp(1i*J2*theta);
Espp3 = 1i^J3*besselj(J3,kspp*R).*exp(1i*J3*theta);
Espp4 = 1i^J4*besselj(J4,kspp*R).*exp(1i*J4*theta);

data_0 = load('J_0.txt');
E_0=data_0(:,3);

data_1 = load('J_1.txt');
E_1=data_1(:,3);

data_2 = load('J_2.txt');
E_2=data_2(:,3);

data_3 = load('J_3.txt');
E_3=data_3(:,3);

data_4 = load('J_4.txt');
E_4=data_4(:,3);

data_5 = load('J_5.txt');
E_5=data_5(:,3);

%figure;
%surf(X,Y,intensity./max(max(intensity)));
%shading interp;
%view(2);
%axis([-1 1 -1 1]);
%colormap('hot');
%xlabel('x (\mum)','FontSize',20,'FontName','Times New Roman');
%ylabel('y (\mum)','FontSize',20,'FontName','Times New Roman');
%set(gca,'FontName','Times New Roman');;
%set(gca,'FontWeight','bold');
%set(gca,'fontsize',16);
%colorbar;

data0=reshape(flipud(Espp0)',1,250000);
data1=reshape(flipud(Espp1)',1,250000);
data2=reshape(flipud(Espp2)',1,250000);
data3=reshape(flipud(Espp3)',1,250000);
data4=reshape(flipud(Espp4)',1,250000);

%LSQ
x0=data0;
x1=data1;
x2=data2;
x3=data3;
x4=data4;
y0=E_0';
y1=E_1';
y2=E_2';
y3=E_3';
y4=E_4';
%a0=[0 0 0 0 0];
comp0=[10 10 10 10 10];
a0(1:5)=real(comp0);
a0(6:10)=imag(comp0);
%options=optimset('lsqnonlin');
options=optimset('algorithm','levenberg-marquardt','TolFun',.000000000000001,'TolX',1e-8,...
    'MaxFunEvals',5000,'MaxIter',5000,'Display','iter');
a_0=lsqnonlin(@EE_new,a0,[],[],options,[x0;x1;x2;x3;x4],y0);
a_1=lsqnonlin(@EE_new,a0,[],[],options,[x0;x1;x2;x3;x4],y1);
a_2=lsqnonlin(@EE_new,a0,[],[],options,[x0;x1;x2;x3;x4],y2);
a_3=lsqnonlin(@EE_new,a0,[],[],options,[x0;x1;x2;x3;x4],y3);
a_4=lsqnonlin(@EE_new,a0,[],[],options,[x0;x1;x2;x3;x4],y4);
k_0=a_0(1:5)+1i*a_0(6:10); k0=k_0(1)
k_1=a_1(1:5)+1i*a_1(6:10); k1=k_1(2)
k_2=a_2(1:5)+1i*a_2(6:10); k2=k_2(3)
k_3=a_3(1:5)+1i*a_3(6:10); k3=k_3(4)
k_4=a_4(1:5)+1i*a_4(6:10); k4=k_4(5)

Espp_0=k0*Espp0;
Espp_1=k1*Espp1;
Espp_2=k2*Espp2;
Espp_3=k3*Espp3;
Espp_4=k4*Espp4;

x_0=reshape(flipud(Espp_0)',1,250000);
x_1=reshape(flipud(Espp_1)',1,250000);
x_2=reshape(flipud(Espp_2)',1,250000);
x_3=reshape(flipud(Espp_3)',1,250000);
x_4=reshape(flipud(Espp_4)',1,250000);
y=E_5';
comp0=[10 10 10 10 10];
a0(1:5)=real(comp0);
a0(6:10)=imag(comp0);
options=optimset('algorithm','levenberg-marquardt','TolFun',.000000000000001,'TolX',1e-8,...
    'MaxFunEvals',5000,'MaxIter',5000,'Display','iter');
a=lsqnonlin(@EE_new,a0,[],[],options,[x_0;x_1;x_2;x_3;x_4],y);
k=a(1:5)+1i*a(6:10)
