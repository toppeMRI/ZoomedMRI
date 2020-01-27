function [b1] = sinc_pulse_gen(tip_angle,slice_thickness,k_max,slew_rate,g_max,dt)

gam=26751;                          % rad/sec/g
gambar = gam/2/pi;                  % Hz/g
g_std_max = 4;

gradient_info.slew_rate = slew_rate;
gradient_info.g_max = g_max;


% a reference RF amplitude in Tesla for a 30 degree flip angle pulse
% when an idea gradient(no rise time) 4g/cm is applied.                                    
b1amp = 3.3368e-005*tip_angle/30;   
                                    
g1 = gen_trapezoid(k_max,gradient_info,dt);
g2 = gen_trapezoid(-k_max/2,gradient_info,dt);

g = [g1;g2];
k = -(sum(g)-cumsum(g))*dt*gambar;
b1 = sinc(k(1:length(g1))*slice_thickness).*abs(g1/g_std_max)*b1amp*slice_thickness;

b1 = b1.*hamming(length(b1));
b1 = [b1(1:end-1);0];

%b1 = [b1;zeros(length(g2),1)];
disp(sprintf('unit pulse length %f ms',length(b1)*dt*10^3));
return;
%{
% Bloch simulation 
z = [-3:0.1:3]; %z samples in cm
Mi = [0;0;1]*ones(1,length(z));
bx = b1*ones(1,length(z));
by = zeros(size(b1))*ones(1,length(z));
bz = g*z*10^(-4); %gauss to Tesla
T1 = 1000*ones(size(z));
T2 = 100*ones(size(z));
[mx,my,mz] = blochsim2(Mi, bx, by, bz, T1, T2, dt*10^(3));
mx_f = mx(end,:);
my_f = my(end,:);
m_trav = mx(end,:)+1i*my(end,:);
figure,plot(z,abs(m_trav));
figure,plot(z,angle(m_trav));
figure,plot(z,mx_f);


return;
gambar = 42.57e3;
dt = 4e-3; %in msec
z = [-3:0.1:3]; %z samples in cm

t_up = 4/15000;   %sec
t_flat = 3/2/4257-t_up; %usec
t_up_duration = floor(t_up*(10^6)/4);
t_flat_duration = floor(t_flat*(10^6)/4);
t = [0:1:3*t_up_duration+3*t_flat_duration-1]*4*10^(-3);
g_up = [0:1:t_up_duration-1]*4*10^(-6)*dgdtmaxtrap;
g_flat = ones(1,t_flat_duration)*4;
g_down = g_up(length(g_up):-1:1);
g_flat_n = -[g_flat];
%g_flat_n = -[g_flat(1:end) 4 4 4 4 4 4 4 4 4 4 4 4];
g = [g_up g_flat g_flat g_down -g_up g_flat_n];
g = g*(10^-4);
k = -(sum(g)-cumsum(g))*dt*gambar;
b1 = (sinc(0.5*k(1:length([g_up g_flat])*2)));
b1 = b1.*hamming(length(b1))';
b1 = b1amp*[b1';zeros(length([g_up g_flat_n]),1)];
%figure;plot(k,b1);
%mask = find(k<-3 | k > 3);
%b1(mask)=0;
%mask = find(-3<=k & k<=3);
%b1(mask)= b1(mask).*abs(g(mask)')/4;
b1 = b1.*abs(g')/4*10^4;
b1 = [0;0;b1(1:end-2)];
%b1 = [0;b1(1:end-2)];

%{
t1 = [0:4:176*2-4]*(10^(-3)); %time samples
t2 = [176*2:4:176*3]*(10^(-3));
t = [t1 t2];
g1 = ones(length(t1),1)*gmax; %4g/cm
g2 = -ones(length(t2),1)*gmax;
g = [g1;g2];
k = -(sum(g)-cumsum(g))*dt*gambar;
b1 = (sinc(k(1:length(t1)))).*hamming(length(t1));
b1 = [b1;zeros(length(g2),1)]*b1amp;
%}

%plot(t,k,t,g,t,b1*10^5)
%b1(mask)=b1(mask)./abs(g(mask)');
%}
%}


%t = [0:1:t_up_duration-1];
%t = [t t(end)+1:1:t(end)+t_flat_duration-1];
%t = [t t(end)+1:1:t(end)+t_flat_duration-1];
%t = [t t(end)+1:1:t(end)+t_up_duration-1];

%{
t1 = [0:4:176*2-4]*(10^(-3)); %time samples
t2 = [176*2:4:176*3]*(10^(-3));
t = [t1 t2];
g1 = ones(length(t1),1)*gmax; %4g/cm
g2 = -ones(length(t2),1)*gmax;
g = [g1;g2];
k = -(sum(g)-cumsum(g))*dt*gambar;
b1 = (sinc(k(1:length(t1)))).*hamming(length(t1));
b1 = [b1;zeros(length(g2),1)]*b1amp;

Mi = [0;0;1]*ones(1,length(z));
bx = b1*ones(1,length(z));
by = zeros(size(b1))*ones(1,length(z));
bz = g*z;
T1 = 1000*ones(size(z));
T2 = 100*ones(size(z));
[mx,my,mz] = blochsim2(Mi, bx, by, bz, T1, T2, dt);
mx_f = mx(end,:);
my_f = my(end,:);
m_trav = mx(end,:)+1i*my(end,:);
figure,plot(z,abs(m_trav));
figure,plot(z,angle(m_trav));
return;
%}
