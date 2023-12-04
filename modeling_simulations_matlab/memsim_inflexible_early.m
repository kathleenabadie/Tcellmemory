clear all;
close all;

% Inflexible decision model (early decision)
% --------------------------------------------
% x - naive cells, cells/uL
% y - memory cells, cells/uL
% z - effector cells, cells/uL
% v - viral load,
% m - inactive macrophage: cells/uL
% n - activated macrophage: cells/uL

tspan = [0 5000]; % hr
y0 = [0.02 0 0 0.001 1 0];     % x y z v m n

params.ax = 5;   % naive cell activation rate / viral load / hr
params.gy = 0.25;   % maximal memory proliferation rate per unit antigen  /hr
params.Ky1 = 0.1;   % viral load required for half-maximal memory proliferation
params.ay = 0.12;  % maximal rate of memory to effector differentiation
params.Ky2 = 0.1;   % viral load required for half-maximal memory differentiation
params.gz = 0.25;   % maximal effector proliferation rate, 4 hrs
params.Kz1 = 0.1;  % viral concentration needed for half-maximal effector proliferation
params.bz = 0;  % maximum rate of effector de-differentiation
params.Kz2 = 0.0025;  % viral concentration needed for half-maximal inhibition of effector de-differentiation
params.dz = .016;  % rate of effector apoptosis
params.gv =  0;  % pathogen replication rate
params.N = 50;    % sharpness of extinction coefficient for pathogen
params.ep = 0.0001;   % epsilon, load of pathogen at which extinction happens, 10^(-4)
params.dv1 = 0.0045;  % T cell dependent viral killing rate
params.dv2 = 0.0015;   % T cell independent macrophage killing rate
params.dv3 = 0.0015;  %10^(-5);   % T cell dependent macrophage killing rate per unit effector T cell
params.am = 2.5;   % rate of virus-dependent macrophage activation / viral titer
params.gn = 0.02;    % this is the maximum rate of macrophage proliferation
params.Kn = 0.01;   % this is the viral load at which macrophage proliferation rate is half-maximal
params.dn = 0.0002;   % this is the macrophage turnover rate


% this is the fraction predicted from the analytical results
predictedfrac = params.bz/(params.bz+params.dz);

scan = -2.5:0.025:-0.35 % differences in viral replication rate
frac = [];  % the fraction of memory cells
peak = [];  % the total number of cells at the peak
mem = [];   % the number of memory cells at the end

examples = [46 82];   % indices for scanned parameters for plotting;

for i = 1:length(scan)

    i
    params.gv = 10.^(scan(i));   % changing the viral replication rate

    [t,data] = ode23t(@(t,y) fmem(t,y,params), tspan, y0);
    % these are the dynamical variables
    x = data(:,1);
    y = data(:,2);
    z = data(:,3);
    v = data(:,4);  
    m = data(:,5);
    n = data(:,6);
    total = x+y+z;

    frac = [frac y(end)/max(total)];
    peak = [peak max(total)];
    mem = [mem y(end)];

    tmax = 800;      
  
  if ismember(i,examples)
        figure(i);
        subplot(3,1,1)
        plot(t,y); hold on;
        plot([0 tmax], [0.05*peak(end) 0.05*peak(end)])
        set(gca,'Xlim',[0 tmax])
        set(gca,'YLim',[0 0.25*peak(end)]);
        xlabel('time')
        title(['memory cell (y) parameter value = ' num2str(scan(i))]);

        subplot(3,1,2)
        plot(t,x+y+z); hold on;
        set(gca,'Xlim',[0 tmax])
        set(gca,'YLim',[0 1.1*peak(end)]);
        xlabel('time')
        title('total T cells (T=x+y+z)')
  
        subplot(3,1,3)
        plot(t,v);
        set(gca,'Xlim',[0 tmax])
        xlabel('time')
        title('pathogen (v)')
        saveas(gcf,['Example' num2str(i)],'epsc');
    end
end

figure(i+1);
subplot(2,1,1);
semilogy(scan,frac,'k'); hold on;
semilogy([scan(1) scan(end)],[frac(examples(1)) frac(examples(1))]);
semilogy([scan(examples(1)), scan(examples(1))], [0 0.6],'k');
semilogy([scan(examples(2)), scan(examples(2))], [0 0.6],'k');
xlabel('replication rate');
ylabel('fraction memory');
title('reversible epigenetic switch')
set(gca,'XLim',[scan(1) scan(end)]);
set(gca,'YLim',[10^(-2.3) 1]);
axis square

subplot(2,1,2);
semilogy(scan,peak,'k'); hold on;
semilogy(scan,mem,'b');   % plot the memory cell number
semilogy([scan(1) scan(end)], [y0(1) y0(1)],'LineStyle','-');   %plot the naive cell number
semilogy([scan(examples(1)), scan(examples(1))], [10^(-2.5) 10^(3)],'k');
semilogy([scan(examples(2)), scan(examples(2))], [10^(-2.5) 10^(3)],'k');
xlabel('replication rate');
ylabel('maximum cells');
set(gca,'XLim',[scan(1) scan(end)]);
set(gca,'YLim',[10^(-2.5) 10^(3)]);
axis square

saveas(gcf,'Scan','epsc');


%% here we plot an example here 


function f = fmem(t,vec,params)

ax = params.ax;   % naive differentiation rate /hr
ay = params.ay;  % maximum rate of effector differentiation
gy = params.gy;   % maximal rate of memory proliferation /hr
Ky1 = params.Ky1;    % viral load at which proliferation is half-maximal
Ky2 = params.Ky2;    % viral load at which differentiation is half-maximal
bz = params.bz;    % maximal rate of effector de-differentiation
Kz1 = params.Kz1;  % viral concentration needed for half-maximal effector proliferation
Kz2 = params.Kz2;  % viral concentration needed for half-maximal effector de-differentiation
gz = params.gz;   % maximal rate of effector proliferation
dz = params.dz;  % rate of effector apoptosis
gv =  params.gv;  % viral replication rate     % also changed to 0.012 for low infection
dv1 = params.dv1;  % viral clearance rate by cytotoxic T cell
ep = params.ep;   % viral load at which extinction happens
N = params.N;    % sharpness of the extinction coefficient at the base

am = params.am;    % rate of macrophage activation
dv2 = params.dv2;   % T cell independent killing rate
dv3 = params.dv3;   % T cell dependent killing rate  per unit effector T cell
Kn = params.Kn;    % the viral load at which proliferation is half-maximal
gn = params.gn;    % the maximal rate of macrophage proliferation
dn = params.dn;    % death rate of activated macrophage

x = vec(1);   % naive
y = vec(2);   % memory
z = vec(3);   % effector
v = vec(4);   % virus
m = vec(5);   % inactive macrophage
n = vec(6);  % activated macrophage
f = zeros(6,1);

f(1) = -ax*v*x     ; %x - naive
f(2) = ax*v*x + gy*v/(v+Ky1)*y - ay*v/(v+Ky2)*y + bz/(v/Kz2+1)*z; %y - memory
f(3) = gz*v/(v+Kz1)*z + ay*v/(v+Ky2)*y - bz/(v/Kz2+1)*z - dz*z ; %z - effector
f(4) = gv*(v^(N)/(ep^(N-1)+v^(N-1))) - (dv1*z + dv2*n + dv3*n*z)*v;     ; %v - virus
f(5) = -am*v*m;    % m - inactive macrophage
f(6) = am*v*m + gn*v/(Kn+v)*n - dn*n; % mp - activated macrophage

end
