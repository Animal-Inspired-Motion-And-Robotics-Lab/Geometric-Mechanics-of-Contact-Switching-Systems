dom_lim = [-1,1];
dom_thresh = 10*diff(dom_lim)/100; 
dom = [dom_lim(1)-dom_thresh, dom_lim(2)+dom_thresh]; dom = linspace(dom(1),dom(2),1000);
dom_mid = dom_lim(1) + diff(dom_lim)/2;
dom_inp = dom;
% make sure the extended the points are flat outside the domain
dom(dom < dom_lim(1)) = dom_lim(1); dom(dom > dom_lim(2)) = dom_lim(2);
f1 = @(d,d_mid) 2./(1 + exp(-5*(d-d_mid))) - 1;
f2 = @(d,d_mid) -(2./(1 + exp(-5*(d-d_mid))) - 1);
p1 = f1(dom,dom_mid); p2 = f2(dom,dom_mid);

figure()
plot(dom_inp,p1,'r','linewidth',1.2,'DisplayName','c_1')
hold on; grid on;
plot(dom_inp,p2,'b','linewidth',1.2,'DisplayName','c_2')
legend('location','north')
xlabel(['Contact Domain [',dom_lim(1),',',dom_lim(2),']'])
xlabel('Contact Codomain (0,1)')
title('Sigmoidal Contact Interpolation map')
% Multiple of 10 is arbitrarily chosen since the difference between 0,1 are
% less than 1e-3

f1 = @(d,d_mid) d;
f2 = @(d,d_mid) -d;
p1 = f1(dom,dom_mid); p2 = f2(dom,dom_mid);

figure()
plot(dom_inp,p1,'r','linewidth',1.2,'DisplayName','c_1')
hold on; grid on;
plot(dom_inp,p2,'b','linewidth',1.2,'DisplayName','c_2')
legend('location','north')
xlabel(['Contact Domain [',dom_lim(1),',',dom_lim(2),']'])
xlabel('Contact Codomain [0,1]')
title('Linear Contact Interpolation map')

f1 = @(d,d_mid) cos(0.5*pi*(d+1));
f2 = @(d,d_mid) -cos(0.5*pi*(d+1));
p1 = f1(dom,dom_mid); p2 = f2(dom,dom_mid);

figure()
plot(dom_inp,p1,'r','linewidth',1.2,'DisplayName','c_1')
hold on; grid on;
plot(dom_inp,p2,'b','linewidth',1.2,'DisplayName','c_2')
legend('location','north')
xlabel(['Contact Domain [',dom_lim(1),',',dom_lim(2),']'])
xlabel('Contact Codomain [0,1]')
title('Cosine Contact Interpolation map')

% %%
% % In this section, we explore a contact conservation law based on the 
% % 2-norm instead of a simple 1-norm,
% 
% dom_lim = [-1,1];
% dom_thresh = 10*diff(dom_lim)/100; 
% dom = [dom_lim(1)-dom_thresh, dom_lim(2)+dom_thresh]; dom = linspace(dom(1),dom(2),1000);
% dom_mid = dom_lim(1) + diff(dom_lim)/2;
% dom_inp = dom;
% dom(dom < dom_lim(1)) = dom_lim(1); dom(dom > dom_lim(2)) = dom_lim(2);
% f1 = @(d,d_mid) 1./(1 + exp(-5*(d))); %-d_mid
% f2 = @(d,d_mid) sqrt(1 - (1./(1 + exp(-5*(d)))).^2); %-d_mid
% p1 = f1(dom,dom_mid); p2 = f2(dom,dom_mid);
% 
% figure()
% plot(dom_inp,p1,'r','linewidth',1.2,'DisplayName','c_1')
% hold on; grid on;
% plot(dom_inp,p2,'b','linewidth',1.2,'DisplayName','c_2')
% legend('location','north')
% xlabel(['Contact Domain [',dom_lim(1),',',dom_lim(2),']'])
% xlabel('Contact Codomain (0,1)')
% title('Sigmoidal Contact Interpolation map')
