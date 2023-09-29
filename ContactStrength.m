figure()
plot(linspace(-1,1,1000),(sin(pi/2*linspace(-1,1,1000))).^100,'r','linewidth',1.2);
grid on;
xlabel('Contact Variable c');
ylabel('Contact Strength');