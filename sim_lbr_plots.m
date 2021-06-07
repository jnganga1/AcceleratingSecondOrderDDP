%%
b=3;a=20;
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');
set(0, 'DefaultAxesFontSize', b);

time = 0:dt:dt*(params.N-1); eTime = ceil(dt*(params.N-1)); 
figure; 
subplot(2,2,1); hold on
for i=1:7 
    s = sprintf('$q_%i$',i);
    plot(time,xbar(i,:),'DisplayName',s,'LineWidth',b);
end
grid on; grid minor
hold off; xlabel('Time (s)'); ylabel('State q (rad)')
xlim([time(1) time(end)]); l = legend; l.Interpreter = 'latex';
set(gca,'FontSize',a)

subplot(2,2,2); hold on
for i=1:7 
    s= sprintf('$\\dot{q}_%i$',i);
    plot(time,xbar(i+7,:),'DisplayName',s,'LineWidth',b);
end
grid on; grid minor
hold off; xlabel('Time (s)'); ylabel('State $\dot{q}$ (rad/s)','interpreter','latex')
xlim([time(1) time(end)]);l = legend; l.Interpreter = 'latex';
set(gca,'FontSize',a)


subplot(2,2,3); hold on
for i=1:7 
    s = sprintf('$\\tau_%i$',i);
    plot(time,ubar(i,:),'DisplayName',s,'LineWidth',b);
end
grid on; grid minor
hold off; xlabel('Time (s)'); ylabel('Torque (N)')
xlim([time(1) time(end)]);l = legend; l.Interpreter = 'latex';
set(gca,'FontSize',a)


subplot(2,2,4)
semilogy(Vstore-Vstore(end),'LineWidth',b);
grid on; grid minor
hold off; ylabel('Cost'); xlabel('Iterations')
set(gca,'FontSize',a)

