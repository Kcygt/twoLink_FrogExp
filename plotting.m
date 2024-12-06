figure(1);hold on, grid off, xlabel('x-Direction [m]'), ylabel('z-Direction [m]'), axis equal
step = 2:1000:17001;
g1=plot(x_d_star(step,1),x_d_star(step,3),'*',LineWidth=2);
g2=plot(trajectory_first(step,3),trajectory_first(step,5),'*','Color',"#7E2F8E",LineWidth=2);
% g3=plot(trajectory_second(2:end,3),trajectory_second(2:end,5),'Color',"#EDB120",LineWidth=2);
% g4=plot(trajectory_third(2:end,3),trajectory_third(2:end,5),'Color',"#4DBEEE",LineWidth=2);
% legend([g1 g2 g3 g4],'The Given Trajectory','The First Followed Trajectory','The Second Followed Trajectory','The Third Followed Trajectory','Location','northeastoutside','Orientation','vertical')
% figure(2);hold on, grid off, xlabel('Time [s]'), ylabel('Force [N]'),%title('Force Measurement')
plot(trajectory_first(2:end,2),trajectory_first(2:end,12), 'Color',"#7E2F8E",'LineWidth',0.5);
plot(trajectory_second(2:end,2),trajectory_second(2:end,12), 'r','LineWidth',0.5);
plot(trajectory_third(2:end,2),trajectory_third(2:end,12), 'r','LineWidth',0.5);

plot(t_vec, f_d_star(:, 3),'r','LineWidth',2);
plot(t_vec+32, f_d_star(:, 3),'b','LineWidth',2);
plot(t_vec+64, f_d_star(:, 3),'b','LineWidth',2)

x = trajectory_first(step,3);
z = trajectory_first(step,5);
fAct = trajectory_first(step,12);
