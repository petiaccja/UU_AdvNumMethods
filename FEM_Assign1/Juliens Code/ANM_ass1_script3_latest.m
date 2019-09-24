%% Problem 1.1
%
% %hmax = 1/8
% h = 1/8;
% g = @circleg;
% [p,e,t] = initmesh(g, 'hmax', h);
% [ U, M ] = ANM_ass1_FEM_unstable( p, e, t );
%
%
% for i=1:size(U,2)
%     pdesurf(p,t,U(:,i));
%     xlim([-1,1])
%     ylim([-1,1])
%     zlim([-0.2,1])
%     drawnow;
% end %for, i

% %hmax = 1/16
% h = 1/16;
% g = @circleg;
% [p,e,t] = initmesh(g, 'hmax', h);
% [ U, M ] = ANM_ass1_FEM_unstable( p, e, t );
%
% for i=1:size(U,2)
%     pdesurf(p,t,U(:,i));
%     xlim([-1,1])
%     ylim([-1,1])
%     zlim([-0.2,1])
%     drawnow;
% end %for, i

%% Problem 1.2
%
% h = [1/4,1/8,1/16,1/32];
%
% L2E = zeros(length(h), 1);
%
% for i=1:length(h)
%     %mesh
%     g = @circleg;
%     [p,e,t] = initmesh(g, 'hmax', h(i));
%     [ U, M ] = ANM_ass1_FEM_unstable( p, e, t );
%     e = U(:,1)-U(:,end);
%     L2E(i) = sqrt(e'*M*e);
% end %for, hmax
%
% ls = polyfit(h', L2E, 2);
% x = linspace(0,1/4,26);
% y = zeros(1, 26);
% for i=1:length(x)
%     y(i) = ls(1)*x(i)^2+ls(2)*x(i)+ls(3);
% end %i
% size(x)
% size(y)
% subplot(1,2,1)
% plot(h, L2E, 'r', 'linewidth', 2)
% hold on
% loglog(h, h.^2, 'b', 'linewidth', 2)
% legend('Error function', 'h_{max} versus h_{max}^{\alpha}')
% xlabel('h')
% ylabel('L^2-error')
% subplot(1,2,2)
% plot(h, L2E, 'r', x, y, 'b', 'linewidth', 2)
% legend('Error function', 'Quadratic interpol.')
% xlabel('h')
% ylabel('L^2-error')
% print('ANM_Ass1_data1_convergence','-djpeg')

%% Problem 1.3
% %hmax = 1/8
% h = 1/8;
% g = @circleg;
% [p,e,t] = initmesh(g, 'hmax', h);
% [ U, M ] = ANM_ass1_FEM_unstable( p, e, t );
%
%
% for i=1:size(U,2)
%     pdesurf(p,t,U(:,i));
%     drawnow;
% end %for, i
%
% %hmax = 1/16
% h = 1/16;
% g = @circleg;
% [p,e,t] = initmesh(g, 'hmax', h);
% [ U, M ] = ANM_ass1_FEM_unstable( p, e, t );
%
% for i=1:size(U,2)
%     pdesurf(p,t,U(:,i));
%     xlim([-1,1])
%     ylim([-1,1])
%     zlim([-0.5,1.5])
%     drawnow;
% end %for, i


%
%
% h = [1/4,1/8,1/16,1/32];
%
% L2E = zeros(length(h), 1);
%
% for i=1:length(h)
%     %mesh
%     g = @circleg;
%     [p,e,t] = initmesh(g, 'hmax', h(i));
%     [ U, M ] = ANM_ass1_FEM_unstable( p, e, t );
%     e = U(:,1)-U(:,end);
%     L2E(i) = sqrt(e'*M*e);
% end %for, hmax
% L2E
% ls = polyfit(h', L2E, 1);
% x = linspace(0,1/4,26);
% y = zeros(1, 26);
% for i=1:length(x)
%     y(i) = ls(1)*x(i)+ls(2);%*x(i)+ls(3);
% end %i
%
% subplot(1,2,1)
% plot(h, L2E, 'r', 'linewidth', 2)
% hold on
% loglog(h, h.^1, 'b', 'linewidth', 2)
% legend('Error function', 'h_{max} versus h_{max}^{\alpha}')
% xlabel('h')
% ylabel('L^2-error')
% subplot(1,2,2)
% plot(h, L2E, 'r', x, y, 'b', 'linewidth', 2)
% legend('Error function', 'Linear interpol.')
% xlabel('h')
% ylabel('L^2-error')
% print('ANM_Ass1_data2_convergence','-djpeg')




