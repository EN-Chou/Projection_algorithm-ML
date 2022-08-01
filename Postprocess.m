%m_grid=xlsread('result.csv', 'A1:AO41');
%m_grid=xlsread('result.csv', 'A1:CC81');
%m_grid=xlsread('result.csv', 'A1:FE161');
%m_grid=xlsread('result.csv', 'A1:IG241');

u_grid=xlsread('u.csv', 'A1:CC81');
v_grid=xlsread('v.csv', 'A1:CC81');
p_grid=xlsread('p.csv', 'A1:CC81');
vel_grid=xlsread('velocity.csv', 'A1:CC81');

figure(1)
imagesc(u_grid);
axis square
figure(2)
imagesc(v_grid);
axis square
figure(3)
imagesc(p_grid);
axis square
figure(4)
imagesc(vel_grid);
axis square