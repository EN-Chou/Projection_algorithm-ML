%m_grid=xlsread('result.csv', 'A1:AO41');
%m_grid=xlsread('result.csv', 'A1:CC81');
%m_grid=xlsread('result.csv', 'A1:FE161');
%m_grid=xlsread('result.csv', 'A1:IG241');

load u.dat
load v.dat
load p.dat
% load pred_p.dat
% u=fliplr(u);
% v=fliplr(v);
u=u.';
v=v.';
p=p.';
pred_p=pred_p.';
velocity=(u.^2+v.^2).^0.5;

figure(1)
contourf(u);
xlabel('i')
ylabel('j')
colorbar
axis square
set(gcf,'position',[800 500 500 500])

figure(2)
contourf(v);
xlabel('i')
ylabel('j')
colorbar
axis square
set(gcf,'position',[1300 500 500 500])

figure(3)
imagesc(velocity);
xlabel('i')
ylabel('j')
colorbar
axis square
set(gca,'YDir','normal') 
set(gcf,'position',[800 0 500 500])

figure(6)
imagesc(p);
xlabel('i')
ylabel('j')
colorbar
axis square
set(gca,'YDir','normal') 
set(gcf,'position',[1300 0 500 500])

% figure(7)
% imagesc(pred_p);
% xlabel('i')
% ylabel('j')
% colorbar
% axis square
% set(gca,'YDir','normal') 
% set(gcf,'position',[1300 0 500 500])

%Benchmark
x_b=[0 0.0625 0.0703 0.0781 0.0938 0.1563 0.2266 0.2344 0.5000 0.8047 0.8594 0.9063 0.9453 0.9531 0.9609 0.9688 1.0];
y_b=[0 0.0547 0.0625 0.0703 0.1016 0.1719 0.2813 0.4531 0.5000 0.6172 0.7344 0.8516 0.9531 0.9609 0.9688 0.9766 1.0];
Re_100_v=[0 0.09233 0.10091 0.10890 0.12317 0.16077 0.17507 0.17527 0.05454 -0.24533 -0.22445 -0.16914 -0.10313 -0.08864 -0.07391 -0.05906 0];
Re_100_u=[0 -0.03717 -0.04192 -0.04775 -0.06434 -0.10150 -0.15662 -0.21090 -0.20581 -0.13641 0.00332 0.23151 0.68717 0.73722 0.78871 0.84123 1];
%grid point
for i=1:81
    pt(i)=(i-1)/80;
end
% plot x-v
figure(4)
hold on
%plot(x_b, Re_100_v, 'o')
plot(pt, v(41, :), 'x')
axis square
hold off
set(gcf,'position',[0 500 500 500])
%plot y-u
figure(5)
hold on
%plot(Re_100_u, y_b'o')
plot(u(:, 41), pt, 'x')
axis square
hold off
set(gcf,'position',[0 0 500 500])
