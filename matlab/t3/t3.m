clc;
clear;
base_count=load('out/case023_base_count.txt');
user_count=load('out/case023_user_count.txt');
deminsion = load('out/case023_dimension.txt');
[base]=load('out/case023_base_location.txt');
[t_t]=load('out/case023_toa.txt');
p1=load('out/case023_slope.txt');
p2=load('out/case023_intercept.txt');
base_location_start_row = 3+1;
base_locationend_row = 3 + base_count;
user_toa_start_row = base_count+3+1;
user_toa_end_row = user_toa_start_row + user_count -1;
[per_length]=3*10^8.*([t_t]*p1+p2);%距离修正

detector_number=6;
[QQ QC]=sort(per_length,2);
d=QQ(:,1:detector_number);
base_N=QC(:,1:detector_number);
co_base_j=base(base_N,:);
for r=1:detector_number
    co_base(:,(2*r-1):2*r)=co_base_j(((r-1)*user_count+1):r*user_count,1:2);
end
combos=combntns(1:detector_number,2)
for n=1:user_count
for m=1:detector_number*(detector_number-1)/2
A(m,1:2)=co_base(n,(combos(m,1)*2-1):combos(m,1)*2)-co_base(n,(combos(m,2)*2-1):combos(m,2)*2);
B(m,1)=(d(n,combos(m,2)))^2-(d(n,combos(m,1)))^2+(co_base(n,(combos(m,1)*2-1):combos(m,1)*2)-co_base(n,(combos(m,2)*2-1):combos(m,2)*2))*(co_base(n,(combos(m,1)*2-1):combos(m,1)*2)+co_base(n,(combos(m,2)*2-1):combos(m,2)*2))';
end
x(n,1:2)=(2.*A(:,1:2))\B;
end

yy1=plot(x(:,1),x(:,2),'.','MarkerSize',6);
hold on;

yy=polyfit(x(:,1),x(:,2),2);
y1=polyval(yy,x(:,1));
yy2=plot(x(:,1),y1,'.'); 
legend('终端坐标','拟合运动轨迹');

title('终端运动轨迹case023','FontSize',15);
xlabel('X坐标／米','FontSize',15);
ylabel('Y坐标／米','FontSize',15);
set(gca,'fontsize',12); 

saveas(1,'终端运动轨迹case023','jpg');
save('output_case_023.txt','x','-ascii');