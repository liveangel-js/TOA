clc;
clear;
%t=load('out/case011_true_toa.txt');
base_count=load('out/case011_base_count.txt');
user_count=load('out/case011_user_count.txt');
deminsion = load('out/case011_dimension.txt');
[base]=load('out/case011_base_location.txt');
[t_t]=load('out/case011_toa.txt');
p1=load('out/case011_slope.txt');
p2=load('out/case011_intercept.txt');
%[x_ans] = load('out/case011_ans.txt');
base_location_start_row = 3+1;
base_location_end_row = 3 + base_count;
user_toa_start_row = base_count+3+1;
user_toa_end_row = user_toa_start_row + user_count -1;
%[t_t_m]=reshape(t_t,base_count*user_count,1);
%[t_m]=reshape(t,base_count*user_count,1);
[per_length]=3*10^8.*([t_t]*p1+p2);%距离修正

for detector_number=4:11
[QQ QC]=sort(per_length,2);
d=QQ(:,1:detector_number);
base_N=QC(:,1:detector_number);
co_base_j=base(base_N,:);
for r=1:detector_number
    co_base(:,(3*r-2):3*r)=co_base_j(((r-1)*user_count+1):r*user_count,1:3);
end
combos=combntns(1:detector_number,2)
for n=1:user_count
    for j=10:20
for m=1:detector_number*(detector_number-1)/2
A(m,1:3)=co_base(n,(combos(m,1)*3-2):combos(m,1)*3)-co_base(n,(combos(m,2)*3-2):combos(m,2)*3);
B(m,1)=(d(n,combos(m,2)))^2-(d(n,combos(m,1)))^2+(co_base(n,(combos(m,1)*3-2):combos(m,1)*3)-co_base(n,(combos(m,2)*3-2):combos(m,2)*3))*(co_base(n,(combos(m,1)*3-2):combos(m,1)*3)+co_base(n,(combos(m,2)*3-2):combos(m,2)*3))';
end
x_n(1,1:2)=(2.*A(:,1:2))\(B-A(:,3).*j*0.1);
x_n(1,3)=j*0.1;
deta(j-9) = sum(abs(sqrt(sum((repmat(x_n,detector_number,1)-base(QC(n,1:detector_number),1:3)).^2,2))-(per_length(n,QC(n,1:detector_number))')));
    end
    deta_0(n,:)=deta;
[PP PC]=min(deta,[],2);
PP_0(n,:)=PP/detector_number;
x(n,3)=(PC-1)/10+1;
x(n,1:2)=(2.*A(:,1:2))\(B-A(:,3).*x(n,3));
end

%dx=x(:,1:3)-x_ans(:,1:3);
%error=sqrt(dx(:,1).^2+dx(:,2).^2);
%error1=sqrt(dx(:,1).^2+dx(:,2).^2+dx(:,3).^2);

figure(1);
[GG GC]=hist(PP_0,50);
GG_cdf=cumsum(GG)/sum(GG);

for k=1:50
    if GG_cdf(k)>0.9
    GC_cdf=GC(k);%90%的可能性
    break
    end
end

GC_cdf_x(detector_number-3)=GC_cdf;
end

%plot([4:1:detector_number],GC_cdf_x);
xx=[4:1:detector_number];

%yy=polyfit(xx,GC_cdf_x,3);
%y1=polyval(yy,xx);
%plot(xx,y1); 
%hold on
%plot(xx,GC_cdf_x,'.'); 

figure(1);
sigma=GC_cdf_x-1./xx;
yy=polyfit(xx,sigma,3);
y1=polyval(yy,xx);
%plot(xx,y1); 
%hold on
%plot(xx,sigma,'.'); 
plot(xx,y1+1./xx); 
hold on
plot(xx,GC_cdf_x,'.'); 
title('误差与基站个数的关系case011','FontSize',15);
xlabel('基站数／个','FontSize',15);
ylabel('CDF','FontSize',15);
set(gca,'fontsize',12); 

saveas(1,'误差与基站个数的关系case011','jpg');
%saveas(2,'终端分布距离误差图case011','jpg');
save('xs_output_case_011.txt','yy','-ascii');
%save('output_case_011_2.txt','xxx','-ascii');
%save('距离cdf_case_011.txt','GC_cdf','-ascii');