clc;
clear;
%t=load('out/case001_true_toa.txt');
base_count=load('out/case001_base_count.txt');
user_count=load('out/case001_user_count.txt');
deminsion = load('out/case001_dimension.txt');
[base]=load('out/case001_base_location.txt');
[t_t]=load('out/case001_toa.txt');
p1=load('out/case001_slope.txt');
p2=load('out/case001_intercept.txt');
%[x_ans] = load('out/case001_ans.txt');
base_location_start_row = 3+1;
base_location_end_row = 3 + base_count;
user_toa_start_row = base_count+3+1;
user_toa_end_row = user_toa_start_row + user_count -1;
%[t_t_m]=reshape(t_t,base_count*user_count,1);
%[t_m]=reshape(t,base_count*user_count,1);
[per_length]=3*10^8.*([t_t]*p1+p2);%距离修正

detector_number=6;
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
[GX G1 G2]=plotyy(GC,GG,GC,GG_cdf,'bar','plot');
set(GX,'Xlim',[0 10],'FontSize',12);
set(GX(1),'ylim',[0 150],'ytick',[0:10:150]);
set(GX(2),'ylim',[0 1],'ytick',[0:0.1:1]);
%set(H1,'linestyly','-','color','b','linewidth',2.5);
set(G2,'linestyle','-','color','r','linewidth',1.5);
set(get(GX(1),'Ylabel'),'string','终端个数/个','FontSize',15);
set(get(GX(2),'Ylabel'),'string','CDF','FontSize',15,'FontName','Times New Roman');
set(get(GX(1),'Xlabel'),'string','拟合终端距离与优化测量距离的误差／米','FontSize',15);
title('距离误差分析sample case001','FontSize',15);
legend([G1 G2],'终端','CDF');
for k=1:50
    if GG_cdf(k)>0.9
    GC_cdf=GC(k);%90%的可能性
    break
    end
end

figure(2);
H3=plot(base(:,1),base(:,2),'.','color','r','MarkerSize',20);
legend(H3,'探测器位置');
hold on;
scatter3(x(:,1),x(:,2),PP_0,[],PP_0);
colorbar;
title('终端分布距离误差图sample case001','FontSize',15);
xlabel('X坐标／米','FontSize',15);
ylabel('Y坐标／米','FontSize',15);
set(gca,'fontsize',12); 

%saveas(1,'距离误差分析case001','jpg');
%saveas(2,'终端分布距离误差图case001','jpg');
save('output_case_001.txt','x','-ascii');
%save('output_case_001_2.txt','xxx','-ascii');
%save('距离cdf_case_001.txt','GC_cdf','-ascii');