t=load('sample_case001_true_time_out');
base_count=load('sample_case001_input_base_count');
user_count=load('sample_case001_input_user_count');
deminsion = load('sample_case001_input_de');
[base]=load('sample_case001_input_base_location');
[t_t]=load('sample_case001_input_toa');
[x_ans] = load('sample_case001_ans.txt');
base_location_start_row = 3+1;
base_location_end_row = 3 + base_count;
user_toa_start_row = base_count+3+1;
user_toa_end_row = user_toa_start_row + user_count -1;
[t_t_m]=reshape(t_t,base_count*user_count,1);
[t_m]=reshape(t,base_count*user_count,1);
p=polyfit(t_t_m,t_m,1);
[per_length]=3*10^8.*([t_t]*p(1)+p(2));%¾àÀëĞŞÕı

detector_number=8;
[QQ QC]=sort(per_length,2);
d=QQ(:,1:detector_number);
base_N=QC(:,1:detector_number);
co_base_j=base(base_N,:);
co_base=cell2mat(mat2cell(co_base_j,[user_count user_count user_count user_count user_count user_count user_count user_count])');
combos=combntns(1:detector_number,2)
for n=1:user_count
for m=1:detector_number*(detector_number-1)/2
A(m,1:3)=co_base(n,(combos(m,1)*3-2):combos(m,1)*3)-co_base(n,(combos(m,2)*3-2):combos(m,2)*3);
B(m,1)=(d(n,combos(m,2)))^2-(d(n,combos(m,1)))^2+(co_base(n,(combos(m,1)*3-2):combos(m,1)*3)-co_base(n,(combos(m,2)*3-2):combos(m,2)*3))*(co_base(n,(combos(m,1)*3-2):combos(m,1)*3)+co_base(n,(combos(m,2)*3-2):combos(m,2)*3))';
end
x(n,1:2)=(2.*A(:,1:2))\(B-A(:,3).*1.5);
end
dx=x-x_ans(:,1:2);
figure(1);
plot(base(:,1),base(:,2),'.');
hold on;
error=sqrt(dx(:,1).^2+dx(:,2).^2);
scatter3(x_ans(:,1),x_ans(:,2),error,[],error);
colorbar;
figure(2);
hist(dx,80);
