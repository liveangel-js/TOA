
%[per_length]=3*10^8.*b(44:1243,1:40);
%[terminal_coordinate]=zeros(1100,3);
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
