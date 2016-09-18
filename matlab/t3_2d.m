clc;
clear;
base_count=load('out/case024_base_count.txt');
user_count=load('out/case024_user_count.txt');
deminsion = load('out/case024_dimension.txt');
[base]=load('out/case024_base_location.txt');
[t_t]=load('out/case024_toa.txt');
base_location_start_row = 3+1;
base_location_end_row = 3 + base_count;
user_toa_start_row = base_count+3+1;
user_toa_end_row = user_toa_start_row + user_count -1;
[p]=[0.749758269 -1.26E-09];
[per_length]=3*10^8.*([t_t]*p(1)+p(2));%æ‡¿Î–ﬁ’˝

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

plot(x(:,1),x(:,2),'.');
