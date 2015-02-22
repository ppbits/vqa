function s = combineSaliency(s1, s2)
len = size(s1,3);

s = zeros(size(s1));
for i = 1:len
    f1 = s1(:,:,i);
    f2 = s2(:,:,i);
    min1 = min(f1(:));
    min2 = min(f2(:));
    max1 = max(f1(:));
    max2 = max(f2(:));
    f1_n = (f1 - min1)/(max1-min1);
    f2_n = (f2 - min2)/(max2-min2);
%     std1 = std(f1(:));
%     std2 = std(f2(:));
    s(:,:,i) = f1_n + f1_n.*f2_n;
end

end