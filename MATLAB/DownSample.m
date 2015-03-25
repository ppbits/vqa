function v2 = DownSample(v1, scale)

[M, N, L] = size(v1);
f = max(1,round(min(M,N)/scale));
%downsampling by f
%use a simple low-pass filter 

if(f>1)
    v2 = zeros( ceil(M/f), ceil(N/f), L);
%     size(v2)
    lpf = ones(f,f);
    lpf = lpf/sum(lpf(:));
    for i = 1:L
        v1(:,:,i) = imfilter(v1(:,:,i),lpf,'symmetric','same');
        v2(:,:,i) = v1(1:f:end,1:f:end, i);
    end
    v2 = uint8(v2);
else
    v2 = v1;
end

end
