function Sp_p = Sp_A(A_tensor, p)
% A_tensor: m × d × V   
Y = shiftdim(A_tensor, 1);   % d × V × m
Yf = fft(Y, [], 3);
Sp_p = 0;
for j = 1:size(Yf,3)         % 每个切片: d × V
    s = svd(Yf(:,:,j), 'econ');
    Sp_p = Sp_p + sum(s.^p);
end
end
