function a=jifft2(Fn)
a=fftshift(ifft2(ifftshift(Fn)));
end