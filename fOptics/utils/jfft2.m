function a=jfft2(Fn)
a=fftshift(fft2(ifftshift(Fn)));
end