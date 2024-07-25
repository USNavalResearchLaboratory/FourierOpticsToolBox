function propped = fresnel_prop(obj, zf, width1)
    %{
    Propagate an image (assumed to start in real space) a distance zf.
    %}
    arguments
        obj %array that we want to propagate
        zf
        width1 = obj.width
    end
    H = fresnel_propagator(zf,  width1, obj.N, ...
                            obj.wvl);

    ft = (jifft2(obj.data));
    proppedFt = ft .* H;
    propped = (jfft2(proppedFt));
end