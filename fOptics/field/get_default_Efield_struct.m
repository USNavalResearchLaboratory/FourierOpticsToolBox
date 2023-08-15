function Efield_default = get_default_Efield_struct()
    target = struct();
    target.wvl = 650e-9;
    target.N = 1024;
    
    target.D = 0.002;
    target.dx = target.D/10;
    target.roi = target.dx * target.N;
    target.position = 0;
    target.targetType = 'gauss'; 
    target.distance = 1;
    target.specks = 1; 
    target = target_class(target); 
    target = target.genTarget_simplegauss(); 
    Efield_default = Efield(target);
    %Efield_default.showme
end
