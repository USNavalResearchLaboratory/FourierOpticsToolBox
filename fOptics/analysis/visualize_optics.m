function visualize_optics(varargin)

%Each argument is a Optic instance or field instance
figure; 
hold on

for i = 1:length(varargin)%for each argument

    obj = varargin{i};
    image = abs(obj.data); 
    image = mat2gray(image); 
    imsurf(image, [0,0,obj.position], [0,0,1], [], [], 'FaceAlpha',0.5);

end

view(3)
 
end