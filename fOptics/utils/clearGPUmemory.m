function varargout = clearGPUmemory(varargin)
% Custom function to clear up GPU and pull data to CPU RAM
for i = 1:length(varargin)
    obj = varargin{i};
    obj.data = gather(obj.data);
    varargout{i} = obj; 
end