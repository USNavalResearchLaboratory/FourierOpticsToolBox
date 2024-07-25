
function fig = figureFormats(fig,type)
    %Fig is figure object.
    % varargin are options for common types of figures like field, target,
    % detector, etc. 

    arguments 
        fig
        type = 'None'
    end
    
    
%    fontsize(fig,18, 'pixels')
%     fig.FontWeight = 'Bold'; %doesn't work. Need to know how to modify
%     fig
%     fig.FontName = 'Arial';
    axis image
    switch type
        case 'None'
            %do nothing
        case 'field'
            %something
        case 'intensity'
            %something
    end


end
