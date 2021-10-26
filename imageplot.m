function h=imageplot(varargin)

imagesc(varargin{:}); set(gca,'Ydir','Normal'); %daspect([1,1,1]);

% set(gca,'Color',0.5*[1,1,1]);

if(nargout>0)
    h=gca;
end

