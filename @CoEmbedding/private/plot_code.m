function [hx,hy]=plot_code(phi_x,psi_y,x_names,y_names,x_clr,y_clr);
%
% PLOT_CODE
% 
% [HX,HY]=PLOT_CODE(PHI_X,PSI_Y);
% [HX,HY]=PLOT_CODE(PHI_X,PSI_Y);
%
% Plot 2D embedding obtained with code
%
% Inputs: 
%   phi_x - embedding of x objects 
%   psi_x - embedding of y objects 
%   x_names - cell aray with string description of x objects
%             (optional)
%   y_names - cell aray with string description of y objects
%             (optional) 
%
% Outputs: 
%   HX,HY - handles for the printed x and y objects. 
% 

if(nargin<5), x_clr = 'b';end;
if(nargin<6), y_clr = 'k';end;

nx= size(phi_x,1);
ny= size(psi_y,1);

hold on ;
hx=[];
hy=[];

% Plot X objects 
if(exist('x_names') & ~isempty(x_names))
  for(ix=1:nx), 
    hx(ix) = text(phi_x(ix,1),phi_x(ix,2),x_names{ix},'Color',x_clr);  
  end
else
  for(ix=1:nx), 
    hx(ix) = plot(phi_x(ix,1),phi_x(ix,2),'.','Color',x_clr); 
  end  
end

% Plot Y objects 
if(exist('y_names')& ~isempty(y_names))
  for(iy=1:ny), 
    hy(iy) = text(psi_y(iy,1),psi_y(iy,2),y_names{iy},'Color',y_clr);
  end
else
  for(iy=1:ny), 
    hy(iy) = plot(psi_y(iy,1),psi_y(iy,2),'.','Color',y_clr); 
  end  
end

return
