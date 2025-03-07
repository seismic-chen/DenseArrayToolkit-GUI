% function[fillhandle,msg]=jbfill_1(xpoints,upper,lower,color,edge,add,transparency, parent)
% %USAGE: [fillhandle,msg]=jbfill(xpoints,upper,lower,color,edge,add,transparency)
% %This function will fill a region with a color between the two vectors provided
% %using the Matlab fill command.
% %
% %fillhandle is the returned handle to the filled region in the plot.
% %xpoints= The horizontal data points (ie frequencies). Note length(Upper)
% %         must equal Length(lower)and must equal length(xpoints)!
% %upper = the upper curve values (data can be less than lower)
% %lower = the lower curve values (data can be more than upper)
% %color = the color of the filled area 
% %edge  = the color around the edge of the filled area
% %add   = a flag to add to the current plot or make a new one.
% %transparency is a value ranging from 1 for opaque to 0 for invisible for
% %the filled color only.
% %
% %John A. Bockstege November 2006;
% %Example:
% %     a=rand(1,20);%Vector of random data
% %     b=a+2*rand(1,20);%2nd vector of data points;
% %     x=1:20;%horizontal vector
% %     [ph,msg]=jbfill(x,a,b,rand(1,3),rand(1,3),0,rand(1,1))
% %     grid on
% %     legend('Datr')
% % 如果没有传入 parent 参数，则默认使用 gca
% if nargin < 8
%     parent = gca;
% end
% if nargin<7;transparency=.5;end %default is to have a transparency of .5
% if nargin<6;add=1;end     %default is to add to current plot
% if nargin<5;edge='k';end  %dfault edge color is black
% if nargin<4;color='b';end %default color is blue
% 
% if length(upper)==length(lower) && length(lower)==length(xpoints)
%     msg='';
%     filled=[upper,fliplr(lower)];
%     xpoints=[xpoints,fliplr(xpoints)];
%     if add
%         hold on
%     end
%     fillhandle = fill(parent, xpoints, filled, color);
%     set(fillhandle, 'EdgeColor', edge, 'FaceAlpha', transparency, 'EdgeAlpha', transparency);
%     if add
%         hold(parent, 'off')
%     end
% else
%     msg='Error: Must use the same number of points in each vector';
% end
function [fillhandle,msg] = jbfill(xpoints, upper, lower, color, edge, add, transparency, parent)
    % 如果没有传入 parent 参数，则默认使用 gca
    if nargin < 8
        parent = gca;
    end
    if nargin < 7; transparency = 0.5; end
    if nargin < 6; add = 1; end
    if nargin < 5; edge = 'k'; end
    if nargin < 4; color = 'b'; end

    if length(upper)==length(lower) && length(lower)==length(xpoints)
        msg = '';
        filled = [upper, fliplr(lower)];
        xpoints = [xpoints, fliplr(xpoints)];
        if add
            hold(parent, 'on')
        end
        % 这里使用 'Parent' 指定父对象
        fillhandle = fill(parent, xpoints, filled, color);
        set(fillhandle, 'EdgeColor', edge, 'FaceAlpha', transparency, 'EdgeAlpha', transparency);
        if add
            hold(parent, 'off')
        end
    else
        msg = 'Error: Must use the same number of points in each vector';
    end
end
