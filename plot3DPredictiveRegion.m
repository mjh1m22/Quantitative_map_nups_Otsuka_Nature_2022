% plot3DSpindle(X, Y, Z, volSpindle, volOocyte, thresh)
% plot two volumes in space set by X, Y, Z and thersholded with thresh.
%   thresh - 2x1 vector
% plot3DSpindle(X, Y, Z, volSpindle, volOocyte, thresh, T)
% Translate according to T (3x1 vector) 
% 
% author: Antonio Politi, mail@apoliti.de
function ifig = plot3DPredictiveRegion(ifig, X, Y, Z, vol1, vol2, vol3, vol4, vol5, thresh, color, alpha,  varargin)
if nargin > 12
    visible = varargin{1};
else
    visible = 'on';
end
smooth = 0;
if nargin > 12
    smooth = varargin{2};
end
try
    set(ifig, 'Visible',visible, 'Position',  [100 100 800 800], 'Renderer', 'OpenGL');
catch
    figure(ifig);
    set(ifig, 'Visible',visible, 'Position',  [100 100 800 800], 'Renderer', 'OpenGL');
end


if ~isempty(vol2)
    assert(all(size(vol2) == size(vol1)), 'Matrix size for two volumes must be identical');
end
%assert(length(X) == size(vol1,1), 'X dimension has not same size as volume');
%assert(length(Y) == size(vol1,2), 'Y dimension has not same size as volume');
%assert(length(Z) == size(vol1,3), 'Z dimension has not same size as volume');

if smooth
    vol1 = smooth3(vol1);
    if ~isempty(vol2)
        vol2 = smooth3(vol2);
    end
end
%patch creates vortexes out of surface. These can then be ligthened
patch1 = patch(isosurface(X, Y, Z, vol1, thresh(1)));
isonormals(X, Y, Z, vol1>= thresh(1), patch1);
set(patch1, 'facecolor', color{1}, 'edgecolor', 'none');
set(patch1, 'FaceAlpha', 'Flat', 'FaceVertexAlphaData', alpha(1));
if ~isempty(vol2)
    patch2 = patch(isosurface(X, Y, Z, vol2>=thresh(2), thresh(2)));
    isonormals(X, Y, Z, vol2, patch2);
    set(patch2, 'facecolor', color{2}, 'edgecolor', 'none');
    set(patch2, 'FaceAlpha', 'Flat', 'FaceVertexAlphaData', alpha(2));
end
if ~isempty(vol3)
    patch3 = patch(isosurface(X, Y, Z, vol3>=thresh(3), thresh(3)));
    isonormals(X, Y, Z, vol3, patch3);
    set(patch3, 'facecolor', color{3}, 'edgecolor', 'none');
    set(patch3, 'FaceAlpha', 'Flat', 'FaceVertexAlphaData', alpha(3));
end
if ~isempty(vol4)
    patch4 = patch(isosurface(X, Y, Z, vol4>=thresh(4), thresh(4)));
    isonormals(X, Y, Z, vol4, patch4);
    set(patch4, 'facecolor', color{4}, 'edgecolor', 'none');
    set(patch4, 'FaceAlpha', 'Flat', 'FaceVertexAlphaData', alpha(4));
end
if ~isempty(vol5)
    patch5 = patch(isosurface(X, Y, Z, vol5>=thresh(5), thresh(5)));
    isonormals(X, Y, Z, vol5, patch5);
    set(patch5, 'facecolor', color{5}, 'edgecolor', 'none');
    set(patch5, 'FaceAlpha', 'Flat', 'FaceVertexAlphaData', alpha(5));
end
daspect([1 1 1]);
%view(3);
axis tight;
camlight; 
%lighting flat;
lighting gouraud;
campos([20,30,25]);
end