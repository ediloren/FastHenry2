function plotfastH(name);
% Load and plot a 3D fasthenry structure produced with "zbuf -m name" 
% where name is without the ".mat"

xt = []; yt = []; zt =  []; 
xq = []; yq = []; zq = []; 

eval(['load ' name]);
fprintf(1, 'loaded %d panels\n', length(xt) + length(xq)); 
hold off; 
if length(xt) > 0, 
  ht = fill3(xt, yt, zt, 'r'); 
  hold on; 
end
X = max([max(xt) max(yt) max(zt)]); 
Y = min([min(xt) min(yt) min(zt)]); 

if length(xq) > 0, 
  hq = fill3(xq, yq, zq, 'y'); 
end; 

X = max([X max(xq) max(yq) max(zq)]); 
Y = min([Y min(xq) min(yq) min(zq)]); 


axis([Y X Y X Y X]); 
fprintf(1, 'finished filling polygons\n'); 

%return; 


axis('square');
%return; 


set(ht, 'FaceColor', 'w')
set(ht, 'EdgeColor', 'k')

set(hq, 'FaceColor', 'w')
set(hq, 'EdgeColor', 'k')
axis('square'); 
axis('off'); 
return; 

%g = get(hq(1), 'LineWidth'); set(hq, 'LineWidth', 2*g); 

f = gcf; 
set(f, 'Color', [1 1 1]); 
set(f, 'InvertHardcopy', 'off'); 

%set(f, 'PaperOrientation', 'landscape'); 
print -deps panels.ps

return; 

fprintf(1, 'printing...\n'); 
print -dps -Plouvre
!lpq -Plouvre

