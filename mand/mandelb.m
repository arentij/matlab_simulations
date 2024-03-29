d = 500;
x0 = -2;%-1.713;%-1.23856;
x1 = 1;%-1.708;%-1.23851;
y0 = -1;%-8e-3;%0.41867;
y1 = 1;%-1*y0;%0.41877;

x = x0:1/d*(x1-x0):x1;
y = y0:1/d*(y1-y0):y1;

z = zeros(length(x),length(y));

for i = 1:length(x)
    for j = 1:length(y)
        c = x(i) + y(j)*1i;
        z(i,j) = length_inside(c).^-1;
    end
end
figure(1);
c = pcolor(x,y,z');
set(c, 'EdgeColor', 'none');
% colorpmap hot

h = zoom;
set(h, 'ActionPreCallback',  @myprecallback);
set(h, 'ActionPostCallback', @mypostcallback);
function myprecallback(h, eventdata)
% disp(eventdata);
end
function mypostcallback(h, eventdata)
disp(eventdata);
ac = gca;
x0 = ac.XLim(1);
x1 = ac.XLim(2);
y0 = ac.YLim(1);
y1 = ac.YLim(2);
[x, y, z] = upd_mand(x0,x1,y0,y1);
c = pcolor(x,y,z');
set(c, 'EdgeColor', 'none');
end
