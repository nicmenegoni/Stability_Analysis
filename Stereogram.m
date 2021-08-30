function Stereogram(dipdir, dip, nplane, figurename, plane0_Intersection1, idx3, xp, yp)
% Stereogram.m takes into account clusters and uses colour coding to
% distinguis them
%
% Note: Original script to draw Schmidt's net is by Gerry Middleton,
% November 1995.
%
% 12.10.2010 M. Markovaara-Koivisto, Aalto University School of Science and
% Technology, Finland

% if plane pole plot --> plane0_Intersection1 = 0
% if line vector plot -> plane0_Intersection1 = 1
close(findobj('type','figure','name','Stereogram coloured'))
figure('name', figurename, 'NumberTitle', 'off')
schmidt
% Draws the oriented data into Schmidt's net (G. Middleton, November 1995)
theta=(90-dipdir); % Poles are at the opposite direction to dip direction
r=sqrt(2)*sind((90-dip-90)/2); % Poles are perpendicular to the dip

% Coordinates on the strereographic projection
m=nplane;;
for i=1:m;
    xp(i) = r(i)*cosd(theta(i));
    yp(i) = r(i)*sind(theta(i));
    
    % Stereographic projection with colourcodes according to index idx3
    if idx3(i)==1
        if plane0_Intersection1 == 0
            plot(xp(i),yp(i),'+y')
        else
            plot(xp(i),yp(i),'oy')
        end
    elseif idx3(i)==2
        if plane0_Intersection1 == 0
            plot(xp(i),yp(i),'+r')
        else
            plot(xp(i),yp(i),'or')
        end
    elseif idx3(i)==3
        if plane0_Intersection1 == 0
            plot(xp(i),yp(i),'+b')
        else
            plot(xp(i),yp(i),'ob')
        end
    elseif idx3(i)==4
        plot(xp(i),yp(i),'vg')
    elseif idx3(i)==5
        plot(xp(i),yp(i),'vm')
    elseif idx3(i)==6
        plot(xp(i),yp(i),'vw')
    else
        if plane0_Intersection1 == 0
            plot(xp(i),yp(i),'vk')
        else
            plot(xp(i),yp(i),'ok')
        end
    end
    
end
