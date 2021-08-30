clear all
close all

%% Define slope orientation dipdir/dip and friction and lateral limit values
%values dipidir/dip of the slope 'parete3 Ormea' are 300/69

slope_min_dip= 75;
slope_max_dip= 75;
dipstep=1;
slope_min_dipdir=280;
slope_max_dipdir=320;
dipdirstep=2;
    %Define friction angle
    friction = 30;
    
    %Define lateral limits
    latlimit = 20;
plot=0; %specify 1 or 0 if you wnt or not the plots

%

%% Select data to load
%plane data with set
[filename, pathname] = uigetfile({'*.xlsx', 'Select the XLSX file of discontinuity plane with set interpretation'},'Select XLSX file of discontinuity plane with set interpretation',...
    'F:\Menegoni\Ormea\Nuovi voli giugno 2017\Albris_SenseFly\parete\obj\measures')
% pathname='F:\Menegoni\Ormea\Nuovi voli giugno 2017\Albris_SenseFly\parete\obj\measures\Fit_Matlab\interpretationset';
% filename='Group.xlsx';
data=readtable(fullfile(pathname,filename));

% % uiwait(msgbox(['Data of discontinuity Intersections have been select' newline 'SELECT INTERSECTION DISCONTINUITY DATA']));
% %Intersection data
% [Int_filename, Int_pathname] = uigetfile({'*.xlsx', 'Select the XLSX file of the discontinutiy intersections'},'Select XLSX file of the discontinutiy intersections',...
%     'F:\Menegoni\Ormea\Nuovi voli giugno 2017\Albris_SenseFly\parete\obj\measures')
% % Int_filename='Intersection_Fit_Group..xlsx';
% % Int_pathname='F:\Menegoni\Ormea\Nuovi voli giugno 2017\Albris_SenseFly\parete\obj\measures\Fit_Matlab\intersections_Group';
% Int_data=readtable(fullfile(Int_pathname,Int_filename));

%Selection of folder in which images will be saved
if plot == 1
uiwait(msgbox('Select the path in which images will be saved'));
pathFig=uigetdir([...
    'F:\Menegoni\Ormea\Nuovi voli giugno 2017\Albris_SenseFly\parete\obj\measures\Fit_Matlab']);
end

%% Read data and write Dip, Dip Direction, Set, Pole Dip Direction and Dip vectors
%Discontinuity planes
num_disc=numel(data)/3;
Dip=data.Dip(:);%read Dip value of discontinuities
DipDir=data.DipDirection(:);%read Dip Direction value of discontinuities
Set=data.Set(:);%read Set value of discontinuities (random values are defined as NaN as default)
Set(isnan(Set))=0;%change Set value of random discontinuities from NaN to 0
Set_name=unique(Set);
nSet=numel(unique(Set));
for i= 1 : num_disc %calculate discontinuity line pole (dip and dip direction)
    %this is an trick to plot plunge and trend line vector
    pole_Dip(i,1) = 90 - Dip(i);%discontinuity pole dip
    if DipDir(i) < 180%discontinuity pole dip direction
        pole_DipDir(i,1) = DipDir(i) + 180;
    else
        pole_DipDir(i,1) = DipDir(i) - 180;
    end
    
end
% %Discontinuity Intersection
% num_Int=numel(Int_data.Trend);
% Trend=Int_data.Trend;
% Plunge=Int_data.Plunge;
% Int_Sets(:,1)=Int_data.Set_i;
% Int_Sets(:,2)=Int_data.Set_j;
% 
% for i= 1 : num_Int %calculate discontinuity line pole (dip and dip direction)
%     %this is an trick to plot plunge and trend line vector
%     pole_Plunge(i,1) = 90 - Plunge(i);%discontinuity pole dip
%     if Trend(i) < 180%discontinuity pole dip direction
%         pole_Trend(i,1) = Trend(i) + 180;
%     else
%         pole_Trend(i,1) = Trend(i) - 180;
%     end
%     
% end


%Plot discontinuity plane poles and discontinuity intersections lower hemisphere projections
if plot==1
    Stereogram(DipDir, Dip, num_disc, 'Discontinutiy plane poles', 0, Set )
%     Stereogram(pole_Trend, pole_Plunge, num_Int, 'Discontinuity intersection', 0, Trend )
end
disp('DipDir Dip %PS %FT ')
%% Kynematic analysis
for dipdirloop = slope_min_dipdir: dipdirstep : slope_max_dipdir
for diploop= slope_min_dip: dipstep : slope_max_dip
    clearvars planeslope polesslope strikeslope sxlatlimit dxlatlimit
    planeslope = zeros(1,2);
    planeslope = [dipdirloop, diploop];
    
    %Setting properly planeslope and poleslope
    if planeslope(1,1) < 180%calculate poles dip and dip direction of the slope
        poleslope = [planeslope(1,1)+180, 90 - planeslope(1,2)];
    else
        poleslope = [planeslope(1,1)-180, 90 - planeslope(1,2)];
    end
    
    if planeslope(1,1)<90%calculate strike of the slope
        strikeslope = 360 + planeslope(1,1) - 90;
    else
        strikeslope = planeslope(1,1) - 90;
    end
    
    if planeslope(1,1) < latlimit
        sxlatlimit = 360 + planeslope(1,1) - latlimit;%sxlatlimit=laterl limit at left of dip direction of slope
        dxlatlimit = planeslope(1,1) + latlimit;%dxlatlimit=laterl limit at rigth of dip direction of slope
    elseif planeslope(1,1) >= 360 - latlimit
        sxlatlimit = planeslope(1,1) - latlimit;
        dxlatlimit = planeslope(1,1) + latlimit - 360;
    else
        sxlatlimit = planeslope(1,1) - latlimit;
        dxlatlimit = planeslope(1,1) + latlimit;
    end
    
%% Planare sliding (PS) kinematic analysis
%A plane could act as a sliding plane if is dip-slope and dipping with an
%angle lower than apparent angle of the slope along dipdirection of the plane
%and higher than frction angle.

app_angle=zeros(num_disc,1);
PlanarSliding=zeros(num_disc, 1);
for i = 1 : num_disc % calculate if discontinuity plane could be a sliding plane
    
    if (DipDir(i) > strikeslope && DipDir(i) < strikeslope + 180)%Calculate apparent angle of the slope along Dip Direction of the discontinuity
        app_angle(i,1) = atand(tand(planeslope(1,2)) * cosd(DipDir(i)- planeslope(1,1)));
    elseif (strikeslope+180> 360 && DipDir(i) < strikeslope - 180)
        app_angle(i,1) = atand(tand(planeslope(1,2)) * cosd(360 + DipDir(i)- planeslope(1,1)));
    else
        app_angle(i,1)=NaN;
    end
    if Dip(i) > friction && Dip(i) < app_angle(i)
        if (DipDir(i)>sxlatlimit && DipDir(i)<dxlatlimit) || (sxlatlimit >dxlatlimit && DipDir(i)<sxlatlimit && DipDir(i)<dxlatlimit)|| (sxlatlimit >dxlatlimit && DipDir(i)>sxlatlimit && DipDir(i)>dxlatlimit)
            %1st case if sx is < dx latlimit, 2nd and 3rd when dx < sx
            %lat limit with DipDir of discontinuity >340 && <360 (3rdcase) and
            %DipDir >0 && <20 (2nd case).
            PlanarSliding(i) = 2;%value assigned to critical discontinuity inside the lateral limits
        else
            PlanarSliding(i) = 1;%value assigned to crtitical discontinuity outside the lateral limits
        end
        
    else
        PlanarSliding(i) = NaN;
    end
end

%Calculate and write statistic
PS_Hcrit=zeros;
    PS_Lcrit=zeros;
PS_Hcrit=sum(PlanarSliding(:) == 2);
PS_Lcrit=sum(PlanarSliding(:) == 1);
if plot==1
    Stereogram(DipDir, Dip, num_disc, 'PlanarSliding', 0, PlanarSliding )
    hold on
    title(['PlanarSliding, N = ', num2str(num_disc) newline 'critical with lateral limits (red) =', num2str(PS_Hcrit) newline 'critical without lateral limits (colored) =', num2str(PS_Hcrit + PS_Lcrit)])
    hold off
end
% saveas(gcf, fullfile(pathFig,'PlanarSliding.fig'))
% saveas(gcf, fullfile(pathFig,'PlanarSliding.png'))
%remember to write a code to export the results.


%% Flexural Toppling (FT) kinematic analysis
slopelimit=[planeslope(1,1),  planeslope(1,2)-friction];
pole_app_angle=zeros(num_disc,1);
FlexuralToppling = zeros(num_disc,1);
for i = 1 : num_disc
    %Calculate apparent angle of the slope along Dip Direction of the discontinuity
    if (pole_DipDir(i)>sxlatlimit && pole_DipDir(i)<dxlatlimit) || (sxlatlimit >dxlatlimit && pole_DipDir(i)<sxlatlimit && pole_DipDir(i)<dxlatlimit)|| (sxlatlimit >dxlatlimit && pole_DipDir(i)>sxlatlimit && pole_DipDir(i)>dxlatlimit)
        %if cylcle for Pole DipDirections included in the lateral limit
        
        if (pole_DipDir(i) > strikeslope && pole_DipDir(i) < strikeslope + 180)
            pole_app_angle(i,1) = atand(tand(slopelimit(1,2)) * cosd(pole_DipDir(i)- slopelimit(1,1)));
        elseif (strikeslope+180> 360 && pole_DipDir(i) < strikeslope - 180)
            pole_app_angle(i,1) = atand(tand(slopelimit(1,2)) * cosd(360 + pole_DipDir(i)- slopelimit(1,1)));
        else
            pole_app_angle(i,1)=NaN;
        end
    end
    
    %Calculate if discontinuity could be affected by flexural toppling
    if pole_Dip(i) < pole_app_angle(i)
        FlexuralToppling(i) = 2;
    else
        FlexuralToppling(i) = NaN;
    end
    
end
%Calculate and write statistic
  FT_crit=zeros;
FT_crit=sum(FlexuralToppling(:) == 2);
if plot==1
    Stereogram(DipDir, Dip, num_disc, 'FlexuralToppling', 0, FlexuralToppling )
    hold on
    title(['FlexuralToppling, N = ', num2str(num_disc) newline 'critical with lateral limits =', num2str(FT_crit)])
    hold off
end
% saveas(gcf, fullfile(pathFig,'FlexuralToppling.fig'))
% saveas(gcf, fullfile(pathFig,'FlexuralToppling.png'))

% %% WedgeSliding (WS) kinematic analysis
% Int_app_angle=zeros(num_Int,1);
% Int_app_angle_friction=zeros(num_Int,1);
% WedgeSliding=zeros(num_Int, 1);
% for i = 1 : num_Int % calculate if discontinuity plane could be a sliding plane
%     
%     if (Trend(i) > strikeslope && Trend(i) < strikeslope + 180)%Calculate apparent angle of the slope along Dip Direction of the discontinuity
%         Int_app_angle(i,1) = atand(tand(planeslope(1,2)) * cosd(Trend(i)- planeslope(1,1)));
%         Int_app_angle_friction(i,1) = atand(tand(friction) * cosd(Trend(i)- planeslope(1,1)));
%         
%     elseif (strikeslope+180> 360 && Trend(i) < strikeslope - 180)
%         
%         Int_app_angle(i,1) = atand(tand(planeslope(1,2)) * cosd(360 + Trend(i)- planeslope(1,1)));
%         Int_app_angle_friction(i,1) = atand(tand(friction) * cosd(Trend(i)- planeslope(1,1)));
%     else
%         Int_app_angle(i,1)=NaN;
%         Int_app_angle_friction(i,1)=NaN;
%     end
%     
%     if Plunge(i) < Int_app_angle(i)
%         
%         if Plunge(i)>friction
%             WedgeSliding(i) = 2;%value assigned to critical discontinuity inside the lateral limits
%         elseif Plunge(i)<friction && Plunge(i)>Int_app_angle_friction(i)
%             WedgeSliding(i) = 1;%value assigned to crtitical discontinuity outside the lateral limits
%             
%         else
%             WedgeSliding(i) = NaN;
%         end
%     else
%         WedgeSliding(i) = NaN;
%     end
%     
%     
% end
% %Calculate and write statistic
% WS_Hcrit=sum(WedgeSliding(:) == 2);
% WS_Lcrit=sum(WedgeSliding(:) == 1);
% if plot==1
%     Stereogram(pole_Trend, pole_Plunge, num_Int, 'WedgeSliding', 1, WedgeSliding )
%     hold on
%     title(['WedgeSliding, N = ', num2str(num_Int) newline ...
%         'wedge that could slide on both the parental discontinuies (red) =', num2str(WS_Hcrit) newline ...
%         'critical without lateral limits (pink) =', num2str(WS_Lcrit)])
%     hold off
% end
% % saveas(gcf, fullfile(pathFig,'WedgeSliding.fig'))
% % saveas(gcf, fullfile(pathFig,'WedgeSliding.png'))
% 
% %% DirectToppling (DT) and ObliqueToppling (OT) (with BasalPlane, BP)  kinematic analysis
% DirectToppling=zeros(num_Int,1);
% ObliqueToppling=zeros(num_Int,1);
% BasalPlane=zeros(num_disc,1);
% upperTopplingLimit = 90 - planeslope(1,2);%angle between slope face and intersection
% %(ho qualche dubbio sul considerare questo limite costante e non variabile
% %in base all'inclinazione apparente della parete lungo la direzione
% %d'immersione dell'intersezione.
% 
% if sxlatlimit(1,1) < 180%calculate the opposite of the lateral limit (sx)
%     oppsxlatlimit = sxlatlimit+180;
% else
%     oppsxlatlimit = sxlatlimit-180;
% end
% if dxlatlimit(1,1) < 180%calculate the opposite of the lateral limits (dx)
%     oppsxlatlimit = dxlatlimit+180;
% else
%     oppdxlatlimit = dxlatlimit-180;
% end
% 
% for i = 1 : num_Int %Discontinuity Intersection
%     if (Trend(i) > oppsxlatlimit && Trend(i) < oppdxlatlimit) || ...
%             (oppdxlatlimit < oppsxlatlimit && Trend(i) > oppsxlatlimit && Trend(i) > oppdxlatlimit) || ...
%             (oppdxlatlimit < oppsxlatlimit && Trend(i) < oppsxlatlimit && Trend(i) < oppdxlatlimit)
%         
%         if Plunge(i) > upperTopplingLimit && Plunge(i) < (90 - friction)
%             DirectToppling(i)=2;
%         elseif Plunge(i) > (90 - friction)
%             DirectToppling(i)=1;
%         else
%             DirectToppling(i)=NaN;
%         end
%         
%     elseif (Plunge(i) > (90 - friction) && strikeslope <= 180 && planeslope(1,1)<= 180 && (Trend(i) < strikeslope || Trend(i)> (strikeslope+180))) || ...
%             (Plunge(i) > (90 - friction) && strikeslope>= 180 && (planeslope(1,1)<= 360 || planeslope(1,1)<= 0) && Trend(i)>(strikeslope-180) && Trend(i)<strikeslope)
%         ObliqueToppling(i)=1;
%         DirectToppling(i)=NaN;
%     else
%         DirectToppling(i)=NaN;
%         ObliqueToppling(i)=NaN;
%         
%     end
%     
% end
% for i = 1 : num_disc %Discontinuity Planes
%     if (pole_DipDir(i) > oppsxlatlimit && pole_DipDir(i) < oppdxlatlimit) || ...
%             (oppdxlatlimit < oppsxlatlimit && pole_DipDir(i) > oppsxlatlimit && pole_DipDir(i) > oppdxlatlimit) || ...
%             (oppdxlatlimit < oppsxlatlimit && pole_DipDir(i) < oppsxlatlimit && pole_DipDir(i) < oppdxlatlimit)
%         if pole_Dip(i) > upperTopplingLimit && pole_Dip(i) < (90 - friction)
%             BasalPlane(i) = 2;
%         elseif pole_Dip(i)>(90-friction)
%             BasalPlane(i) = 1;
%         else
%             BasalPlane(i) = NaN;
%         end
%     elseif (pole_Dip(i) > (90 - friction) && strikeslope <= 180 && planeslope(1,1)<= 180 && (pole_DipDir(i) < strikeslope || pole_DipDir(i)> (strikeslope+180))) || ...
%             (pole_Dip(i) > (90 - friction) && strikeslope>= 180 && (planeslope(1,1)<= 360 || planeslope(1,1)<= 0) && pole_DipDir(i)>(strikeslope-180) && pole_DipDir(i)<strikeslope)
%         BasalPlane(i) = 3;
%     else
%         BasalPlane(i)=NaN;
%     end
% end
% %DirectToppling
% DT_Hcrit=sum(DirectToppling(:) == 2);
% DT_Lcrit=sum(DirectToppling(:) == 1);
% if plot==1
%     Stereogram(pole_Trend, pole_Plunge, num_Int, 'DirectToppling', 1, DirectToppling )
%     hold on
%     title(['DirectToppling, N = ', num2str(num_Int) newline ...
%         'wedge that could be affected by direct toppling and with a sliding release mode (red) =', num2str(DT_Hcrit) newline ...
%         'wedge that could be affected by direct toppling and without a sliding release mode (pink) =', num2str(DT_Lcrit)])
%     hold off
% end
% % saveas(gcf, fullfile(pathFig,'DirectToppling.fig'))
% % saveas(gcf, fullfile(pathFig,'DirectToppling.png'))
% 
% %ObliqueToppling
% OT_crit=sum(ObliqueToppling(:) == 1);
% if plot==1
%     Stereogram(pole_Trend, pole_Plunge, num_Int, 'ObliqueToppling', 1, ObliqueToppling )
%     hold on
%     title(['ObliqueToppling, N = ', num2str(num_Int) newline ...
%         'wedge that could be affected by oblique toppling and without a sliding basal plane (pink) =', num2str(OT_crit)])
%     hold off
% end
% % saveas(gcf, fullfile(pathFig,'ObliqueToppling.fig'))
% % saveas(gcf, fullfile(pathFig,'ObliqueToppling.png'))
% 
% %BasalPlane
% BP_Hcrit=sum(DirectToppling(:) == 2);
% BP_Mcrit=sum(DirectToppling(:) == 1);
% BP_Lcrit=sum(BasalPlane(:) == 3);
% if plot==1
%     Stereogram(DipDir, Dip, num_disc, 'BasalPlane', 0, BasalPlane )
%     hold on
%     title(['BasalPlane, N = ', num2str(num_Int) newline ...
%         'BasalPlane that could be need to activate direct toppling with sliding release mode(red) =', num2str(BP_Hcrit) newline...
%         'BasalPlane that could be need to activate direct toppling without sliding release mode(yellow) =', num2str(BP_Mcrit) newline...
%         'BasalPlane that could be need to activate oblique toppling without sliding release mode(pink) =', num2str(BP_Lcrit)])
%     hold off
% end
% % saveas(gcf, fullfile(pathFig,'BasalPlane.fig'))
% % saveas(gcf, fullfile(pathFig,'BasalPlane.png'))

%% Calculation statistics and display them
%PlanarSliding-Statistic
Set_PS=zeros(nSet+1,5);
critic_PS_Hc=zeros;
critic_PS_Lc=zeros;
for k = 1 : nSet
    Hc=0;
    Lc=0;
    Set_PS(k,1) = Set_name(k);%set number value
    Set_PS(k,2) = sum(Set(:) == (Set_name(k)));%total number of discontinuity of set
    
    for i = 1 : num_disc%value(%) of critical discontinuity
        if Set(i)== Set_name(k) && PlanarSliding(i) == 2
            Hc = Hc + 1;
            critic_PS_Hc=critic_PS_Hc+1;
        elseif Set(i)== Set_name(k) == 1 && PlanarSliding(i) == 1
            Lc=Lc+1;
            critic_PS_Lc=critic_PS_Lc+1;
        end
    end
    Set_PS(k,4) =Hc ;
    Set_PS(k,5) =Lc ;
    Set_PS(k,3) =Set_PS(k,4) + Set_PS(k,5);
end

for i = 1 : numel(Set_PS(1,:))
    if i == 1
        Set_PS(nSet+1,i)=NaN;
    elseif i==2
        Set_PS(nSet+1,i)=sum(Set_PS(:,i));
    elseif i==3
        Set_PS(nSet+1,i)=(critic_PS_Hc+critic_PS_Lc);
    elseif i==4
        Set_PS(nSet+1,i)=(critic_PS_Hc);
    elseif i==5
        Set_PS(nSet+1,i)=(critic_PS_Lc);
    end
end
clearvars PS_T
PS_T=table(Set_PS(:,1),Set_PS(:,2),Set_PS(:,3),Set_PS(:,4),Set_PS(:,5));
PS_T.Properties.VariableNames{'Var1'} = 'Set';
PS_T.Properties.VariableNames{'Var2'} = 'TotDiscontinuities';
PS_T.Properties.VariableNames{'Var3'} = 'TotCritical';
PS_T.Properties.VariableNames{'Var4'} = 'HighCritical';
PS_T.Properties.VariableNames{'Var5'} = 'LowCritical';

% %WedgeSliding-Statistic
% Set_WS=zeros(1,6);
% index_Set_WS=0;
% tot_int_sets=0;
% critic_WS_Hc=0;
% critic_WS_Lc=0;
% for k = 1 : nSet
%     for l = k : nSet
%         Set_WS(index_Set_WS+1,1) = Set_name(k);
%         Set_WS(index_Set_WS+1,2) = Set_name(l);
%         index_Set_WS = numel(Set_WS(:,1));
%         tot_int_sets = 0;
%         Hc = 0;
%         Lc = 0;
%         for i = 1 : num_Int
%             if (Int_Sets(i,1) == Set_name(k) && Int_Sets(i,2) == Set_name(l) )||(Int_Sets(i,2) == Set_name(k) && Int_Sets(i,1) == Set_name(l))
%                 tot_int_sets = tot_int_sets + 1;
%                 if WedgeSliding(i) == 2
%                     Hc = Hc + 1;
%                     critic_WS_Hc=critic_WS_Hc+1;
%                 elseif WedgeSliding(i) == 1
%                     Lc = Lc +1;
%                     critic_WS_Lc=critic_WS_Lc+1;
%                 end
%             end
%         end
%         Set_WS(index_Set_WS,3) = tot_int_sets;
%         Set_WS(index_Set_WS,4) = (Hc + Lc);
%         Set_WS(index_Set_WS,5) = Hc;
%         Set_WS(index_Set_WS,6) = Lc;
%     end
% end
% for i = 1 : numel(Set_WS(1,:))
%     if i ==1 || i==2
%         Set_WS(index_Set_WS+1,i)=NaN;
%     elseif i==3
%         Set_WS(index_Set_WS+1,i)=sum(Set_WS(:,i));
%     elseif i==4
%         Set_WS(index_Set_WS+1,i)=(critic_WS_Hc+critic_WS_Lc);
%     elseif i==5
%         Set_WS(index_Set_WS+1,i)=(critic_WS_Hc);
%     elseif i==6
%         Set_WS(index_Set_WS+1,i)=(critic_WS_Lc);
%     end
% end
% WS_T=table(Set_WS(:,1),Set_WS(:,2),Set_WS(:,3),Set_WS(:,4),Set_WS(:,5),Set_WS(:,6));
% WS_T.Properties.VariableNames{'Var1'} = 'Set_i';
% WS_T.Properties.VariableNames{'Var2'} = 'Set_j';
% WS_T.Properties.VariableNames{'Var3'} = 'TotIntersections';
% WS_T.Properties.VariableNames{'Var4'} = 'TotCritical';
% WS_T.Properties.VariableNames{'Var5'} = 'HighCritical';
% WS_T.Properties.VariableNames{'Var6'} = 'LowCritical';

%FlexuraToppling-Statistic
Set_FT=zeros(nSet+1,5);
critic_FT=zeros;
for k = 1 : nSet
    Hc=0;
    Set_FT(k,1) = Set_name(k);%set number value
    Set_FT(k,2) = sum(Set(:) == (Set_name(k)));%total number of discontinuity of set
    
    for i = 1 : num_disc%value(%) of critical discontinuity
        if Set(i)== Set_name(k) && FlexuralToppling(i) == 2
            Hc = Hc + 1;
            critic_FT=critic_FT+1;
        end
    end
    Set_FT(k,3) =Hc ;
    
end

for i = 1 : numel(Set_FT(1,:))
    if i == 1
        Set_FT(nSet+1,i)=NaN;
    elseif i==2
        Set_FT(nSet+1,i)=sum(Set_FT(:,i));
    else
        Set_FT(nSet+1,i)=critic_FT;
    end
end
clearvars FT_T
FT_T=table(Set_FT(:,1),Set_FT(:,2),Set_FT(:,3));
FT_T.Properties.VariableNames{'Var1'} = 'Set';
FT_T.Properties.VariableNames{'Var2'} = 'TotDiscontinuities';
FT_T.Properties.VariableNames{'Var3'} = 'TotCritical';

% %DirectToppling-Statistic
% Set_DT=zeros(1,6);
% index_Set_DT=0;
% tot_int_sets=0;
% critic_DT=0;
% critic_DT_Hc=0;
% critic_DT_Lc=0;
% for k = 1 : nSet
%     for l = k : nSet
%         Set_DT(index_Set_DT+1,1) = Set_name(k);
%         Set_DT(index_Set_DT+1,2) = Set_name(l);
%         index_Set_DT = numel(Set_DT(:,1));
%         tot_int_sets = 0;
%         Hc = 0;
%         Lc = 0;
%         for i = 1 : num_Int
%             if (Int_Sets(i,1) == Set_name(k) && Int_Sets(i,2) == Set_name(l) )||(Int_Sets(i,2) == Set_name(k) && Int_Sets(i,1) == Set_name(l))
%                 tot_int_sets = tot_int_sets + 1;
%                 if DirectToppling(i) == 2
%                     Hc = Hc + 1;
%                     critic_DT_Hc=critic_DT_Hc+1;
%                 elseif DirectToppling(i) == 1
%                     Lc = Lc +1;
%                     critic_DT_Lc=critic_DT_Lc+1;
%                 end
%             end
%         end
%         Set_DT(index_Set_DT,3) = tot_int_sets;
%         Set_DT(index_Set_DT,4) = (Hc + Lc);
%         Set_DT(index_Set_DT,5) = Hc;
%         Set_DT(index_Set_DT,6) = Lc;
%     end
% end
% for i = 1 : numel(Set_DT(1,:))
%     if i ==1 || i==2
%         Set_DT(index_Set_DT+1,i)=NaN;
%     elseif i==3
%         Set_DT(index_Set_DT+1,i)=sum(Set_DT(:,i));
%     elseif i==4
%         Set_DT(index_Set_DT+1,i)=(critic_DT_Hc+critic_DT_Lc);
%     elseif i==5
%         Set_DT(index_Set_DT+1,i)=critic_DT_Hc;
%     elseif i==6
%         Set_DT(index_Set_DT+1,i)=critic_DT_Lc;
%     end
% end
% DT_T=table(Set_DT(:,1),Set_DT(:,2),Set_DT(:,3),Set_DT(:,4),Set_DT(:,5),Set_DT(:,6));
% DT_T.Properties.VariableNames{'Var1'} = 'Set_i';
% DT_T.Properties.VariableNames{'Var2'} = 'Set_j';
% DT_T.Properties.VariableNames{'Var3'} = 'TotIntersections';
% DT_T.Properties.VariableNames{'Var4'} = 'TotCritical';
% DT_T.Properties.VariableNames{'Var5'} = 'SlidingFailiureMode';
% DT_T.Properties.VariableNames{'Var6'} = 'NoSlidingFailureMode';
% 
% %ObliqueToppling-Statistic
% Set_OT=zeros(1,4);
% index_Set_OT=0;
% tot_int_sets=0;
% critic_OT=0;
% for k = 1 : nSet
%     for l = k : nSet
%         Set_OT(index_Set_OT+1,1) = Set_name(k);
%         Set_OT(index_Set_OT+1,2) = Set_name(l);
%         index_Set_OT = numel(Set_OT(:,1));
%         tot_int_sets = 0;
%         Hc = 0;
%         for i = 1 : num_Int
%             if (Int_Sets(i,1) == Set_name(k) && Int_Sets(i,2) == Set_name(l) )||(Int_Sets(i,2) == Set_name(k) && Int_Sets(i,1) == Set_name(l))
%                 tot_int_sets = tot_int_sets + 1;
%                 if ObliqueToppling(i) == 1
%                     Hc = Hc + 1;
%                     critic_OT=critic_OT+1;
%                     
%                 end
%             end
%         end
%         Set_OT(index_Set_OT,3) = tot_int_sets;
%         Set_OT(index_Set_OT,4) = Hc;
%     end
% end
% for i = 1 : numel(Set_OT(1,:))
%     if i ==1 || i==2
%         Set_OT(index_Set_OT+1,i)=NaN;
%     elseif i==3
%         Set_OT(index_Set_OT+1,i)=sum(Set_OT(:,i));
%     elseif i==4
%         Set_OT(index_Set_OT+1,i)=critic_OT;
%         
%     end
% end
% OT_T=table(Set_OT(:,1),Set_OT(:,2),Set_OT(:,3),Set_OT(:,4));
% OT_T.Properties.VariableNames{'Var1'} = 'Set_i';
% OT_T.Properties.VariableNames{'Var2'} = 'Set_j';
% OT_T.Properties.VariableNames{'Var3'} = 'TotIntersections';
% OT_T.Properties.VariableNames{'Var4'} = 'TotCritical';
% 
% %BasalPlane-Statistic
% Set_BP=zeros(nSet+1,5);
% critic_BP_Hc=0;
% critic_BP_Lc=0;
% for k = 1 : nSet
%     Hc=0;
%     Lc=0;
%     Set_BP(k,1) = Set_name(k);%set number value
%     Set_BP(k,2) = sum(Set(:) == (Set_name(k)));%total number of discontinuity of set
%     
%     for i = 1 : num_disc%value(%) of critical discontinuity
%         if Set(i)== Set_name(k) && BasalPlane(i) == 2
%             Hc = Hc + 1;
%             critic_BP_Hc=critic_BP_Hc+1;
%         elseif Set(i)== Set_name(k) == 1 && (BasalPlane(i) == 1 ||BasalPlane(i) ==3)
%             Lc=Lc+1;
%             critic_BP_Lc=critic_BP_Lc+1;
%         end
%     end
%     Set_BP(k,4) =Hc ;
%     Set_BP(k,5) =Lc ;
%     Set_BP(k,3) =Set_BP(k,4) + Set_BP(k,3);
% end
% 
% for i = 1 : numel(Set_BP(1,:))
%     if i == 1
%         Set_BP(nSet+1,i)=NaN;
%     elseif i==2
%         Set_BP(nSet+1,i)=sum(Set_BP(:,i));
%     elseif i==3
%         Set_BP(nSet+1,i)=(critic_BP_Hc+critic_BP_Lc);
%     elseif i==4
%         Set_BP(nSet+1,i)=critic_BP_Hc;
%     elseif i==5
%         Set_BP(nSet+1,i)=critic_BP_Lc;
%     end
% end
% 
% BP_T=table(Set_BP(:,1),Set_BP(:,2),Set_BP(:,3),Set_BP(:,4),Set_BP(:,5));
% BP_T.Properties.VariableNames{'Var1'} = 'Set';
% BP_T.Properties.VariableNames{'Var2'} = 'TotDiscontinuities';
% BP_T.Properties.VariableNames{'Var3'} = 'TotBasalPlane';
% BP_T.Properties.VariableNames{'Var4'} = 'SlidingBasalPlane';
% BP_T.Properties.VariableNames{'Var5'} = 'NoSlidingBasalPlane';

% % % % % % disp('PlanarSlidingKynematicAnalysis')
% % % % % % disp(PS_T)
% disp(' ')
% disp('WedgeSlidingKynematicAnalysis')
% disp(WS_T)
% disp(' ')
% % % % % % disp('FlexuralTopplingKynematicAnalysis')
% % % % % % disp(FT_T)
% disp(' ')
% disp('DirectTopplingKynematicAnalysis')
% disp(DT_T)
% disp(' ')
% disp('ObliqueTopplingKynematicAnalysis')
% disp(OT_T)
% disp(' ')
% disp('BasalPlaneOfTopplingKynematicAnalysis')
% disp(BP_T)
% % disp('CriticalPercentual')
% % disp('%PlanarSliding %FlexuralToppling')
disp([num2str(planeslope(1,1)),' ',num2str(planeslope(1,2)),' ',num2str(PS_T.TotCritical(end)/PS_T.TotDiscontinuities(end)),' ', num2str(FT_T.TotCritical(end)/FT_T.TotDiscontinuities(end))])
% 
% disp('%PlanarSliding %WedgeSliding %FlexuralToppling %DirectToppling %ObliqueToppling')
% disp([num2str(PS_T.TotCritical(end)/PS_T.TotDiscontinuities(end)),' ',...
%     num2str(WS_T.TotCritical(end)/WS_T.TotIntersections(end)),' ',...
%     num2str(FT_T.TotCritical(end)/FT_T.TotDiscontinuities(end)),' ',...
%     num2str(DT_T.TotCritical(end)/DT_T.TotIntersections(end)),' ',...
%     num2str(OT_T.TotCritical(end)/OT_T.TotIntersections(end))])
end
end
