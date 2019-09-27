%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   function index = SOFAfindPosition(Obj,azimuth_target,elevation_target, distance_target)
%
%   This function outputs the index of a HRIR in a SOFA file (of type
%   'SimpleFreeFieldHRIR') that is the closest to a desired incoming
%   direction described in spherical coordinates.
%
%   Note: if azimuth_target is an empty array, SOFAfindPosition() returns
%   the index of all the HRIRs whose elevation angle is the closest
%   available to the desired elevation angle (iso-elevation set of
%   positions, no matter of the azimuth angle), and vice-versa with
%   elevation_traget being an empty array.
%
%   Julien De Muynke, 05/02/2019
%   julien.demuynke@eurecat.org
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function index = SOFAfindPosition(Obj,azimuth_target, elevation_target, distance_target)

index = [];

if strcmp(Obj.GLOBAL_SOFAConventions,'AmbisonicsDRIR')
    position = Obj.EmitterPosition;
    position_Type = Obj.EmitterPosition_Type;
    position_Units = Obj.EmitterPosition_Units;
else
    position = Obj.SourcePosition;
    position_Type = Obj.SourcePosition_Type;
    position_Units = Obj.SourcePosition_Units;
end

if nargin < 4
    
    % Force coordinates system to spherical, force angles to degree
    
    if strcmp(position_Type, 'cartesian')
        coord_sph = SOFAconvertCoordinates(position,'cartesian','spherical');
    else
        coord_sph = position;
        if contains(position_Units,'radian')
            coord_sph(:,1) = coord_sph(:,1)*180/pi;
            coord_sph(:,2) = coord_sph(:,2)*180/pi;
        end
    end
    
    % Force azimuth angle between -179º and 180º
    
    coord_sph(find(coord_sph(:,1)>180),1) = coord_sph(find(coord_sph(:,1)>180),1)-360;
    
    
    % For Sfear 2.x, only the 'if' section has to be implemented.
    % The closest position is searched for recursively, consequently the
    % distance to the desired position is initialized as the diameter of the
    % sphere that crosses the furthest measured position throughout the SOFA
    % file.
    % Note that the distance between the current and the desired positions is
    % calculated as the length of the chord of the sphere that crosses the
    % current position (since the desired position is also forced to belong to
    % the sphere's surface)
    
    if ~isempty(azimuth_target)
        
        if ~isempty(elevation_target)
            
            sphere_radius = max(unique(coord_sph(:,3)));
            dist_min = 2*sphere_radius;
            
            for k = 1:length(azimuth_target)
                
                dist = zeros(1,Obj.API.M);
                for i = 1:Obj.API.M
                    
                    az_diff = abs(coord_sph(i,1) - azimuth_target(k));
                    if az_diff > 180
                        az_diff = 360 - az_diff;
                    end
                    elev_diff = abs(coord_sph(i,2) - elevation_target);
                    
                    dist(i) = 2*sphere_radius*sqrt(1-(cosd(az_diff)+cosd(elev_diff))/2);
                    
                    %                 if dist < dist_min
                    %                     dist_min = dist;
                    %                     target_index = i;
                    %                 end
                    
                end
                [~,target_index] = find(dist == min(dist));
                index = [index target_index];
            end
            
        else
            
            azim_diff = abs(round(coord_sph(:,1),1) - azimuth_target);
            [mini,first_index] = min(azim_diff);
            for i = 1:Obj.API.M
                if(azim_diff(i)) == mini
                    index = [index i];
                end
            end
            
        end
        
    else
        
        elev_diff = abs(round(coord_sph(:,2),1) - elevation_target);
        [mini,first_index] = min(elev_diff);
        for i = 1:Obj.API.M
            if(elev_diff(i)) == mini
                index = [index i];
            end
        end
        
        
    end
    
else
    
    [x_target, y_target, z_target] = sph2cart(azimuth_target*pi/180, elevation_target*pi/180, distance_target);
    
    if strcmp(position_Type, 'spherical')
        dist_min = max(position(:,3));
        if contains(position_Units,'radian')
            position_degree = [position(:,1)*180/pi position(:,2)*180/pi position(:,3)];
        else
            position_degree = position;
        end
        coord_cartesian = SOFAconvertCoordinates(position_degree,'spherical','cartesian');
    else
        dist_min = max(sqrt(position(:,1).^2+position(:,2).^2+position(:,3).^2));
        coord_cartesian = position;
    end
    
    for i = 1:Obj.API.M
        
        dist = sqrt((x_target-coord_cartesian(i,1))^2+(y_target-coord_cartesian(i,2))^2+(z_target-coord_cartesian(i,3))^2);
        
        if dist < dist_min
            index = i;
            dist_min = dist;
        elseif abs(dist - dist_min) <= 1e-3
            index = [index i];
        end
        
    end
    
end