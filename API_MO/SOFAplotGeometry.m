function SOFAplotGeometry(Obj, index)
% SOFAplotGeometry(Obj) plots the geometry found in the Obj.
%
% SOFAplotGeometry(Obj, index) plots the geometry for the measurements
% given in the index.

% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.1 or ? as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License.

global legendEntries
global legendStrings

switch Obj.GLOBAL_SOFAConventions
    %%
    case {'SimpleFreeFieldHRIR','SingleRoomDRIR','SimpleFreeFieldTF'}
        
        legendStrings = {'ListenerPosition','ListenerView','Receivers','SourcePosition'};
        
        if ~exist('index','var')
            index=1:Obj.API.M;
        end
        
        % Expand entries to the same number of measurement points
        Obj = SOFAexpand(Obj);
        % See if the room geometry is specified
        if strcmp(Obj.GLOBAL_RoomType,'shoebox')
            figure('Position',[1 1 (Obj.RoomCornerB(1)-Obj.RoomCornerA(1))*1.2 Obj.RoomCornerB(2)-Obj.RoomCornerA(2)]*100);
            box on; hold on; h=[];
            % plot the room
            rectangle('Position',[Obj.RoomCornerA(1) ...
                Obj.RoomCornerA(2) ...
                Obj.RoomCornerB(1)-Obj.RoomCornerA(1) ...
                Obj.RoomCornerB(2)-Obj.RoomCornerA(2)]);
        else
            figure; hold on;
        end
        legendEntries = [];
        title(sprintf('%s, %s',Obj.GLOBAL_SOFAConventions,Obj.GLOBAL_RoomType));
        % Get ListenerPosition, ListenerView, ReceiverPosition, and SourcePosition
        % NOTE: ListenerPosition is set to [0 0 0] for SimpleFreeFieldHRIR
        LP = SOFAconvertCoordinates(Obj.ListenerPosition(index,:),Obj.ListenerPosition_Type,'cartesian');
        LV = SOFAconvertCoordinates(Obj.ListenerView(index,:),Obj.ListenerView_Type,'cartesian');
        RP = SOFAconvertCoordinates(Obj.ReceiverPosition(:,:,index),Obj.ReceiverPosition_Type,'cartesian');
        S  = SOFAconvertCoordinates(Obj.SourcePosition(index,:),Obj.SourcePosition_Type,'cartesian');
        % Use only unique listener and source positions
        uniquePoints = unique([LP LV S],'rows','stable');
        LP = uniquePoints(:,1:3);
        LV = uniquePoints(:,4:6);
        S  = uniquePoints(:,7:9);
        % Plot ListenerPosition
        legendEntries(end+1) = plot3(LP(:,1),LP(:,2),LP(:,3),'ro','MarkerFaceColor',[1 0 0]);
        % Plot ListenerView
        for ii=1:size(LV,1)
            % Scale size of ListenerView vector smaller
            LV(ii,:) = 0.2*LV(ii,:)./norm(LV(ii,:));
            % Plot line for ListenerView vector
            line([LP(ii,1), LV(ii,1)+LP(ii,1)], [LP(ii,2) LV(ii,2)+LP(ii,2)], 'Color',[1 0 0]);
        end
        legendEntries(end+1) = plot3(LV(:,1),LV(:,2),LV(:,3),'ro','MarkerFaceColor',[1 1 1]);
        % Plot ReceiverPositon (this is plotted only for the first ListenerPosition)
        if ndims(RP)>2
            % If ReceiverPositon has more than two dimesnions reduce it to the first
            % ListenerPosition
            RP = shiftdim(RP,2);
            RP = squeeze(RP(1,:,:));
        end
        legendEntries(end+1) = plot3(LP(1,1)+RP(1,1), LP(1,2)+RP(1,2), LP(1,3)+RP(1,3),'rx');
        for ii=2:size(RP,1)
            plot3(LP(1,1)+RP(ii,1), LP(1,2)+RP(ii,2), LP(1,3)+RP(ii,3),'rx');
        end
        % Plot SourcePosition
        %if ~isfield(Obj,'GLOBAL_OriginalSOFAlabels')
        legendEntries(end+1)=plot3(S(:,1),S(:,2),S(:,3),'k.');
        %     else
        %         labels = Obj.GLOBAL_OriginalSOFAlabels;
        %         legendStrings = legendStrings(1:end-1);
        %         markers = {'k.','m.','g.','r.','c.','y.','k+','ko','k>','k<','kp','kh'};
        %         list = flipud(unique(labels(index,:),'rows','stable'));
        %         labels_cell = cellstr(labels);
        %         for k = 1:size(list,1)
        %             source_position_index = all(ismember(labels_cell,strtrim(list(k,:))),2);lidx = find(source_position_index);
        %             legendEntries(end+1)=plot3(S(lidx,1),S(lidx,2),S(lidx,3),markers{k});
        %             legendStrings = [legendStrings, ['SourcePosition ' strtrim(list(k,:))]];
        %         end
        %     end
        
        legend(legendEntries,legendStrings,'Location','NorthEastOutside');
        xlabel(['X / ' Obj.ListenerPosition_Units]);
        ylabel(['Y / ' Obj.ListenerPosition_Units]);
        zlabel(['Z / ' Obj.ListenerPosition_Units]);
        % Set fixed aspect ratio
        axis equal;
        
    case {'AmbisonicsDRIR'}
        
        legendStrings = {'ListenerPosition','ListenerView','EmitterPosition'};
        
        if ~exist('index','var')
            index=1:Obj.API.E;
        end
        
        
        % Expand entries to the same number of measurement points
        %Obj = SOFAexpand(Obj);
        % See if the room geometry is specified
        if strcmp(Obj.GLOBAL_RoomType,'shoebox')
            figure('Position',[1 1 (Obj.RoomCornerB(1)-Obj.RoomCornerA(1))*1.2 Obj.RoomCornerB(2)-Obj.RoomCornerA(2)]*100);
            box on; hold on; h=[];
            % plot the room
            rectangle('Position',[Obj.RoomCornerA(1) ...
                Obj.RoomCornerA(2) ...
                Obj.RoomCornerB(1)-Obj.RoomCornerA(1) ...
                Obj.RoomCornerB(2)-Obj.RoomCornerA(2)]);
        else
            figure; hold on;
        end
        legendEntries = [];
        title(sprintf('%s, %s',Obj.GLOBAL_SOFAConventions,Obj.GLOBAL_RoomType));
        
        LP = SOFAconvertCoordinates(Obj.ListenerPosition,Obj.ListenerPosition_Type,'cartesian');
        LV = SOFAconvertCoordinates(Obj.ListenerView,Obj.ListenerView_Type,'cartesian');
        SP = SOFAconvertCoordinates(Obj.SourcePosition,Obj.SourcePosition_Type,'cartesian');
        
        if size(LP,1) == 1
            E = SOFAconvertCoordinates(Obj.EmitterPosition(index,:),Obj.EmitterPosition_Type,'cartesian');
        else
            if ndims(squeeze(Obj.EmitterPosition)) == 3
                E  = SOFAconvertCoordinates(transpose(Obj.EmitterPosition(:,index)),Obj.EmitterPosition_Type,'cartesian');
            else
                E  = SOFAconvertCoordinates(Obj.EmitterPosition(:,:,index),Obj.EmitterPosition_Type,'cartesian');
            end
        end
        
        % Plot ListenerPosition
        legendEntries(end+1) = plot3(LP(:,1),LP(:,2),LP(:,3),'ro','MarkerFaceColor',[1 0 0]);
        % Plot ListenerView
        for ii=1:size(LV,1)
            % Scale size of ListenerView vector smaller
            LV(ii,:) = 1*LV(ii,:)./norm(LV(ii,:));
            % Plot line for ListenerView vector
            line([LP(ii,1), LV(ii,1)+LP(ii,1)], [LP(ii,2) LV(ii,2)+LP(ii,2)], [LP(ii,3), LV(ii,3)+LP(ii,3)], 'Color',[1 0 0]);
        end
        legendEntries(end+1) = plot3(LV(:,1)+LP(:,1),LV(:,2)+LP(:,2),LV(:,3)+LP(:,3),'ro','MarkerFaceColor',[1 1 1]);
        
        % Plot EmitterPosition
        legendEntries(end+1)=scatter3(E(:,1)+SP(:,1),E(:,2)+SP(:,2),E(:,3)+SP(:,3),'filled');
        
        %Plot EmitterView if any
        if isfield(Obj,'EmitterView')
            EV = SOFAconvertCoordinates(Obj.EmitterView,Obj.EmitterView_Type,'cartesian');
            EP = SOFAconvertCoordinates(Obj.EmitterPosition,Obj.EmitterPosition_Type,'cartesian');
            for ii=1:size(EV,1)
                % Scale size of EmitterView vector smaller
                EV(ii,:) = 1*EV(ii,:)./norm(EV(ii,:));
                % Plot line for EmitterView vector
                line([EP(ii,1), EV(ii,1)+EP(ii,1)], [EP(ii,2) EV(ii,2)+EP(ii,2)], [EP(ii,3), EV(ii,3)+EP(ii,3)]);
            end
            legendStrings = [legendStrings,{'EmitterView'}];
            legendEntries(end+1) = plot3(EV(:,1)+EP(:,1),EV(:,2)+EP(:,2),EV(:,3)+EP(:,3),'bo','MarkerFaceColor',[1 1 1]);
        end
        
        
        legend(legendEntries,legendStrings,'Location','NorthEastOutside');
        xlabel(['X / ' Obj.ListenerPosition_Units]);
        ylabel(['Y / ' Obj.ListenerPosition_Units]);
        zlabel(['Z / ' Obj.ListenerPosition_Units]);
        
        view(-30,25);
        % Set fixed aspect ratio
        axis equal;
        
        if strcmp(Obj.GLOBAL_RoomType,'shoebox')            
            xlim([0 Obj.RoomCornerB(1)]);
            ylim([0 Obj.RoomCornerB(2)]);
            zlim([0 Obj.RoomCornerB(3)]);   
        end
        
        
    otherwise
        error('This SOFAConventions is not supported for plotting');
end


% Add a little bit extra space at the axis
if ~strcmp(Obj.GLOBAL_RoomType,'shoebox')
    
    axisLimits = axis();
    paddingSpace = 0.2 * max(abs(axisLimits(:)));
    axisLimits([1 3]) = axisLimits([1 3]) - paddingSpace;
    axisLimits([2 4]) = axisLimits([2 4]) + paddingSpace;
    axis(axisLimits);
    
end
grid on;
