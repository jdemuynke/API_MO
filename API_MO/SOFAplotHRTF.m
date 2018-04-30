function M=SOFAplotHRTF(Obj,type,ch)
% SOFAplotHRTF(OBJ, TYPE, CH) plots the CH channel of HRTFs given in OBJ.
%  The following TYPEs are supported:
%  'EtcHorizontal'  energy-time curve in the horizontal plane (+/- 5 deg)
%  'EtcMedian'      energy-time curve in the median plane (+/- 2 deg)
%  'MagMedian'      magnitude spectrum in the median plane (+/- 2 deg)
%
%  OBJ must be in SimpleFreeFieldHRIR or SimpleFreeFieldSOS.
%
% M=SOFAplotHRTF... returns the matrix M displayed in the figure.
%

% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.1 or ? as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License.


if ~exist('ch','var')
    ch=1;
end

%% Convert data to FIR
switch Obj.GLOBAL_SOFAConventions
    case 'SimpleFreeFieldSOS'
        N=512;
        impulse=[1; zeros(N-1,1)];
        T=SOFAgetConventions('SimpleFreeFieldHRIR');
        Obj.GLOBAL_SOFAConventions=T.GLOBAL_SOFAConventions;
        Obj.GLOBAL_SOFAConventionsVersion=T.GLOBAL_SOFAConventionsVersion;
        Obj.GLOBAL_DataType=T.GLOBAL_DataType;
        Obj.API.Dimensions.Data.IR=Obj.API.Dimensions.Data.SOS;
        Obj.Data.IR=zeros(Obj.API.M, Obj.API.R, N);
        for ii=1:Obj.API.M
            for jj=1:Obj.API.R
                Obj.Data.IR(ii,jj,:)=sosfilt(reshape(squeeze(Obj.Data.SOS(ii,jj,:)),6,[])',impulse);
            end
        end
        Obj.Data=rmfield(Obj.Data,'SOS');
        Obj.API.Dimensions.Data=rmfield(Obj.API.Dimensions.Data,'SOS');
        Obj=SOFAupdateDimensions(Obj);
    case {'SimpleFreeFieldTF' 'GeneralTF'}
        fs=max(Obj.N)*2;
        N=fs/min([min(diff(Obj.N)) Obj.N(1)]);
        N=2*(round(N/2+1)-1);
        T=SOFAgetConventions('SimpleFreeFieldHRIR');
        Obj.GLOBAL_SOFAConventions=T.GLOBAL_SOFAConventions;
        Obj.GLOBAL_SOFAConventionsVersion=T.GLOBAL_SOFAConventionsVersion;
        Obj.GLOBAL_DataType=T.GLOBAL_DataType;
        Obj.API.Dimensions.Data.IR=Obj.API.Dimensions.Data.Real;
        Obj.Data.SamplingRate=fs;
        Obj.Data.SamplingRate_Units='Hertz';
        Obj.Data.IR_LongName=Obj.Data.Real_LongName;
        Obj.Data.IR_Units=Obj.Data.Real_Units;
        Obj.Data.IR=zeros(Obj.API.M, Obj.API.R, N);
        for ii=1:Obj.API.M
            for jj=1:Obj.API.R
                s=zeros(N/2+1,1);
                s(Obj.N*N/fs+1)=squeeze(Obj.Data.Real(ii,jj,:))+1i*squeeze(Obj.Data.Imag(ii,jj,:));
                Obj.Data.IR(ii,jj,:)=myifftreal(s,N);
            end
            Obj.SourcePosition(ii,:)=SOFAconvertCoordinates(Obj.SourcePosition(ii,:),Obj.SourcePosition_Type,T.SourcePosition_Type,Obj.SourcePosition_Units,T.SourcePosition_Units);
        end
        Obj.Data.Delay=zeros(1,Obj.API.R);
        Obj.SourcePosition_Type=T.SourcePosition_Type;
        Obj.SourcePosition_Units=T.SourcePosition_Units;
        Obj=rmfield(Obj,{'N','N_LongName','N_Units'});
        Obj.Data=rmfield(Obj.Data,{'Real','Imag','Real_LongName','Imag_LongName','Real_Units','Imag_Units'});
        Obj.API.Dimensions.Data=rmfield(Obj.API.Dimensions.Data,{'Real','Imag'});
        Obj=SOFAupdateDimensions(Obj);
    case 'SimpleFreeFieldHRIR'
        T=SOFAgetConventions('SimpleFreeFieldHRIR');
        Obj.GLOBAL_SOFAConventions=T.GLOBAL_SOFAConventions;
        Obj.GLOBAL_SOFAConventionsVersion=T.GLOBAL_SOFAConventionsVersion;
        Obj.GLOBAL_DataType=T.GLOBAL_DataType;
        %Obj.API.Dimensions.Data.IR=Obj.API.Dimensions.Data.Real;
        %Obj.Data.SamplingRate=fs;
        %Obj.Data.SamplingRate_Units='Hertz';
        %Obj.Data.IR_LongName=Obj.Data.Real_LongName;
        %Obj.Data.IR_Units=Obj.Data.Real_Units;
        %     Obj.Data.IR=zeros(Obj.API.M, Obj.API.R, N);
        %     for ii=1:Obj.API.M
        %       for jj=1:Obj.API.R
        %         s=zeros(N/2+1,1);
        %         s(Obj.N*N/fs+1)=squeeze(Obj.Data.Real(ii,jj,:))+1i*squeeze(Obj.Data.Imag(ii,jj,:));
        %         Obj.Data.IR(ii,jj,:)=myifftreal(s,N);
        %       end
        %       Obj.SourcePosition(ii,:)=SOFAconvertCoordinates(Obj.SourcePosition(ii,:),Obj.SourcePosition_Type,T.SourcePosition_Type,Obj.SourcePosition_Units,T.SourcePosition_Units);
        %     end
        %     Obj.Data.Delay=zeros(1,Obj.API.R);
        if strcmp(Obj.SourcePosition_Type,'cartesian')
            for i = 1:Obj.API.M
                [azim_rad,elev_rad,dist] = cart2sph(Obj.SourcePosition(i,1),Obj.SourcePosition(i,2),Obj.SourcePosition(i,3));
                azim = rad2deg(azim_rad);
                elev= rad2deg(elev_rad);
                Obj.SourcePosition(i,1) = azim;
                Obj.SourcePosition(i,2) = elev;
                Obj.SourcePosition(i,3) = dist;
            end
        elseif ~strcmp(Obj.SourcePosition_Type,'spherical')
            error('Unknown coordinates system!');
        end
        Obj.SourcePosition(:,1) = mod(Obj.SourcePosition(:,1)+360,360);
        
        Obj.SourcePosition_Type=T.SourcePosition_Type;
        Obj.SourcePosition_Units=T.SourcePosition_Units;
        
        %     Obj=rmfield(Obj,{'N','N_LongName','N_Units'});
        %     Obj.Data=rmfield(Obj.Data,{'Real','Imag','Real_LongName','Imag_LongName','Real_Units','Imag_Units'});
        %     Obj.API.Dimensions.Data=rmfield(Obj.API.Dimensions.Data,{'Real','Imag'});
        %     Obj=SOFAupdateDimensions(Obj);
    otherwise
        error('Conventions not supported');
end

fs = Obj.Data.SamplingRate;

%% Plot according to the type
switch lower(type)
    % Energy-time curve (ETC) in the horizontal plane
    case 'etchorizontal'
        noisefloor=-50;
        ele=0;
        thr=5;
        Obj=SOFAexpand(Obj,'Data.Delay');
        hM=double(squeeze(Obj.Data.IR(:,ch,:)));
        pos=Obj.SourcePosition;
        pos(pos(:,1)>180,1)=pos(pos(:,1)>180,1)-360;
        idx=find(pos(:,2)<(ele+thr) & pos(:,2)>(ele-thr));
        M=(20*log10(abs(hM(idx,:))));
        pos=pos(idx,:);
        del=round(Obj.Data.Delay(idx,ch));
        M2=noisefloor*ones(size(M)+[0 max(del)]);
        for ii=1:size(M,1)
            M2(ii,del(ii)+(1:Obj.API.N))=M(ii,:);
        end
        [azi,i]=sort(pos(:,1));
        M=M2(i,:);
        M=M-max(max(M));
        M(M<=noisefloor)=noisefloor;
        surface(0:1/fs*1000:(size(M,2)-1)/fs*1000,azi,M(:,:));
        set(gca,'FontName','Arial','FontSize',10);
        set(gca, 'TickLength', [0.02 0.05]);
        set(gca,'LineWidth',1);
        cmap=colormap(hot);
        cmap=flipud(cmap);
        shading flat
        colormap(cmap);
        box on;
        colorbar;
        xlabel('Time (ms)');
        ylabel('Azimuth (deg)');
        title([Obj.GLOBAL_Title '; channel: ' num2str(ch)],'Interpreter','none');
        
        % Magnitude spectrum in the median plane
    case 'magmedian'
        noisefloor=-50;
        azi=0;
        thr=2;
        hM=double(squeeze(Obj.Data.IR(:,ch,:)));
        pos=Obj.SourcePosition;
        idx=find(abs(pos(:,1))>90);
        pos(idx,2)=180-pos(idx,2);
        pos(idx,1)=180-pos(idx,1);
        idx=find(pos(:,1)<(azi+thr) & pos(:,1)>(azi-thr));
        M=(20*log10(abs(fft(hM(idx,:)')')));
        M=M(:,1:floor(size(M,2)/2));  % only positive frequencies
        pos=pos(idx,:);
        M=M-max(max(M));
        M(M<noisefloor)=noisefloor;
        [ele,i]=sort(pos(:,2));
        M=M(i,:);
        surface(0:fs/size(hM,2):(size(M,2)-1)*fs/size(hM,2),ele,M(:,:));
        shading flat
        xlabel('Frequency (Hz)');
        ylabel('Elevation (deg)');
        title([Obj.GLOBAL_Title '; channel: ' num2str(ch)],'Interpreter','none');
        
        % Magnitude spectrum in the median plane
    case 'maghorizontal'
        
        if ~strcmp(Obj.GLOBAL_SOFAConventions,'SimpleFreeFieldHRIR')
            error('This SOFA Convention is currently not supported for this type of plot!\n');
        end
        
        noisefloor = -50;
        elev = 0;
        thr = 4;
        
        while(1)
            
            hM1 = double(squeeze(Obj.Data.IR(:,1,:)));
            pos = Obj.SourcePosition;
            
            idx = find(pos(:,2)<(elev+thr) & pos(:,2)>(elev-thr));
            
            % Check if there are several positions with the same azimuth
            if length(unique(Obj.SourcePosition(idx,1))) ~= length(Obj.SourcePosition(idx,1))
                fprintf('Doublons in azimuth have been detected for an elevation threshold of %.2f degrees.\nTrying again with a threshold of %.2f degrees\n',thr,thr/2);
                thr = thr/2;
                continue;
            else
                hM2 = double(squeeze(Obj.Data.IR(:,2,:)));
                break;
            end
            
        end
        
        M1 = (20*log10(abs(fft(hM1(idx,:)')')));
        M1 = M1(:,1:floor(size(M1,2)/2));  % only positive frequencies
        M2 = (20*log10(abs(fft(hM2(idx,:)')')));
        M2 = M2(:,1:floor(size(M2,2)/2));  % only positive frequencies
        freq = (1/size(M1,2):1/size(M1,2):1).*fs/2;
        
        maxi = max(max(max(M1)),max(max(M2)));
        
        M1 = M1 - maxi;
        M1(M1 < noisefloor) = noisefloor;
        M2 = M2 - maxi;
        M2(M2 < noisefloor) = noisefloor;
        
        pos = pos(idx,:);
        [azi,i] = sort(pos(:,1));
        
        M1 = M1(i,:);
        M2 = M2(i,:);
        subplot(211),
        surface(azi,freq,transpose(M1(:,:)));
        shading interp;
        ylabel('Frequency (Hz)');
        ylim([20 16000]);
        xlabel('Azimuth (deg)');
        xlim([0 360]);
        view(-25,45);
        grid on;
        title([Obj.GLOBAL_ListenerShortName '; channel: 1'],'Interpreter','none');
        subplot(212),
        surface(azi,freq,transpose(M2(:,:)));
        shading interp;
        ylabel('Frequency (Hz)');
        ylim([20 16000]);
        xlabel('Azimuth (deg)');
        xlim([0 360]);
        view(25,45);
        grid on;
        title([Obj.GLOBAL_ListenerShortName '; channel: 2'],'Interpreter','none');
        
    case 'isd'
        
        if ~strcmp(Obj.GLOBAL_SOFAConventions,'SimpleFreeFieldHRIR')
            error('This SOFA Convention is currently not supported for this type of plot!\n');
        end
        noisefloor = -50;
        elev = 0;
        thr = 4;
        while(1)
            hM1=double(squeeze(Obj.Data.IR(:,1,:)));
            pos=Obj.SourcePosition;
            
            idx=find(pos(:,2)<(elev+thr) & pos(:,2)>(elev-thr));
            
            % Check if there are several positions with the same azimuth
            if length(unique(Obj.SourcePosition(idx,1))) ~= length(Obj.SourcePosition(idx,1))
                fprintf('Doublons in azimuth have been detected for an elevation threshold of %.2f degrees.\nTrying again with a threshold of %.2f degrees\n',thr,thr/2);
                thr = thr/2;
                continue;
            else
                hM2=double(squeeze(Obj.Data.IR(:,2,:)));
                break;
            end
        end
        
        M1=(20*log10(abs(fft(hM1(idx,:)')')));
        M1=M1(:,1:floor(size(M1,2)/2));  % only positive frequencies
        M2=(20*log10(abs(fft(hM2(idx,:)')')));
        M2=M2(:,1:floor(size(M2,2)/2));  % only positive frequencies
        freq = (1/size(M1,2):1/size(M1,2):1).*fs/2;
        
        maxi = max(max(max(M1)),max(max(M2)));
        
        M1=M1-maxi;
        M1(M1<noisefloor)=noisefloor;
        M2=M2-maxi;
        M2(M2<noisefloor)=noisefloor;
        
        pos=pos(idx,:);
        [azi,i]=sort(pos(:,1));
        
        M1=M1(i,:);
        M2=M2(i,:);
        
        subplot(211),
        semilogx(freq,transpose(M1(find(azi == 90),:)-M1(find(azi == 270),:)));
        xlim([freq(1) 20000]);
        grid on;
        title('Interaural Spectral Difference, left ear, source at 90deg azimuth, 0deg elevation');
        xlabel('Frequency (Hz)');
        ylabel('Magnitude (dB)');
        ylim([-20 50]);
        
        
        subplot(212),
        semilogx(freq,transpose(M2(find(azi == 270),:)-M2(find(azi == 90),:)));
        xlim([freq(1) 20000]);
        grid on;
        title('Interaural Spectral Difference, right ear, source at 270deg azimuth, 0deg elevation');
        xlabel('Frequency (Hz)');
        ylabel('Magnitude (dB)');
        ylim([-20 50]);
        
    case 'ildhorizontal'
        
        if ~strcmp(Obj.GLOBAL_SOFAConventions,'SimpleFreeFieldHRIR')
            error('This SOFA Convention is currently not supported for this type of plot!\n');
        end
        noisefloor = -50;
        elev = 0;
        thr = 4;
        
        while(1)
            
            hM1 = double(squeeze(Obj.Data.IR(:,1,:)));
            pos = Obj.SourcePosition;
            
            idx = find(pos(:,2)<(elev+thr) & pos(:,2)>(elev-thr));
            
            % Check if there are several positions with the same azimuth
            if length(unique(Obj.SourcePosition(idx,1))) ~= length(Obj.SourcePosition(idx,1))
                fprintf('Doublons in azimuth have been detected for an elevation threshold of %.2f degrees.\nTrying again with a threshold of %.2f degrees\n',thr,thr/2);
                thr = thr/2;
                continue;
            else
                hM2 = double(squeeze(Obj.Data.IR(:,2,:)));
                break;
            end
            
        end
        
        M1 = (20*log10(abs(fft(hM1(idx,:)')')));
        M1 = M1(:,1:floor(size(M1,2)/2));  % only positive frequencies
        M2 = (20*log10(abs(fft(hM2(idx,:)')')));
        M2 = M2(:,1:floor(size(M2,2)/2));  % only positive frequencies
        freq = (1/size(M1,2):1/size(M1,2):1).*fs/2;
        
        maxi = max(max(max(M1)),max(max(M2)));
        
        M1 = M1 - maxi;
        M1(M1 < noisefloor) = noisefloor;
        M2 = M2 - maxi;
        M2(M2 < noisefloor) = noisefloor;
        
        pos = pos(idx,:);
        [azi,i] = sort(pos(:,1));
        
        M1 = M1(i,:);
        M2 = M2(i,:);
        
        ildhorizontal = zeros(size(azi,1),size(M1,2));
        
        for k = 1:size(azi,1)
            
            %I = find(azi == mod(360-azi(k),360));
            
            if azi(k) <= 180
            
                ildhorizontal(k,:) = M1(k,:) - M2(k,:);
                
            else
                
                ildhorizontal(k,:) = M2(k,:) - M1(k,:);
                
            end

        end
        
        freq_set = [1000 2000 5000 10000];
        plotColor = {'b','k','r','g'};
        polar_plot = zeros(size(azi,1)+1,length(freq_set));
        
        for k = length(freq_set):-1:1
            
            [mini,ind] = min(abs(freq - freq_set(k)));
            
            polar_plot(:,k) = [abs(ildhorizontal(:,ind));abs(ildhorizontal(1,ind))];
            
            polar([azi*pi/180;2*pi], [abs(ildhorizontal(:,ind));abs(ildhorizontal(1,ind))], plotColor{k});hold on;
            
            legendItem{length(freq_set)-k+1} = sprintf('%.1f Hz', freq(ind));
            
        end
        
        legend(legendItem);
        view(-90,90);
        title(Obj.GLOBAL_ListenerShortName);
        
        % ETC in the median plane
    case 'etcmedian'
        noisefloor = -50;
        azi = 0;
        thr = 2;
        Obj = SOFAexpand(Obj,'Data.Delay');
        hM = double(squeeze(Obj.Data.IR(:,ch,:)));
        pos = Obj.SourcePosition;
        idx = find(abs(pos(:,1))>90);
        pos(idx,2) = 180-pos(idx,2);
        pos(idx,1) = 180-pos(idx,1);
        idx = find(pos(:,1)<(azi+thr) & pos(:,1)>(azi-thr));
        M = (20*log10(abs(hM(idx,:))));
        pos = pos(idx,:);
        del = round(Obj.Data.Delay(idx,ch));
        M2 = zeros(size(M)+[0 max(del)]);
        for ii = 1:size(M,1)
            M2(ii,del(ii)+(1:Obj.API.N)) = M(ii,:);
        end
        M = M2-max(max(M2));
        M(M<noisefloor) = noisefloor;
        [ele,i] = sort(pos(:,2));
        M = M(i,:);
        surface(0:1/fs*1000:(size(M,2)-1)/fs*1000,ele,M(:,:));
        set(gca,'FontName','Arial','FontSize',10);
        set(gca, 'TickLength', [0.02 0.05]);
        set(gca,'LineWidth',1);
        cmap = colormap(hot);
        cmap = flipud(cmap);
        shading flat
        colormap(cmap);
        box on;
        colorbar;
        xlabel('Time (ms)');
        ylabel('Elevation (deg)');
        title([Obj.GLOBAL_Title '; channel: ' num2str(ch)],'Interpreter','none');
end


function f = myifftreal(c,N) % thanks goto the LTFAT <http://ltfat.sf.net>
if rem(N,2) == 0
    f = [c; flipud(conj(c(2:end-1,:)))];
else
    f = [c; flipud(conj(c(2:end,:)))];
end;
f = real(ifft(f,N,1));
