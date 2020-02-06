function IR = SOFAgetIR(Obj,mic_pos_idx,lpk_pos_idx)

if ~strcmp(Obj.GLOBAL_SOFAConventions,'AmbisonicsDRIR') || ~strcmp(Obj.GLOBAL_SOFAConventionsVersion,'2.0')
    
    error('Only AmbisonicsDRIR v2.0 convention is supported!');
    
end

nof_mic_pos = length(nonzeros(any(Obj.ListenerSourcePairIndex,2)));

if nargin < 3
    
    if nargin < 2
        
        mic_pos_idx = input(sprintf('Enter the desired microphone position (1->%i) [1]:\n',nof_mic_pos));
        
        if isempty(mic_pos_idx)            
            mic_pos_idx = 1;            
        elseif mic_pos_idx < 1 || mic_pos_idx > nof_mic_pos            
            error('The desired microphone position must be an integer between 1 and %i!',nof_mic_pos);
        end
    
        lpk_pos_idx = input(sprintf('Enter the desired loudspeaker position (1->%i) [all positions associated to microphone position %i]: \n',length(find(~isnan(Obj.ListenerSourcePairIndex(mic_pos_idx,:)))),mic_pos_idx));
    
        if isempty(lpk_pos_idx)
             lpk_pos_idx = Obj.ListenerSourcePairIndex(mic_pos_idx,find(~isnan(Obj.ListenerSourcePairIndex(mic_pos_idx,:)))); 
        elseif lpk_pos_idx < 1 || lpk_pos_idx > length(find(~isnan(Obj.ListenerSourcePairIndex(mic_pos_idx,:))))
            error('The desired loudspeaker position must be an integer between 1 and %i!',length(find(~isnan(Obj.ListenerSourcePairIndex(mic_pos_idx,:)))));
        end
        
    end
    
end

fprintf('Desired microphone position: %i\n\tDesired loudspeaker position: %s\n',mic_pos_idx,sprintf('%i ',lpk_pos_idx));

[~,ind] = find(~isnan(Obj.ListenerSourcePairIndex(mic_pos_idx,:)));

IR = Obj.Data.IR(Obj.ListenerSourcePairIndex(mic_pos_idx,ind(lpk_pos_idx)),:,:,:);