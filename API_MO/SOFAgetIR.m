%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   function IR = SOFAgetIR(Obj,mic_pos_idx,lpk_pos_idx)
%
%   This function extracts the IR data based on the indices of both the
%   microphone position and loudspeaker position. 
%
%   mic_pos_idx and lpk_pos_idx are optional arguments. If not specified,
%   the user will be prompted for choosing them. If so, mic_pos_idx can
%   still be left empty: in that case, IR is an array containing all 
%   valid IR data, sorted by increasing index of microphone positions (this
%   is equivalent to removing the all-zeros data from Obj.Data.IR along the
%   1rst dimension). If mic_pos_idx is specified, lpk_pos_idx can also be
%   left empty: in that case, IR is an array containing all valid IR data
%   for the specified microphone position.
%
%   The index of the loudspeaker position is global (global indexing 
%   across the whole measurements dataset).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function IR = SOFAgetIR(Obj,mic_pos_idx,lpk_pos_idx)

if ~strcmp(Obj.GLOBAL_SOFAConventions,'AmbisonicsDRIR') || ~strcmp(Obj.GLOBAL_SOFAConventionsVersion,'2.0')
    
    error('Only AmbisonicsDRIR v2.0 convention is supported!');
    
end

nof_mic_pos = length(nonzeros(any(Obj.ListenerSourcePairIndex,2)));

if nargin < 3
    
    if nargin < 2
        
        mic_pos_idx = input(sprintf('Enter the desired microphone position (1->%i) [1]:\n',nof_mic_pos));
        
        if isempty(mic_pos_idx)
            ListenerSourcePairIndex_tranposed = transpose(Obj.ListenerSourcePairIndex);
            mic_pos_idx = ListenerSourcePairIndex_tranposed(find(transpose(~isnan(Obj.ListenerSourcePairIndex))));
            IR = Obj.Data.IR(mic_pos_idx,:,:,:);
            fprintf('Removing the all-zeros data from Obj.Data.IR along the 1rst dimension.\n');
            return;
        elseif mic_pos_idx < 1 || mic_pos_idx > nof_mic_pos
            error('The desired microphone position must be an integer between 1 and %i!',nof_mic_pos);
        end
        
    end    
    
    lpk_pos_idx = input(sprintf('Enter the desired loudspeaker position from the list: %s [all positions]: \n',num2str(find(~isnan(Obj.ListenerSourcePairIndex(mic_pos_idx,:))))));
    
    if isempty(lpk_pos_idx)
        lpk_pos_idx = find(~isnan(Obj.ListenerSourcePairIndex(mic_pos_idx,:)));
    elseif ~find(~isnan(Obj.ListenerSourcePairIndex(mic_pos_idx,:))==lpk_pos_idx)
        error('The desired loudspeaker position must be an integer from the list %s!',num2str(find(~isnan(Obj.ListenerSourcePairIndex(mic_pos_idx,:)))));
    end 
    
end

fprintf('Desired microphone position: %i\n\tDesired loudspeaker position: %s\n',mic_pos_idx,sprintf('%i ',lpk_pos_idx));
IR = Obj.Data.IR(Obj.ListenerSourcePairIndex(mic_pos_idx,lpk_pos_idx),:,:,:);