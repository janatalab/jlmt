% Copied from https://www.mathworks.com/help/audio/ug/convert-midi-files-into-midi-messages.html#ConvertMIDIFilesIntoMIDIMessagesExample-7

function [valueOut,byteLength] = findVariableLength(lengthIndex,readOut)

    byteStream = zeros(4,1);
    
    for i = 1:4
        valCheck = readOut(lengthIndex+i);
        byteStream(i) = bitand(valCheck,127);   % Mask MSB for value
        if ~bitand(valCheck,uint32(128))        % If MSB is 0, no need to append further
            break
        end
    end
    
    valueOut = polyval(byteStream(1:i),128);    % Base is 128 because 7 bits are used for value
    byteLength = i;

end