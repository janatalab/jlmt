function msgArray = readmidi(fname, BPM)
    msgArray = [];

    if nargin < 2
        BPM = 120; % Assume a tempo of 120 beats-per-minute by default
    end

    % Get a file handle
    fh = fopen(fname);

    if fh < 0
        error("Could not open file %s", fname)
    end

    % Read the file
    [readOut, byteCount] = fread(readme);

    % Close the file
    fclose(readme);

    % Extract the tick count
    ticksPerQNote = polyval(readOut(13:14),256);

    % Initialize timestamp value
    ts = 0;

    % Parse track chunks in outer loop
    chunkIndex = 14;
    while chunkIndex < byteCount
        
        % Read header of track chunk, find chunk length   
        % Add 8 to chunk length to account for track chunk header length
        chunkLength = polyval(readOut(chunkIndex+(5:8)),256)+8;
        
        ptr = 8+chunkIndex;             % Determine start for MIDI event parsing
        statusByte = -1;                % Initialize statusByte. Used for running status support
        
        % Parse MIDI track events in inner loop
        while ptr < chunkIndex+chunkLength
            % Read delta-time
            [deltaTime,deltaLen] = findVariableLength(ptr,readOut);  
            % Push pointer to beginning of MIDI message
            ptr = ptr+deltaLen;
            
            % Read MIDI message
            [statusByte,messageLen,message] = interpretMessage(statusByte,ptr,readOut);
            % Extract relevant data - Create midimsg object
            [ts,msg] = createMessage(message,ts,deltaTime,ticksPerQNote,BPM);
            
            % Add midimsg to msgArray
            msgArray = [msgArray;msg];
            % Push pointer to next MIDI message
            ptr = ptr+messageLen;
        end
        
        % Push chunkIndex to next track chunk
        chunkIndex = chunkIndex+chunkLength;
    end

end