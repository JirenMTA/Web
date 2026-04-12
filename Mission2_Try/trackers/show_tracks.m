function show_tracks(tracks)
    % Hien thi danh sach track con song

    if ~isfield(tracks, 'list') || isempty(tracks.list)
        fprintf('[show_tracks] Khong co track nao.\n');
        return;
    end

    nTracks = numel(tracks.list);
    nShown  = 0;

    for i = 1:nTracks
        tr = tracks.list(i);

        % Bo qua track da bi xoa
        if isfield(tr, 'status') && (strcmpi(tr.status, 'deleted') || strcmpi(tr.status, "deleted"))
            continue;
        end

        % Lay id
        if isfield(tr, 'id')
            id = tr.id;
        else
            id = i;
        end

        % Lay state = [range; velocity]
        if isfield(tr, 'state') && numel(tr.state) >= 2
            rangeVal = tr.state(1);
            velVal   = tr.state(2);
        else
            rangeVal = NaN;
            velVal   = NaN;
        end

        % In status neu co
        if isfield(tr, 'status')
            statusStr = string(tr.status);
        else
            statusStr = "unknown";
        end

        fprintf('[Track %d] Status: %s | Range: %.3f m | Velocity: %.3f m/s\n', ...
            id, statusStr, rangeVal, velVal);
        nShown = nShown + 1;
    end

    if nShown == 0
        fprintf('[show_tracks] Khong co track hop le de hien thi.\n');
    end
end