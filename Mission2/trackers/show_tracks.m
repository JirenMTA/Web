function show_tracks(tracks)
    for i=1:numel(tracks.list) 
        fprintf("[Tracked id %d]: Range: %f - Velocity: %f\n", ...
            tracks.list(i).id, tracks.list(i).state(1), tracks.list(i).state(2));
    end
end
