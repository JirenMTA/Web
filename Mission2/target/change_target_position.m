function target = change_target_position(target, time)
   target.Range = target.Range + time*target.Velocity;
end
