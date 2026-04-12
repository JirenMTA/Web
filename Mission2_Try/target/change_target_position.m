function target = change_target_position(target, time)
   target.Range(1) = target.Range(1) + time*target.Velocity(1);
end
