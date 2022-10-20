function idof = global_to_active(iglob, active)
idof = find(~(active-iglob));
end