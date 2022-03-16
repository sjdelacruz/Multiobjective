function gd = GD(PF_current,PF_true)
Distance = min(pdist2(PF_current,PF_true),[],2);
gd    = norm(Distance)/length(Distance);
end

