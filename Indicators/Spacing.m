function sp = Spacing(PF_current)
Distance = pdist2(PF_current,PF_current,'cityblock');
Distance(logical(eye(size(Distance,1)))) = inf;
sp = std(min(Distance,[],2));
end

