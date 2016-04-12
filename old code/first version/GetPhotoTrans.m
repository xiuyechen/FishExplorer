function phototrans = GetPhotoTrans(photostate)
% 4 photostates: 0 = all black; 1 = black/white; 2 = white/black; 3 = all white;
phototrans = zeros(size(photostate));
sw_ = find([1 diff(photostate)]);
interval = mode(diff(sw_));
sw = sw_(2):interval:length(photostate);
if sw_(2)-sw_(1)>interval,
    sw = [1, sw_(2)-interval, sw];
else
    sw = [1, sw];
end
sw_last = [length(photostate)-interval+1,sw]; % fake-circular
sw_end = [sw(2:end),length(photostate)+1]; % fake-circular
for k = 1:length(sw),
    i = photostate(sw_last(k)); % last digit
    j = photostate(sw(k)); % current digit
    phototrans(sw(k):sw_end(k)-1) = i*4 + j;
end

end