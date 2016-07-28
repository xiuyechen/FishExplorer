function fpsec = GetFishFreq(i_fish)

M_fpsec = {1.97,1.97,1.97,1.97,1.97,1.97,1.97,...% in Hz % Fish 1-7
    1.97,1.93,1.97,1.84,... % Fish 8-11 % Actually Fish 10 freq unrecorded
    2.56,2.33,2.27,... % Fish 12-14
    2.33,2.33,2.38,2.38}; % Fish 15-18

fpsec = M_fpsec{i_fish};
end