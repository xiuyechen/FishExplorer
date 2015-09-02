% reformatting frame_turn
figure;plot(frame_turn(:,17))
ylim([-10,200])
hold on;plot(frame_turn(:,18),'r')
%%
st2 = frame_turn(:,18);
IX_start = find(diff([1;st2])>10)
IX_stop = find(diff([1;st2])<-10)
% IX_stop = [IX_stop,8248];
%%
for i = 1:2,
    IX = IX_start(i):IX_stop(i)-1;
    frame_turn(IX,17) = frame_turn(IX,17)+40;
end