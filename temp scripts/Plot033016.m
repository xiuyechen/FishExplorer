
%%
range_stim = 1:12;
n = round(length(range_stim)*1.1);
clrmap = hsv(max(1,n));

anat_yx_norm = getappdata(hfig,'anat_yx_norm');

figure;
hold on;
image(anat_yx_norm)
axis equal
axis ij
axis off

range_fish = 8:11;
for i_fish = range_fish,
    IX = ~isnan(Fish{i_fish}.xyz_norm_avr(:,1));
    scatter(Fish{i_fish}.xyz_norm_avr(IX,2),Fish{i_fish}.xyz_norm_avr(IX,1),20,clrmap(IX,:))
end

%%
range_stim = 1:12;
n = round(length(range_stim)*1.1);
clrmap = hsv(max(1,n));

anat_yx_norm = getappdata(hfig,'anat_yx_norm');

figure;
hold on;
image(anat_yx_norm)
axis equal
axis ij
axis off

% range_fish = 9;
% for i_fish = range_fish,
%     gIX = Fish{i_fish}.gIX;
%     IX = find(gIX == 9);
%     scatter(Fish{i_fish}.xyz_norm(IX,2),Fish{i_fish}.xyz_norm(IX,1),1);%,20,clrmap(IX,:))
% end
i_fish = 9;
IX = find(VAR(i_fish).ClusGroup{3}.gIX == 37);
cIX_abs = VAR(i_fish).ClusGroup{3}.cIX_abs(IX);
f.LoadFullFish(hfig,i_fish);
CellXYZ_norm = getappdata(hfig,'CellXYZ_norm');
xyz_norm = CellXYZ_norm(cIX_abs,:);    
scatter(xyz_norm(:,2),xyz_norm(:,1),1);%,20,clrmap(IX,:))

i_fish = 10;
IX = find(VAR(i_fish).ClusGroup{3}.gIX == 87);
cIX_abs = VAR(i_fish).ClusGroup{3}.cIX_abs(IX);
f.LoadFullFish(hfig,i_fish);
CellXYZ_norm = getappdata(hfig,'CellXYZ_norm');
xyz_norm = CellXYZ_norm(cIX_abs,:);    
scatter(xyz_norm(:,2),xyz_norm(:,1),1);%,20,clrmap(IX,:))
%%
range_stim = 1:12;
n = round(length(range_stim)*1.1);
clrmap = hsv(max(1,n));

anat_yx_norm = getappdata(hfig,'anat_yx_norm');

figure;
hold on;
image(anat_yx_norm)
axis equal
axis ij
axis off


i_fish = 8;
    gIX = FC{i_fish}.gIX;
    IX = find(gIX <100);
    scatter(FC{i_fish}.xyz_norm(IX,2),FC{i_fish}.xyz_norm(IX,1),1,'r');%,20,clrmap(IX,:))

i_fish = 9;
    scatter(FC{i_fish}.xyz_norm(:,2),FC{i_fish}.xyz_norm(:,1),1,'g');%,20,clrmap(IX,:))
    
    i_fish = 10;
    gIX = FC{i_fish}.gIX;
    IX = find(gIX <100);
    scatter(FC{i_fish}.xyz_norm(IX,2),FC{i_fish}.xyz_norm(IX,1),1,'b');%,20,clrmap(IX,:))
    
%% 
anat_yx_norm = getappdata(hfig,'anat_yx_norm');

figure;
hold on;
image(anat_yx_norm)
axis equal
axis ij
axis off

clrmap = hsv(round(max(bins)*1.1));

for i_node = 1:length(Nodes),
   i_fish = Nodes(i_node,1);
   clusIX = Nodes(i_node,2);
   
   temp = VAR(i_fish).ClusGroup{3}.gIX;
   IX = find(temp == clusIX);
   temp = VAR(i_fish).ClusGroup{3}.cIX_abs;
   cIX_abs = temp(IX);
   LoadFishDataWithoutTS(hfig,i_fish);
   CellXYZ_norm = getappdata(hfig,'CellXYZ_norm');
   xyz_norm = CellXYZ_norm(cIX_abs,:);
%    cmap = repmat(clrmap(bins(i_node),:)',1,length(cIX_abs));
   scatter(xyz_norm(:,2),xyz_norm(:,1),2,clrmap(bins(i_node),:));%,20,clrmap(IX,:))
end