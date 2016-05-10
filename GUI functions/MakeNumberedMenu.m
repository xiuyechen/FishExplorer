function menu = MakeNumberedMenu(name) % e.g. name = {Cluster.name} (note {})
menu = [{'(choose)'},name];
for j=2:length(menu),menu(j)={[num2str(j-1) ': ' menu{j}]};end
end
