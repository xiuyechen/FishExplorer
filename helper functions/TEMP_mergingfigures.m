    % Prepare subplots
newfig = figure('Position',[800,200,1400,1000]);
h(1)=subplot_tight(1,2,1);%axis ij;axis off
h(2)=subplot_tight(1,2,2);%axis ij;axis off
f(1) = fig1;
f(2) = fig3;
% Paste figures on the subplots
for j = 1:2
    ax(j) = get(f(j),'CurrentAxes');
    
    % copy axes properties
    % axis ij
    str_prop = {'XDir','YDir','Visible'};
    str_flag{1} = get(ax(j),'XDir');
    str_flag{2} = get(ax(j),'YDir');
    % axis off
    str_flag{3} = get(ax(j),'Visible');
%     str_flag{4} = get(get(ax(j),'Title'),'Visible');
    
    
    copyobj(allchild(ax(j)),h(j));
    for i = 1:length(str_prop)
        set(h(j),str_prop{i},str_flag{i});
    end
end