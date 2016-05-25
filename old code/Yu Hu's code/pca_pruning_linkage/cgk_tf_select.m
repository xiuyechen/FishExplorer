function cgk_new = cgk_tf_select(cgk, tf_select, hfig_number)
if ~exist('hfig_number', 'var')
    hfig_number = 1;
end

cgk_new = cgk;
cgk_new{1} = cgk{1}(tf_select);
cgk_new{2} = cgk{2}(tf_select);

if hfig_number > 0
    if ~isempty(cgk_new{1})
        push_cgk(figure(hfig_number), cgk_new);
    else
        'empty set!'
    end       
end

end