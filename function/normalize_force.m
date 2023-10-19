function force_norm=normalize_force(force,mvc)

[row_num,column_num]=size(force);

for i=1:row_num
    for j=1:column_num
        force_tmp=force{i,j};
        for u=1:size(force_tmp,1)
            for v=1:size(force_tmp,2)
                if(force_tmp(u,v)<0)
                    force_tmp(u,v)=force_tmp(u,v)/mvc(v,1);
                else
                    force_tmp(u,v)=force_tmp(u,v)/mvc(v,2);
                end
            end
        end
        force_norm{i,j}=force_tmp;
    end
end