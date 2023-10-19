function force_cal=calibrate_force(force,cal)

[row_num,column_num]=size(force);

for i=1:row_num
    for j=1:column_num
        force_tmp=force{i,j};
        cal_tmp=cal{i,j};
        force_cal{i,j}=force_tmp(:,2:6)-mean(cal_tmp(:,2:6));
    end
end