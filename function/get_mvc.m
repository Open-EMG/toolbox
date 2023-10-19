function mvc=get_mvc(path,subject,session)

%path: the location of dataset

%subject_idx: the index of a specific subject (1-20)

%session_idx: the index of a specific session (1 or 2)

%mvc: the mvc value (5 x 2 matrix in double). The 1st and 2nd column 
%represents the flexion and extension of a specific finger, respectively. 
%The 1st to 5th row represents thumb, index, middle, ring and little finger,
%respectively. Data in mvc are absolute force value (constantly positive).

force_data=load_mvc(path,subject,session,'force');

mvc=zeros(5,2);

for i=1:10
    finger=ceil(i/2);
    direction=i-(finger-1)*2;
    force_tmp=abs(force_data{i,1}(:,finger));
    force_sort=sort(force_tmp,'descend');
    mvc(finger,direction)=mean(force_sort(1:200));
end