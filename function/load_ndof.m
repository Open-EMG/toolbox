function data=load_ndof(path,subject,session,sig_type)

if(strcmp(sig_type,'force'))
    M=5;
else
    M=256;
end

for i=1:15
    for j=1:2
        Fname = [path,'/ndof_dataset/subject',subject,'_session',num2str(session),'/ndof_',sig_type,'_combination',num2str(i),'_sample',num2str(j),'.dat'];         % Name with path.
        fid = fopen(Fname,'r','n');            % Open for reading.
        if fid<0, error(['Failed to open: ' Fname]); end
        data_tmp = fread(fid, [M,inf], 'int16');  % Read.
        fclose(fid);
        data_tmp=data_tmp';
        
        Fname = [path,'/ndof_dataset/subject',subject,'_session',num2str(session),'/ndof_',sig_type,'_combination',num2str(i),'_sample',num2str(j),'.hea'];         % Name with path.
        head_info=textread(Fname,'%s');
        idx = find( strcmp( head_info , ['ndof_',sig_type,'_combination',num2str(i),'_sample',num2str(j),'.dat'] ));
    
        for u=1:length(idx)
            str_tmp=head_info(idx(u)+2);
            idx2=strfind( str_tmp , '(' );
            gain=str2num(str_tmp{1,1}(1:(idx2{1,1}-1)));
            idx3=strfind( str_tmp , ')' );
            baseline=str2num(str_tmp{1,1}((idx2{1,1}+1):(idx3{1,1}-1)));
        
            data_tmp(:,u)=(data_tmp(:,u)-baseline)/gain;
        end
        data{i,j}=data_tmp;
    end
end
