clear all
close all

medpc_start_dir=uigetdir('');

output_dir=('');
subject_list=dir(medpc_start_dir);
subject_list = subject_list(~ismember({subject_list.name}, {'.', '..','.DS_Store'}))

index = 1;

for subj = 1:length(subject_list)
    subject = subject_list(subj).name;
    subj_dir = sprintf('%s/%s',medpc_start_dir,subject);
    subj_odors = dir(subj_dir);
    subj_odors = subj_odors(~ismember({subj_odors.name}, {'.', '..','.DS_Store'}))
    %subj_odors(1:3) = [];
    output_dir_sub = sprintf('%s%s', output_dir, subject);

    for odor = 1:length(subj_odors)
        odor_dir = sprintf('%s/%s', subj_dir, subj_odors(odor).name)
        subj_files = dir(odor_dir);
        subj_files = subj_files(~ismember({subj_files.name}, {'.', '..','.DS_Store'}))
        %subj_files(1:2) = [];
    
        ConcatenatedSessions = [];
        for medpc_file=1:length(subj_files)
            medpc_file_name = sprintf('%s/%s',subj_files(medpc_file).folder,subj_files(medpc_file).name);
            medpcfn = subj_files(medpc_file).name
            medpc_date = subj_files(medpc_file).date(1:11);

            tmp = MEDPCparser(medpc_file_name);
            rwidx = find(nansum(tmp,2)==0);
            tmp(rwidx,:) = [];
        
            ParsedMEDPCoutput{medpc_file} = tmp;
            ConcatenatedSessions = cat(1,ConcatenatedSessions,tmp);
            clear tmp rwidx
        end


        %%%%%% RUN FUNCTION(S)
        
        [SDTout, Dout, f2] = SDTcalculator2(ConcatenatedSessions);
        
        data(index).file = medpcfn(1:end-4);
        data(index).date = medpc_date;
        data(index).subject = subject;
        data(index).odor = subj_odors(odor).name;
        
        
        data(index).trials = SDTout(1);
        data(index).dprime = SDTout(2);
        data(index).criterion = SDTout(3);
        
        data(index).dprimeRollingWindow = Dout';
        

        output_dir_odor = sprintf('%s%s_%s', output_dir, subject, subj_odors(odor).name);
        fig_name1 = sprintf('%s_%s_FIG1.jpg',output_dir_odor,medpcfn(1:end-4))
        print (f2,'-djpeg',fig_name1)
        
        index = index+1;
        
        close all
        clear SDTout Dout f2 fig_name1
        
    end
end
 
   
    
