% exam_rabbit_blood_pressure
pathold = path; path(pathold,'.\MatlabFuns');
close all; clear; 
format compact

% %% load data
% X = cell(1,3);
% spmd (3)
%     X{labindex} = getBPtxtData(['.\All-Rabbits-Data-WY\XR00',num2str(labindex),'\XR00',num2str(labindex),'.txt']);
% end
% 
% data = X{1}; save('.\All-Rabbits-Data-WY\XR001_BP.mat','data');
% data = X{2}; save('.\All-Rabbits-Data-WY\XR002_BP.mat','data'); 
% data = X{3}; save('.\All-Rabbits-Data-WY\XR003_BP.mat','data');
% 
% return

%% analysis
base_path = '.\All-Rabbits-Data-WY\';
sub_path = '';

caseid = 1;

switch caseid
    case 1
        casename = 'XR001'; sub_path = [casename,'\'];
        desc = {'XR001'};
    case 2
        casename = 'XR002'; sub_path = [casename,'\'];
        desc = {'XR002'};
    case 3
        casename = 'XR003'; sub_path = [casename,'\'];
        desc = {'XR003'};
    otherwise
        disp('Incorrect caseId');
end

if ~exist([base_path,casename,'_BP.mat'],'file')
    filename = [base_path, sub_path, casename, '.txt'];
    disp(filename);
    data = getBPtxtData(filename);
    save([base_path,casename,'_BP.mat'],'data');
else
    load([base_path,casename,'_BP.mat']);
end
