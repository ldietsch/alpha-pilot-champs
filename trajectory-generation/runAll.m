function [simRuns] = runAll(~)
    %RUNALL Summary of this function goes here
    %   Detailed explanation goes here
    tic
    input_file = '../input_files/Inputs.xlsx';
    output_file = '../results/output.xlsx';
    
    dateFormat = 'yyyy-mm-dd.HH-MM-SS';
    output_file = replace(output_file,'.xlsx', ...
        strcat('.',datestr(datetime,dateFormat),'.xlsx'));

    % parse run data
    [~,~,runs] = xlsread(input_file,'runs');

    startRun = runs{1,2};
    stopRun = runs{2,2};
    [~,~,params] = xlsread(input_file,'params');
    % Error Checking
    if isnumeric(startRun) && isnumeric(stopRun)
        if ~isnumeric(params{startRun+1,1}) || ~isnumeric(params{stopRun+1,1})
            disp('Start and Stop Run IDs do not match sheet "params"');
            return;
        end
    end
    
    numRuns = stopRun-startRun+1;
    Nsteps = params{numRuns+1,2};
    results = cell(numRuns+1,4*Nsteps + 1);
    try
        % test to see if we can write to the file
        xlswrite(output_file,[1]);
    catch ME
        disp('File is locked, please close.');
        disp('Press F5 to continue.');
        keyboard
    end
    
    % First row of results is header row
    results(1,1:2) = {'Run_ID','uMPC'};
    
    % Check for Parallel Processing
    runParallel = license('test','Distrib_Computing_Toolbox');

    % Comment out this line to run in parallel mode

%      runParallel = false;

    if runParallel
        for i = 1:numRuns
            F(i) = parfeval(@funcTrajectoryALM_SUMT_FixedTime,1, ...
                input_file,startRun+i-1);
        end
        
        % Build a waitbar to track progress
        h = waitbar(0,['Waiting for ', num2str(numRuns), ' runs to complete...']);        
        
        for i = 1:numRuns
            [completedIdx, r] = fetchNext(F);
            results(completedIdx+1,:) = r;
            waitbar(i/numRuns,h,sprintf('%d / %d runs completed.',i,numRuns));
        end
        delete(h);
    else
        for i=1:numRuns
            results(i+1,:) = funcTrajectoryALM_SUMT_FixedTime(input_file,...
                startRun+i-1);
        end
    end
    
    % Write outputs to file
    xlswrite(output_file,results);
   toc 
end