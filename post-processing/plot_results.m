function X = plot_results(input_file, results, startRun)
close all
addpath('../discretization');
addpath('../motion-model');
addpath('../post-processing');
addpath('../standard-problem');
addpath('../test-functions');


X = processData(results, input_file, startRun);

x = zeros(size(X,1),size(X,2)/12);
y = zeros(size(X,1),size(X,2)/12);
z = zeros(size(X,1),size(X,2)/12);

Gates = xlsread('../trajectory-generation/Gate Locations.xlsx');
gates = Gates(1:6,1);
gatex = Gates(1:6,2);
gatey = Gates(1:6,3);
gatez = Gates(1:6,4);
run = startRun;
for i = 1:size(X,1)
    
   [x(i,:), y(i,:), z(i,:)] = extractPos(X(i,:)', size(X,2)/12); 
   plot3(x(i,:),y(i,:),-z(i,:));
   hold on
   ID_labels{i} = "ID = " + run;
   run = run + 1;
end
ID_labels{end+1} = "Gates";
title("Calculated trajectories for N = 5")
xlabel("X [m]")
ylabel("Y [m]")
zlabel("Z [m]")
scatter3(gatex, gatey, gatez)
legend(ID_labels)
grid
setFont();
end