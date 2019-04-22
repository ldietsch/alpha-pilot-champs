clear all
close all
addpath('../discretization');
addpath('../motion-model');
addpath('../post-processing');
addpath('../standard-problem');
addpath('../test-functions');

input_file = "../input_files/Inputs.xlsx";
results = "../results/output.2019-04-22.00-00-12.xlsx";

X = processData(results, input_file);