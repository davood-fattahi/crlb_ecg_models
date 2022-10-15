clear all
close all
clc

%% load the case names in the directory
DA='..\mcode\database\ptb-gaussian-parameters';
files=dir([DA '\*.mat']); 
NumCases=length(files);

for i=1:NumCases
    GausParam=load([DA '\' files(i).name]);
    for j=1:12
        pmean=[1e-3*GausParam.params{j}.a; GausParam.params{j}.b; GausParam.params{j}.theta];
        tm=1/1000:1/1000:4;
        rng(i) % fixing the seed for random variation in beats
        [signal(j,:), ~]=ECGResembler(tm, pmean, 0, 70/60, 0);       
    end
plot(signal(:));

end
