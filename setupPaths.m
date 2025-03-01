function setupPaths()
    addpath ./processRFmatlab-master/
    processRFmatlab_startup;
    addpath(genpath('./data_io/'));
    addpath(genpath('./preprocessing/'));
    addpath(genpath('./deconvolution/'));
    addpath(genpath('./array_processing/'));
    addpath(genpath('./imaging/'));
    addpath(genpath('./visualization/'));
    addpath(genpath('./utilities/'));
    javaaddpath('./utilities/FMI/lib/FMI.jar');
end