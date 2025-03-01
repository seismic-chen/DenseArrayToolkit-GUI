function write_MigResult(outFile, MigResult)
% WRITE_MIGRESULT  Saves the migration result structure to a .mat file.
%
% Usage:
%   write_MigResult('migResult.mat', MigResult)
%
% Inputs:
%   outFile   : Name of the output file (e.g., 'migResult.mat'). If the path
%               includes subfolders, those folders will be created if needed.
%   MigResult : The result structure to save, typically containing fields like:
%                 .x, .z, .mig, .migls, .param, ...
%
% Example:
%   MigResult.x = x;
%   MigResult.z = z;
%   MigResult.mig = mig;
%   MigResult.migls = migls;
%   write_MigResult('results/migResult.mat', MigResult);

    if nargin < 1 || isempty(outFile)
        outFile = 'migResult.mat';  % default filename if not provided
    end

    % Ensure the output directory exists (create if necessary)
    outDir = fileparts(outFile);
    if ~isempty(outDir) && ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    % Save the MigResult structure
    save(outFile, 'MigResult', '-v7.3');
    fprintf('[%s] Migration results saved to "%s"\n', datestr(now, 31), outFile);
end
