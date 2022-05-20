% getPhaseDiffPDFs.m
%
% Function that takes phase differences (from getLegPhaseDiffs()) and
%  converts them into probability density functions (PDFs)
%
% INPUTS:
%   phaseDiffs - struct of phase differences from getLegPhaseDiffs()
%   numBins - number of bins to divide PDF into
%   units - 'radians' or 'degrees' for units of leg phase
%
% OUTPUTS:
%   phasePDFs - struct with same fields as phaseDiffs, each field being PDF
%   phaseBinMids - midpoints of bins
%   phaseBinEdges - edges of bins
%
% CREATED: 10/26/21 - HHY
%
% UPDATED:
%   10/26/21 - HHY
%
function [phasePDFs, phaseBinMids, phaseBinEdges] = getPhaseDiffPDFs(...
    phaseDiffs, numBins, units)

    % appropriate scaling depending on units
    if (strcmpi(units, 'radians'))
        binMaxVal = 2 * pi;
    elseif (strcmpi(units, 'degrees'))
        binMaxVal = 360;
    else
        disp('units must be radians or degrees');
        phasePDFs = [];
        return;
    end
    
    % start value for bins
    binMinVal = 0;
    
    % bin edges
    phaseBinEdges = linspace(binMinVal,binMaxVal, numBins + 1);

    % get midpoints of bins
    phaseBinMids = zeros(1,numBins);
    for i = 1:numBins
        phaseBinMids(i) = (phaseBinEdges(i+1) - phaseBinEdges(i))/2 + ...
            phaseBinEdges(i);
    end

    % names of all fields of phaseDiffs (which leg relations to use)
    phaseFieldnames = fieldnames(phaseDiffs);

    % get pdfs, all time points
    for i = 1:length(phaseFieldnames)
        [phasePDFs.(phaseFieldnames{i}), ~] = histcounts(...
            phaseDiffs.(phaseFieldnames{i}), phaseBinEdges, ...
            'Normalization', 'pdf');
    end
end