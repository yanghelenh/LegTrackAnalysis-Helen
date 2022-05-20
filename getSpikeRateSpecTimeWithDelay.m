% getSpikeRateSpecTimeWithDelay.m
%
% Function that takes in times of spikes and converts it to a spike rate
%  for a given time vector. Introduces a specified delay between the time
%  vector and the spike rate. 
% Computes spike rate by finding inter-spike-interval and getting inverse
%
% INPUTS:
%   spikeTimes - times (in sec) at which spikes occured
%   t - timing vector to get spike rate at (in sec)
%   delay - offset (in sec) b/w actual spike times and times in t; negative
%       delays are spike rate before and positive delays are spike rates
%       after
%
% OUTPUTS:
%   spikeRate - spike rate at each time in t, with offset
%
% CREATED: 10/27/21 - HHY
%
% UPDATED:
%   10/27/21 - HHY
%
function spikeRate = getSpikeRateSpecTimeWithDelay(spikeTimes, t, delay)

    % inter-spike-interval
    intSpikeInt = gradient(spikeTimes);
    
    % inverse of the inter-spike-interval
    invIntSpikeInt = ones(size(intSpikeInt))./intSpikeInt;
    
    % interpolate vector to same size as time vector, instead of same size
    %  as spikeInd vector
    % introduce delay
    spikeRate = interp1(spikeTimes, invIntSpikeInt, t + delay);
end