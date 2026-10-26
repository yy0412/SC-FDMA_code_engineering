function [waveform,nominalStartSample,appliedTimingOffsetSamples] = ...
        mynrApplyTimeFrequencyOffset(txWaveform,timingOffsetSamples, ...
        frequencyOffsetHz,sampleRate,guardSamples)
%mynrApplyTimeFrequencyOffset Apply signed timing and carrier offsets.
%   The returned capture contains guard samples before and after the
%   nominal TB position. This keeps all subsequent indexing valid even
%   when the configured timing offset is outside the receiver search range.
%
%   NOMINALSTARTSAMPLE is the one-based first sample of the zero-offset TB.
%   A compensated TB beginning at timing estimate D therefore starts at
%   NOMINALSTARTSAMPLE + D.

if ~isscalar(timingOffsetSamples) || ~isfinite(timingOffsetSamples)
    error('timingOffsetSamples must be a finite scalar.');
end
if ~isscalar(frequencyOffsetHz) || ~isfinite(frequencyOffsetHz)
    error('frequencyOffsetHz must be a finite scalar.');
end
if ~isscalar(sampleRate) || ~isfinite(sampleRate) || sampleRate <= 0
    error('sampleRate must be a positive finite scalar.');
end
if nargin < 5 || isempty(guardSamples) || ~isfinite(guardSamples)
    guardSamples = 0;
end

appliedTimingOffsetSamples = round(double(timingOffsetSamples));
guardSamples = max([0 round(abs(double(guardSamples))) abs(appliedTimingOffsetSamples)]);
nominalStartSample = guardSamples+1;
txStartSample = nominalStartSample+appliedTimingOffsetSamples;

numberOfOutputSamples = size(txWaveform,1)+2*guardSamples;
waveform = complex(zeros(numberOfOutputSamples,size(txWaveform,2), 'like',txWaveform));
txRange = txStartSample+(0:size(txWaveform,1)-1);
waveform(txRange,:) = txWaveform; % Add time delay

sampleIndex = (0:numberOfOutputSamples-1).';
waveform = waveform .* exp(1i*2*pi*double(frequencyOffsetHz) * sampleIndex/double(sampleRate)); % Add frequency offset
end
