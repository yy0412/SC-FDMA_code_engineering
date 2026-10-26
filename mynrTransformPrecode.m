function out = mynrTransformPrecode(in,msc)
%mynrTransformPrecode nrTransformPrecode with explicit subcarrier count.

if isempty(in)
    out = in;
    return;
end
msc = double(msc);
nSymbols = size(in,1)/msc;
if nSymbols ~= fix(nSymbols)
    error('mynrTransformPrecode:InvalidLength', ...
        'Input rows (%d) must be an integer multiple of Msc (%d).',size(in,1),msc);
end
xtilde = reshape(in,msc,numel(in)/msc);
out = reshape(fft(xtilde)/sqrt(msc),size(in));
end
