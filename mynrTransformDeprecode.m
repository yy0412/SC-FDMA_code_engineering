function out = mynrTransformDeprecode(in,msc)
%mynrTransformDeprecode nrTransformDeprecode with explicit subcarrier count.

if isempty(in)
    out = in;
    return;
end
msc = double(msc);
nSymbols = size(in,1)/msc;
if nSymbols ~= fix(nSymbols)
    error('mynrTransformDeprecode:InvalidLength', ...
        'Input rows (%d) must be an integer multiple of Msc (%d).',size(in,1),msc);
end
ytilde = reshape(in,msc,numel(in)/msc);
out = reshape(ifft(ytilde)*sqrt(msc),size(in));
end
