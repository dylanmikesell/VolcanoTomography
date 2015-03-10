function swddat = ftanFilter(symc, fmin, fmax, df, dt, wpm, alpha )

fArray = fmin : df : fmax;
nfreq = numel(fArray);
freq = 1/dt;

% frequency filter to see dispersion
swddat = zeros(nfreq,wpm);

% assuming wpm is even
omgas = 2 * pi * (freq / wpm) * [ -(wpm/2) : ( (wpm/2) - 1) ];

% frequency filter for FTAN
for ii = 1 : nfreq
    
    omga0 = 2 *pi * fArray(ii);
    
    gfl = exp( -alpha * ( ( ( omgas - omga0 ) ./ omga0 ).^2 ) );
    %gfl = exp(-alpha*(((abs(omgas)-omga0)./omga0).^2));
    
    swddat(ii,:) =...
        ifft( ifftshift( fftshift( fft( hilbert( symc ) ) ) .* gfl ) );
end

return