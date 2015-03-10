function [qq1, qq2, dtims] = pickFTANvelocities( fmin, fmax, df,...
    vmin, vmax, dv, vUnits, dt, offset, wpm, ftanMat )


% check input units on velocity and modify distance units to match
switch vUnits
    case 'm'
        offset = offset * 1000; % convert interstation distance to meters
    case 'km'
        offset = offset; % interstation distance is already in kilometers
end
% setup velocity vector
vArray = vmax : -dv : vmin;
nVel   = numel(vArray);

% new time axis for velocity vector we specify
tArray2 = offset ./ vArray; 
tArray = (1:wpm) .* dt; % original time axis

% setup the frequency vector
fArray = fmin : df : fmax;
nfreq  = numel(fArray);

% allocate the interpolated FTAN matrix
ftanMatInt = zeros( nfreq, nVel );

% interpolate to velocity using a spline at this frequency
for ff = 1 : nfreq
    ftanMatInt( ff, : ) = spline( tArray, abs( ftanMat( ff, : ) ), tArray2 );
end

% pick maximum at each frequency
[qq1, qq2] = max( abs( ftanMatInt ), [], 2 );

% Matt's version before generalizing the vmin,vmax,dv parameters
% dtims = offset ./ (vmax-(dv*qq2));

% this is different than Matt's version
dtims = tArray2(qq2);

return