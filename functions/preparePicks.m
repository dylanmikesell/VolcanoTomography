function [ an2, ttp2, ntrac, pairn5 ] = preparePicks( fcmp, dtims, dtimsin, latlon, pairn3, nsrc)


npair = size( dtims, 1); % number of traces (i.e. correlation pairs)

npicks = sum( dtimsin( fcmp,:) ); % number of picks at this frequency (1=good, 0=bad from checkPicks.m)

% load up all times and positions
an  = zeros( npair, 2 ); % alllocate latlon of good picks to write
ttp = zeros( npair, 1 ); % allocate travel times of good picks to write

for ii = 1 : npair
    an( ii, : ) = latlon( pairn3( ii, 2 ), : );
    ttp( ii )   = dtims( ii, fcmp );
end

% allocate quality controlled picks 
an2   = zeros( npicks, 2 );
ttp2  = zeros( npicks, 1 );

% first make pairn4, a continuous record of source number
jmpd = 0;

pairn4 = zeros( npair - 1 );

pairn4(1) = pairn3( 1, 1 );

for ii = 1 : ( npair - 1 )
    
    if ( pairn3( ii+1, 1 ) > pairn3( ii, 1 ) + 1 )
        jmpd = jmpd + 1;
    end
    
    pairn4( ii+1 ) =  pairn3( ii+1, 1 ) - jmpd;

end

% make pairn5, a record of original source index
jmpd = 0;

pairn5(1) = pairn3(1,1);

for ii = 1 : ( npair - 1 )

    if ( pairn3( ii+1, 1 ) > pairn3( ii, 1 ) )
        jmpd             = jmpd + 1;
        pairn5( jmpd+1 ) =  pairn3( ii+1, 1 ) ;
    end
    
end

% pick qualifiers
cnt = 0;

ntrac = zeros( 1, nsrc ); % allocate -- keep track of good source-receiver pair data

for ii = 1 : npair % loop through each trace
    
    if ( dtimsin( fcmp, ii ) == 1 ) % if travel time good then save to matrix
        cnt                 = cnt + 1; % update local counter of good picks
        an2( cnt, : )       = an( ii, : ); % save latlon
        ttp2( cnt )         = ttp( ii ); % save travel time
        ntrac( pairn4(ii) ) = ntrac( pairn4(ii) ) + 1; % add a good receiver number to this source
    end
    
end

return