function [vtracks,ntracks,meanlength,rmslength] = PredictiveTracker(inputnames,threshold,max_disp,bground_name)
% Usage: [vtracks,ntracks,meanlength,rmslength] = PredictiveTracker(inputnames,threshold,max_disp,[bground_name])
% Given a movie of particle motions saved as a series of image files
% specified by "inputnames" (e.g., 'movie*.tif'), PredictiveTracker 
% produces Lagrangian particle tracks using a predictive three-frame 
% best-estimate algorithm. A steady background, read from the file
% "bground_name", is subtracted from each frame so that objects that do not
% move are not tracked; see BackgroundImage.m. Particles are identified 
% wherever a pixel has intensity larger than its neighbors and larger than 
% "threshold". Each particle is tracked using a kinematic prediction, and a
% track is broken when no particle lies within "max_disp" pixels of the
% predicted location. The results are returned in the structure "vtracks", 
% whose fields "len", "X", "Y", "T", "U", and "V" contain the length, 
% horizontal coordinates, vertical coordinates, times, horizontal 
% velocities, and vertical velocities of each track, respectively. The 
% total number of tracks is returned as "ntracks"; the mean and root-mean-
% square track lengths are returned in "meanlength" and "rmslength", 
% respectively. This file can be downloaded from 
% http://leviathan.eng.yale.edu/software.

% Written by Nicholas T. Ouellette September 2010. 

  bground_name_default = 'background.tif';

  if nargin<1
      error(['Usage: [vtracks,ntracks,meanlength,rmslength] = ' ...
          mfilename '(inputnames,threshold,max_disp,[bground_name])'])
  end
  if ~exist('bground_name','var') || isempty('bground_name')
      bground_name=bground_name_default;
  end

  % set up a few constants and parameters
  
  % minimum track length to store
  minlength = 10;
  
  % create the file names of the images
  names = dir(inputnames);
  filepath = inputnames(1:max(findstr(inputnames,'/')));
  if ~isempty(filepath) && filepath(end)~='/'
      filepath(end+1)='/';
  end

  % pre-compute logarithms used for locating particle centers
  d=imfinfo([filepath names(1).name]);
  color_depth=2^(d.BitDepth);
  logs = 1:color_depth;
  logs = [log(0.0001) log(logs)];

  % read in the background image
  if exist(bground_name,'file')==2
      background = imread(bground_name);
  else
      warning('MATLAB:PredictiveTracker:noBackgroundFile', ...
        ['Cannot find background file ' bground_name ...
        '. Using zeros instead.'])
      background = imread([filepath names(1).name])*0;
  end
  if ndims(background)==3 % check for RGB instead of grayscale
      background=mean(background,3);
  end
  
  % process the first frame to get things moving
  im = imread([filepath names(1).name]);
  if ndims(im)==3
      im=mean(im,3); % convert to grayscale if necessary
  end
  im = im - background;
  im(im<0) = 0;

  fr0 = FindParticles(im, threshold, logs);
  nparticles = numel(fr0(:,1));
  % each of these particles will start a new track. 
  % to hold them, we'll use a struct with fields 'len','X','Y','T'.
  % preallocate the array:
  tracks = repmat(struct('len',[],'X',[],'Y',[],'T',[]),nparticles,1);
  for ii = 1:nparticles
    tracks(ii) = struct('len',1,'X',fr0(ii,1),'Y',fr0(ii,2),'T',1);
  end
  % we'll also keep a record of which of these tracks are 'active';
  % all of these initial tracks are active.
  active = 1:nparticles;
  n_active = numel(active);
  disp(['Processed ' names(1).name])
  disp(['    Number of particles found: ' num2str(nparticles,'%.0f')])
  disp(['    Number of active tracks: ' num2str(n_active,'%.0f')])
  disp(['    Total number of tracks: ' num2str(numel(tracks),'%.0f')])
  
  % loop over frames
  for t = 2:numel(names)
    % read the next image
    im = imread([filepath names(t).name]);
    if ndims(im)==3
        im=mean(im,3); % convert to grayscale if necessary
    end
    im = im - background;
    im(im<0) = 0;
    % and find the particles
    fr1 = FindParticles(im, threshold, logs);
    nfr1 = numel(fr1(:,1));

    % now we want to figure out how to match the tracks. we do this by
    % estimating the particle velocities and using simple kinematics to
    % project these positions into fr1.

    % for convenience, we'll grab the relevant positions from the tracks
    now = zeros(n_active,2);
    prior = zeros(n_active,2);
    for ii = 1:n_active
      tr = tracks(active(ii));
      now(ii,1) = tr.X(end);
      now(ii,2) = tr.Y(end);
      if tr.len > 1
        prior(ii,1) = tr.X(end-1);
        prior(ii,2) = tr.Y(end-1);
      else
        prior(ii,:) = now(ii,:);
      end
    end
    
    % estimate a velocity for each particle in fr0
    velocity = now - prior;
    % and use kinematics to estimate a future position
    estimate = now + velocity;
    
    % define cost and link arrays
    costs = zeros(n_active,1);
    links = zeros(n_active,1);
    
    % loop over active tracks
    for ii = 1:n_active
      % now, compare this estimated positions with particles in fr1
      dist_fr1 = (estimate(ii,1)-fr1(:,1)).^2 + (estimate(ii,2)-fr1(:,2)).^2;
      % save its cost and best match
      costs(ii) = min(dist_fr1);
      if costs(ii) > max_disp^2
        continue;
      end
      bestmatch = find(dist_fr1 == costs(ii));
      % if there is more than one best match, we are confused; stop
      if numel(bestmatch) ~= 1
        continue;
      end
      % has another track already matched to this particle?
      ind = links == bestmatch;
      if sum(ind) ~= 0
        if costs(ind) > costs(ii)
          % this match is better
          links(ind) = 0;
        else
          continue;
        end
      end
      links(ii) = bestmatch;
    end
    
    % now attach the matched particles to their tracks
    matched = zeros(nfr1,1);
    inactive = find(links == 0);
    for ii = 1:n_active
      if links(ii) ~= 0 
        % this track found a match
        tracks(active(ii)).X(end+1) = fr1(links(ii),1);
        tracks(active(ii)).Y(end+1) = fr1(links(ii),2);
        tracks(active(ii)).len = tracks(active(ii)).len + 1;
        tracks(active(ii)).T(end+1) = t;
        matched(links(ii)) = 1;
      end
    end
    active = setdiff(active, inactive);
    % and start new tracks with the particles in fr1 that found no match
    unmatched = find(matched == 0);
    newtracks = repmat(struct('len',[],'X',[],'Y',[],'T',[]),numel(unmatched),1);
    for ii = 1:numel(unmatched)
      newtracks(ii) = struct('len',1,'X',fr1(unmatched(ii),1),...
                             'Y',fr1(unmatched(ii),2),'T',t);
    end
    active = [active numel(tracks):numel(tracks)+numel(newtracks)];
    tracks = [tracks ; newtracks];
    n_active = numel(active);
    
    disp(['Processed ' names(t).name])
    disp(['    Number of particles found: ' num2str(nfr1,'%.0f')])
    disp(['    Number of active tracks: ' num2str(n_active,'%.0f')])
    disp(['    Number of new tracks started here: ' ...
        num2str(numel(unmatched),'%.0f')])
    disp(['    Number of tracks that found no match: ' ...
        num2str(sum(links==0),'%.0f')])
    disp(['    Total number of tracks: ' num2str(numel(tracks),'%.0f')])

  end

  disp('Pruning...');
  % prune off tracks that are too short
  tracks = tracks([tracks.len] >= minlength);
  ntracks = numel(tracks);
  meanlength = mean([tracks.len]);
  rmslength = sqrt(mean([tracks.len].^2));
  
  % now compute the velocities
  disp('Differentiating...');
  
  % define the convolution kernel
  w = 1;
  L = 3*w;
  Av = 1.0/(0.5*w^2 * (sqrt(pi)*w*erf(L/w) - 2.0*L*exp(-L^2/w^2)));
  vkernel = -L:L;
  vkernel = Av.*vkernel.*exp(-vkernel.^2./w^2);

  % loop over tracks
  vtracks = repmat(struct('len',[],'X',[],'Y',[],'T',[],'U',[],'V',[]),ntracks,1);
  for ii = 1:ntracks
    u = -conv(tracks(ii).X,vkernel,'valid');
    v = -conv(tracks(ii).Y,vkernel,'valid');
    vtracks(ii) = struct('len',tracks(ii).len - 2*L,...
                         'X',tracks(ii).X(L+1:end-L),...
                         'Y',tracks(ii).Y(L+1:end-L),...
                         'T',tracks(ii).T(L+1:end-L),...
                         'U',u,...
                         'V',v);
  end
  disp('Done.')

end

function pos = FindParticles(image, threshold, logs)

  s = size(image);
  
  % identify the local maxima that are above threshold  
  maxes = find(image >= threshold & ...
               image > circshift(image,[0 1]) & ...
               image > circshift(image,[0 -1]) & ...
               image > circshift(image,[1 0]) & ...
               image > circshift(image,[-1 0]));
             
  % now turn these into subscripts
  [x,y] = ind2sub(s, maxes);
  
  % throw out unreliable maxes in the outer ring
  good = find(x~=1 & y~=1 & x~=s(1) & y~=s(2));
  x = x(good);
  y = y(good);
  
  % find the horizontal positions
  
  % look up the logarithms of the relevant image intensities
  z1 = logs(image(sub2ind(s,x-1,y)) + 1)';
  z2 = logs(image(sub2ind(s,x,y)) + 1)';
  z3 = logs(image(sub2ind(s,x+1,y)) + 1)';
  
  % compute the centers
  xcenters = -0.5 * (z1.*(-2*x-1) + z2.*(4*x) + z3.*(-2*x+1)) ./ ...
                    (z1 + z3 - 2*z2);
                  
  % do the same for the vertical position
  z1 = logs(image(sub2ind(s,x,y-1)) + 1)';
  z3 = logs(image(sub2ind(s,x,y+1)) + 1)';
  ycenters = -0.5 * (z1.*(-2*y-1) + z2.*(4*y) + z3.*(-2*y+1)) ./ ...
                    (z1 + z3 - 2*z2);
                 
  % make sure we have no bad points
  good = find(isfinite(xcenters) & isfinite(ycenters));
                  
  % fix up the funny coordinate system used by matlab
  pos = [ycenters(good), xcenters(good)];

end