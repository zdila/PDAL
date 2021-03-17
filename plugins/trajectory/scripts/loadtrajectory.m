function dat = loadtrajectory(esttraj, groundtruth)
% Read in ground truth (sbet) data (t, x, y, z) in file groundtruth and
% estimated trajectory (t, x, y, z) in file esttraj.  Interpolate ground
% truth to times for estimated trajectory return an array
% [t, xest, yest, zest, xgt, ygt, zgt]

% After
%
%   dat = loadtrajectory('traj0.txt', 'sbet.txt');
%
% compute the RMS/max error as follows
%
%   err = dat(:,2:4)-dat(:,5:7);
%   verr = err(:,3);
%   herr = hypot(err(:,1), err(:,2));
%   maxerr = max(abs([herr, verr]))
%   rmserr = sqrt(mean([herr, verr].^2))
%
% plot the error with, e.g.,
%
%   t = dat(:,1);
%   plot(t, err(:,1), t, err(:,2), t, err(:,3));
%
% N.B. ";" suppress output from a command.  Without ";" the results are
% printed.

  t = load(esttraj);
  g = load(groundtruth);
  % linear interpolation is OK in sbet data
  dat = [t(:,1:4), interp1(g(:,1),g(:,2:4),t(:,1))];
  err = dat(:,2:4)-dat(:,5:7);
  verr = err(:,3);
  herr = hypot(err(:,1), err(:,2));
  maxerr = max(abs([herr, verr]));
  rmserr = sqrt(mean([herr, verr].^2));
  printf(['max errs: horiz = %.3f vert = %.3f\n', ...
          'rms errs: horiz = %.3f vert = %.3f\n'], ...
         maxerr,rmserr);
end
