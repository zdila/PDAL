function res = trajresidue(traj,pulses)
  sensor = interp1(traj(:,1),traj(:,2:4),pulses(:,1));
  r = pulses(:,2:4);
  d = pulses(:,5);
  v = pulses(:,6:8);
  p = sensor-r;
  v2 = v(:,1) .* v(:,1) + v(:,2) .* v(:,2);
  % See proj.mac for derivation
  M00 = (v(:,1).*v(:,1) .* v(:,3) + v(:,2).*v(:,2)) ./ v2;
  M11 = (v(:,2).*v(:,2) .* v(:,3) + v(:,1).*v(:,1)) ./ v2;
  M01 = -(1 - v(:,3)) .* v(:,1) .* v(:,2) ./ v2;
  M10 = M01;
  M02 = -v(:,1);
  M12 = -v(:,2);
  M20 = v(:,1);
  M21 = v(:,2);
  M22 = v(:,3);
  % Scale first two rows to get projection error at v
  M00 = d .* M00;
  M01 = d .* M01;
  M02 = d .* M02;
  M10 = d .* M10;
  M11 = d .* M11;
  M12 = d .* M12;
  xt = M00 .* p(:,1) + M01 .* p(:,2) + M02 .* p(:,3);
  yt = M10 .* p(:,1) + M11 .* p(:,2) + M12 .* p(:,3);
  zt = M20 .* p(:,1) + M21 .* p(:,2) + M22 .* p(:,3);
  res = hypot(xt./zt, yt./zt);
end
