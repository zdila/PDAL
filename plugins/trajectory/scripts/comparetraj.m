function comparetraj(dat1, dat2)
% Plot errors for trajectories compute two different ways.  dat[12] is
% the output of loadtrajectory
  legend1 = 'method 1';
  legend2 = 'method 2';
  % legend1 = 'ceres  ';
  % legend2 = 'G+McG  ';
  mult = 1000;
  if mult == 1000
    unit = 'mm';
  else
    mult = 1;
    unit = 'm';
  end
  t1 = dat1(:,1);
  err1 = dat1(:,2:4) - dat1(:,5:7);
  t2 = dat2(:,1);
  err2 = dat2(:,2:4) - dat2(:,5:7);
  torg = floor(min([t1;t2]));
  t1 = t1-torg;
  t2 = t2-torg;
  for i = 1:3
    subplot(3,1,i);
    plot(t1,err1(:,i)*mult,'b-', t2,err2(:,i)*mult,'b--');
    xlabel('t (s)')
    if i == 1
      % axis([0,70,-400,400]);
      ylabel(['x error (' unit ')'])
    elseif i == 2
      % axis([0,70,-400,400]);
      ylabel(['y error (' unit ')'])
    else
      % axis([0,70,-1000,1000]);
      ylabel(['z error (' unit ')']);
    end
    legend(legend1, legend2);
  end
end

