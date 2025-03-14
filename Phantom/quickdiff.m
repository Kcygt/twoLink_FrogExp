function dvar=quickdiff(tt,var,method)
% diff_of_var=QUICKDIFF(timevec,var,method)
% A quick way to calculate the differential. 
% Estimated as $dx/dt \approx \frac{x_{n+1}-x_n}{t_{n+1}-t_n)$ and
% this is interpolated back on to the time vector. Thus the time vector
% does not need to have equal increments.
%
% Best to convert the time vector to numbers before the call
%
% Most of the code is handling time structures, the actual work is
% done with midpts=... and dvar=... This could be simplified!
%
% 'method' is not implimented for now, but could just be the method
% passed to interp1
%
% can var be a matrix?

    if isduration(tt) || isdatetime(tt)
        warning('Best to convert the time vector to numbers');
    end

    deltaTimes=diff(tt);
    if isduration(deltaTimes)
        if ~strcmp(deltaTimes.Format,'s') % will assume seconds
            warning(sprintf('Time vector format is %s',deltaTimes.Format))
        end
        if isdatetime(tt)
            tt=timeofday(tt);
        end
        tt=seconds(tt); % finally a number?        
        midpts=(tt(1:end-1)+tt(2:end))/2; % find midpoints in time vector
        dvar=interp1(midpts,diff(var)./seconds(deltaTimes),tt); % interpolate diffs to the original time vector
    else % assume time is just numbers
        midpts=(tt(1:end-1)+tt(2:end))/2; % find midpoints in time vector

        dvar=interp1(midpts,diff(var)./diff(tt),tt); % interpolate diffs to the original time vector
    end
end

%% test code
%w=1.5;tt=0:.01:20;xx=sin(w*tt);
%ttasdatetime=datetime('now')+seconds(tt);
%ttasduration=seconds(tt);
%
%figure(1);plot(tt,w*sin(w*tt),tt,quickdiff(tt,xx))
%figure(2);plot(tt,w*sin(w*tt),tt,quickdiff(ttasduration,xx))
%figure(3);plot(tt,w*sin(w*tt),tt,quickdiff(ttasdatetime,xx))