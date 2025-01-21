function [values,isterminal,direction] = myevent(t,y,tstart)
 %  Don't let t cross zero...a dummy "event" to illustrate how 
 %  one might handle other events in conjunction with the time
 %  constraint.  Not necessary, but I put it in just in case.
 values(1) = t;
 %  Don't let integration go for more than 90 seconds.
 values(2) = toc(tstart) < 180;
 isterminal = true(size(values));
 direction = zeros(size(values));
end