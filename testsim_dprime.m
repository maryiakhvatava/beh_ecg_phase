function [d,beta,c] = testsim_dprime(pHit,pFA)
% pHit	- hit probability
% pFA	- probability of False Alarms
% http://en.wikipedia.org/wiki/D%27

% outputs

% d'    = z(pHit) - z(pFA)

% ratio of the height of the signal distribution to the noise distribution for the value of the threshold (criterion)
% log(beta) = c * d'

% criterion for yes/no: if positive - conservative (more no responses), if negative - liberal (more yes responses)
% then the ideal observer
% criterion c =  -(zHit + zFA) / 2


% for special cases when rates are 0 or 100 see:
% https://stats.stackexchange.com/questions/134779/d-prime-with-100-hit-rate-probability-and-0-false-alarm-probability

% Convert to Z scores
zHit = norminv(pHit);
zFA  = norminv(pFA);
% Calculate d-prime
d = zHit - zFA;

% Calculate BETA (the criterion value)
% if (nargout > 1)
	yHit = normpdf(zHit);
	yFA  = normpdf(zFA);
	beta = yHit ./ yFA;
    
    % T. DeCarlo / Journal of Mathematical Psychology 56 (2012) 196�207, page 197 - WHY? - is this the special case of yes/no?
    % c = -norminv(pFA); % which is same as c = -zFA
    
    % standard formula
    c = -0.5*(zHit + zFA);
    disp(['d-prime: ', num2str(d)]);
    disp(['criterion: ', num2str(c)]);
end

% Note that d' derived from  T. DeCarlo / Journal of Mathematical Psychology 56 (2012) 196�207, Eq. (4) results in the following:
% b (bias) = sqrt(2)/2*(zHit + zFA)
% d' = sqrt(2)/2*(zHit - zFA)