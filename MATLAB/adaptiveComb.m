function score = adaptiveComb(score1, score2, beta1, beta2)
% beta1 = 3;
% beta2 = 5;

% B =10;
% A = 15;
% alpha = 1./(1+exp(-(A*score1-B)));

alpha = 1./(1+ beta1*score1.^beta2);
score = score1.^(1-alpha) .* score2.^alpha;


end
