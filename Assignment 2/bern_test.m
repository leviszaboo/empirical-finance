% Bernoulli coverage test in Matlab
function res=bern_test(p,v)
  a=p^(sum(v))*(1-p)^(length(v)-sum(v));
  b=(sum(v)/length(v))^(sum(v))*(1-(sum(v)/length(v)))^(length(v)-sum(v));
  res=-2*log(a/b);
end