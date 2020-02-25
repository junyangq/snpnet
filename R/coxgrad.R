### Computes the gradient of Cox partial likelihood
### f is fitted function from glmnet at a particular lambda
### time is death or censoring time
### d is death indicator; d=0 means censored, d=1 means death
### w is a weight vector of non-negative weights, which will be normalized to sum to 1
coxgrad <- function(f, time, d, w, eps=0.00001){
  if(missing(w))w=rep(1,length(f))
  time = time - d*eps
  d = d * w
  f=scale(f,TRUE,FALSE)
  o = order(time)
  ordered_time = time[o]
  r = rank(ordered_time, ties.method="min")
  rm = rank(ordered_time, ties.method="max")
  ef=(w * exp(f))[o]
  rskden=rev(cumsum(rev(ef)))[r]
  grad =  cumsum(d[o]/rskden)[rm]
  grad[o] = d[o] -  ef * grad
  grad
}


  