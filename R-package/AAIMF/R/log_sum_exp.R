
# log(sum(exp()) is computed in a numerically stable way.
# It relies on the identity:
# log Sum_k exp (a_k) = m + log Sum_k exp(a_k - m), where m = max_k a_k

log_sum_exp <- function(a)
{
  m <- max(a)
  m + log(sum(exp(a - m)))
}
