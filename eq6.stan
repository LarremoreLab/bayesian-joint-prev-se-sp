data {
  int field_pos;
  int field_n;
  int true_neg;
  int total_neg;
  int true_pos;
  int total_pos;
}
parameters {
  real<lower = 0, upper = 1> p;
  real<lower = 0, upper = 1> spec;
  real<lower = 0, upper = 1> sens;
}
model {
  real p_sample;
  p_sample = p * sens + (1 - p) * (1 - spec);
  field_pos ~ binomial(field_n, p_sample);
  true_neg ~ binomial(total_neg, spec);
  true_pos ~ binomial(total_pos, sens);
}