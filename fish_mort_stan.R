# Full model
stancode1 <- "data{
int<lower=1> N;
int<lower=1> N_AGD;
int<lower=1> N_Des;
int<lower=1> N_cyst;
int<lower=1> N_PRV;
int<lower=1> N_pox;
int<lower=1> N_Ten;
int<lower=1> N_sample;
int<lower=1> N_cage;
int mort[N];
real AGD[N_AGD];
real Des[N_Des];
real cyst[N_cyst];
real PRV[N_PRV];
real pox[N_pox];
real Ten[N_Ten];
real logpop[N];
int sample[N];
int cage[N];
real Temp[N];
real time_w[N];
matrix[N_sample,N_sample] Dmat;
real Tense[N];
real poxse[N];
real PRVse[N];
real cystse[N];
real Desse[N];
real AGDse[N];
}



parameters{
vector[N_AGD] AGD_est;
vector[N_AGD] Des_est;
vector[N_AGD] cyst_est;
vector[N_AGD] PRV_est;
vector[N_AGD] pox_est;
vector[N_AGD] Ten_est;
real a;
real bagd;
real bdes;
real bcyst;
real bprv;
real bpox;
real bten;
real bagd_des;
real bagd_cyst;
real bagd_prv;
real bagd_pox;
real bagd_ten;
real btemp;
real btime_w;
real<lower=0> theta;
vector[N_sample] a_sample;
real<lower=0> etasq;
real<lower=0> rhosq;
real<lower=0> sigma;
}

model{
matrix[N_sample,N_sample] SIGMA_Dmat;
vector[N] lambda;
rhosq ~ cauchy( 0 , 1 );
etasq ~ cauchy( 0 , 1 );
sigma ~ cauchy( 0 , 1 );
for ( i in 1:(N_sample-1) )
for ( j in (i+1):N_sample ) {
SIGMA_Dmat[i,j] = etasq*exp(-rhosq*pow(Dmat[i,j],2));
SIGMA_Dmat[j,i] = SIGMA_Dmat[i,j];
}
for ( k in 1:N_sample )
SIGMA_Dmat[k,k] = etasq + sigma;
a_sample ~ multi_normal( rep_vector(0,N_sample) , SIGMA_Dmat );
theta ~ exponential( 1 );
btime_w ~ normal( 0 , 0.5 );
btemp ~ normal( 0 , 0.5 );
bagd_ten ~ normal( 0 , 0.5 );
bagd_pox ~ normal( 0 , 0.5 );
bagd_prv ~ normal( 0 , 0.5 );
bagd_cyst ~ normal( 0 , 0.5 );
bagd_des ~ normal( 0 , 0.5 );
bten ~ normal( 0 , 0.5 );
bpox ~ normal( 0 , 0.5 );
bprv ~ normal( 0 , 0.5 );
bcyst ~ normal( 0 , 0.5 );
bdes ~ normal( 0 , 0.5 );
bagd ~ normal( 0 , 0.5 );
a ~ normal( 0 , 5 );
for ( i in 1:N ) {
lambda[i] = logpop[i] + a + a_sample[sample[i]] + bagd * AGD_est[i] + bdes * Des_est[i] 
+ bcyst * cyst_est[i] + bprv * PRV_est[i] + bpox * pox_est[i] + bten * Ten_est[i] + bagd_des * AGD_est[i] * Des_est[i] 
+ bagd_cyst * AGD_est[i] * cyst_est[i] + bagd_prv * AGD_est[i] * PRV_est[i] + bagd_pox * AGD_est[i] * pox_est[i] 
+ bagd_ten * AGD_est[i] * Ten_est[i] + btemp * Temp[i] + btime_w * time_w[i];
lambda[i] = exp(lambda[i]);
}
Ten ~ normal( Ten_est , Tense );
pox ~ normal( pox_est , poxse );
PRV ~ normal( PRV_est , PRVse );
cyst ~ normal( cyst_est , cystse );
Des ~ normal( Des_est , Desse );
AGD ~ normal( AGD_est , AGDse );
mort ~ neg_binomial_2( lambda , theta );
}

generated quantities{
matrix[N_sample,N_sample] SIGMA_Dmat;
vector[N] lambda;
vector[N] log_lik;
vector[N] y_pred;

for ( i in 1:(N_sample-1) ) 
for ( j in (i+1):N_sample ) {
SIGMA_Dmat[i,j] = etasq*exp(-rhosq*pow(Dmat[i,j],2));
SIGMA_Dmat[j,i] = SIGMA_Dmat[i,j];
}

for ( k in 1:N_sample )
SIGMA_Dmat[k,k] = etasq + sigma;

for ( i in 1:N ) {
lambda[i] = logpop[i] + a + a_sample[sample[i]] + bagd * AGD_est[i] + bdes * Des_est[i] 
+ bcyst * cyst_est[i] + bprv * PRV_est[i] + bpox * pox_est[i] + bten * Ten_est[i] + bagd_des * AGD_est[i] * Des_est[i] 
+ bagd_cyst * AGD_est[i] * cyst_est[i] + bagd_prv * AGD_est[i] * PRV_est[i] + bagd_pox * AGD_est[i] * pox_est[i] 
+ bagd_ten * AGD_est[i] * Ten_est[i] + btemp * Temp[i] + btime_w * time_w[i];
lambda[i] = exp(lambda[i]);

log_lik[i] = neg_binomial_2_lpmf( mort[i] | lambda[i] , theta );
y_pred[i] = neg_binomial_2_rng(lambda[i], theta);

}

}


"
# ZINB model...comverges like shit. make it run well
stancode2= "data {
int<lower=1> N;
int<lower=1> N_AGD;
int<lower=1> N_Des;
int<lower=1> N_cyst;
int<lower=1> N_PRV;
int<lower=1> N_pox;
int<lower=1> N_Ten;
int<lower=1> N_sample;
int<lower=1> N_cage;
int mort[N];
real AGD[N_AGD];
real Des[N_Des];
real cyst[N_cyst];
real PRV[N_PRV];
real pox[N_pox];
real Ten[N_Ten];
real logpop[N];
int sample[N];
int cage[N];
real Temp[N];
real time_w[N];
matrix[N_sample,N_sample] Dmat;
real Tense[N];
real poxse[N];
real PRVse[N];
real cystse[N];
real Desse[N];
real AGDse[N];
}

parameters {
real a;
vector[N_sample] a_sample;
real b;
vector[N_sample] b_sample;
real<lower=0> sigma_a;
real<lower=0> sigma_b;
real<lower=0> theta;
real<lower=0> etasq_a;
real<lower=0> rhosq_a;
real<lower=0> etasq_b;
real<lower=0> rhosq_b;
}


model {
matrix[N_sample,N_sample] SIGMA_Dmat_a;
matrix[N_sample,N_sample] SIGMA_Dmat_b;
vector[N] p;
vector[N] lambda;

sigma_a ~ cauchy( 0 , 1);
sigma_b ~ cauchy( 0 , 1);
rhosq_a ~ cauchy( 0 , 1 );
etasq_a ~ cauchy( 0 , 1 );
rhosq_b ~ cauchy( 0 , 1 );
etasq_b ~ cauchy( 0 , 1 );

for ( i in 1:(N_sample-1) )
for ( j in (i+1):N_sample ) {
SIGMA_Dmat_a[i,j] = etasq_a*exp(-rhosq_a*pow(Dmat[i,j],2));
SIGMA_Dmat_a[j,i] = SIGMA_Dmat_a[i,j];
}
for ( k in 1:N_sample )
SIGMA_Dmat_a[k,k] = etasq_a + sigma_a;
a_sample ~ multi_normal( rep_vector(0,N_sample) , SIGMA_Dmat_a );

for ( i in 1:(N_sample-1) )
for ( j in (i+1):N_sample ) {
SIGMA_Dmat_b[i,j] = etasq_b*exp(-rhosq_b*pow(Dmat[i,j],2));
SIGMA_Dmat_b[j,i] = SIGMA_Dmat_b[i,j];
}
for ( k in 1:N_sample )
SIGMA_Dmat_b[k,k] = etasq_b + sigma_b;
b_sample ~ multi_normal( rep_vector(0,N_sample) , SIGMA_Dmat_b );

a ~ normal( 0 , 5);
b ~ normal( 0 , 5);
theta ~ exponential( 1 );

for ( i in 1:N ) {
lambda[i] = logpop[i] + a + a_sample[sample[i]];
lambda[i] = exp(lambda[i]);
}

for ( i in 1:N ) {
p[i] = inv_logit(b + b_sample[sample[i]]);
}
for (n in 1:N) {
if (mort[n] == 0){
target += log_sum_exp(bernoulli_lpmf(1 | p[n]),
bernoulli_lpmf(0 | p[n])
+ neg_binomial_2_lpmf(mort[n] | lambda[n], theta));
}else{
target += bernoulli_lpmf(0 | p[n])
+ neg_binomial_2_lpmf(mort[n] | lambda[n], theta);
}
}
}

generated quantities {
real lambda[N];
real p[N];
real pred[N];
real resid[N];
real log_lik[N];

for ( i in 1:N ) {
lambda[i] = logpop[i] + a + a_sample[sample[i]];
lambda[i] = exp(lambda[i]);
}

for ( i in 1:N ) {
p[i] = inv_logit(b + b_sample[sample[i]]);
}

for (n in 1:N) {
if (mort[n] == 0){
log_lik[n] = log_sum_exp(bernoulli_lpmf(1 | p[n]),
bernoulli_lpmf(0 | p[n])
+ neg_binomial_2_lpmf(mort[n] | lambda[n], theta));
}else{
log_lik[n] = bernoulli_lpmf(0 | p[n])
+ neg_binomial_2_lpmf(mort[n] | lambda[n], theta);
}
pred[n] = (1-p[n]) * lambda[n];
resid[n] = (mort[n] - pred[n])/sqrt(pred[n] + theta * pred[n]^2);
}
}
"

# ZINB model no Gaussian effects
stancode3= "data {
int<lower=1> N;
int<lower=1> N_AGD;
int<lower=1> N_Des;
int<lower=1> N_cyst;
int<lower=1> N_PRV;
int<lower=1> N_pox;
int<lower=1> N_Ten;
int<lower=1> N_sample;
int<lower=1> N_cage;
int mort[N];
real AGD[N_AGD];
real Des[N_Des];
real cyst[N_cyst];
real PRV[N_PRV];
real pox[N_pox];
real Ten[N_Ten];
real logpop[N];
int sample[N];
int cage[N];
real Temp[N];
real time_w[N];
matrix[N_sample,N_sample] Dmat;
real Tense[N];
real poxse[N];
real PRVse[N];
real cystse[N];
real Desse[N];
real AGDse[N];
}

parameters {
real a;
vector[N_sample] a_sample_raw;
real b;
vector[N_sample] b_sample_raw;
real<lower=0> sigma_a;
real<lower=0> sigma_b;
real<lower=0> theta;
}

transformed parameters{
vector[N_sample] a_sample;
vector[N_sample] b_sample;

a_sample = a_sample_raw*sigma_a;
b_sample = b_sample_raw*sigma_b;
}

model {
vector[N] p;
vector[N] lambda;
sigma_a ~ cauchy( 0 , 1);
sigma_b ~ cauchy( 0 , 1);
a_sample_raw ~ normal( 0 , 1 );
b_sample_raw ~ normal( 0 , 1 );
a ~ normal( 0 , 5);
b ~ normal( 0 , 5);
theta ~ exponential( 1 );

for ( i in 1:N ) {
lambda[i] = logpop[i] + a + a_sample[sample[i]];
lambda[i] = exp(lambda[i]);
}

for ( i in 1:N ) {
p[i] = inv_logit(b + b_sample[sample[i]]);
}
for (n in 1:N) {
if (mort[n] == 0){
target += log_sum_exp(bernoulli_lpmf(1 | p[n]),
bernoulli_lpmf(0 | p[n])
+ neg_binomial_2_lpmf(mort[n] | lambda[n], theta));
}else{
target += bernoulli_lpmf(0 | p[n])
+ neg_binomial_2_lpmf(mort[n] | lambda[n], theta);
}
}
}

generated quantities {
real lambda[N];
real p[N];
real pred[N];
real resid[N];
real log_lik[N];

for ( i in 1:N ) {
lambda[i] = logpop[i] + a + a_sample[sample[i]];
lambda[i] = exp(lambda[i]);
}

for ( i in 1:N ) {
p[i] = inv_logit(b + b_sample[sample[i]]);
}

for (n in 1:N) {
if (mort[n] == 0){
log_lik[n] = log_sum_exp(bernoulli_lpmf(1 | p[n]),
bernoulli_lpmf(0 | p[n])
+ neg_binomial_2_lpmf(mort[n] | lambda[n], theta));
}else{
log_lik[n] = bernoulli_lpmf(0 | p[n])
+ neg_binomial_2_lpmf(mort[n] | lambda[n], theta);
}
pred[n] = (1-p[n]) * lambda[n];
resid[n] = (mort[n] - pred[n])/sqrt(pred[n] + theta * pred[n]^2);
}
}
"