# Full model
stancode1 <-"
data{
int<lower=1> N;
int<lower=1> N_sample;
int<lower=1> N_cage;
int histo[N];
real AGD[N];
real Des[N];
real cyst[N];
real PRV[N];
real pox[N];
real Ten[N];
int sample[N];
int cage[N];
real Temp[N];
real time_w[N];
matrix[N_sample,N_sample] Dmat;
}
parameters{
ordered[3] cutpoints;
real bt;
real bwater_t;
real bagd;
real bdes;
real bcyst;
real bpox;
real bten;
real bagd_des;
real bagd_cyst;
real bagd_pox;
real bagd_ten;
vector[N_sample] a_sample;
real<lower=0> etasq;
real<lower=0> rhosq;
real<lower=0> sigma;
}

model{
matrix[N_sample,N_sample] SIGMA_Dmat;
vector[N] phi;
sigma ~ cauchy( 0 , 1 );
rhosq ~ cauchy( 0 , 1 );
etasq ~ cauchy( 0 , 1 );
for ( i in 1:(N_sample-1) )
for ( j in (i+1):N_sample ) {
SIGMA_Dmat[i,j] = etasq*exp(-rhosq*pow(Dmat[i,j],2));
SIGMA_Dmat[j,i] = SIGMA_Dmat[i,j];
}
for ( k in 1:N_sample )
SIGMA_Dmat[k,k] = etasq + sigma;
a_sample ~ multi_normal( rep_vector(0,N_sample) , SIGMA_Dmat );
cutpoints ~ normal( 0 , 5 );
bagd_ten ~ normal( 0 , 0.5 );
bagd_pox ~ normal( 0 , 0.5 );
bagd_cyst ~ normal( 0 , 0.5 );
bagd_des ~ normal( 0 , 0.5 );
bten ~ normal( 0 , 0.5 );
bpox ~ normal( 0 , 0.5 );
bcyst ~ normal( 0 , 0.5 );
bdes ~ normal( 0 , 0.5 );
bagd ~ normal( 0 , 0.5 );
bwater_t ~ normal( 0 , 0.5 );
bt ~ normal( 0 , 0.5 );
for ( i in 1:N ) {
phi[i] = bt * Temp[i] + bwater_t * time_w[i] + bagd * AGD[i] 
+ bdes * Des[i] + bcyst * cyst[i] + bpox * pox[i] + bten * Ten[i] 
+ bagd_des * AGD[i] * Des[i] + bagd_cyst * AGD[i] * cyst[i] 
+ bagd_pox * AGD[i] * pox[i] + bagd_ten * AGD[i] * Ten[i] 
+ a_sample[sample[i]];
}
for ( i in 1:N )
histo[i] ~ ordered_logistic( phi[i] , cutpoints );
}
generated quantities{
matrix[N_sample,N_sample] SIGMA_Dmat;
vector[N] phi;
real log_lik[N];

for ( i in 1:(N_sample-1) )
for ( j in (i+1):N_sample ) {
SIGMA_Dmat[i,j] = etasq*exp(-rhosq*pow(Dmat[i,j],2));
SIGMA_Dmat[j,i] = SIGMA_Dmat[i,j];
}
for ( k in 1:N_sample )
SIGMA_Dmat[k,k] = etasq + sigma;

for ( i in 1:N ) {
phi[i] = bt * Temp[i] + bwater_t * time_w[i] + bagd * AGD[i] 
+ bdes * Des[i] + bcyst * cyst[i] + bpox * pox[i] + bten * Ten[i] 
+ bagd_des * AGD[i] * Des[i] + bagd_cyst * AGD[i] * cyst[i] 
+ bagd_pox * AGD[i] * pox[i] + bagd_ten * AGD[i] * Ten[i] 
+ a_sample[sample[i]];
}


for ( i in 1:N )
log_lik[i] = ordered_logistic_lpmf( histo[i] | phi[i] , cutpoints );
}
"
# Null model
stancode2 <-"
data{
int<lower=1> N;
int<lower=1> N_sample;
int<lower=1> N_cage;
int histo[N];
real AGD[N];
real Des[N];
real cyst[N];
real PRV[N];
real pox[N];
real Ten[N];
int sample[N];
int cage[N];
real Temp[N];
real time_w[N];
matrix[N_sample,N_sample] Dmat;
}
parameters{
ordered[3] cutpoints;
vector[N_sample] a_sample;
real<lower=0> etasq;
real<lower=0> rhosq;
real<lower=0> sigma;
}

model{
matrix[N_sample,N_sample] SIGMA_Dmat;
vector[N] phi;
sigma ~ cauchy( 0 , 1 );
rhosq ~ cauchy( 0 , 1 );
etasq ~ cauchy( 0 , 1 );
for ( i in 1:(N_sample-1) )
for ( j in (i+1):N_sample ) {
SIGMA_Dmat[i,j] = etasq*exp(-rhosq*pow(Dmat[i,j],2));
SIGMA_Dmat[j,i] = SIGMA_Dmat[i,j];
}
for ( k in 1:N_sample )
SIGMA_Dmat[k,k] = etasq + sigma;
a_sample ~ multi_normal( rep_vector(0,N_sample) , SIGMA_Dmat );
cutpoints ~ normal( 0 , 5 );
for ( i in 1:N ) {
phi[i] = a_sample[sample[i]];
}
for ( i in 1:N )
histo[i] ~ ordered_logistic( phi[i] , cutpoints );
}
generated quantities{
matrix[N_sample,N_sample] SIGMA_Dmat;
vector[N] phi;
real log_lik[N];

for ( i in 1:(N_sample-1) )
for ( j in (i+1):N_sample ) {
SIGMA_Dmat[i,j] = etasq*exp(-rhosq*pow(Dmat[i,j],2));
SIGMA_Dmat[j,i] = SIGMA_Dmat[i,j];
}
for ( k in 1:N_sample )
SIGMA_Dmat[k,k] = etasq + sigma;

for ( i in 1:N ) {
phi[i] = a_sample[sample[i]];
}

for ( i in 1:N )
log_lik[i] = ordered_logistic_lpmf( histo[i] | phi[i] , cutpoints );
}
" 
# temp, time_fw, and agd
stancode3 <-"
data{
int<lower=1> N;
int<lower=1> N_sample;
int<lower=1> N_cage;
int histo[N];
real AGD[N];
real Des[N];
real cyst[N];
real PRV[N];
real pox[N];
real Ten[N];
int sample[N];
int cage[N];
real Temp[N];
real time_w[N];
matrix[N_sample,N_sample] Dmat;
}
parameters{
ordered[3] cutpoints;
real bt;
real bwater_t;
real bagd;
vector[N_sample] a_sample;
real<lower=0> etasq;
real<lower=0> rhosq;
real<lower=0> sigma;
}

model{
matrix[N_sample,N_sample] SIGMA_Dmat;
vector[N] phi;
sigma ~ cauchy( 0 , 1 );
rhosq ~ cauchy( 0 , 1 );
etasq ~ cauchy( 0 , 1 );
for ( i in 1:(N_sample-1) )
for ( j in (i+1):N_sample ) {
SIGMA_Dmat[i,j] = etasq*exp(-rhosq*pow(Dmat[i,j],2));
SIGMA_Dmat[j,i] = SIGMA_Dmat[i,j];
}
for ( k in 1:N_sample )
SIGMA_Dmat[k,k] = etasq + sigma;
a_sample ~ multi_normal( rep_vector(0,N_sample) , SIGMA_Dmat );
cutpoints ~ normal( 0 , 5 );
bagd ~ normal( 0 , 0.5 );
bwater_t ~ normal( 0 , 0.5 );
bt ~ normal( 0 , 0.5 );
for ( i in 1:N ) {
phi[i] = bt * Temp[i] + bwater_t * time_w[i] + bagd * AGD[i] + a_sample[sample[i]];
}
for ( i in 1:N )
histo[i] ~ ordered_logistic( phi[i] , cutpoints );
}
generated quantities{
matrix[N_sample,N_sample] SIGMA_Dmat;
vector[N] phi;
real log_lik[N];

for ( i in 1:(N_sample-1) )
for ( j in (i+1):N_sample ) {
SIGMA_Dmat[i,j] = etasq*exp(-rhosq*pow(Dmat[i,j],2));
SIGMA_Dmat[j,i] = SIGMA_Dmat[i,j];
}
for ( k in 1:N_sample )
SIGMA_Dmat[k,k] = etasq + sigma;

for ( i in 1:N ) {
phi[i] = bt * Temp[i] + bwater_t * time_w[i] + bagd * AGD[i] + a_sample[sample[i]];
}

for ( i in 1:N )
log_lik[i] = ordered_logistic_lpmf( histo[i] | phi[i] , cutpoints );
}
"

# Desmozoan, agd, temp, time fw

stancode4 <-"
data{
int<lower=1> N;
int<lower=1> N_sample;
int<lower=1> N_cage;
int histo[N];
real AGD[N];
real Des[N];
real cyst[N];
real PRV[N];
real pox[N];
real Ten[N];
int sample[N];
int cage[N];
real Temp[N];
real time_w[N];
matrix[N_sample,N_sample] Dmat;
}
parameters{
ordered[3] cutpoints;
real bagd;
real bt;
real bwater_t;
real bdes;
vector[N_sample] a_sample;
real<lower=0> etasq;
real<lower=0> rhosq;
real<lower=0> sigma;
}

model{
matrix[N_sample,N_sample] SIGMA_Dmat;
vector[N] phi;
sigma ~ cauchy( 0 , 1 );
rhosq ~ cauchy( 0 , 1 );
etasq ~ cauchy( 0 , 1 );
for ( i in 1:(N_sample-1) )
for ( j in (i+1):N_sample ) {
SIGMA_Dmat[i,j] = etasq*exp(-rhosq*pow(Dmat[i,j],2));
SIGMA_Dmat[j,i] = SIGMA_Dmat[i,j];
}
for ( k in 1:N_sample )
SIGMA_Dmat[k,k] = etasq + sigma;
a_sample ~ multi_normal( rep_vector(0,N_sample) , SIGMA_Dmat );
cutpoints ~ normal( 0 , 5 );
bt ~ normal( 0 , 0.5 );
bwater_t ~ normal( 0 , 0.5 );
bagd ~ normal( 0 , 0.5 );
bdes ~ normal( 0 , 0.5 );
for ( i in 1:N ) {
phi[i] = bt*Temp[i] + bwater_t * time_w[i] + bagd * AGD[i] + bdes * Des[i] + a_sample[sample[i]];
}
for ( i in 1:N )
histo[i] ~ ordered_logistic( phi[i] , cutpoints );
}
generated quantities{
matrix[N_sample,N_sample] SIGMA_Dmat;
vector[N] phi;
real log_lik[N];

for ( i in 1:(N_sample-1) )
for ( j in (i+1):N_sample ) {
SIGMA_Dmat[i,j] = etasq*exp(-rhosq*pow(Dmat[i,j],2));
SIGMA_Dmat[j,i] = SIGMA_Dmat[i,j];
}
for ( k in 1:N_sample )
SIGMA_Dmat[k,k] = etasq + sigma;

for ( i in 1:N ) {
phi[i] = bt*Temp[i] + bwater_t * time_w[i] + bagd * AGD[i] + bdes * Des[i] + a_sample[sample[i]];
}

for ( i in 1:N )
log_lik[i] = ordered_logistic_lpmf( histo[i] | phi[i] , cutpoints );
}
"

# temp and FW treatment

stancode5 <-"
data{
int<lower=1> N;
int<lower=1> N_sample;
int<lower=1> N_cage;
int histo[N];
real AGD[N];
real Des[N];
real cyst[N];
real PRV[N];
real pox[N];
real Ten[N];
int sample[N];
int cage[N];
real Temp[N];
real time_w[N];
matrix[N_sample,N_sample] Dmat;
}
parameters{
ordered[3] cutpoints;
real bt;
real bwater_t;
vector[N_sample] a_sample;
real<lower=0> etasq;
real<lower=0> rhosq;
real<lower=0> sigma;
}

model{
matrix[N_sample,N_sample] SIGMA_Dmat;
vector[N] phi;
sigma ~ cauchy( 0 , 1 );
rhosq ~ cauchy( 0 , 1 );
etasq ~ cauchy( 0 , 1 );
for ( i in 1:(N_sample-1) )
for ( j in (i+1):N_sample ) {
SIGMA_Dmat[i,j] = etasq*exp(-rhosq*pow(Dmat[i,j],2));
SIGMA_Dmat[j,i] = SIGMA_Dmat[i,j];
}
for ( k in 1:N_sample )
SIGMA_Dmat[k,k] = etasq + sigma;
a_sample ~ multi_normal( rep_vector(0,N_sample) , SIGMA_Dmat );
cutpoints ~ normal( 0 , 5 );
bt ~ normal( 0 , 0.5 );
bwater_t ~ normal( 0 , 0.5 );
for ( i in 1:N ) {
phi[i] = bt * Temp[i] + bwater_t * time_w[i] + a_sample[sample[i]];
}
for ( i in 1:N )
histo[i] ~ ordered_logistic( phi[i] , cutpoints );
}
generated quantities{
matrix[N_sample,N_sample] SIGMA_Dmat;
vector[N] phi;
real log_lik[N];

for ( i in 1:(N_sample-1) )
for ( j in (i+1):N_sample ) {
SIGMA_Dmat[i,j] = etasq*exp(-rhosq*pow(Dmat[i,j],2));
SIGMA_Dmat[j,i] = SIGMA_Dmat[i,j];
}
for ( k in 1:N_sample )
SIGMA_Dmat[k,k] = etasq + sigma;

for ( i in 1:N ) {
phi[i] = bt * Temp[i] + bwater_t * time_w[i] + a_sample[sample[i]];
}

for ( i in 1:N )
log_lik[i] = ordered_logistic_lpmf( histo[i] | phi[i] , cutpoints );
}
"

# b. c. cysticola, agd, temp, time fw

stancode6 <-"
data{
int<lower=1> N;
int<lower=1> N_sample;
int<lower=1> N_cage;
int histo[N];
real AGD[N];
real Des[N];
real cyst[N];
real PRV[N];
real pox[N];
real Ten[N];
int sample[N];
int cage[N];
real Temp[N];
real time_w[N];
matrix[N_sample,N_sample] Dmat;
}
parameters{
ordered[3] cutpoints;
real bt;
real bwater_t;
real bagd;
real bcyst;
vector[N_sample] a_sample;
real<lower=0> etasq;
real<lower=0> rhosq;
real<lower=0> sigma;
}

model{
matrix[N_sample,N_sample] SIGMA_Dmat;
vector[N] phi;
sigma ~ cauchy( 0 , 1 );
rhosq ~ cauchy( 0 , 1 );
etasq ~ cauchy( 0 , 1 );
for ( i in 1:(N_sample-1) )
for ( j in (i+1):N_sample ) {
SIGMA_Dmat[i,j] = etasq*exp(-rhosq * pow(Dmat[i,j],2));
SIGMA_Dmat[j,i] = SIGMA_Dmat[i,j];
}
for ( k in 1:N_sample )
SIGMA_Dmat[k,k] = etasq + sigma;
a_sample ~ multi_normal( rep_vector(0,N_sample) , SIGMA_Dmat );
cutpoints ~ normal( 0 , 5 );
bt ~ normal( 0 , 0.5 );
bwater_t ~ normal( 0 , 0.5 );
bagd ~ normal( 0 , 0.5 );
bcyst ~ normal( 0 , 0.5 );
for ( i in 1:N ) {
phi[i] = bt * Temp[i] + bwater_t * time_w[i] + bagd * AGD[i] + bcyst * cyst[i] 
+ a_sample[sample[i]];
}
for ( i in 1:N )
histo[i] ~ ordered_logistic( phi[i] , cutpoints );
}
generated quantities{
matrix[N_sample,N_sample] SIGMA_Dmat;
vector[N] phi;
real log_lik[N];

for ( i in 1:(N_sample-1) )
for ( j in (i+1):N_sample ) {
SIGMA_Dmat[i,j] = etasq*exp(-rhosq * pow(Dmat[i,j],2));
SIGMA_Dmat[j,i] = SIGMA_Dmat[i,j];
}
for ( k in 1:N_sample )
SIGMA_Dmat[k,k] = etasq + sigma;

for ( i in 1:N ) {
phi[i] = bt * Temp[i] + bwater_t * time_w[i] + bagd * AGD[i] + bcyst * cyst[i] 
+ a_sample[sample[i]];
}

for ( i in 1:N )
log_lik[i] = ordered_logistic_lpmf( histo[i] | phi[i] , cutpoints );
}
"

# Main effects model

stancode7 <-"
data{
int<lower=1> N;
int<lower=1> N_sample;
int<lower=1> N_cage;
int histo[N];
real AGD[N];
real Des[N];
real cyst[N];
real PRV[N];
real pox[N];
real Ten[N];
int sample[N];
int cage[N];
real Temp[N];
real time_w[N];
matrix[N_sample,N_sample] Dmat;
}
parameters{
ordered[3] cutpoints;
real bt;
real bwater_t;
real bagd;
real bdes;
real bcyst;
real bpox;
real bten;
vector[N_sample] a_sample;
real<lower=0> etasq;
real<lower=0> rhosq;
real<lower=0> sigma;
}

model{
matrix[N_sample,N_sample] SIGMA_Dmat;
vector[N] phi;
sigma ~ cauchy( 0 , 1 );
rhosq ~ cauchy( 0 , 1 );
etasq ~ cauchy( 0 , 1 );
for ( i in 1:(N_sample-1) )
for ( j in (i+1):N_sample ) {
SIGMA_Dmat[i,j] = etasq*exp(-rhosq*pow(Dmat[i,j],2));
SIGMA_Dmat[j,i] = SIGMA_Dmat[i,j];
}
for ( k in 1:N_sample )
SIGMA_Dmat[k,k] = etasq + sigma;
a_sample ~ multi_normal( rep_vector(0,N_sample) , SIGMA_Dmat );
cutpoints ~ normal( 0 , 5 );
bten ~ normal( 0 , 0.5 );
bpox ~ normal( 0 , 0.5 );
bcyst ~ normal( 0 , 0.5 );
bdes ~ normal( 0 , 0.5 );
bagd ~ normal( 0 , 0.5 );
bwater_t ~ normal( 0 , 0.5 );
bt ~ normal( 0 , 0.5 );
for ( i in 1:N ) {
phi[i] = bt * Temp[i] + bwater_t * time_w[i] + bagd * AGD[i] 
+ bdes * Des[i] + bcyst * cyst[i] + bpox * pox[i] + bten * Ten[i] 
+ a_sample[sample[i]];
}
for ( i in 1:N )
histo[i] ~ ordered_logistic( phi[i] , cutpoints );
}
generated quantities{
matrix[N_sample,N_sample] SIGMA_Dmat;
vector[N] phi;
real log_lik[N];

for ( i in 1:(N_sample-1) )
for ( j in (i+1):N_sample ) {
SIGMA_Dmat[i,j] = etasq*exp(-rhosq*pow(Dmat[i,j],2));
SIGMA_Dmat[j,i] = SIGMA_Dmat[i,j];
}
for ( k in 1:N_sample )
SIGMA_Dmat[k,k] = etasq + sigma;

for ( i in 1:N ) {
phi[i] = bt * Temp[i] + bwater_t * time_w[i] + bagd * AGD[i] 
+ bdes * Des[i] + bcyst * cyst[i] + bpox * pox[i] + bten * Ten[i] 
+ a_sample[sample[i]];
}

for ( i in 1:N )
log_lik[i] = ordered_logistic_lpmf( histo[i] | phi[i] , cutpoints );
}
"

# pox, agd, temp, time fw

stancode9 <-"
data{
int<lower=1> N;
int<lower=1> N_sample;
int<lower=1> N_cage;
int histo[N];
real AGD[N];
real Des[N];
real cyst[N];
real PRV[N];
real pox[N];
real Ten[N];
int sample[N];
int cage[N];
real Temp[N];
real time_w[N];
matrix[N_sample,N_sample] Dmat;
}
parameters{
ordered[3] cutpoints;
real bt;
real bwater_t;
real bagd;
real bpox;
vector[N_sample] a_sample;
real<lower=0> etasq;
real<lower=0> rhosq;
real<lower=0> sigma;
}

model{
matrix[N_sample,N_sample] SIGMA_Dmat;
vector[N] phi;
sigma ~ cauchy( 0 , 1 );
rhosq ~ cauchy( 0 , 1 );
etasq ~ cauchy( 0 , 1 );
for ( i in 1:(N_sample-1) )
for ( j in (i+1):N_sample ) {
SIGMA_Dmat[i,j] = etasq*exp(-rhosq*pow(Dmat[i,j],2));
SIGMA_Dmat[j,i] = SIGMA_Dmat[i,j];
}
for ( k in 1:N_sample )
SIGMA_Dmat[k,k] = etasq + sigma;
a_sample ~ multi_normal( rep_vector(0,N_sample) , SIGMA_Dmat );
cutpoints ~ normal( 0 , 5 );
bpox ~ normal( 0 , 0.5 );
bagd ~ normal( 0 , 0.5 );
bwater_t ~ normal( 0 , 0.5 );
bt ~ normal( 0 , 0.5 );
for ( i in 1:N ) {
phi[i] = bt * Temp[i] + bwater_t * time_w[i] + bagd * AGD[i] 
+ bpox * pox[i] + a_sample[sample[i]];
}
for ( i in 1:N )
histo[i] ~ ordered_logistic( phi[i] , cutpoints );
}
generated quantities{
matrix[N_sample,N_sample] SIGMA_Dmat;
vector[N] phi;
real log_lik[N];

for ( i in 1:(N_sample-1) )
for ( j in (i+1):N_sample ) {
SIGMA_Dmat[i,j] = etasq*exp(-rhosq*pow(Dmat[i,j],2));
SIGMA_Dmat[j,i] = SIGMA_Dmat[i,j];
}
for ( k in 1:N_sample )
SIGMA_Dmat[k,k] = etasq + sigma;

for ( i in 1:N ) {
phi[i] = bt * Temp[i] + bwater_t * time_w[i] + bagd * AGD[i] 
+ bpox * pox[i] + a_sample[sample[i]];
}

for ( i in 1:N )
log_lik[i] = ordered_logistic_lpmf( histo[i] | phi[i] , cutpoints );
}
"

# ten, agd, temp, time fw

stancode10 <-"
data{
int<lower=1> N;
int<lower=1> N_sample;
int<lower=1> N_cage;
int histo[N];
real AGD[N];
real Des[N];
real cyst[N];
real PRV[N];
real pox[N];
real Ten[N];
int sample[N];
int cage[N];
real Temp[N];
real time_w[N];
matrix[N_sample,N_sample] Dmat;
}
parameters{
ordered[3] cutpoints;
real bt;
real bwater_t;
real bagd;
real bten;
vector[N_sample] a_sample;
real<lower=0> etasq;
real<lower=0> rhosq;
real<lower=0> sigma;
}

model{
matrix[N_sample,N_sample] SIGMA_Dmat;
vector[N] phi;
sigma ~ cauchy( 0 , 1 );
rhosq ~ cauchy( 0 , 1 );
etasq ~ cauchy( 0 , 1 );
for ( i in 1:(N_sample-1) )
for ( j in (i+1):N_sample ) {
SIGMA_Dmat[i,j] = etasq*exp(-rhosq*pow(Dmat[i,j],2));
SIGMA_Dmat[j,i] = SIGMA_Dmat[i,j];
}
for ( k in 1:N_sample )
SIGMA_Dmat[k,k] = etasq + sigma;
a_sample ~ multi_normal( rep_vector(0,N_sample) , SIGMA_Dmat );
cutpoints ~ normal( 0 , 5 );
bten ~ normal( 0 , 0.5 );
bagd ~ normal( 0 , 0.5 );
bwater_t ~ normal( 0 , 0.5 );
bt ~ normal( 0 , 0.5 );
for ( i in 1:N ) {
phi[i] = bt * Temp[i] + bwater_t * time_w[i] + bagd * AGD[i] 
+ bten * Ten[i] + a_sample[sample[i]];
}
for ( i in 1:N )
histo[i] ~ ordered_logistic( phi[i] , cutpoints );
}
generated quantities{
matrix[N_sample,N_sample] SIGMA_Dmat;
vector[N] phi;
real log_lik[N];

for ( i in 1:(N_sample-1) )
for ( j in (i+1):N_sample ) {
SIGMA_Dmat[i,j] = etasq*exp(-rhosq*pow(Dmat[i,j],2));
SIGMA_Dmat[j,i] = SIGMA_Dmat[i,j];
}
for ( k in 1:N_sample )
SIGMA_Dmat[k,k] = etasq + sigma;

for ( i in 1:N ) {
phi[i] = bt * Temp[i] + bwater_t * time_w[i] + bagd * AGD[i] 
+ bten * Ten[i] + a_sample[sample[i]];
}

for ( i in 1:N )
log_lik[i] = ordered_logistic_lpmf( histo[i] | phi[i] , cutpoints );
}
"
