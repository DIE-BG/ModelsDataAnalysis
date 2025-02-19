# Block Bootstrap Simulation Exercise

## Bootstraping methods for dependent data
In the literature there are some methods for dependent data. Some of them are named as "first generation methods" and take as a basis a block that is of fixed length or changes at each iteration and then resampling theese bloks until they replicates the length of the original sample.

## Moving Block Bootstrap (MBB)
We start from a given sample size $X_n \equiv \{X_1,...,X_n\}$. And then we define a block of length $l$, with $l$ being an even number that satisfies $1\leq l \leq n$. Once we have defined the block length we form different blocks of length $l$ contained in $X_n$ as follows.

$B_1 = (X_1, X_2, ... ,X_l)$
$B_2 = (X_2, X_3, ... ,X_{l+1})$
$... ... ...$
$B_N = (X_{n-l+1}, ... ,X_n)$

$N = n-l+1$ denotes the number of total blocks of length $l$ that will be formed and on which the bootstrap will be applied. Allowing the overlapping of the blocks in each iteration until the original sample length is reached. After the n bootstrap observations are formed, we proceed to estimate the statistic of interest (mean, variance, autocorrelation function or some other).

## Stationary Block Bootstrap (SBB)
As Moving Block Bootstrap we have a sample of $N$ observations of which we assume that they are stationary and time dependent. This method has the property of producing stationary pseudo time series.

The SBB method shares with the MBB the characteristic of working with blocks of time series that when joined together generate a “single block” of length $N$ as that of the original sample. The difference between SBB and MBB is that the block length in SBB changes at each iteration and in MBB the length is fixed. The algorithm as follows:

Let $X_1^*$ be an observation drawn at random from the original $N$ observations. $X_1^*$ is the first observation of the first block, so $X_1^*$ is defined as $X_1^* = X_{I_1}$. The next one has probability $p$ of being randomly drawn from the original $N$ observations and probability $1-p$ of being $X_2^* = X_{I_1+1}$ which be the next one in the original sample. The process ends when the $X_j^*$ is to be drawn at random from the original $N$ observations, this being the first observation of the next block. The iteration ends when we have a number of blocks that being joined together formed a "single block" of length $N$. So, the average length of the blocks is $\frac{1}{p}$.

## Basic simulation exercise
In our exercise, the objetive is to generate pseu time series using a bootstrap method suitable for our sample series. The pseudo time series must adequately replicate the statistical characteristics of the original series. All series are year-on-year rates of change and these are:

1. real GDP of USA
2. PCE core inflation
3. Effective Federal Funds Rate
4. Domestic real GDP
5. Total domestic inflation
6. Domestic core inflation
7. exchange rate
8. Monetary base
9. Monetary policy rate

In order to generate robust results in this exercise, we generated 10,000 pseudo time series across a window of all possible block lengths. To check the correct replication of the statistical characteristics we calculated some statistics that are important for our purposes. These are:

1. mean
2. standartd deviation
3. Autocorrelation fuction
4. Variance-covariance matrix

For each iteration (pseudo time series for given block length), we calculate the above statistics and compare them to the original sample statistics. We calculate a modified mean square error to measure the error between the sample statistic and the original sample statistic. The modified MSE is definen as follows:

$MSE=\frac{1}{B}\sum_{i=1}^{B}\frac{\hat{\theta}_i-\theta}{\sigma_l}$

where:
- $i \text{ Number of the iteration}$
- $B \text{ Total numer of iterations}$
- $l \text{ Length of the block}$
- $\hat{\theta} \text{ Sample statistic}$
- $\theta \text{ Original sample statistic}$
- $\sigma_l \text{ Standar deviation from de sample distribution of } \hat{\theta}$

The modified MSE takes values $\geq 1$ for each statistic so it is very convenient as it allows us to aggregate each one and obtain a single measure that captures the error for each serie and statistic. Finally the bootstrap method that obtains the minimum MSE for all possible blocks of length $l$ will be the best bootstrap method for the evaluation exercise.

## Results for the sample mean
The excecise campares the performance of two methods Moving and Stationary Block Bootstrap for replicate the sample mean for all possible block lengths which, given our sample size (91 observations) $l$ should be a number that satisfies $1\leq l \leq 91$. For each $l$ we generate 10,000 pseudo time series with the Moving and Stationary methods respectively and estimate the mean for each and then compare these means to the historical sample mean. The process described above will be de same for standard deviation, autocorrelation fuction and variance-covariance matrix.

The original MSE[^1] measures the error for one series and probably the error in all series if we decide to average all series, but it is not the best way to measure the error if we have different statistics (we have mean, standard deviation and variance-covariance matrix), since we want to generate a measure that embeds the error in all statistics, we implemented the modified MSE which consists of a normalized MSE by the standard deviation of the historical estimator of the statistic for each iteration. This modification allows us to aggregate the error for all the statistics. The modified MSE for the historical estimator of the mean for all series is the average of the modified MSE of all series for all possible block lengths. The original MSE of all series is also obtained in the same way.

![](images/simulation_study/mean/Original_MSE.png)
![](images/simulation_study/mean/Modified_MSE.png)

The first graphic shows the average of the original MSE of all series and the second one the modified MSE. We note that in both cases the best method is the Stationary, however for all possible block lengths the historical mean estimator of Stationary Block Bootstrap, for modified MSE, is significantly more unbiased than the estimator of the historical mean of the Moving Block Bootstrap, so the Stationary method is the best method to replicate the historical mean in each series.

As we have said, we prefer the modified MSE since it will allow us to aggregate the error across the various statistics.

## Results for the sample standard deviation





[^1]: $MSE=\frac{1}{B}\sum_{i=1}^{B}\hat{\theta}_i-\theta$
<!-- Use the folder images/simulation_study for the images of the presentation -->