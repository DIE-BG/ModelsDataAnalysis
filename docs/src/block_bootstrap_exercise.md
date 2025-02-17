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

Let $X_1^*$ be an observation drawn at random from the original $N$ observations. $X_1^*$ is the first observation of the first block, so $X_1^*$ is defined as $X_1^* = X_{I_1}$. The next one has probability $p$ of being randomly drawn from the original $N$ observations and probability $1-p$ of being $X_2^* = X_{I_1+1}$ which be the next one in the original sample. The process ends when the $X_j^*$ is to be drawn at random from the original $N$ observations, this being the first observation of the next block. The iteration ends when we have a number of blocks that being joined together formed a "single block" of length $N$. So, the average length of the blocks is $\frac{1}{p}$,

<!-- Use the folder images/simulation_study for the images of the presentation -->