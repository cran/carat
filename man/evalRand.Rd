\name{evalRand}
\alias{evalRand}
\alias{evalRand.careval}
\title{Evaluation of Randomization Procedures}
\usage{
\method{evalRand}{careval}(data, method = "HuHuCAR", N = 500, ...)
}
\description{
Evaluates a specific randomization procedure based on several different quantities of imbalances. 
}
\arguments{
  \item{data}{a dataframe. A row of the dataframe contains the covariate profile of a patient.}
  \item{N}{iteration number.}
  \item{method}{the randomization to be used in allocating patients. The default randomization \dQuote{\code{HuHuCAR}} uses Hu and Hu's general covariate-adaptive randomization; the alternatives are \dQuote{\code{PocSimMIN}}, \dQuote{\code{StrBCD}}, \dQuote{\code{StrPBR}}, \dQuote{\code{DoptBCD}} as well as  \dQuote{\code{AdjBCD}}.}
  \item{\dots}{arguments to be passed to methods. It depends on the \code{method} and the following arguments are accepted:
  \describe{
  \item{omega}{vector of weights at the overall, within-stratum as well as maginal levels. It is required that at least one element is larger than 0. Notice that \code{omega} is only needed when \code{HuHuCAR} is to be assessed.}
  \item{weight}{vector of weights for all involved margins. It is required that at least one element is NOT 0 and \code{length(weight) = cov_num}. Notice that \code{weight} is only needed when \code{PocSimMIN} is to be assessed.}
  \item{p}{probabilities of assigning one patinet to treatment 1. \code{p} should be larger than \code{1/2} to obtain balance. Notice that \code{p} is only needed when \code{"HuHuCAR"},  \code{"PocSimMIN"} and \code{"StrBCD"} are to be assessed.}
  \item{a}{a design parameter. As \code{a} goes to \eqn{\infty}, the design becomes more deteministic. Notice that \code{a} is only needed when \code{"AdjBCD"} is to be assessed.}
  \item{bsize}{block size for stratified permuted block randomization. It should be multiple of \code{2}. Notice that \code{bsize} is only needed when \code{"StrPBR"} is to be assessed.}
  }}
}
\details{
    \code{data} is designed for \code{N} times using \code{method}.
}
\value{
It returns an object of \code{\link{class}} \code{"careval"}. 

The functions \code{\link{print}} is used to obtain results. The generic accessor functions \code{Assig}, \code{Diff}, \code{data}, \code{All strata} and so on extract various useful features of the value returned by that function. 

An object of class \code{"careval"} is a list containing at least the following components: 

\item{N}{number of patients.}
\item{Assig}{a \code{n*N} matrix containing assignments for each patient for N iterationsd. }
\item{Imb}{matrix containing maximum, 95\%-quantile, median, loss and mean of absolute imbalances at overall, each within-stratum and each marginal levels.}
\item{Within-strt. by num of pats}{a vector containing mean absolute imbalances of the strata with \eqn{i} patients falling in, where \eqn{i = 1, \dots, bsize}.}
\item{Data Type}{data type. \code{Real} or \code{Simulated}. }
}
\seealso{
  See \code{\link{evalRand.sim}} to evaluate a randomization procedure by generating simulated data.
}
\references{
Atkinson A C. \emph{Optimum biased coin designs for sequential clinical trials with prognostic factors}[J]. Biometrika, 1982, 69(1): 61-67.

Baldi Antognini A, Zagoraiou M. \emph{The covariate-adaptive biased coin design for balancing clinical trials in the presence of prognostic factors}[J]. Biometrika, 2011, 98(3): 519-535.

Hu Y, Hu F. \emph{Asymptotic properties of covariate-adaptive randomization}[J]. The Annals of Statistics, 2012, 40(3): 1794-1815.

Pocock S J, Simon R. \emph{Sequential treatment assignment with balancing for prognostic factors in the controlled clinical trial}[J]. Biometrics, 1975: 103-115.

Shao J, Yu X, Zhong B. \emph{A theory for testing hypotheses under covariate-adaptive randomization}[J]. Biometrika, 2010, 97(2): 347-360. 

Zelen M. \emph{The randomization and stratification of patients to clinical trials}[J]. Journal of chronic diseases, 1974, 27(7): 365-375.
}
\examples{
# a simple use
## Access by real data
## create a dataframe
df <- data.frame("gender" = sample(c("female", "male"), 1000, TRUE, c(1 / 3, 2 / 3)), 
                 "age" = sample(c("0-30", "30-50", ">50"), 1000, TRUE), 
                 "jobs" = sample(c("stu.", "teac.", "others"), 1000, TRUE))
Res <- evalRand(data = df, method = "HuHuCAR", N = 500, 
                omega = c(1, 2, rep(1, ncol(df))), p = 0.85)
## view the output
Res
\donttest{
## view all patients' assignments
Res$Assig}

## Assess by simulated data
cov_num <- 3
level_num <- c(2, 3, 5)
pr <- c(0.35, 0.65, 0.25, 0.35, 0.4, 0.25, 0.15, 0.2, 0.15, 0.25)
n <- 1000
N <- 50
omega = c(1, 2, 1, 1, 2)
# assess Hu and Hu's procedure with the same group of patients
Res.sim <- evalRand.sim(n = n, N = N, Replace = FALSE, cov_num = cov_num, 
                        level_num = level_num, pr = pr, method = "HuHuCAR", 
                        omega, p = 0.85)
\donttest{ 
## Compare four procedures
cov_num <- 3
level_num <- c(2, 10, 2)
pr <- c(rep(0.5, times = 2), rep(0.1, times = 10), rep(0.5, times = 2))
n <- 100
N <- 200 # <<adjust according to CPU
bsize <- 4
## set weights for HuHuCAR
omega <- c(1, 2, rep(1, cov_num)); 
## set weights for PocSimMIN
weight = rep(1, cov_num); 
## set biased probability
p = 0.80
# assess Hu and Hu's procedure
RH <- evalRand.sim(n = n, N = N, Replace = FALSE, cov_num = cov_num, 
                   level_num = level_num, pr = pr, method = "HuHuCAR", 
                   omega = omega, p = p)
# assess Pocock and Simon's method
RPS <- evalRand.sim(n = n, N = N, Replace = FALSE, cov_num = cov_num, 
                    level_num = level_num, pr = pr, method = "PocSimMIN", 
                    weight, p = p)
# assess Shao's procedure
RS <- evalRand.sim(n = n, N = N, Replace = FALSE, cov_num = cov_num, 
                   level_num = level_num, pr = pr, method = "StrBCD", 
                   p = p)
# assess stratified randomization
RSR <- evalRand.sim(n = n, N = N, Replace = FALSE, cov_num = cov_num, 
                    level_num = level_num, pr = pr, method = "StrPBR", 
                    bsize)

# create containers
C_All = C_M = C_O = matrix(NA, nrow = 4, ncol = 5)
colnames(C_All) = colnames(C_M) = colnames(C_O) =
  c("max", "95\%quan", "med", "mean", "loss")
C_WS = matrix(NA, nrow = 4, ncol = bsize + 5)
colnames(C_WS) = c("max", "95\%quan", "med", "mean", "loss",
                   "num = 1", "num = 2", "num = 3", "num = 4")
rownames(C_All) = rownames(C_M) = rownames(C_O) = rownames(C_WS) = 
  c("HH", "PocSim", "Shao", "StraRand")

# access the imbalaces at overall, within-stratum, and marginal levels
C_All[1, ] = RH$`All Imbalances`
C_All[2, ] = RPS$`All Imbalances`
C_All[3, ] = RS$`All Imbalances`
C_All[4, ] = RSR$`All Imbalances`
# view the result
C_All

# assess the overall imbalance
C_O[1, ] = RH$`Overall Imbalances`
C_O[2, ] = RPS$`Overall Imbalances`
C_O[3, ] = RS$`Overall Imbalances`
C_O[4, ] = RSR$`Overall Imbalances`
# view the result
C_O

# assess the marginal imbalances
C_M[1, ] = RH$`Marginal Imbalances`
C_M[2, ] = RPS$`Marginal Imbalances`
C_M[3, ] = RS$`Marginal Imbalances`
C_M[4, ] = RSR$`Marginal Imbalances`
# view the result
C_M

# assess the within-stratum imbalances
C_WS[1, ] = RH$`Within-strt. Imbalances`
C_WS[2, ] = RPS$`Within-strt. Imbalances`
C_WS[3, ] = RS$`Within-strt. Imbalances`
C_WS[4, ] = RSR$`Within-strt. Imbalances`
# view the result
C_WS

# Compare the four procedures through plots
meth = rep(c("Hu", "PS", "Shao", "STR"), times = 4)
shape <- rep(1 : 4, times = 4)
crt <- rep(1 : 4, each = 4)
crt_c <- rep(c("All", "O", "M", "WS"), each = 4)
mean <- c(C_All[, 4], C_O[, 4], C_M[, 4], C_WS[, 4])
loss <- c(C_All[, 5], C_O[, 5], C_M[, 5], C_WS[, 5])
df_1 <- data.frame(meth, shape, crt, crt_c, mean, loss)

require(ggplot2)
require(gridExtra)
p1 <- ggplot(df_1, aes(x = meth, y = mean, color = crt_c, group = crt,
                       linetype = crt_c, shape = crt_c)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  xlab("method") +
  ylab("absolute mean") +
  theme(plot.title = element_text(hjust = 0.5))

# analyze within-stratum imbalances especially
df_2 <- as.data.frame(c(C_WS[, 6], C_WS[, 7], C_WS[, 8], C_WS[, 9]))
colnames(df_2) = "val"
df_2$num <- rep(c("num = 1", "num = 2", "num = 3", "num = 4"), each = 4)
df_2$meth <- meth

p2 <- ggplot(df_2, aes(x = num, y = val, color = meth, group = meth,
                       linetype = meth, shape = meth)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  xlab("numbers of patients for each strata") +
  ylab("absolute within-stratum mean") +
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(p1, p2, ncol = 1)
}
}
