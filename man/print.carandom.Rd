\name{print.carandom}
\alias{print.carandom}
\alias{print.careval}
\alias{print.carcomp}
\alias{print.carseq}
\alias{AdBCDOne}
\alias{AdjBCD_BT}
\alias{AdjBCD_BT_In}
\alias{AdjBCD_BT_power}
\alias{AdjBCD_CT_power}
\alias{AdjBCD_In}
\alias{AdjBCD_RT}
\alias{AdjBCD_RT_In}
\alias{AdjBCD_RT_power}
\alias{AdjBCD_getData}
\alias{AtBCDOne}
\alias{BBCDname}
\alias{Bpert}
\alias{CTT}
\alias{CTT_In}
\alias{C_AdjustBCD}
\alias{C_AtkinBCD}
\alias{C_HPS}
\alias{C_StrR}
\alias{C_RAdjustBCD}
\alias{C_RAtkinBCD}
\alias{C_RStrR}
\alias{C_RSummarize}
\alias{C_Summarize}
\alias{DoptBCD_BT}
\alias{DoptBCD_BT_In}
\alias{DoptBCD_BT_power}
\alias{DoptBCD_CT_power}
\alias{DoptBCD_In}
\alias{DoptBCD_RT}
\alias{DoptBCD_RT_In}
\alias{DoptBCD_RT_power}
\alias{DoptBCD_getData}
\alias{EnterCov}
\alias{EnterCovPrf}
\alias{EnterLev}
\alias{EnterSigCov}
\alias{EnterSigLev}
\alias{EnterWeig}
\alias{HPSOne}
\alias{StrROne}
\alias{HuHuCAR_BT}
\alias{HuHuCAR_BT_In}
\alias{HuHuCAR_BT_power}
\alias{HuHuCAR_CT_power}
\alias{HuHuCAR_RT}
\alias{HuHuCAR_RT_In}
\alias{HuHuCAR_RT_power}
\alias{HuHuCAR_getData}
\alias{PStrGen}
\alias{PocSimMIN_BT}
\alias{PocSimMIN_BT_In}
\alias{PocSimMIN_BT_power}
\alias{PocSimMIN_CT_power}
\alias{PocSimMIN_RT}
\alias{PocSimMIN_RT_In}
\alias{PocSimMIN_RT_power}
\alias{PocSimMIN_getData}
\alias{Preprocess}
\alias{Prob_S}
\alias{Rprod}
\alias{StrBCD_BT}
\alias{StrBCD_BT_In}
\alias{StrBCD_BT_power}
\alias{StrBCD_CT_power}
\alias{StrBCD_RT}
\alias{StrBCD_RT_In}
\alias{StrBCD_RT_power}
\alias{StrBCD_getData}
\alias{StrPBR_BT}
\alias{StrPBR_BT_In}
\alias{StrPBR_BT_power}
\alias{StrPBR_CT_power}
\alias{StrPBR_RT}
\alias{StrPBR_RT_In}
\alias{StrPBR_RT_power}
\alias{StrPBR_getData}
\alias{genData_sim}
\alias{nameString}
\alias{pathout}
\alias{C_RHPS}
\title{
Print Methods for Sequences, Evaluation and Compariason for Covariate-adaptive Randomization Procedures
}
\description{
Printing objects of class \code{"carandom"}, \code{"careval"}, \code{"carcomp"} or \code{carseq} respectively for different functions, by simple \code{\link{print}} method.
}
\usage{
\method{print}{carandom}(x, digits = getOption("digits"), prefix = "\t", ...)

\method{print}{careval}(x, digits = getOption("digits"), prefix = "\t", ...)

\method{print}{carcomp}(x, digits = getOption("digits"), prefix = "\t", ...)

\method{print}{carseq}(x, digits = getOption("digits"), prefix = "\t", ...)
}
\arguments{
  \item{x}{objects of class \code{"carandom"}, \code{"careval"}, \code{"carcomp"}, or \code{carseq}. }
  \item{digits}{number of significant digits to be used.}
  \item{prefix}{string, passed to \code{\link{strwrap}} for displaying
    the \code{method} component of the \code{carandom} object, \code{careval} object, \code{carcomp} object and \code{carseq} object.}
    \item{\dots}{further arguments to be passed to or from methods.}
}
\value{
the argument \code{x}, invisibly, as for all \code{\link{print}}
  methods.
}
\seealso{
\code{\link{HuHuCAR}}, \code{\link{evalRand}}, \code{\link{compRand}}, 
\code{\link{HuHuCAR.ui}}
}
\keyword{ carandom }
\keyword{ careval }
\keyword{ carcomp }
\keyword{ carseq }
