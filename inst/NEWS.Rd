\name{NEWS}
\title{News for Package \pkg{carat}}
\newcommand{\cpkg}{\href{https://CRAN.R-project.org/package=#1}{\pkg{#1}}}
\newcommand{\issue}{\href{https://github.com/yexiaoqingruc/#1}{(#1)}}

\section{version 1.0 (May 4th, 2020)}{
  \itemize{
    \item This version is the released version of \pkg{carat}.
    \item Due to the limitation of R on \pkg{OpenMP}, the \pkg{OpenMP}-supported version, named \pkg{caratOMP}, can be loaded via \issue{carat}. 
  }
}

\section{version 1.1 (May 15th, 2020)}{
  \itemize{
    \item A small modification was added to the function \code{getData}. 
  }
}

\section{version 1.2 (August 31th, 2020)}{
  \itemize{
     \item A our own defined class \code{plotrandtest} was added and applied to the function \code{rand.test}. 
     \item The option \code{plot} was delected in \code{rand.test}. 
     \item Character strings were directly used as arguments values for \code{method} and \code{test} of the testing functions and \code{getData}. 
     \item We standardized the usual function to check if an object has from a specific class in the function \code{compRand}, using \code{is()}. 
     \item We corrected the man pages of the corresponding randomization functions, such as \code{PocSimMIN}. 
     \item The underlying S3 methods are not exported to users. 
     \item The generation of names vector in \code{getData} was altered. 
     \item All of the modifications are consistent with the \pkg{OpenMP}-supported version in \issue{carat}. 
  }
}

\section{version 1.3 (September 10th, 2020)}{
\itemize{
    \item A small modification was added to the function \code{evalRand}. 
    \item The modification is consistent with the \pkg{OpenMP}-supported version in \issue{carat}. 
  }
}

\section{version 1.4 (September 26th, 2020)}{
\itemize{
    \item References were added to the man files. 
    \item The modification is consistent with the \pkg{OpenMP}-supported version in \issue{carat}. 
  }
}
