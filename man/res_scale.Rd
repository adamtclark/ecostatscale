% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ecostatscale.R
\name{res_scale}
\alias{res_scale}
\title{Resilience scaling function}
\usage{
res_scale(
  mvar_b,
  murho_b_abundance,
  mucov_b_abundance = NULL,
  msd_b,
  murho_b_disturbance,
  mucov_b_disturbance = NULL,
  b = 1,
  B,
  lambda
)
}
\arguments{
\item{mvar_b}{Mean abundance variance observed at scale b}

\item{murho_b_abundance}{Mean Pearson correlation coefficient of abundance values observed at scale b}

\item{mucov_b_abundance}{Mean covariance of abundance values observed at scale b. Ignored unless murho_b_abundance is NULL Defaults to NULL.}

\item{msd_b}{Mean disturbance standard deviation observed at scale b}

\item{murho_b_disturbance}{Mean Pearson correlation coefficient of disturbances observed at scale b}

\item{mucov_b_disturbance}{Mean covariance of disturbances observed at scale b. Ignored unless murho_b_abundance is NULL Defaults to NULL.}

\item{b}{Size of observed scale. Defaults to 1.}

\item{B}{Larger scale being extrapolated to (e.g. total number of species, or size of patch B relative to b)}

\item{lambda}{Mean disturbance frequency.}
}
\value{
Extrapolated median resilience at scale of M species.
}
\description{
Extrapolate resilience observed at the scale of a single spatial or ecological scale b (e.g. a patch or species)
to a larger scale, B (e.g. functional group or landscape).
}
\examples{
# extrapolate from scale of 1 species to 10 species
res_scale(mvar_b = 0.25, murho_b_abundance = -0.034, msd_b = sqrt(0.1),
           murho_b_disturbance = 0, B = 30, lambda=1)

# plot relationship for groups of 1 to 30 species
plot(1:30, res_scale(mvar_b = 0.25, murho_b_abundance = -0.034, msd_b = sqrt(0.1),
      murho_b_disturbance = 0, B = 1:30, lambda=1),
      xlab="ecological scale", ylab="resilience, r", type="b")
}
