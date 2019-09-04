#' @template survival_doc_template
#' @param formula a \code{formula}; the response 
#' must be the time scale to compute survival time function estimates
#' over, e.g. \code{fot ~ sex}. Variables on the right-hand side of the formula
#' separated by \code{+} are considered stratifying variables, for which 
#' estimates are computed separately. May contain usage of \code{adjust()} 
#' --- see Details and Examples.
#' @param data since popEpi 0.4.0, a \code{data.frame}
#' containing variables used in \code{formula} and other arguments.
#' \code{aggre} objects are recommended as they contain information on any
#' time scales and are therefore safer; for creating \code{aggre} objects see
#' \code{\link{as.aggre}} when your data is already aggregated and \code{aggre}
#' for aggregating split \code{Lexis} objects.
#' 
#' @param surv.breaks a vector of breaks on the 
#' survival time scale. Optional if \code{data} is an \code{aggre} object
#' and mandatory otherwise. Must define each intended interval;
#' e.g. \code{surv.breaks = 0:5} when data has intervals defined by 
#' breaks \code{seq(0, 5, 1/12)} will aggregate to wider intervals first.
#' It is generally recommended (and sufficient; 
#' see Seppa, Dyban and Hakulinen (2015)) to use monthly
#' intervals where applicable.
#' 
#' @param n variable containing counts of subjects at-risk at the start of a 
#' time interval; e.g. \code{n = "at.risk"}. 
#' Required when \code{surv.method = "lifetable"}.
#' \link[=flexible_argument]{Flexible input}.
#' 
#' @param d variable(s) containing counts of subjects experiencing an event. 
#' With only one type of event, e.g. \code{d = "deaths"}. With multiple types of 
#' events (for CIF or cause-specific survival estimation), supply e.g.
#' \code{d = c("canD", "othD")}. If the survival time function to be estimated
#' does not use multiple types of events, supplying more than one variable
#' to \code{d} simply causes the variables to be added together. 
#' Always required. \link[=flexible_argument]{Flexible input}.
#' 
#' @param n.cens variable containing counts of subjects censored during a 
#' survival time interval; E.g. \code{n.cens = "alive"}.
#' Required when \code{surv.method = "lifetable"}. 
#' \link[=flexible_argument]{Flexible input}.

#' @param pyrs variable containing total subject-time accumulated within a 
#' survival time interval; E.g. \code{pyrs = "pyrs"}. 
#' Required when \code{surv.method = "hazard"}. Flexible input.

#' @param d.exp variable denoting total "expected numbers of events" 
#' (typically computed \code{pyrs * pop.haz}, where 
#' \code{pop.haz} is the expected hazard level) 
#' accumulated within a survival time interval; E.g. \code{pyrs = "pyrs"}.
#' Required when computing EdererII relative survivals or 
#' CIFs based on excess counts of events. Flexible input.

#' @param n.pp variable containing total Pohar-Perme weighted counts of
#' subjects at risk in an interval,
#' supplied as argument \code{n} is supplied. 
#' Computed originally on the subject
#' level as analogous to \code{pp * as.integer(status == "at-risk")}.
#' Required when \code{relsurv.method = "pp"}. Flexible input.
#' 
#' @param d.pp variable(s) containing Pohar-Perme weighted counts of events,
#' supplied as argument \code{d} is supplied. Computed originally on the subject
#' level as analogous to \code{pp * as.integer(status == some_event)}.
#' Required when \code{relsurv.method = "pp"}. Flexible input.

#' @param d.pp.2 variable(s) containing total Pohar-Perme 
#' "double-weighted" counts of events,
#' supplied as argument \code{d} is supplied. Computed originally on the subject
#' level as analogous to \code{pp * pp * as.integer(status == some_event)}.
#' Required when \code{relsurv.method = "pp"}. Flexible input.

#' @param n.cens.pp variable containing total Pohar-Perme weighted counts 
#' censorings,
#' supplied as argument \code{n.cens} is supplied. 
#' Computed originally on the subject
#' level as analogous to \code{pp * as.integer(status == "censored")}.
#' Required when \code{relsurv.method = "pp"}. Flexible input.

#' @param pyrs.pp variable containing total Pohar-Perme weighted subject-times,
#' supplied as argument \code{pyrs} is supplied. 
#' Computed originally on the subject
#' level as analogous to \code{pp * pyrs}.
#' Required when \code{relsurv.method = "pp"}. Flexible input.

#' @param d.exp.pp variable containing total Pohar-Perme weighted counts 
#' of excess events,
#' supplied as argument \code{pyrs} is supplied. 
#' Computed originally on the subject
#' level as analogous to \code{pp * d.exp}.
#' Required when \code{relsurv.method = "pp"}. Flexible input.
#' 
#' 
#' @section Data requirements:
#' 
#' \code{survtab_ag} computes estimates of survival time functions using 
#' pre-aggregated data. For using subject-level data directly, use 
#' \code{\link{survtab}}. For aggregating data, see \code{\link{lexpand}}
#' and \code{\link{aggre}}. 
#' 
#' By default, and if data is an \code{aggre} object (not mandatory), 
#' \code{survtab_ag} makes use of the exact same breaks that were used in 
#' splitting the original data (with e.g. \code{lexpand}), so it is not 
#' necessary to specify any \code{surv.breaks}. If specified, the 
#' \code{surv.breaks} must be a subset of the pertinent 
#' pre-existing breaks. When data is not an \code{aggre} object, breaks
#' must always be specified. Interval lengths (\code{delta} in output) are 
#' also calculated based on whichever breaks are used, 
#' so the upper limit of the breaks should
#' therefore be meaningful and never e.g. \code{Inf}. 
#' 
#' 
#' @examples
#' ## see more examples with explanations in vignette("survtab_examples")
#' 
#' #### survtab_ag usage
#' 
#' data("sire", package = "popEpi")
#' ## prepare data for e.g. 5-year "period analysis" for 2008-2012
#' ## note: sire is a simulated cohort integrated into popEpi.
#' BL <- list(fot=seq(0, 5, by = 1/12),
#'            per = c("2008-01-01", "2013-01-01"))
#' x <- lexpand(sire, birth = bi_date, entry = dg_date, exit = ex_date,
#'              status = status %in% 1:2,
#'              breaks = BL,
#'              pophaz = popmort,
#'              aggre = list(fot))
#'              
#' ## calculate relative EdererII period method
#' ## NOTE: x is an aggre object here, so surv.breaks are deduced
#' ## automatically
#' st <- survtab_ag(fot ~ 1, data = x)
#' 
#' summary(st, t = 1:5) ## annual estimates
#' summary(st, q = list(r.e2 = 0.75)) ## 1st interval where r.e2 < 0.75 at end
#' \dontrun{
#' plot(st)
#' 
#' 
#' ## non-aggre data: first call to survtab_ag would fail
#' df <- data.frame(x)
#' # st <- survtab_ag(fot ~ 1, data = x)
#' st <- survtab_ag(fot ~ 1, data = x, surv.breaks = BL$fot)
#' 
#' ## calculate age-standardised 5-year relative survival ratio using 
#' ## Ederer II method and period approach 
#' 
#' sire$agegr <- cut(sire$dg_age,c(0,45,55,65,75,Inf),right=F)
#' BL <- list(fot=seq(0, 5, by = 1/12),
#'            per = c("2008-01-01", "2013-01-01"))
#' x <- lexpand(sire, birth = bi_date, entry = dg_date, exit = ex_date,
#'              status = status %in% 1:2,
#'              breaks = BL,
#'              pophaz = popmort,
#'              aggre = list(agegr, fot))
#' 
#' ## age standardisation using internal weights (age distribution of 
#' ## patients diagnosed within the period window)
#' ## (NOTE: what is done here is equivalent to using weights = "internal")
#' w <- aggregate(at.risk ~ agegr, data = x[x$fot == 0], FUN = sum)
#' names(w) <- c("agegr", "weights")
#' 
#' st <- survtab_ag(fot ~ adjust(agegr), data = x, weights = w)
#' plot(st, y = "r.e2.as", col = c("blue"))
#' 
#' ## age standardisation using ICSS1 weights
#' data(ICSS)
#' cut <- c(0, 45, 55, 65, 75, Inf)
#' agegr <- cut(ICSS$age, cut, right = FALSE)
#' w <- aggregate(ICSS1~agegr, data = ICSS, FUN = sum)
#' names(w) <- c("agegr", "weights")
#'
#' st <- survtab_ag(fot ~ adjust(agegr), data = x, weights = w)
#' lines(st, y = "r.e2.as", col = c("red"))
#' 
#' 
#' ## cause-specific survival
#' sire$stat <- factor(sire$status, 0:2, c("alive", "canD", "othD"))
#' x <- lexpand(sire, birth = bi_date, entry = dg_date, exit = ex_date,
#'              status = stat,
#'              breaks = BL,
#'              pophaz = popmort,
#'              aggre = list(agegr, fot))
#' st <- survtab_ag(fot ~ adjust(agegr), data = x, weights = w,
#'                  d = c("fromalivetocanD", "fromalivetoothD"),
#'                  surv.type = "surv.cause")
#' plot(st, y = "surv.obs.fromalivetocanD.as")
#' lines(st, y = "surv.obs.fromalivetoothD.as", col = "red")
#' 
#' 
#' }
#' @export
survtab_ag <- function(
  formula = NULL,
  data, 
  adjust = NULL,
  weights = NULL,
  surv.breaks = NULL, 
  n = "at.risk",
  d = "from0to1",
  n.cens = "from0to0",
  pyrs = "pyrs",
  d.exp = "d.exp",
  n.pp = NULL,
  d.pp = "d.pp",
  d.pp.2 = "d.pp.2",
  n.cens.pp = "n.cens.pp",
  pyrs.pp = "pyrs.pp",
  d.exp.pp = "d.exp.pp",
  surv.type = "surv.rel", 
  surv.method = "hazard", 
  relsurv.method = "e2",  
  subset = NULL,
  conf.level = 0.95, 
  conf.type = "log-log",
  verbose = FALSE
) {
  
  parent_frame <- parent.frame(1)
  user_call <- match.call()
  full_call_args <- full_call_args()
  
  # nse ------------------------------------------------------------------------
  nse_arg_nms <- c(
    "formula", "adjust", "subset", "d",
    switch(
      surv.method,
      hazard = "pyrs",
      lifetable = c("n", "n.cens")
    )
  )
  if (surv.type == "surv.rel") {
    nse_arg_nms <- c(
      nse_arg_nms,
      switch(
        relsurv.method,
        e2 = "d.exp",
        pp = c(
          "d.exp.pp", "d.pp", "d.pp.2",
          switch(
            surv.method,
            hazard = c("pyrs.pp"),
            lifetable = c("n.pp", "n.cens.pp")
          )
        )
      )
    )
  }
  
  nse_args <- full_call_args[nse_arg_nms]
  nse_args["subset"] <- list(nse_eval(
    dt = data, 
    j = nse_args[["subset"]], 
    enclos = parent_frame
  ))
  nse_args <- nse_args[nse_arg_nms]
  nonsubset_nse_arg_nms <- setdiff(names(nse_args), "subset")
  nse_args[nonsubset_nse_arg_nms] <- lapply(
    nse_args[nonsubset_nse_arg_nms], 
    function(expr) {
      subset <- nse_args[["subset"]]
      arg_list <- list(dt = data, j = expr, enclos = parent_frame)
      if (!is.null(subset)) {
        arg_list[["i"]] <- subset
      } 
      do.call(nse_eval, arg_list)
    }
  )
  
  nse_args
}




