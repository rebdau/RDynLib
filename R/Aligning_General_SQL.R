#' @title
#'
#' @description
#'
#' The `Aligning_General_SQL()` function ...
#'
#' @param FT_con
#'
#' @param QTOF_con
#'
#' @param FT_expnr
#'
#' @param QTOF_expnr
#'
#' @param err
#'
#' @param t.ini
#'
#' @return
#'
#' @author Ahlam Mentag
#'
#' @export
Aligning_General_SQL <- function(FT_con, QTOF_con, FT_expnr, QTOF_expnr, err,
                                 t.ini) {

    ft_sql  <- dbGetQuery(FT_con, paste0("SELECT * FROM ms_compound where ",
                                         "expid = ", FT_expnr))
    syn_sql <- dbGetQuery(QTOF_con, paste0("SELECT * FROM ms_compound where ",
                                           "expid = ", QTOF_expnr))

    ft_mode  <- dbGetQuery(FT_con, "SELECT mode FROM experiment")$mode[1]
    syn_mode <- dbGetQuery(QTOF_con, "SELECT mode FROM experiment")$mode[1]

    ft.exp <- ft_sql[, c("retention_time", "mass_measured", "compound_id",
                         "expid", "name", "smiles", "nodename")]
    syn.exp <- syn_sql[, c("retention_time", "mass_measured", "compound_id",
                           "expid", "name", "smiles", "nodename")]
    ft.exp.o <- ft.exp[order(ft.exp$retention_time), ]

    reg <- nrow(ft.exp.o)
    if (reg %% 2 != 0) reg <- reg - 1
    COMPID <- ft.exp.o[reg/2, "compound_id"]

    ft.sh <- Make_loc_shortlist(COMPID, reg, ft.exp)

    if (ft_mode != syn_mode) {
        syn.exp$mass_measured <- syn.exp$mass_measured - 2.01456
    }

    ft.sh <- Remov_empties(err, syn.exp, ft.sh, reg)

    LCal <- LCalign(err, t.ini, syn.exp, ft.sh)

    if (ft_mode != syn_mode) {
        LCal[,5] <- as.numeric(LCal[,5]) + 2.01456
    }

    return(LCal)
}
