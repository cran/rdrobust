## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  fig.width = 7,
  fig.height = 4.5
)

## ----data---------------------------------------------------------------------
library(rdrobust)
data(rdrobust_RDsenate)
df <- rdrobust_RDsenate
head(df)

## ----data-arg-basic-----------------------------------------------------------
r <- rdrobust(y = vote, x = margin, data = df)
round(r$coef, 3)

## ----data-arg-expr------------------------------------------------------------
r_cov <- rdrobust(y = vote, x = margin,
                  covs = cbind(class, termshouse, termssenate),
                  data = df)
round(r_cov$coef, 3)

# subset via an expression evaluated inside `data`
r_sub <- rdrobust(y = vote, x = margin,
                  subset = year >= 1950, data = df)
round(r_sub$coef, 3)

## ----covs-formula-factor------------------------------------------------------
df$classf <- factor(df$class, labels = c("A", "B", "C"))

r_fm <- rdrobust(y = vote, x = margin,
                 covs = ~ classf + termshouse + I(termssenate^2),
                 data = df)
round(r_fm$coef, 3)

## ----covs-char----------------------------------------------------------------
covs_names <- c("class", "termshouse", "termssenate")
r_ch <- rdrobust(y = vote, x = margin, covs = covs_names, data = df)
round(r_ch$coef, 3)

## ----covs-matrix--------------------------------------------------------------
Z <- with(df, cbind(class, termshouse, termssenate))
r_mat <- rdrobust(y = vote, x = margin, covs = Z, data = df)
round(r_mat$coef, 3)

## ----rdplot, fig.cap="RD plot: U.S. Senate incumbency advantage."-------------
rdplot(y = vote, x = margin, data = df,
       title   = "Senate: incumbency advantage",
       x.label = "Vote margin (t)",
       y.label = "Vote share (t+2)")

## ----hc-variants--------------------------------------------------------------
hc_fits <- lapply(c("nn", "hc0", "hc1", "hc2", "hc3"), function(v) {
  rdrobust(y = vote, x = margin, vce = v, data = df)
})
names(hc_fits) <- c("NN", "HC0", "HC1", "HC2", "HC3")

hc_tbl <- do.call(rbind, lapply(names(hc_fits), function(nm) {
  r <- hc_fits[[nm]]
  data.frame(vce    = nm,
             coef   = round(r$coef[2], 3),   # bias-corrected
             se_rb  = round(r$se[3],   3),   # robust SE
             pval   = round(r$pv[3],   3))
}))
print(hc_tbl, row.names = FALSE)

## ----cluster-counts-----------------------------------------------------------
cat("Number of clusters:", length(unique(df$state)), "\n")

## ----crv-variants-------------------------------------------------------------
r_cr1 <- rdrobust(y = vote, x = margin, cluster = state, vce = "cr1", data = df)
r_cr2 <- rdrobust(y = vote, x = margin, cluster = state, vce = "cr2", data = df)
r_cr3 <- rdrobust(y = vote, x = margin, cluster = state, vce = "cr3", data = df)

cl_tbl <- data.frame(
  vce    = c("CR1", "CR2", "CR3"),
  coef   = c(r_cr1$coef[2],  r_cr2$coef[2],  r_cr3$coef[2]),
  se_rb  = c(r_cr1$se[3],    r_cr2$se[3],    r_cr3$se[3]),
  pval   = c(r_cr1$pv[3],    r_cr2$pv[3],    r_cr3$pv[3])
)
cl_tbl[, 2:4] <- round(cl_tbl[, 2:4], 3)
print(cl_tbl, row.names = FALSE)

## ----cr3-summary--------------------------------------------------------------
summary(r_cr3)

## ----plot-basic, fig.cap="Binscatter with local polynomial fits and 95% CI bands."----
plot(r_cr1, df$vote, df$margin,
     title   = "Senate: incumbency advantage",
     x.label = "Vote margin (t)",
     y.label = "Vote share (t+2)")

## ----plot-effect, fig.height=6, fig.cap="Main plot with RD effect panel."-----
plot(r_cr1, df$vote, df$margin,
     show_effect = TRUE,
     title       = "Senate: incumbency advantage",
     x.label     = "Vote margin (t)",
     y.label     = "Vote share (t+2)")

## ----plot-cr3, fig.cap="CR3 clustered variance, with effect panel."-----------
plot(r_cr3, df$vote, df$margin,
     show_effect = TRUE,
     title       = "Senate: CR3 clustered SE",
     x.label     = "Vote margin (t)",
     y.label     = "Vote share (t+2)")

## ----rdplot-hide--------------------------------------------------------------
out <- rdplot(y = vote, x = margin, c = 0,
              hide = TRUE, ci = 95, data = df)
str(out$vars_bins, give.attr = FALSE)
str(out$vars_poly, give.attr = FALSE)
out$coef

## ----rdplot-manual, fig.cap="Custom RD plot built from rdplot() output."------
library(ggplot2)
poly_l <- subset(out$vars_poly, rdplot_x <  0)
poly_r <- subset(out$vars_poly, rdplot_x >= 0)

ggplot() +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_errorbar(data = out$vars_bins,
                aes(x = rdplot_mean_bin,
                    ymin = rdplot_ci_l, ymax = rdplot_ci_r),
                color = "grey60", width = 0) +
  geom_point(data = out$vars_bins,
             aes(x = rdplot_mean_bin, y = rdplot_mean_y),
             color = "grey40", size = 1.5) +
  geom_line(data = poly_l, aes(x = rdplot_x, y = rdplot_y)) +
  geom_line(data = poly_r, aes(x = rdplot_x, y = rdplot_y)) +
  theme_bw() +
  labs(x = "Vote margin (t)", y = "Vote share (t+2)",
       title = "Senate: regression function fit")

## ----rdplot-from-coef---------------------------------------------------------
eval_fit <- function(xs, coef, c)
  sapply(xs, function(xi) sum(coef * (xi - c) ^ (0:(length(coef) - 1))))

xs <- seq(-50, 50, length.out = 5)
data.frame(
  x       = xs,
  y_left  = ifelse(xs <  0, eval_fit(xs, out$coef[, "Left"],  0), NA),
  y_right = ifelse(xs >= 0, eval_fit(xs, out$coef[, "Right"], 0), NA)
)

## ----broom-tidy, eval=requireNamespace("broom", quietly=TRUE)-----------------
library(broom)
tidy(r_cr3)

## ----broom-glance, eval=requireNamespace("broom", quietly=TRUE)---------------
glance(r_cr3)

## ----broom-collect, eval=requireNamespace("broom", quietly=TRUE)--------------
specs <- list(CR1 = r_cr1, CR2 = r_cr2, CR3 = r_cr3)

# Robust BC estimates side-by-side
res <- do.call(rbind, lapply(names(specs), function(nm) {
  cbind(spec = nm, tidy(specs[[nm]])[3, ])  # row 3 = robust
}))
res[, c("spec", "estimate", "std.error", "conf.low", "conf.high")]

## ----broom-glance-compare, eval=requireNamespace("broom", quietly=TRUE)-------
gl <- do.call(rbind, lapply(names(specs), function(nm) {
  cbind(spec = nm, glance(specs[[nm]]))
}))
gl[, c("spec", "nobs.effective.left", "nobs.effective.right",
       "h.left", "vce", "clusters.left", "clusters.right")]

## ----broom-hc-vs-cl, eval=requireNamespace("broom", quietly=TRUE)-------------
r_hc3 <- rdrobust(y = vote, x = margin, vce = "hc3", data = df)

all_specs <- list("HC3"  = r_hc3,
                  "CR1"  = r_cr1,
                  "CR2"  = r_cr2,
                  "CR3"  = r_cr3)

comparison <- do.call(rbind, lapply(names(all_specs), function(nm) {
  g <- glance(all_specs[[nm]])
  t <- tidy(all_specs[[nm]])[3, ]   # robust row
  data.frame(spec      = nm,
             estimate  = round(t$estimate, 3),
             std.error = round(t$std.error, 3),
             ci        = paste0("[", round(t$conf.low, 3), ", ", round(t$conf.high, 3), "]"),
             h         = round(g$h.left, 3),
             clusters  = g$clusters.left)
}))
print(comparison, row.names = FALSE)

## ----table-fits---------------------------------------------------------------
covs_form <- ~ class + termshouse + termssenate

m1 <- rdrobust(y = vote, x = margin, data = df)
m2 <- rdrobust(y = vote, x = margin, covs = covs_form, data = df)
m3 <- rdrobust(y = vote, x = margin, vce = "cr1", cluster = state, data = df)
m4 <- rdrobust(y = vote, x = margin, vce = "cr1", cluster = state,
               covs = covs_form, data = df)

## ----table-row-builder--------------------------------------------------------
rdr_row <- function(fit, label) {
  data.frame(
    spec      = label,
    estimate  = fit$coef[1, 1],                                 # tau_conv
    rbc_ci    = sprintf("[%.3f, %.3f]", fit$ci[3, 1], fit$ci[3, 2]),
    rbc_pval  = fit$pv[3, 1],
    h         = fit$bws[1, 1],                                  # h_left
    N_eff     = fit$N_h[1] + fit$N_h[2],
    clusters  = if (!is.null(fit$n_clust)) sum(fit$n_clust) else NA_integer_
  )
}

tab <- rbind(
  rdr_row(m1, "sharp"),
  rdr_row(m2, "sharp + covs"),
  rdr_row(m3, "cluster"),
  rdr_row(m4, "cluster + covs")
)

## ----table-kable--------------------------------------------------------------
knitr::kable(
  tab,
  digits   = 3,
  col.names = c("Specification", "Estimate", "Robust 95% CI",
                "RBC p-value", "h", "Eff. N", "Clusters"),
  caption  = "RD estimates: conventional point estimate with RBC confidence interval."
)

## ----table-paper--------------------------------------------------------------
paper <- t(tab[, c("estimate", "rbc_ci", "h", "N_eff")])
colnames(paper) <- tab$spec
rownames(paper) <- c("RD effect (conv.)", "Robust 95% CI", "h", "Eff. N")

knitr::kable(
  paper,
  caption = "Senate incumbency: RD estimates across specifications."
)

## ----table-latex, eval=FALSE--------------------------------------------------
# tex <- knitr::kable(
#   paper,
#   format    = "latex",
#   booktabs  = TRUE,
#   caption   = "RD estimates for Senate incumbency.",
#   label     = "tab:rdrobust"
# )
# writeLines(tex, "rdrobust_table.tex")

## ----session------------------------------------------------------------------
sessionInfo()

