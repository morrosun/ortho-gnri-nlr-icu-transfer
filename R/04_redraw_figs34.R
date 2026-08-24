# ============================================================
# Reproducibility header (added for public release)
#   Set BASE to the repository root before running, e.g.:
#     BASE <- "D:/path/to/ortho-gnri-nlr-icu-transfer"
#   or use the 'here' package:  BASE <- here::here()
#   Input data = data/deriv/cohort_analysis.csv (de-identified, derived
#   from MIMIC-IV / INSPIRE / eICU under PhysioNet credentialed access).
#   Local 395-case validation data is NOT included; supply your own
#   CSV and point LOCAL_DATA at it.
# ============================================================
BASE <- "."
LOCAL_DATA <- file.path(BASE, "data/local_validation.csv")  # user-provided, not in repo

# -*- coding: utf-8 -*-
# ============================================================
# Figure 3 (ZH, v7): 模型性能三合一（共同模型 5 变量）
#   ROC / 校准 / DCA — 开发拟合 → 固定系数外部验证
# Figure 4 (ZH, v7): 验证集分层剂量反应
# ============================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(patchwork)
  library(pROC); library(showtext); library(sysfonts)
})
font_add("SimHei", "C:/Windows/Fonts/simhei.ttf")
showtext_auto(); showtext_opts(dpi = 300)
OUT <- "BASE/投稿图表集"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

col_dev <- "#3182BD"; col_val <- "#D24B40"; grey_ref <- "#8C8C8C"

pub <- read.csv("file.path(BASE, "data/deriv/cohort_analysis.csv")",
                na.strings = c("", "NA", "nan"))
pub <- pub[pub$db != "eicu", ]
pub$icu <- as.integer(pub$icu_admit == 1); pub$sex_m <- as.integer(pub$sex == "M")
vars5 <- c("gnri", "nlr", "age", "sex_m", "bmi")
comp_pub <- pub[complete.cases(pub[, vars5]), ]

loc <- read.csv("LOCAL_DATA  # not included in repo; provide your local validation data",
                na.strings = c("", "NA", "nan"))
loc$case_icu <- as.integer(loc$icu_transfer == "是"); loc$sex_m <- as.integer(loc$gender == "男")
main_loc <- loc[loc$icu_transfer == "是" | loc$group == "对照组", ]
comp_loc <- main_loc[complete.cases(main_loc[, vars5]), ]

m_dev5 <- glm(icu ~ gnri + nlr + age + sex_m + bmi, data = comp_pub, family = binomial)
p_dev <- predict(m_dev5, type = "response")
X_val <- model.matrix(~ gnri + nlr + age + sex_m + bmi, data = comp_loc)
p_val_fixed <- plogis(as.numeric(X_val %*% coef(m_dev5)))

roc_dev_m <- roc(comp_pub$icu, p_dev, quiet = TRUE)
roc_val_m <- roc(comp_loc$case_icu, p_val_fixed, quiet = TRUE)
roc_dev_g <- roc(comp_pub$icu, comp_pub$gnri, quiet = TRUE)
roc_dev_n <- roc(comp_pub$icu, comp_pub$nlr, quiet = TRUE)

roc_to_df <- function(r, series) {
  data.frame(fpr = 1 - r$specificities, tpr = r$sensitivities, series = series)
}
roc_df <- bind_rows(
  roc_to_df(roc_dev_m, "开发 共同模型 (AUC 0.712)"),
  roc_to_df(roc_dev_g, "开发 GNRI (AUC 0.651)"),
  roc_to_df(roc_dev_n, "开发 NLR (AUC 0.675)"),
  roc_to_df(roc_val_m, "验证 固定系数 (AUC 0.726)"))
roc_df$series <- factor(roc_df$series, levels = c(
  "开发 共同模型 (AUC 0.712)", "开发 GNRI (AUC 0.651)",
  "开发 NLR (AUC 0.675)", "验证 固定系数 (AUC 0.726)"))

calib <- function(y, p) {
  ord <- order(p); y <- y[ord]; p <- p[ord]
  k <- 10; n <- length(y)
  idx <- ceiling(seq(1, n, length.out = k + 1))
  obs <- pred <- rep(NA, k)
  for (i in 1:k) {
    sel <- idx[i]:(idx[i + 1] - 1); if (i == k) sel <- idx[i]:n
    if (length(sel) > 0) { pred[i] <- mean(p[sel]); obs[i] <- mean(y[sel]) }
  }
  data.frame(pred = pred, obs = obs)
}
cal_df <- bind_rows(
  transform(calib(comp_pub$icu, p_dev), cohort = "开发队列"),
  transform(calib(comp_loc$case_icu, p_val_fixed), cohort = "验证队列"))

dca <- function(y, p, pt_seq) {
  nb_model <- nb_all <- nb_none <- numeric(length(pt_seq))
  prev <- mean(y)
  for (i in seq_along(pt_seq)) {
    pt <- pt_seq[i]
    tp <- sum(p >= pt & y == 1); fp <- sum(p >= pt & y == 0)
    n <- length(y)
    nb_model[i] <- tp / n - fp / n * pt / (1 - pt)
    nb_all[i] <- prev - (1 - prev) * pt / (1 - pt)
    nb_none[i] <- 0
  }
  data.frame(pt = pt_seq, model = nb_model, all = nb_all, none = nb_none)
}
pt_seq <- seq(0.01, 0.6, by = 0.01)
dca_df <- bind_rows(
  transform(dca(comp_pub$icu, p_dev, pt_seq), cohort = "开发 共同模型"),
  transform(dca(comp_loc$case_icu, p_val_fixed, pt_seq), cohort = "验证 固定系数"))
ref_df <- data.frame(pt = pt_seq, all = dca(comp_pub$icu, p_dev, pt_seq)$all)

theme_nf <- function(base = 8) {
  theme_classic(base_size = base, base_family = "SimHei") +
    theme(axis.line = element_line(linewidth = 0.35, colour = "black"),
          axis.ticks = element_line(linewidth = 0.35, colour = "black"),
          axis.title = element_text(size = base),
          axis.text = element_text(size = base - 0.5, colour = "black"),
          panel.grid = element_blank(),
          plot.margin = margin(2, 2, 2, 2, unit = "mm"))
}

p_roc <- ggplot(roc_df, aes(fpr, tpr, colour = series, linetype = series)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = grey_ref, linewidth = 0.4) +
  geom_line(linewidth = 0.85) +
  scale_colour_manual(values = c(col_dev, col_dev, col_dev, col_val)) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "solid")) +
  labs(x = "1 - 特异度", y = "敏感度") +
  theme_nf() +
  theme(legend.position = c(0.99, 0.02), legend.justification = c(1, 0),
        legend.title = element_blank(), legend.text = element_text(size = 5.5),
        legend.key.size = unit(0.30, "cm"), legend.spacing.y = unit(0.02, "cm"),
        legend.background = element_blank())

p_cal <- ggplot(cal_df, aes(pred, obs, colour = cohort)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = grey_ref, linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = cohort), size = 1.6, fill = "white", stroke = 0.6) +
  scale_colour_manual(values = c(col_dev, col_val)) +
  scale_shape_manual(values = c(16, 17)) +
  coord_equal(xlim = c(0, 0.6), ylim = c(0, 0.6)) +
  labs(x = "预测概率", y = "实际发生率") +
  theme_nf() +
  theme(legend.position = c(0.02, 0.99), legend.justification = c(0, 1),
        legend.title = element_blank(), legend.text = element_text(size = 6),
        legend.key.size = unit(0.30, "cm"), legend.background = element_blank())

p_dca <- ggplot(dca_df, aes(pt, model, colour = cohort)) +
  geom_hline(yintercept = 0, colour = grey_ref, linewidth = 0.35) +
  geom_line(data = ref_df, aes(pt, all), linetype = "dashed", colour = grey_ref, linewidth = 0.5) +
  geom_line(linewidth = 0.85) +
  scale_colour_manual(values = c(col_dev, col_val)) +
  coord_cartesian(xlim = c(0, 0.6), ylim = c(-0.05, 0.35)) +
  labs(x = "阈值概率", y = "净获益") +
  theme_nf() +
  theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1),
        legend.title = element_blank(), legend.text = element_text(size = 6),
        legend.key.size = unit(0.30, "cm"), legend.background = element_blank())

fig3 <- p_roc + p_cal + p_dca + plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

W3 <- 190; H3 <- 62
pdf(file.path(OUT, "Figure2_模型性能三合一.pdf"), width = W3 / 25.4, height = H3 / 25.4)
print(fig3)
dev.off()
png(file.path(OUT, "Figure2_模型性能三合一.png"), width = W3, height = H3, units = "mm", res = 300)
print(fig3)
dev.off()
cat("Figure 3 (ZH v7) done\n")

# ── Figure 4 (ZH): 分层剂量反应 ──
g_grp <- cut(comp_loc$gnri, breaks = c(-Inf, 82, 92, 98, Inf), right = FALSE,
             labels = c("严重 (<82)", "中度 (82-91)", "轻度 (92-97)", "正常 (≥98)"))
n_grp <- ifelse(comp_loc$nlr >= 2.96, "高NLR", "低NLR")
g_rate <- tapply(comp_loc$case_icu, g_grp, mean) * 100
n_rate <- tapply(comp_loc$case_icu, n_grp, mean) * 100

g_df <- data.frame(strata = names(g_rate), rate = as.numeric(g_rate))
g_df$strata <- factor(g_df$strata,
  levels = c("严重 (<82)", "中度 (82-91)", "轻度 (92-97)", "正常 (≥98)"))
n_df <- data.frame(strata = names(n_rate), rate = as.numeric(n_rate))
n_df$strata <- factor(n_df$strata, levels = c("高NLR", "低NLR"))

p_g <- ggplot(g_df, aes(strata, rate, fill = strata)) +
  geom_col(width = 0.62, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%", rate)), vjust = -0.6, size = 2.4) +
  scale_fill_manual(values = c("#D24B40", "#E28E2C", "#F2C14E", "#66A85E")) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.06))) +
  labs(x = NULL, y = "ICU 转入率 (%)") +
  theme_nf() + theme(axis.text.x = element_text(size = 7))

p_n <- ggplot(n_df, aes(strata, rate, fill = strata)) +
  geom_col(width = 0.5, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%", rate)), vjust = -0.6, size = 2.4) +
  scale_fill_manual(values = c("#D24B40", "#3182BD")) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.06))) +
  labs(x = NULL, y = "ICU 转入率 (%)") +
  theme_nf() + theme(axis.text.x = element_text(size = 7))

fig4 <- p_g + p_n + plot_layout(ncol = 2, widths = c(1.5, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

W4 <- 150; H4 <- 68
pdf(file.path(OUT, "Figure4_分层剂量反应_验证集.pdf"), width = W4 / 25.4, height = H4 / 25.4)
print(fig4)
dev.off()
png(file.path(OUT, "Figure4_分层剂量反应_验证集.png"), width = W4, height = H4, units = "mm", res = 300)
print(fig4)
dev.off()
cat("Figure 4 (ZH) done\n")
