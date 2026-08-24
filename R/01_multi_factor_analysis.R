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

# ============================================================
# 公共数据库多因素全流程分析（ICU转入）
# 数据：MIMIC-IV + INSPIRE（剔除eICU，6,430例）
# 流程：Table 1 → 单/多因素OR → 多因素ROC → 校准曲线 → DCA
# 模型：icu ~ gnri + nlr + age + sex + bmi + hb + crp + hip_fracture
#       （albumin/中性粒/淋巴 为 GNRI/NLR 组分，不重复纳入；charlson缺失88%剔除）
# ============================================================
suppressPackageStartupMessages({
  library(pROC)
  library(rms)
  library(showtext)
  library(sysfonts)
})
font_add("SimHei", "C:/Windows/Fonts/simhei.ttf")
showtext_auto()

OUT <- "file.path(BASE, "results")"
FIG <- "file.path(BASE, "figures")"

df <- read.csv("file.path(BASE, "data/deriv/cohort_analysis.csv")",
               na.strings = c("", "NA", "nan"))
df$icu <- df$icu_admit
df$sex_m <- ifelse(df$sex == "M", 1, 0)
main <- df[df$db != "eicu", ]

cat("═══════════════════════════════════════════════\n")
cat("公共库多因素分析 | 样本:", nrow(main), "| ICU转入:", sum(main$icu), "\n")
cat("═══════════════════════════════════════════════\n")

# ── 1. Table 1（ICU转入 vs 未转入）───────────────────────
cat("\n── 1. Table 1 ──\n")
vars <- c("age", "sex_m", "bmi", "albumin", "neutrophil_abs", "lymph_abs",
          "hb", "crp", "hip_fracture", "charlson", "gnri", "nlr", "los_days")
t1 <- data.frame()
for (v in vars) {
  c1 <- main[[v]][main$icu == 1]; c0 <- main[[v]][main$icu == 0]
  if (length(unique(main[[v]])) <= 3 & !v %in% c("charlson")) {
    # 分类变量
    r1 <- mean(c1, na.rm = T) * 100; r0 <- mean(c0, na.rm = T) * 100
    pv <- tryCatch(fisher.test(main[[v]], main$icu)$p.value, error = function(e) NA)
    t1 <- rbind(t1, data.frame(变量 = v, 病例组 = sprintf("%.1f%%", r1),
                               对照组 = sprintf("%.1f%%", r0), P = signif(pv, 3)))
  } else {
    if (v == "nlr") {  # 非正态用中位数
      r1 <- sprintf("%.2f [%.2f, %.2f]", median(c1, na.rm=T), quantile(c1, .25, na.rm=T), quantile(c1, .75, na.rm=T))
      r0 <- sprintf("%.2f [%.2f, %.2f]", median(c0, na.rm=T), quantile(c0, .25, na.rm=T), quantile(c0, .75, na.rm=T))
    } else {
      r1 <- sprintf("%.2f±%.2f", mean(c1, na.rm=T), sd(c1, na.rm=T))
      r0 <- sprintf("%.2f±%.2f", mean(c0, na.rm=T), sd(c0, na.rm=T))
    }
    pv <- tryCatch(wilcox.test(c1, c0)$p.value, error = function(e) NA)
    t1 <- rbind(t1, data.frame(变量 = v, 病例组 = r1, 对照组 = r0, P = signif(pv, 3)))
  }
}
print(t1, row.names = FALSE)
write.csv(t1, file.path(OUT, "table1_多因素_ICU转入.csv"), row.names = FALSE,
          fileEncoding = "UTF-8")

# ── 2. 单因素 / 多因素 Logistic ──────────────────────────
cat("\n── 2. 单因素 / 多因素 Logistic 回归 ──\n")
preds <- c("gnri", "nlr", "age", "sex_m", "bmi", "hb", "crp", "hip_fracture")
labels <- c("GNRI(每+1)", "NLR(每+1)", "年龄(每+1岁)", "男性(vs女性)",
            "BMI(每+1)", "血红蛋白(每+1)", "CRP(每+1)", "髋部骨折(vs其他)")

# 用完整病例（多因素分析人群）
comp <- main[complete.cases(main[c("icu", preds)]), ]
cat("多因素分析有效样本:", nrow(comp), "| ICU转入:", sum(comp$icu), "\n")

or_table <- data.frame()
for (i in seq_along(preds)) {
  f <- as.formula(paste("icu ~", preds[i]))
  m <- glm(f, data = comp, family = binomial)
  s <- summary(m); ci <- confint(m)
  or_table <- rbind(or_table, data.frame(
    变量 = labels[i], 分析 = "单因素",
    OR = exp(coef(m)[2]), 低CI = exp(ci[2, 1]), 高CI = exp(ci[2, 2]),
    P = formatC(s$coef[2, 4], format = "e", digits = 2)))
}
# 多因素
m_multi <- glm(icu ~ gnri + nlr + age + sex_m + bmi + hb + crp + hip_fracture,
               data = comp, family = binomial)
s <- summary(m_multi); ci <- confint(m_multi)
for (i in 2:length(coef(m_multi))) {
  nm <- names(coef(m_multi))[i]
  idx <- which(preds == nm)
  or_table <- rbind(or_table, data.frame(
    变量 = labels[idx], 分析 = "多因素",
    OR = exp(coef(m_multi)[i]), 低CI = exp(ci[i, 1]), 高CI = exp(ci[i, 2]),
    P = formatC(s$coef[i, 4], format = "e", digits = 2)))
}
print(or_table, row.names = FALSE)
write.csv(or_table, file.path(OUT, "logistic_多因素_ICU转入.csv"), row.names = FALSE,
          fileEncoding = "UTF-8")

# ── 3. ROC：GNRI / NLR / 多因素模型 ───────────────────────
cat("\n── 3. ROC 对比 ──\n")
roc_g <- roc(comp$icu, comp$gnri, quiet = TRUE)
roc_n <- roc(comp$icu, comp$nlr, quiet = TRUE)
prob_m <- predict(m_multi, type = "response")
roc_m <- roc(comp$icu, prob_m, quiet = TRUE)

cat(sprintf("GNRI      AUC=%.3f (%.3f-%.3f)\n", auc(roc_g), ci.auc(roc_g)[1], ci.auc(roc_g)[3]))
cat(sprintf("NLR       AUC=%.3f (%.3f-%.3f)\n", auc(roc_n), ci.auc(roc_n)[1], ci.auc(roc_n)[3]))
cat(sprintf("多因素模型 AUC=%.3f (%.3f-%.3f)\n", auc(roc_m), ci.auc(roc_m)[1], ci.auc(roc_m)[3]))
# DeLong 检验
cat("DeLong: GNRI vs 多因素 P =", formatC(roc.test(roc_g, roc_m)$p.value, digits = 3), "\n")
cat("DeLong: NLR  vs 多因素 P =", formatC(roc.test(roc_n, roc_m)$p.value, digits = 3), "\n")

# ── 4. 校准曲线（多因素模型）──────────────────────────────
cat("\n── 4. 校准曲线 ──\n")
val_res <- val.prob(prob_m, comp$icu, pl = FALSE)
cat("截距:", round(val_res["Intercept"], 3), "| 斜率:", round(val_res["Slope"], 3), "\n")

# ── 5. DCA（手写）────────────────────────────────────────
cat("\n── 5. 决策曲线分析（DCA）──\n")
dca_calc <- function(y, p) {
  n <- length(y)
  pts <- seq(0.01, 0.99, by = 0.01)
  nb <- sapply(pts, function(pt) {
    tp <- sum(p >= pt & y == 1); fp <- sum(p >= pt & y == 0)
    tp/n - fp/n * pt/(1 - pt)
  })
  nb_all <- mean(y) - (1 - mean(y)) * pts/(1 - pts)
  data.frame(pt = pts, model = nb, all = nb_all, none = 0)
}
dca_m <- dca_calc(comp$icu, prob_m)
dca_g <- dca_calc(comp$icu, predict(glm(icu ~ gnri, data = comp, family = binomial), type = "response"))
dca_n <- dca_calc(comp$icu, predict(glm(icu ~ nlr, data = comp, family = binomial), type = "response"))
cat("DCA 净收益（阈值0.1/0.2/0.3）:\n")
for (pt in c(0.1, 0.2, 0.3)) {
  cat(sprintf("  pt=%.1f: 多因素%.4f | GNRI %.4f | NLR %.4f | 全转 %.4f | 不转 0\n",
              pt,
              dca_m$model[which.min(abs(dca_m$pt - pt))],
              dca_g$model[which.min(abs(dca_g$pt - pt))],
              dca_n$model[which.min(abs(dca_n$pt - pt))],
              dca_m$all[which.min(abs(dca_m$pt - pt))]))
}

# ═══════════════════════════════════════════════
# 图 A：多因素 ROC（含GNRI/NLR对比）
# ═══════════════════════════════════════════════
pdf(file.path(FIG, "Fig1_ROC_多因素模型.pdf"), width = 7, height = 6.5)
par(mar = c(5, 5, 3, 2), family = "SimHei")
plot(roc_g, col = "#2F5597", lwd = 2.5, main = "预测术后ICU转入：多因素模型 vs 单指标",
     xlab = "1 - 特异性", ylab = "敏感性", cex.main = 1.15)
plot(roc_n, col = "#C00000", lwd = 2.5, add = TRUE)
plot(roc_m, col = "#548235", lwd = 3, lty = 2, add = TRUE)
legend("bottomright",
       legend = c(sprintf("GNRI (AUC=%.3f)", auc(roc_g)),
                  sprintf("NLR (AUC=%.3f)", auc(roc_n)),
                  sprintf("多因素模型 (AUC=%.3f)", auc(roc_m))),
       col = c("#2F5597", "#C00000", "#548235"), lwd = c(2.5, 2.5, 3),
       lty = c(1, 1, 2), bty = "n", cex = 0.95)
dev.off()

# ═══════════════════════════════════════════════
# 图 B：校准曲线（手动十分位分箱，更可控）
# ═══════════════════════════════════════════════
pdf(file.path(FIG, "Fig9_校准曲线_多因素.pdf"), width = 6.5, height = 6.5)
par(mar = c(5, 5, 3, 2), family = "SimHei")
# 十分位分箱
q10 <- quantile(prob_m, probs = seq(0, 1, 0.1))
q10[1] <- 0; q10[11] <- 1
bin_idx <- cut(prob_m, breaks = unique(q10), include.lowest = TRUE)
bin_pred <- tapply(prob_m, bin_idx, mean)
bin_obs <- tapply(comp$icu, bin_idx, mean)
plot(bin_pred, bin_obs, xlim = c(0, 0.5), ylim = c(0, 0.5),
     pch = 19, cex = 1.3, col = "#2F5597",
     xlab = "预测概率（十分位）", ylab = "实际ICU转入比例",
     main = "多因素模型校准曲线（预测术后ICU转入）",
     cex.main = 1.1, cex.lab = 1.0)
abline(0, 1, lty = 2, col = "gray60", lwd = 1.5)
lines(bin_pred, bin_obs, col = "#2F5597", lwd = 1.5)
legend("topleft", c("理想校准线", "模型校准"), col = c("gray60", "#2F5597"),
       lty = c(2, 1), lwd = 1.5, bty = "n", cex = 0.9)
text(0.42, 0.04, sprintf("截距=%.3f\n斜率=%.3f", val_res["Intercept"], val_res["Slope"]),
     cex = 0.85, adj = 1)
dev.off()

# ═══════════════════════════════════════════════
# 图 C：DCA
# ═══════════════════════════════════════════════
pdf(file.path(FIG, "Fig10_DCA_多因素.pdf"), width = 7, height = 6)
par(mar = c(5, 5, 3, 2), family = "SimHei")
plot(dca_m$pt, dca_m$model, type = "l", lwd = 2.5, col = "#548235",
     xlim = c(0, 1), ylim = c(-0.05, max(dca_m$model, dca_m$all) * 1.05),
     xlab = "阈值概率", ylab = "净收益", main = "决策曲线分析（DCA）",
     cex.main = 1.1)
lines(dca_g$pt, dca_g$model, lwd = 2, col = "#2F5597")
lines(dca_n$pt, dca_n$model, lwd = 2, col = "#C00000")
lines(dca_m$pt, dca_m$all, lwd = 1.5, lty = 2, col = "gray50")
abline(h = 0, lty = 2, col = "gray80")
legend("topright",
       legend = c("多因素模型", "GNRI", "NLR", "全部转入", "不转入"),
       col = c("#548235", "#2F5597", "#C00000", "gray50", "gray80"),
       lwd = c(2.5, 2, 2, 1.5, 1), lty = c(1, 1, 1, 2, 2), bty = "n", cex = 0.9)
dev.off()

cat("\n✅ 公共库多因素全流程完成\n")
cat("产出: table1_多因素_ICU转入.csv / logistic_多因素_ICU转入.csv\n")
cat("图表: Fig1_ROC_多因素模型 / Fig9_校准曲线 / Fig10_DCA\n")
