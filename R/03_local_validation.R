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
# 病例对照研究统计演练 — 多因素版（与公共库全流程一致）
# ------------------------------------------------------------
# 主分析   ：术后ICU转入（77例） vs 对照组（313例）
#   流程   ：Table 1(SMD) → 单/多因素OR → ROC → 校准曲线 → DCA
# 敏感性分析：院内死亡（15例） vs 对照（仅单因素OR）
# 模型     ：icu ~ gnri + nlr + age + gender + bmi + charlson + hip_fracture
# 数据     ：模拟数据_病例对照_400例.csv
# ============================================================
suppressPackageStartupMessages({
  library(pROC)
  library(rms)
  library(showtext)
  library(sysfonts)
})
font_add("SimHei", "C:/Windows/Fonts/simhei.ttf")
showtext_auto()

OUT <- "BASE"
df <- read.csv(file.path(OUT, "模拟数据_病例对照_400例.csv"), fileEncoding = "UTF-8")

# ── 变量准备 ──────────────────────────────────────────────
df$case_icu <- ifelse(df$icu_transfer == "是", 1, 0)
df$case_death <- ifelse(df$death == "是", 1, 0)
df$gender_m <- ifelse(df$gender == "男", 1, 0)
# 派生 hip_fracture（演练用：骨折内固定/髋关节置换近似髋部手术）
df$hip_fracture <- ifelse(df$surgery_type %in% c("骨折内固定", "髋关节置换"), 1, 0)
df$asa_high <- ifelse(df$asa %in% c("Ⅲ", "Ⅳ"), 1, 0)

main <- df[df$case_icu == 1 | df$group == "对照组", ]
sens <- df[df$case_death == 1 | df$group == "对照组", ]

cat("═══════════════════════════════════════════════\n")
cat("多因素演练 | 主分析:", nrow(main), "例（ICU转入", sum(main$case_icu),
    "+ 对照", sum(main$group == "对照组"), "）\n")
cat("═══════════════════════════════════════════════\n")

# ── 0. Table 1（病例 vs 对照，SMD 检验匹配质量）────────────
cat("\n── 0. Table 1（SMD 检验匹配质量）──\n")
smd_num <- function(x, g) {
  m1 <- mean(x[g == 1], na.rm = T); m0 <- mean(x[g == 0], na.rm = T)
  s1 <- sd(x[g == 1], na.rm = T); s0 <- sd(x[g == 0], na.rm = T)
  abs(m1 - m0) / sqrt((s1^2 + s0^2) / 2)
}
t1_rows <- list(
  c("年龄(岁)", sprintf("%.1f±%.1f", mean(main$age[main$case_icu==1]), sd(main$age[main$case_icu==1])),
    sprintf("%.1f±%.1f", mean(main$age[main$case_icu==0]), sd(main$age[main$case_icu==0])),
    sprintf("%.3f", smd_num(main$age, main$case_icu))),
  c("男性(%)", sprintf("%.1f", mean(main$gender_m[main$case_icu==1])*100),
    sprintf("%.1f", mean(main$gender_m[main$case_icu==0])*100),
    sprintf("%.3f", abs(mean(main$gender_m[main$case_icu==1]) - mean(main$gender_m[main$case_icu==0])))),
  c("BMI(kg/m²)", sprintf("%.1f±%.1f", mean(main$bmi[main$case_icu==1]), sd(main$bmi[main$case_icu==1])),
    sprintf("%.1f±%.1f", mean(main$bmi[main$case_icu==0]), sd(main$bmi[main$case_icu==0])),
    sprintf("%.3f", smd_num(main$bmi, main$case_icu))),
  c("白蛋白(g/dL)", sprintf("%.2f±%.2f", mean(main$albumin[main$case_icu==1]), sd(main$albumin[main$case_icu==1])),
    sprintf("%.2f±%.2f", mean(main$albumin[main$case_icu==0]), sd(main$albumin[main$case_icu==0])),
    sprintf("%.3f", smd_num(main$albumin, main$case_icu))),
  c("GNRI", sprintf("%.1f±%.1f", mean(main$gnri[main$case_icu==1]), sd(main$gnri[main$case_icu==1])),
    sprintf("%.1f±%.1f", mean(main$gnri[main$case_icu==0]), sd(main$gnri[main$case_icu==0])),
    sprintf("%.3f", smd_num(main$gnri, main$case_icu))),
  c("NLR", sprintf("%.2f±%.2f", mean(main$nlr[main$case_icu==1]), sd(main$nlr[main$case_icu==1])),
    sprintf("%.2f±%.2f", mean(main$nlr[main$case_icu==0]), sd(main$nlr[main$case_icu==0])),
    sprintf("%.3f", smd_num(main$nlr, main$case_icu))),
  c("Charlson", sprintf("%.1f±%.1f", mean(main$charlson[main$case_icu==1]), sd(main$charlson[main$case_icu==1])),
    sprintf("%.1f±%.1f", mean(main$charlson[main$case_icu==0]), sd(main$charlson[main$case_icu==0])),
    sprintf("%.3f", smd_num(main$charlson, main$case_icu))),
  c("ASA≥Ⅲ(%)", sprintf("%.1f", mean(main$asa_high[main$case_icu==1])*100),
    sprintf("%.1f", mean(main$asa_high[main$case_icu==0])*100),
    sprintf("%.3f", abs(mean(main$asa_high[main$case_icu==1]) - mean(main$asa_high[main$case_icu==0])))),
  c("手术时间(分)", sprintf("%.0f±%.0f", mean(main$op_time[main$case_icu==1]), sd(main$op_time[main$case_icu==1])),
    sprintf("%.0f±%.0f", mean(main$op_time[main$case_icu==0]), sd(main$op_time[main$case_icu==0])),
    sprintf("%.3f", smd_num(main$op_time, main$case_icu))),
  c("髋部手术(%)", sprintf("%.1f", mean(main$hip_fracture[main$case_icu==1])*100),
    sprintf("%.1f", mean(main$hip_fracture[main$case_icu==0])*100),
    sprintf("%.3f", abs(mean(main$hip_fracture[main$case_icu==1]) - mean(main$hip_fracture[main$case_icu==0]))))
)
t1 <- do.call(rbind, lapply(t1_rows, function(r) data.frame(变量 = r[1], 病例组 = r[2], 对照组 = r[3], SMD = r[4])))
print(t1, row.names = FALSE)
write.csv(t1, file.path(OUT, "结果_Table1_演练.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat("（匹配变量 年龄/性别 SMD<0.1 为合格；暴露变量 SMD 大是预期，不代表匹配失败）\n")

# ── 1. 单因素 / 多因素 OR ─────────────────────────────────
preds <- c("gnri", "nlr", "age", "gender_m", "bmi", "charlson", "hip_fracture")
labels <- c("GNRI(每+1)", "NLR(每+1)", "年龄(每+1岁)", "男性(vs女性)",
            "BMI(每+1)", "Charlson(每+1)", "髋部手术(vs其他)")

or_table <- data.frame()
for (i in seq_along(preds)) {
  f <- as.formula(paste("case_icu ~", preds[i]))
  m <- glm(f, data = main, family = binomial)
  s <- summary(m); ci <- confint(m)
  or_table <- rbind(or_table, data.frame(
    变量 = labels[i], 分析 = "单因素",
    OR = exp(coef(m)[2]), 低CI = exp(ci[2, 1]), 高CI = exp(ci[2, 2]),
    P = formatC(s$coef[2, 4], format = "e", digits = 2)))
}
m_multi <- glm(case_icu ~ gnri + nlr + age + gender_m + bmi + charlson + hip_fracture,
               data = main, family = binomial)
s <- summary(m_multi); ci <- confint(m_multi)
for (i in 2:length(coef(m_multi))) {
  nm <- names(coef(m_multi))[i]
  idx <- which(preds == nm)
  or_table <- rbind(or_table, data.frame(
    变量 = labels[idx], 分析 = "多因素",
    OR = exp(coef(m_multi)[i]), 低CI = exp(ci[i, 1]), 高CI = exp(ci[i, 2]),
    P = formatC(s$coef[i, 4], format = "e", digits = 2)))
}
cat("\n── 单因素 / 多因素 OR ──\n")
print(or_table, row.names = FALSE, digits = 3)
write.csv(or_table, file.path(OUT, "结果_多因素OR_演练.csv"), row.names = FALSE,
          fileEncoding = "UTF-8")

# ── 2. ROC ────────────────────────────────────────────────
roc_g <- roc(main$case_icu, main$gnri, quiet = TRUE)
roc_n <- roc(main$case_icu, main$nlr, quiet = TRUE)
prob_m <- predict(m_multi, type = "response")
roc_m <- roc(main$case_icu, prob_m, quiet = TRUE)

cat("\n── ROC 对比 ──\n")
cat(sprintf("GNRI      AUC=%.3f (%.3f-%.3f)\n", auc(roc_g), ci.auc(roc_g)[1], ci.auc(roc_g)[3]))
cat(sprintf("NLR       AUC=%.3f (%.3f-%.3f)\n", auc(roc_n), ci.auc(roc_n)[1], ci.auc(roc_n)[3]))
cat(sprintf("多因素模型 AUC=%.3f (%.3f-%.3f)\n", auc(roc_m), ci.auc(roc_m)[1], ci.auc(roc_m)[3]))
cat("DeLong: GNRI vs 多因素 P =", formatC(roc.test(roc_g, roc_m)$p.value, digits = 3), "\n")
cat("DeLong: NLR  vs 多因素 P =", formatC(roc.test(roc_n, roc_m)$p.value, digits = 3), "\n")

# ── 3. 校准曲线（十分位分箱）──────────────────────────────
val_res <- val.prob(prob_m, main$case_icu, pl = FALSE)
cat(sprintf("\n校准: 截距=%.3f 斜率=%.3f\n", val_res["Intercept"], val_res["Slope"]))

# ── 4. DCA ────────────────────────────────────────────────
dca_calc <- function(y, p) {
  n <- length(y); pts <- seq(0.01, 0.99, by = 0.01)
  nb <- sapply(pts, function(pt) {
    tp <- sum(p >= pt & y == 1); fp <- sum(p >= pt & y == 0)
    tp/n - fp/n * pt/(1 - pt)
  })
  nb_all <- mean(y) - (1 - mean(y)) * pts/(1 - pts)
  data.frame(pt = pts, model = nb, all = nb_all, none = 0)
}
dca_m <- dca_calc(main$case_icu, prob_m)
dca_g <- dca_calc(main$case_icu, predict(glm(case_icu ~ gnri, data = main, family = binomial), type = "response"))
dca_n <- dca_calc(main$case_icu, predict(glm(case_icu ~ nlr, data = main, family = binomial), type = "response"))
cat("── DCA 净收益（pt=0.1/0.2/0.3）──\n")
for (pt in c(0.1, 0.2, 0.3)) {
  cat(sprintf("  pt=%.1f: 多因素%.4f | GNRI %.4f | NLR %.4f | 全转 %.4f\n", pt,
              dca_m$model[which.min(abs(dca_m$pt - pt))],
              dca_g$model[which.min(abs(dca_g$pt - pt))],
              dca_n$model[which.min(abs(dca_n$pt - pt))],
              dca_m$all[which.min(abs(dca_m$pt - pt))]))
}

# ── 5. 森林图（多因素 OR）────────────────────────────────
or_m <- or_table[or_table$分析 == "多因素", ]
pdf(file.path(OUT, "森林图_多因素_演练.pdf"), width = 8, height = 5)
par(mar = c(4, 7, 3, 6), family = "SimHei")
plot(NA, xlim = c(0.2, max(or_m$高CI) * 1.15), ylim = c(7.5, 0.5),
     xlab = "OR (95% CI) — 多因素模型预测术后ICU转入", ylab = "", yaxt = "n",
     main = "多因素 Logistic 回归森林图（演练）", cex.main = 1.1)
axis(2, at = 1:7, labels = or_m$变量, las = 1, cex.axis = 0.9)
abline(v = 1, lty = 2, col = "gray60")
for (i in 1:7) {
  col_i <- ifelse(or_m$OR[i] > 1, "#C00000", "#2F5597")
  points(or_m$OR[i], i, pch = 15, cex = 1.3, col = col_i)
  segments(or_m$低CI[i], i, or_m$高CI[i], i, lwd = 2.2, col = col_i)
  text(max(or_m$高CI) * 1.14, i, sprintf("%.2f (%.2f-%.2f)", or_m$OR[i], or_m$低CI[i], or_m$高CI[i]),
       cex = 0.78, adj = 1, col = col_i)
}
dev.off()

# ── 6. ROC 图（多因素 vs 单指标）──────────────────────────
pdf(file.path(OUT, "ROC_多因素_演练.pdf"), width = 7, height = 6.5)
par(mar = c(5, 5, 3, 2), family = "SimHei")
plot(roc_g, col = "#2F5597", lwd = 2.5, main = "预测术后ICU转入：多因素模型 vs 单指标（演练）",
     xlab = "1 - 特异性", ylab = "敏感性", cex.main = 1.1)
plot(roc_n, col = "#C00000", lwd = 2.5, add = TRUE)
plot(roc_m, col = "#548235", lwd = 3, lty = 2, add = TRUE)
legend("bottomright",
       legend = c(sprintf("GNRI (AUC=%.3f)", auc(roc_g)),
                  sprintf("NLR (AUC=%.3f)", auc(roc_n)),
                  sprintf("多因素 (AUC=%.3f)", auc(roc_m))),
       col = c("#2F5597", "#C00000", "#548235"), lwd = c(2.5, 2.5, 3),
       lty = c(1, 1, 2), bty = "n", cex = 0.95)
dev.off()

# ── 7. 校准曲线 ──────────────────────────────────────────
pdf(file.path(OUT, "校准曲线_多因素_演练.pdf"), width = 6.5, height = 6.5)
par(mar = c(5, 5, 3, 2), family = "SimHei")
q10 <- quantile(prob_m, probs = seq(0, 1, 0.1)); q10[1] <- 0; q10[11] <- 1
bin_idx <- cut(prob_m, breaks = unique(q10), include.lowest = TRUE)
bin_pred <- tapply(prob_m, bin_idx, mean); bin_obs <- tapply(main$case_icu, bin_idx, mean)
plot(bin_pred, bin_obs, xlim = c(0, 1), ylim = c(0, 1),
     pch = 19, cex = 1.3, col = "#2F5597",
     xlab = "预测概率（十分位）", ylab = "实际病例比例",
     main = "多因素模型校准曲线（演练）", cex.main = 1.1)
abline(0, 1, lty = 2, col = "gray60")
lines(bin_pred, bin_obs, col = "#2F5597")
legend("topleft", c("理想校准线", "模型校准"), col = c("gray60", "#2F5597"),
       lty = c(2, 1), lwd = 1.5, bty = "n", cex = 0.9)
dev.off()

# ── 8. DCA 图 ────────────────────────────────────────────
pdf(file.path(OUT, "DCA_多因素_演练.pdf"), width = 7, height = 6)
par(mar = c(5, 5, 3, 2), family = "SimHei")
plot(dca_m$pt, dca_m$model, type = "l", lwd = 2.5, col = "#548235",
     xlim = c(0, 1), ylim = c(-0.05, max(dca_m$model, dca_m$all) * 1.1),
     xlab = "阈值概率", ylab = "净收益", main = "决策曲线分析（演练）", cex.main = 1.1)
lines(dca_g$pt, dca_g$model, lwd = 2, col = "#2F5597")
lines(dca_n$pt, dca_n$model, lwd = 2, col = "#C00000")
lines(dca_m$pt, dca_m$all, lwd = 1.5, lty = 2, col = "gray50")
abline(h = 0, lty = 2, col = "gray80")
legend("topright",
       legend = c("多因素模型", "GNRI", "NLR", "全部转入", "不转入"),
       col = c("#548235", "#2F5597", "#C00000", "gray50", "gray80"),
       lwd = c(2.5, 2, 2, 1.5, 1), lty = c(1, 1, 1, 2, 2), bty = "n", cex = 0.9)
dev.off()

# ── 9. 死亡敏感性（单因素OR）──────────────────────────────
cat("\n── 死亡敏感性分析（15例 vs 313对照）──\n")
sens_res <- data.frame()
for (v in c("gnri", "nlr", "age")) {
  f <- as.formula(paste("case_death ~", v))
  m <- glm(f, data = sens, family = binomial)
  s <- summary(m); ci <- confint(m)
  sens_res <- rbind(sens_res, data.frame(
    变量 = c("GNRI(每+1)", "NLR(每+1)", "年龄(每+1)")[which(c("gnri", "nlr", "age") == v)],
    OR = exp(coef(m)[2]), 低CI = exp(ci[2, 1]), 高CI = exp(ci[2, 2]),
    P = formatC(s$coef[2, 4], format = "e", digits = 2)))
}
print(sens_res, row.names = FALSE)
write.csv(sens_res, file.path(OUT, "结果_死亡敏感性_演练.csv"), row.names = FALSE,
          fileEncoding = "UTF-8")

cat("\n✅ 多因素演练完成！产出: 结果_多因素OR_演练.csv + 4张图 + 死亡敏感性\n")
