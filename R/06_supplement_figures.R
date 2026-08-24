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
# 投稿级图表补齐（Figure 1/4/5/6 + Table 1/3）
# Figure 1 : STROBE 流程图（研究设计与患者筛选）
# Figure 4 : 校准曲线配对（开发 vs 验证）
# Figure 5 : DCA 配对（开发 vs 验证）
# Figure 6 : GNRI/NLR 分层剂量反应（开发 vs 验证，本地数据）
# Table 1  : 双库基线特征表
# Table 3  : 模型性能指标（截断值/Se/Sp/PPV/NPV/ACC）
# ============================================================
suppressPackageStartupMessages({
  library(pROC)
  library(showtext)
  library(sysfonts)
})
font_add("SimHei", "C:/Windows/Fonts/simhei.ttf")
showtext_auto()

OUT <- "BASE/投稿图表集"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ═══════════════════════════════════════════════
# 数据准备：公共库（开发） + 本地（验证）
# ═══════════════════════════════════════════════
pub <- read.csv("file.path(BASE, "data/deriv/cohort_analysis.csv")",
                na.strings = c("", "NA", "nan"))
pub$icu <- pub$icu_admit
pub$sex_m <- ifelse(pub$sex == "M", 1, 0)
pub <- pub[pub$db != "eicu", ]
pub <- pub[complete.cases(pub[c("icu", "gnri", "nlr", "age", "sex_m", "bmi", "hb", "crp", "hip_fracture")]), ]
m_pub <- glm(icu ~ gnri + nlr + age + sex_m + bmi + hb + crp + hip_fracture,
             data = pub, family = binomial)
prob_pub <- predict(m_pub, type = "response")

loc <- read.csv("LOCAL_DATA  # not included in repo; provide your local validation data",
                fileEncoding = "UTF-8")
loc$case_icu <- ifelse(loc$icu_transfer == "是", 1, 0)
loc$gender_m <- ifelse(loc$gender == "男", 1, 0)
loc$hip_fracture <- ifelse(loc$surgery_type %in% c("骨折内固定", "髋关节置换"), 1, 0)
loc <- loc[loc$case_icu == 1 | loc$group == "对照组", ]
m_loc <- glm(case_icu ~ gnri + nlr + age + gender_m + bmi + charlson + hip_fracture,
             data = loc, family = binomial)
prob_loc <- predict(m_loc, type = "response")

cat("开发集:", nrow(pub), "例 | ICU转入", sum(pub$icu), "\n")
cat("验证集:", nrow(loc), "例 | ICU转入", sum(loc$case_icu), "\n")

# ═══════════════════════════════════════════════
# Figure 4：校准曲线配对（开发 vs 验证）
# ═══════════════════════════════════════════════
cal_plot <- function(y, p, main, col) {
  q <- quantile(p, probs = seq(0, 1, 0.1)); q[1] <- 0; q[11] <- 1
  b <- cut(p, breaks = unique(q), include.lowest = TRUE)
  pred <- tapply(p, b, mean); obs <- tapply(y, b, mean)
  plot(pred, obs, pch = 19, cex = 1.2, col = col,
       xlim = c(0, max(p, 0.1)), ylim = c(0, max(p, 0.1)),
       xlab = "预测概率", ylab = "实际发生比例", main = main,
       cex.main = 1.1, cex.lab = 1.0)
  abline(0, 1, lty = 2, col = "gray60")
  lines(pred, obs, col = col, lwd = 1.5)
  invisible(c(intercept = 0, slope = 1))
}
pdf(file.path(OUT, "Figure4_校准曲线_配对.pdf"), width = 11, height = 5.2)
par(mfrow = c(1, 2), family = "SimHei")
cal_plot(pub$icu, prob_pub, "A. 开发集（MIMIC+INSPIRE）", "#2F5597")
cal_plot(loc$case_icu, prob_loc, "B. 验证集（本地400例）", "#C00000")
dev.off()

# ═══════════════════════════════════════════════
# Figure 5：DCA 配对（开发 vs 验证）
# ═══════════════════════════════════════════════
dca_calc <- function(y, p) {
  n <- length(y); pts <- seq(0.01, 0.99, by = 0.01)
  nb <- sapply(pts, function(pt) {
    tp <- sum(p >= pt & y == 1); fp <- sum(p >= pt & y == 0)
    tp/n - fp/n * pt/(1 - pt)
  })
  nb_all <- mean(y) - (1 - mean(y)) * pts/(1 - pts)
  data.frame(pt = pts, model = nb, all = nb_all, none = 0)
}
dca_pub <- dca_calc(pub$icu, prob_pub)
dca_loc <- dca_calc(loc$case_icu, prob_loc)

pdf(file.path(OUT, "Figure5_DCA_配对.pdf"), width = 11, height = 5.2)
par(mfrow = c(1, 2), family = "SimHei")
for (i in 1:2) {
  d <- if (i == 1) dca_pub else dca_loc
  ylab2 <- if (i == 1) pub$icu else loc$case_icu
  main <- if (i == 1) "A. 开发集（MIMIC+INSPIRE）" else "B. 验证集（本地400例）"
  plot(d$pt, d$model, type = "l", lwd = 2.5, col = "#548235",
       xlim = c(0, 1), ylim = c(-0.05, max(d$model, d$all) * 1.1),
       xlab = "阈值概率", ylab = "净收益", main = main, cex.main = 1.1)
  lines(d$pt, d$all, lwd = 1.5, lty = 2, col = "gray50")
  abline(h = 0, lty = 2, col = "gray80")
  if (i == 1) legend("topright", c("多因素模型", "全部转入", "不转入"),
                     col = c("#548235", "gray50", "gray80"), lwd = c(2.5, 1.5, 1),
                     lty = c(1, 2, 2), bty = "n", cex = 0.85)
}
dev.off()

# ═══════════════════════════════════════════════
# Figure 6：GNRI/NLR 分层剂量反应（验证集，本地）
# ═══════════════════════════════════════════════
loc$gnri_grp <- cut(loc$gnri, breaks = c(0, 82, 92, 98, Inf),
                    labels = c("严重(<82)", "中度(82-92)", "轻度(92-98)", "正常(≥98)"))
loc$nlr_grp <- ifelse(loc$nlr >= 2.96, "高NLR(≥2.96)", "低NLR(<2.96)")

pdf(file.path(OUT, "Figure6_分层剂量反应_验证集.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2), family = "SimHei")
# GNRI 分层
sg <- aggregate(case_icu ~ gnri_grp, data = loc, FUN = function(x) c(n = length(x), rate = mean(x) * 100))
sg <- do.call(rbind, lapply(split(loc, loc$gnri_grp), function(d) data.frame(
  grp = d$gnri_grp[1], n = nrow(d), rate = mean(d$case_icu) * 100)))
sg <- sg[order(match(sg$grp, c("严重(<82)", "中度(82-92)", "轻度(92-98)", "正常(≥98)"))), ]
bp <- barplot(sg$rate, names.arg = sg$grp, col = c("#C00000", "#ED7D31", "#FFC000", "#70AD47"),
              ylim = c(0, 100), border = NA, main = "A. GNRI 营养风险分层",
              xlab = "GNRI 分层", ylab = "ICU转入率（%）", cex.main = 1.1, cex.names = 0.9)
text(bp, sg$rate + 3, sprintf("%.1f%%", sg$rate), cex = 0.9)
text(bp, 8, sprintf("n=%d", sg$n), cex = 0.75, col = "white", font = 2)
# NLR 分层
sn <- aggregate(case_icu ~ nlr_grp, data = loc, FUN = mean)
sn$n <- table(loc$nlr_grp)
bp2 <- barplot(sn$case_icu * 100, names.arg = c("低NLR(<2.96)", "高NLR(≥2.96)"),
               col = c("#70AD47", "#C00000"), ylim = c(0, 100), border = NA,
               main = "B. NLR 分层", xlab = "NLR 分层", ylab = "ICU转入率（%）",
               cex.main = 1.1, cex.names = 0.9)
text(bp2, sn$case_icu * 100 + 3, sprintf("%.1f%%", sn$case_icu * 100), cex = 1.0)
text(bp2, 8, sprintf("n=%d", sn$n), cex = 0.8, col = "white", font = 2)
dev.off()

# ═══════════════════════════════════════════════
# Table 1：双库基线特征
# ═══════════════════════════════════════════════
mk_row <- function(label, p1, l1, p2, l2) data.frame(变量 = label, 公共库_病例 = p1, 公共库_对照 = l1, 本地_病例 = p2, 本地_对照 = l2)
t1 <- rbind(
  mk_row("n", sum(pub$icu), nrow(pub) - sum(pub$icu), sum(loc$case_icu), sum(loc$case_icu == 0)),
  mk_row("年龄(岁)", sprintf("%.1f±%.1f", mean(pub$age[pub$icu==1]), sd(pub$age[pub$icu==1])),
         sprintf("%.1f±%.1f", mean(pub$age[pub$icu==0]), sd(pub$age[pub$icu==0])),
         sprintf("%.1f±%.1f", mean(loc$age[loc$case_icu==1]), sd(loc$age[loc$case_icu==1])),
         sprintf("%.1f±%.1f", mean(loc$age[loc$case_icu==0]), sd(loc$age[loc$case_icu==0]))),
  mk_row("男性(%)", sprintf("%.1f", mean(pub$sex_m[pub$icu==1])*100),
         sprintf("%.1f", mean(pub$sex_m[pub$icu==0])*100),
         sprintf("%.1f", mean(loc$gender_m[loc$case_icu==1])*100),
         sprintf("%.1f", mean(loc$gender_m[loc$case_icu==0])*100)),
  mk_row("BMI(kg/m²)", sprintf("%.1f±%.1f", mean(pub$bmi[pub$icu==1]), sd(pub$bmi[pub$icu==1])),
         sprintf("%.1f±%.1f", mean(pub$bmi[pub$icu==0]), sd(pub$bmi[pub$icu==0])),
         sprintf("%.1f±%.1f", mean(loc$bmi[loc$case_icu==1]), sd(loc$bmi[loc$case_icu==1])),
         sprintf("%.1f±%.1f", mean(loc$bmi[loc$case_icu==0]), sd(loc$bmi[loc$case_icu==0]))),
  mk_row("白蛋白(g/dL)", sprintf("%.2f±%.2f", mean(pub$albumin[pub$icu==1]), sd(pub$albumin[pub$icu==1])),
         sprintf("%.2f±%.2f", mean(pub$albumin[pub$icu==0]), sd(pub$albumin[pub$icu==0])),
         sprintf("%.2f±%.2f", mean(loc$albumin[loc$case_icu==1]), sd(loc$albumin[loc$case_icu==1])),
         sprintf("%.2f±%.2f", mean(loc$albumin[loc$case_icu==0]), sd(loc$albumin[loc$case_icu==0]))),
  mk_row("GNRI", sprintf("%.1f±%.1f", mean(pub$gnri[pub$icu==1]), sd(pub$gnri[pub$icu==1])),
         sprintf("%.1f±%.1f", mean(pub$gnri[pub$icu==0]), sd(pub$gnri[pub$icu==0])),
         sprintf("%.1f±%.1f", mean(loc$gnri[loc$case_icu==1]), sd(loc$gnri[loc$case_icu==1])),
         sprintf("%.1f±%.1f", mean(loc$gnri[loc$case_icu==0]), sd(loc$gnri[loc$case_icu==0]))),
  mk_row("NLR(中位数[IQR])", sprintf("%.2f [%.2f, %.2f]", median(pub$nlr[pub$icu==1]), quantile(pub$nlr[pub$icu==1], .25), quantile(pub$nlr[pub$icu==1], .75)),
         sprintf("%.2f [%.2f, %.2f]", median(pub$nlr[pub$icu==0]), quantile(pub$nlr[pub$icu==0], .25), quantile(pub$nlr[pub$icu==0], .75)),
         sprintf("%.2f [%.2f, %.2f]", median(loc$nlr[loc$case_icu==1]), quantile(loc$nlr[loc$case_icu==1], .25), quantile(loc$nlr[loc$case_icu==1], .75)),
         sprintf("%.2f [%.2f, %.2f]", median(loc$nlr[loc$case_icu==0]), quantile(loc$nlr[loc$case_icu==0], .25), quantile(loc$nlr[loc$case_icu==0], .75))),
  mk_row("Charlson", "—（缺失88.4%，未纳入）", "—",
         sprintf("%.1f±%.1f", mean(loc$charlson[loc$case_icu==1]), sd(loc$charlson[loc$case_icu==1])),
         sprintf("%.1f±%.1f", mean(loc$charlson[loc$case_icu==0]), sd(loc$charlson[loc$case_icu==0]))),
  mk_row("髋部骨折(%)", sprintf("%.1f", mean(pub$hip_fracture[pub$icu==1])*100),
         sprintf("%.1f", mean(pub$hip_fracture[pub$icu==0])*100),
         sprintf("%.1f", mean(loc$hip_fracture[loc$case_icu==1])*100),
         sprintf("%.1f", mean(loc$hip_fracture[loc$case_icu==0])*100))
)
print(t1, row.names = FALSE)
write.csv(t1, file.path(OUT, "Table1_双库基线特征.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# ═══════════════════════════════════════════════
# Table 3：模型性能指标（截断值/Se/Sp/PPV/NPV/ACC）
# ═══════════════════════════════════════════════
perf_row <- function(y, p, label) {
  ro <- roc(y, p, quiet = TRUE)
  ct <- coords(ro, "best", ret = c("threshold", "sensitivity", "specificity",
                                   "ppv", "npv", "accuracy"))
  data.frame(数据集 = label, AUC = round(auc(ro), 3),
             截断值 = round(ct$threshold, 3), 敏感度 = round(ct$sensitivity, 3),
             特异度 = round(ct$specificity, 3), PPV = round(ct$ppv, 3),
             NPV = round(ct$npv, 3), 准确率 = round(ct$accuracy, 3))
}
t3 <- rbind(
  perf_row(pub$icu, prob_pub, "开发集-多因素模型"),
  perf_row(loc$case_icu, prob_loc, "验证集-多因素模型"),
  perf_row(pub$icu, pub$gnri, "开发集-GNRI"),
  perf_row(loc$case_icu, loc$gnri, "验证集-GNRI"),
  perf_row(pub$icu, pub$nlr, "开发集-NLR"),
  perf_row(loc$case_icu, loc$nlr, "验证集-NLR")
)
print(t3, row.names = FALSE)
write.csv(t3, file.path(OUT, "Table3_模型性能指标.csv"), row.names = FALSE, fileEncoding = "UTF-8")

cat("\n✅ 补齐图表完成: Figure4/5/6 + Table1/3\n")
