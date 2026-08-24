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
# 论文级整合图表：公共库（MIMIC+INSPIRE） vs 本地（400例病例对照）
# 图1：跨库 ROC 叠加（多因素模型 + GNRI/NLR 单指标）
# 图2：跨库多因素 OR 森林图（分组对比）
# 图3：AUC 对比条形图（模型 × 库）
# 表 ：跨库结果汇总 CSV
# ============================================================
suppressPackageStartupMessages({
  library(pROC)
  library(showtext)
  library(sysfonts)
})
font_add("SimHei", "C:/Windows/Fonts/simhei.ttf")
showtext_auto()

OUT <- "BASE"

# ── 公共库模型（完整病例 5,324）──────────────────────────
pub <- read.csv("file.path(BASE, "data/deriv/cohort_analysis.csv")",
                na.strings = c("", "NA", "nan"))
pub$icu <- pub$icu_admit
pub$sex_m <- ifelse(pub$sex == "M", 1, 0)
pub <- pub[pub$db != "eicu", ]
pub <- pub[complete.cases(pub[c("icu", "gnri", "nlr", "age", "sex_m", "bmi", "hb", "crp", "hip_fracture")]), ]
m_pub <- glm(icu ~ gnri + nlr + age + sex_m + bmi + hb + crp + hip_fracture,
             data = pub, family = binomial)
roc_pub_m <- roc(pub$icu, predict(m_pub, type = "response"), quiet = TRUE)
roc_pub_g <- roc(pub$icu, pub$gnri, quiet = TRUE)
roc_pub_n <- roc(pub$icu, pub$nlr, quiet = TRUE)

# ── 本地模型（395 例：77 ICU转入 + 318 对照）─────────────
loc <- read.csv(file.path(OUT, "本地数据_utf8.csv"), fileEncoding = "UTF-8")
loc$case_icu <- ifelse(loc$icu_transfer == "是", 1, 0)
loc$gender_m <- ifelse(loc$gender == "男", 1, 0)
loc$hip_fracture <- ifelse(loc$surgery_type %in% c("骨折内固定", "髋关节置换"), 1, 0)
loc <- loc[loc$case_icu == 1 | loc$group == "对照组", ]
m_loc <- glm(case_icu ~ gnri + nlr + age + gender_m + bmi + charlson + hip_fracture,
             data = loc, family = binomial)
roc_loc_m <- roc(loc$case_icu, predict(m_loc, type = "response"), quiet = TRUE)
roc_loc_g <- roc(loc$case_icu, loc$gnri, quiet = TRUE)
roc_loc_n <- roc(loc$case_icu, loc$nlr, quiet = TRUE)

cat("公共库: n=", nrow(pub), " 多因素AUC=", round(auc(roc_pub_m), 3), "\n", sep = "")
cat("本地:   n=", nrow(loc), " 多因素AUC=", round(auc(roc_loc_m), 3), "\n", sep = "")

# ═══════════════════════════════════════════════════════
# 图1：跨库 ROC 叠加
# ═══════════════════════════════════════════════════════
pdf(file.path(OUT, "整合_图1_跨库ROC.pdf"), width = 7.5, height = 6.8)
par(mar = c(5, 5, 3, 2), family = "SimHei")
plot(roc_pub_m, col = "#548235", lwd = 3,
     main = "预测术后ICU转入：公共库 vs 本地验证（多因素模型）",
     xlab = "1 - 特异性", ylab = "敏感性", cex.main = 1.15)
plot(roc_loc_m, col = "#548235", lwd = 3, lty = 2, add = TRUE)
plot(roc_pub_g, col = "#2F5597", lwd = 1.8, add = TRUE)
plot(roc_loc_g, col = "#2F5597", lwd = 1.8, lty = 2, add = TRUE)
plot(roc_pub_n, col = "#C00000", lwd = 1.8, add = TRUE)
plot(roc_loc_n, col = "#C00000", lwd = 1.8, lty = 2, add = TRUE)
legend("bottomright",
       legend = c(sprintf("公共库 多因素 (AUC=%.3f)", auc(roc_pub_m)),
                  sprintf("本地   多因素 (AUC=%.3f)", auc(roc_loc_m)),
                  sprintf("公共库 GNRI (AUC=%.3f)", auc(roc_pub_g)),
                  sprintf("本地   GNRI (AUC=%.3f)", auc(roc_loc_g)),
                  sprintf("公共库 NLR (AUC=%.3f)", auc(roc_pub_n)),
                  sprintf("本地   NLR (AUC=%.3f)", auc(roc_loc_n))),
       col = rep(c("#548235", "#2F5597", "#C00000"), each = 2),
       lwd = c(3, 3, 1.8, 1.8, 1.8, 1.8), lty = c(1, 2, 1, 2, 1, 2),
       bty = "n", cex = 0.85)
dev.off()

# ═══════════════════════════════════════════════════════
# 图2：跨库多因素 OR 森林图
# ═══════════════════════════════════════════════════════
# 共同变量：GNRI/NLR/BMI/年龄/性别；各自特征：公共库HB/CRP/髋部骨折，本地Charlson/髋部手术
or_pub <- read.csv("file.path(BASE, "results")/logistic_多因素_ICU转入.csv",
                   fileEncoding = "UTF-8")
or_loc <- read.csv(file.path(OUT, "结果_多因素OR_本地.csv"), fileEncoding = "UTF-8")
or_pub_m <- or_pub[or_pub$分析 == "多因素", ]
or_loc_m <- or_loc[or_loc$分析 == "多因素", ]

common <- c("GNRI(每+1)", "NLR(每+1)", "BMI(每+1)", "年龄(每+1岁)")
pub_extra <- c("男性(vs女性)", "血红蛋白(每+1)", "CRP(每+1)", "髋部骨折(vs其他)")
loc_extra <- c("男性(vs女性)", "Charlson(每+1)", "髋部手术(vs其他)")

# 合并行
rows <- list()
for (v in common) {
  rp <- or_pub_m[or_pub_m$变量 == v, ]; rl <- or_loc_m[or_loc_m$变量 == v, ]
  rows[[length(rows) + 1]] <- data.frame(变量 = v, 库 = "公共库", OR = rp$OR, 低 = rp$低CI, 高 = rp$高CI)
  rows[[length(rows) + 1]] <- data.frame(变量 = v, 库 = "本地", OR = rl$OR, 低 = rl$低CI, 高 = rl$高CI)
}
for (v in pub_extra) {
  rp <- or_pub_m[or_pub_m$变量 == v, ]
  rows[[length(rows) + 1]] <- data.frame(变量 = v, 库 = "公共库", OR = rp$OR, 低 = rp$低CI, 高 = rp$高CI)
}
for (v in loc_extra) {
  rl <- or_loc_m[or_loc_m$变量 == v, ]
  rows[[length(rows) + 1]] <- data.frame(变量 = v, 库 = "本地", OR = rl$OR, 低 = rl$低CI, 高 = rl$高CI)
}
df_or <- do.call(rbind, rows)

# 变量顺序（公共在前）
var_order <- c(common, "男性(vs女性)", "血红蛋白(每+1)", "CRP(每+1)", "髋部骨折(vs其他)", "Charlson(每+1)", "髋部手术(vs其他)")
df_or$变量 <- factor(df_or$变量, levels = rev(var_order))

# 每个变量组两行（公共库/本地）
ypos <- numeric(nrow(df_or))
uniq_vars <- rev(var_order)
for (v in uniq_vars) {
  idx <- which(as.character(df_or$变量) == v)
  ypos[idx] <- which(uniq_vars == v) * 2 - ifelse(df_or$库[idx] == "公共库", 0.35, -0.35)
}

pdf(file.path(OUT, "整合_图2_跨库OR森林图.pdf"), width = 9, height = 7)
par(mar = c(4, 8, 3, 7), family = "SimHei")
xmax <- max(df_or$高) * 1.25
plot(NA, xlim = c(0.3, xmax), ylim = c(0.5, length(uniq_vars) * 2 + 0.5),
     xlab = "多因素 OR (95% CI) — 预测术后ICU转入", ylab = "", yaxt = "n",
     main = "跨库多因素 Logistic 回归：公共库 vs 本地验证", cex.main = 1.15)
# 变量标签（放在两行的中间）
axis(2, at = seq(1.5, length(uniq_vars) * 2 - 0.5, 2),
     labels = uniq_vars, las = 1, cex.axis = 0.85)
abline(v = 1, lty = 2, col = "gray60")
cols <- ifelse(df_or$库 == "公共库", "#2F5597", "#C00000")
for (i in 1:nrow(df_or)) {
  points(df_or$OR[i], ypos[i], pch = 15, cex = 1.2, col = cols[i])
  segments(df_or$低[i], ypos[i], df_or$高[i], ypos[i], lwd = 2, col = cols[i])
  text(xmax * 1.02, ypos[i], sprintf("%.2f (%.2f-%.2f)", df_or$OR[i], df_or$低[i], df_or$高[i]),
       cex = 0.7, adj = 0, col = cols[i])
}
legend("bottomleft", c("公共库（MIMIC+INSPIRE）", "本地（400例病例对照）"),
       col = c("#2F5597", "#C00000"), pch = 15, bty = "n", cex = 0.9)
dev.off()

# ═══════════════════════════════════════════════════════
# 图3：AUC 对比条形图
# ═══════════════════════════════════════════════════════
auc_df <- data.frame(
  库 = rep(c("公共库", "本地"), each = 3),
  模型 = rep(c("GNRI", "NLR", "多因素"), 2),
  AUC = c(auc(roc_pub_g), auc(roc_pub_n), auc(roc_pub_m),
          auc(roc_loc_g), auc(roc_loc_n), auc(roc_loc_m))
)
auc_df$AUC <- round(auc_df$AUC, 3)

pdf(file.path(OUT, "整合_图3_AUC对比.pdf"), width = 7, height = 5.5)
par(mar = c(5, 5, 4, 2), family = "SimHei")
bp <- barplot(matrix(auc_df$AUC, nrow = 2, byrow = TRUE),
              beside = TRUE, names.arg = c("GNRI", "NLR", "多因素模型"),
              col = c("#2F5597", "#C00000"), ylim = c(0, 1), border = NA,
              main = "预测术后ICU转入：AUC 跨库对比", ylab = "AUC",
              xlab = "模型", cex.main = 1.15, cex.names = 1.0)
vals <- matrix(auc_df$AUC, nrow = 2, byrow = TRUE)
for (i in 1:2) for (j in 1:3) {
  text(bp[i, j], vals[i, j] + 0.03, sprintf("%.3f", vals[i, j]), cex = 0.9)
}
legend("topleft", c("公共库（MIMIC+INSPIRE）", "本地（400例）"),
       fill = c("#2F5597", "#C00000"), bty = "n", cex = 0.9)
dev.off()

# ═══════════════════════════════════════════════════════
# 汇总表 CSV
# ═══════════════════════════════════════════════════════
summary_tab <- data.frame(
  项目 = c("样本量(n)", "ICU转入(n)", "事件率(%)",
           "GNRI AUC (95%CI)", "NLR AUC (95%CI)", "多因素 AUC (95%CI)",
           "GNRI 多因素OR", "NLR 多因素OR", "BMI 多因素OR"),
  公共库 = c(nrow(pub), sum(pub$icu), sprintf("%.1f", mean(pub$icu) * 100),
            sprintf("%.3f (%.3f-%.3f)", auc(roc_pub_g), ci.auc(roc_pub_g)[1], ci.auc(roc_pub_g)[3]),
            sprintf("%.3f (%.3f-%.3f)", auc(roc_pub_n), ci.auc(roc_pub_n)[1], ci.auc(roc_pub_n)[3]),
            sprintf("%.3f (%.3f-%.3f)", auc(roc_pub_m), ci.auc(roc_pub_m)[1], ci.auc(roc_pub_m)[3]),
            "0.925 (0.911-0.940)", "1.024 (1.007-1.040)", "1.171 (1.126-1.218)"),
  本地 = c(nrow(loc), sum(loc$case_icu), sprintf("%.1f", mean(loc$case_icu) * 100),
           sprintf("%.3f (%.3f-%.3f)", auc(roc_loc_g), ci.auc(roc_loc_g)[1], ci.auc(roc_loc_g)[3]),
           sprintf("%.3f (%.3f-%.3f)", auc(roc_loc_n), ci.auc(roc_loc_n)[1], ci.auc(roc_loc_n)[3]),
           sprintf("%.3f (%.3f-%.3f)", auc(roc_loc_m), ci.auc(roc_loc_m)[1], ci.auc(roc_loc_m)[3]),
           "0.889 (0.847-0.930)", "1.447 (1.228-1.735)", "1.278 (1.144-1.439)")
)
print(summary_tab, row.names = FALSE)
write.csv(summary_tab, file.path(OUT, "整合_跨库汇总表.csv"), row.names = FALSE, fileEncoding = "UTF-8")

cat("\n✅ 整合图表完成：图1跨库ROC / 图2跨库森林图 / 图3 AUC对比 + 汇总表\n")
