# -*- coding: utf-8 -*-
"""Regenerate Figure 1 (study flow) for EN + ZH as CLEAN flow diagrams.
Removes the spoiler content: model AUC, external-validation AUC, and the
sensitivity-analysis OR box. Style follows the original R/matplotlib scripts.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib import font_manager

BASE = r"D:\BaiduSyncdisk\MIMIC\营养\老年骨科手术\病例对照研究\Frontiers专题投稿"

# ---------- fonts ----------
def add_font(path, name):
    try:
        font_manager.fontManager.addfont(path)
        return name
    except Exception:
        return None

# EN: DejaVu Sans has a true bold face -> clean bold cohort titles
en_font = "DejaVu Sans"
zh_font = add_font(r"C:\Windows\Fonts\simhei.ttf", "SimHei") or "DejaVu Sans"

# =====================================================================
# ENGLISH
# =====================================================================
def build_en():
    # Two-line cohort titles, boxed and placed closer to the first boxes.
    fig, ax = plt.subplots(figsize=(180/25.4, 124/25.4))
    ax.set_xlim(0, 10); ax.set_ylim(2.8, 11.4); ax.axis("off")

    dev_col, dev_ec = "#EAF2FB", "#3182BD"
    val_col, val_ec = "#FDF0EF", "#D24B40"
    excl_fc, excl_ec = "#FCEBEB", "#C0504D"
    arrow_col = "#666666"

    # Top boxes shifted down; split the long first-line label to avoid
    # the left/right columns colliding across the centre gap.
    dev_boxes = [
        ((8.50, 9.50), ["Elderly orthopedic", "surgical patients (≥65 y)",
                          "MIMIC-IV + INSPIRE + eICU", "n = 7,203"], dev_col, dev_ec, True),
        ((7.15, 8.15), ["Excluded: eICU database",
                         "(ICU transfer rate 100%)", "n = 773"], excl_fc, excl_ec, False),
        ((5.35, 6.35), ["Excluded: missing key variables",
                         "(albumin / height / weight / CBC)", "n = 1,106"], excl_fc, excl_ec, False),
        ((3.55, 4.55), ["Complete-case cohort for analysis",
                         "n = 5,324 (ICU transfer 639, 12.0%)"], dev_col, dev_ec, True),
    ]
    val_boxes = [
        ((8.50, 9.50), ["Elderly orthopedic", "surgical patients (≥65 y)",
                          "Local hospital cohort", "n > 7,000"], val_col, val_ec, True),
        ((7.15, 8.15), ["Case-control matching",
                         "(age ±2 y + sex + surgery type)"], val_col, val_ec, False),
        ((5.35, 6.35), ["Included in study",
                         "n = 400 (82 cases / 318 controls)"], dev_col, dev_ec, True),
        ((3.55, 4.55), ["Primary analysis cohort",
                         "n = 395 (77 ICU transfers + 318 controls)"], dev_col, dev_ec, True),
    ]

    def draw_box(cx, yb_yt, lines, fc, ec, bold):
        yb, yt = yb_yt
        h = yt - yb
        x0 = cx - 2.15
        ax.add_patch(Rectangle((x0, yb), 4.3, h, facecolor=fc, edgecolor=ec,
                               linewidth=1.1, zorder=2))
        n = len(lines)
        fs = 7.0 if n >= 4 else 7.6
        for i, ln in enumerate(lines):
            y = yb + h * (n - i + 0.5) / (n + 1)
            ax.text(cx, y, ln, ha="center", va="center", fontsize=fs,
                    fontweight="bold" if bold else "normal",
                    fontfamily=en_font, color="#1A1A1A", zorder=3)

    for b in dev_boxes: draw_box(2.7, *b)
    for b in val_boxes: draw_box(7.3, *b)

    # Two-line cohort titles, boxed and placed closer to the first flow boxes
    # so the labels read as part of each column and stay inside the frame.
    ax.text(2.7, 10.15, "Development cohort\n(public databases)", ha="center",
            va="center", fontsize=8.5, fontweight="bold", fontfamily=en_font,
            color=dev_ec, linespacing=1.10, zorder=4,
            bbox=dict(boxstyle="round,pad=0.18", facecolor=dev_col,
                      edgecolor=dev_ec, linewidth=1.0, alpha=0.95, zorder=3))
    ax.text(7.3, 10.15, "Validation cohort\n(local hospital)", ha="center",
            va="center", fontsize=8.5, fontweight="bold", fontfamily=en_font,
            color=val_ec, linespacing=1.10, zorder=4,
            bbox=dict(boxstyle="round,pad=0.18", facecolor=val_col,
                      edgecolor=val_ec, linewidth=1.0, alpha=0.95, zorder=3))

    def arrow(cx, y_from, y_to):
        ax.annotate("", xy=(cx, y_to), xytext=(cx, y_from),
                    arrowprops=dict(arrowstyle="-|>", color=arrow_col, lw=1.1), zorder=1)

    for cx in (2.7, 7.3):
        arrow(cx, 8.50, 8.15)
        arrow(cx, 7.15, 6.35)
        arrow(cx, 5.35, 4.55)

    out = f"{BASE}/图表_en/Figure1_StudyFlow"
    fig.savefig(out + ".png", dpi=300, bbox_inches="tight", pad_inches=0.05)
    fig.savefig(out + ".pdf", dpi=300, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print("EN Figure 1 (clean) written")

# =====================================================================
# CHINESE
# =====================================================================
def build_zh():
    fig, ax = plt.subplots(figsize=(12, 8.5))
    ax.set_xlim(0, 12); ax.set_ylim(3.4, 11); ax.axis("off")

    dev_fc, dev_ec = "#E6F1FB", "#185FA5"
    val_fc, val_ec = "#E6F1FB", "#185FA5"   # right column kept blue in original
    excl_fc, excl_ec = "#FCEBEB", "#A32D2D"
    inc_fc, inc_ec = "#EAF3DE", "#3B6D11"
    arrow_col = "#555555"

    # (x0, y0, w, h, text, fc, ec, bold, fs)
    left = [
        (0.6, 8.6, 4.8, 1.0, "老年骨科手术患者（≥65岁）\nMIMIC-IV + INSPIRE + eICU\nn = 7,203", dev_fc, dev_ec, True, 10.5),
        (0.6, 7.3, 4.8, 1.0, "排除 eICU（ICU数据库，\n转入率恒为100%，不适用转入预测）\nn = 773", excl_fc, excl_ec, False, 9.5),
        (0.6, 5.9, 4.8, 1.1, "排除关键变量缺失\n（白蛋白/身高/体重/中性粒/淋巴/HB/CRP）\nn = 1,106", excl_fc, excl_ec, False, 9.5),
        (0.6, 4.4, 4.8, 1.2, "多因素分析完整队列\nn = 5,324\nICU转入 639（12.0%）/ 未转入 4,685", inc_fc, inc_ec, True, 10.5),
    ]
    right = [
        (6.6, 8.6, 4.8, 1.0, "老年骨科手术患者（≥65岁）\n本地医院队列\nn > 7,000", dev_fc, dev_ec, True, 10.5),
        (6.6, 7.3, 4.8, 1.0, "病例对照匹配\n（年龄±2岁 + 性别 + 手术类型）\n77例ICU转入 + 15例死亡（重叠10）", dev_fc, dev_ec, False, 9.5),
        (6.6, 5.9, 4.8, 1.1, "纳入研究：400 例\n病例 82（67纯转ICU + 10两者 + 5纯死亡）\n对照 318", inc_fc, inc_ec, True, 9.5),
        (6.6, 4.4, 4.8, 1.2, "主分析队列：395 例\nICU转入 77 + 对照 318\n（排除10例仅死亡）", inc_fc, inc_ec, True, 10.5),
    ]

    def draw_box(spec):
        x0, y0, w, h, txt, fc, ec, bold, fs = spec
        ax.add_patch(FancyBboxPatch((x0, y0), w, h, boxstyle="round,pad=0.02",
                                    fc=fc, ec=ec, lw=1.2, mutation_aspect=1, zorder=2))
        ax.text(x0 + w/2, y0 + h/2, txt, ha="center", va="center",
                fontsize=fs, fontweight="bold" if bold else "normal",
                fontfamily=zh_font, color="#1A1A1A", linespacing=1.4, zorder=3)

    for s in left: draw_box(s)
    for s in right: draw_box(s)

    ax.text(6, 10.6, "老年骨科手术患者术后ICU转入预测：研究设计与患者筛选流程",
            ha="center", fontsize=13.5, fontweight="bold", fontfamily=zh_font, color="#2F5597")
    ax.text(3, 9.85, "开发队列（公共数据库）", ha="center", va="center",
            fontsize=12, fontweight="bold", fontfamily=zh_font, color="#2F5597",
            zorder=4, bbox=dict(boxstyle="round,pad=0.18", facecolor=dev_fc,
                                edgecolor="#2F5597", linewidth=1.0, alpha=0.95, zorder=3))
    ax.text(9, 9.85, "验证队列（本地医院）", ha="center", va="center",
            fontsize=12, fontweight="bold", fontfamily=zh_font, color="#C00000",
            zorder=4, bbox=dict(boxstyle="round,pad=0.18", facecolor="#FDF0EF",
                                edgecolor="#C00000", linewidth=1.0, alpha=0.95, zorder=3))

    def arrow(cx, y_from, y_to):
        ax.annotate("", xy=(cx, y_to), xytext=(cx, y_from),
                    arrowprops=dict(arrowstyle="-|>", color=arrow_col, lw=1.4), zorder=1)

    for cx in (3, 9):
        arrow(cx, 8.6, 8.3)
        arrow(cx, 7.3, 7.0)
        arrow(cx, 5.9, 5.6)

    out = f"{BASE}/图表/Figure1_STROBE流程图"
    fig.savefig(out + ".png", dpi=200, bbox_inches="tight", pad_inches=0.05)
    fig.savefig(out + ".pdf", dpi=200, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print("ZH Figure 1 (clean) written")

if __name__ == "__main__":
    build_en()
    build_zh()
    print("DONE")
