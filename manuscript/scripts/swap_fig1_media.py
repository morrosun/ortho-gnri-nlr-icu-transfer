# -*- coding: utf-8 -*-
"""Swap the embedded Figure 1 image inside an existing .docx with a new PNG,
and resize the drawing frame to match the new PNG's aspect ratio.
Identifies Figure 1 by locating its caption text, then the nearest <a:blip>
before that caption.
"""
import os, re, shutil, zipfile
from PIL import Image

BASE = r"D:\BaiduSyncdisk\MIMIC\营养\老年骨科手术\病例对照研究\Frontiers专题投稿"

def swap_clean(docx_path, new_png, marker):
    bk = docx_path + ".pre_figfix.bak"
    if not os.path.exists(bk):
        shutil.copy(docx_path, bk)
    with zipfile.ZipFile(docx_path, "r") as z:
        names = z.namelist()
        doc_xml = z.read("word/document.xml").decode("utf-8")
        rels_xml = z.read("word/_rels/document.xml.rels").decode("utf-8")
        data = {n: z.read(n) for n in names}

    idx = doc_xml.find(marker)
    if idx < 0:
        raise SystemExit(f"marker {marker!r} not found in {docx_path}")

    # Figure drawing appears BEFORE its caption -> take the last blip before marker.
    last_blip = None
    for mm in re.finditer(r'<a:blip[^>]*r:embed="([^"]+)"[^>]*>', doc_xml[:idx]):
        last_blip = mm
    if not last_blip:
        raise SystemExit(f"no blip before marker in {docx_path}")
    rid = last_blip.group(1)

    rel = re.search(r'Id="%s"[^>]*Target="([^"]+)"' % re.escape(rid), rels_xml)
    if not rel:
        raise SystemExit(f"rId {rid} not mapped in rels of {docx_path}")
    target = rel.group(1)
    if not target.startswith("word/media/"):
        target = "word/" + target.lstrip("/")
    if target not in data:
        raise SystemExit(f"media {target} not in docx")

    data[target] = open(new_png, "rb").read()

    # new aspect ratio
    with Image.open(new_png) as im:
        new_aspect = im.width / im.height

    # locate the <w:drawing> block containing this blip and update its extents
    start = doc_xml.rfind("<w:drawing>", 0, last_blip.start())
    end = doc_xml.find("</w:drawing>", last_blip.end()) + len("</w:drawing>")
    if start < 0 or end < len("</w:drawing>"):
        raise SystemExit(f"drawing block for blip {rid} not found")
    drawing = doc_xml[start:end]

    def update_extent(s, tag):
        pat = re.compile(r'<%s cx="(\d+)" cy="(\d+)"' % tag)
        m = pat.search(s)
        if m:
            cx = int(m.group(1))
            cy = int(round(cx / new_aspect))
            s = s[:m.start()] + f'<{tag} cx="{cx}" cy="{cy}"' + s[m.end():]
        return s

    drawing = update_extent(drawing, "wp:extent")
    drawing = update_extent(drawing, "a:ext")
    doc_xml = doc_xml[:start] + drawing + doc_xml[end:]
    data["word/document.xml"] = doc_xml.encode("utf-8")

    with zipfile.ZipFile(docx_path, "w", zipfile.ZIP_DEFLATED) as zo:
        for n in names:
            zo.writestr(n, data[n])
    print(f"swapped Figure 1 in {os.path.basename(docx_path)} -> {target} ({len(data[target])} bytes), aspect={new_aspect:.3f}")

if __name__ == "__main__":
    swap_clean(
        os.path.join(BASE, "paper_en_ICUtransfer_v13.docx"),
        os.path.join(BASE, "图表_en/Figure1_StudyFlow.png"),
        "Figure 1.")
    swap_clean(
        os.path.join(BASE, "paper_zh_ICUtransfer_v13.docx"),
        os.path.join(BASE, "图表/Figure1_STROBE流程图.png"),
        "图 1")
    swap_clean(
        os.path.join(BASE, "paper_en_ICUtransfer_v14.docx"),
        os.path.join(BASE, "图表_en/Figure1_StudyFlow.png"),
        "Figure 1.")
    swap_clean(
        os.path.join(BASE, "paper_zh_ICUtransfer_v14.docx"),
        os.path.join(BASE, "图表/Figure1_STROBE流程图.png"),
        "图 1")
    print("ALL SWAPPED")
