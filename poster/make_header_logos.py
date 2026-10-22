"""Build the logo assets for the poster header.

Requires reportlab, svglib and Pillow. The Nathoo Group wordmark uses the same
Crimson Pro Bold font as the live website. The UVic shield and university
wordmark come from the university's official white vector asset.
"""
from pathlib import Path
from PIL import Image

from reportlab.graphics import renderPDF
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas
from reportlab.lib.colors import white
from svglib.svglib import svg2rlg

ROOT = Path(__file__).resolve().parent
FIGURES = ROOT / "figures"

# Lossless format conversion for pdfLaTeX; retain the website's boxed artwork.
Image.open(FIGURES / "nathoo_group_site_original.webp").save(
    FIGURES / "nathoo_group_site_original.png"
)

crest = svg2rlg(str(FIGURES / "uvic_crest_white.svg"))
c = canvas.Canvas(str(FIGURES / "uvic_math_stats_white.pdf"), pagesize=(410, 82))
c.setTitle("University of Victoria Mathematics and Statistics")
c.setAuthor("University of Victoria; department label composed for poster")
c.saveState()
c.translate(0, 4.2)
c.scale(1.15, 1.15)
renderPDF.draw(crest, c, 0, 0)
c.restoreState()
c.setStrokeColor(white)
c.setLineWidth(1.2)
c.line(228, 6, 228, 76)
c.setFillColor(white)
c.setFont("Helvetica", 25)
c.drawString(241, 46, "Mathematics &")
c.drawString(241, 16, "Statistics")
c.showPage()
c.save()

font_name = "CrimsonPro-Bold"
pdfmetrics.registerFont(TTFont(font_name, str(ROOT / "fonts/CrimsonPro-Bold.ttf")))
wordmark = "Nathoo Group"
font_size = 44
width = pdfmetrics.stringWidth(wordmark, font_name, font_size) + 4
c = canvas.Canvas(str(FIGURES / "nathoo_wordmark_white.pdf"), pagesize=(width, 55))
c.setTitle("Nathoo Group wordmark")
c.setAuthor("Nathoo Group")
c.setFillColor(white)
c.setFont(font_name, font_size)
c.drawString(2, 12, wordmark)
c.showPage()
c.save()
print("Built UVic department and Nathoo Group vector wordmarks.")
