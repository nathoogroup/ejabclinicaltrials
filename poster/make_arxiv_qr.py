"""Build a separate, directly scannable link to the original eJAB paper."""
from pathlib import Path

import qrcode
from qrcode.image.svg import SvgPathImage

ROOT = Path(__file__).resolve().parent
URL = "https://arxiv.org/abs/2510.10358"
qr = qrcode.QRCode(error_correction=qrcode.constants.ERROR_CORRECT_M,
                   box_size=12, border=4)
qr.add_data(URL)
qr.make(fit=True)
qr.make_image(fill_color="black", back_color="white").save(
    ROOT / "figures" / "arxiv_qr.png")
qr.make_image(image_factory=SvgPathImage).save(ROOT / "figures" / "arxiv_qr.svg")
print(f"arXiv QR: {URL}; version {qr.version}; {qr.modules_count} modules")
