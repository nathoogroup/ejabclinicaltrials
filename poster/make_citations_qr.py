"""Create a self-contained citations QR in PNG and vector SVG formats."""
from pathlib import Path
import re

import qrcode
from qrcode.image.svg import SvgPathFillImage

ROOT = Path(__file__).resolve().parent
payload = (ROOT / "citations.txt").read_text(encoding="utf-8").strip()
qr = qrcode.QRCode(
    error_correction=qrcode.constants.ERROR_CORRECT_M,
    box_size=12,
    border=4,
)
qr.add_data(payload)
qr.make(fit=True)
qr.make_image(fill_color="black", back_color="white").save(
    ROOT / "figures" / "citations_qr.png"
)
qr.make_image(image_factory=SvgPathFillImage).save(
    ROOT / "figures" / "citations_qr.svg"
)
print(f"QR version {qr.version}; {qr.modules_count} modules per side")
reference_count = len(re.findall(r"^\d+\. ", payload, flags=re.MULTILINE))
print(f"Payload: {len(payload.encode('utf-8'))} bytes; {reference_count} references")
