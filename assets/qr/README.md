# QR code

**Status: no QR files are present in this repository, and none was generated.**

The poster brief states that a print-quality QR code "has already been
generated and is being supplied by the author" — version 6, high error
correction, SVG plus a 1170 × 1170 px PNG, in violet `#6A3D9A` and near-black
`#262626`. Those files were not found anywhere on this machine.

The brief's fallback is to "generate an equivalent using a library already
present in the environment; do not install a new dependency without saying
so." **No QR library is present**: `qrcode`, `segno` and `pyqrcode` are all
absent from the `bc` environment. Nothing was installed, and no QR code was
hand-rolled — an encoder written from scratch and not verified against a real
scanner is worse than no QR code at all, because a poster carrying an
unscannable code fails silently in front of an audience.

So this directory holds the specification and nothing else.

## What the author should supply

| Property | Value |
|---|---|
| Target URL | `https://github.com/almokhtar8-stack/breast_cancer` |
| Version | 6 |
| Error correction | H (high, 30%) |
| Formats | SVG (vector, preferred for print) and PNG at 1170 × 1170 px |
| Foreground, dark background | violet `#6A3D9A` |
| Foreground, light background | near-black `#262626` |
| Quiet zone | 4 modules on all sides, unmodified |
| Printed size | **≥ 40 × 40 mm**; the layout allocates 40 mm |

Expected filenames: `qr_violet.svg`, `qr_violet.png`, `qr_black.svg`,
`qr_black.png`.

## If they need regenerating

With a QR library available (installing one is a decision for the author, not
for this branch):

```bash
python -m pip install segno          # NOT installed by this branch
python - <<'PY'
import segno
qr = segno.make("https://github.com/almokhtar8-stack/breast_cancer",
                error="h", version=6)
qr.save("assets/qr/qr_violet.svg", scale=10, dark="#6A3D9A", light=None)
qr.save("assets/qr/qr_violet.png", scale=30, dark="#6A3D9A", light="#FFFFFF")
qr.save("assets/qr/qr_black.svg",  scale=10, dark="#262626", light=None)
qr.save("assets/qr/qr_black.png",  scale=30, dark="#262626", light="#FFFFFF")
PY
```

**Scan the printed proof before the poster is printed at full size.** A QR
code that resolves on screen can still fail on paper at the wrong contrast or
size.

## Placement

At least 40 × 40 mm at final print size, with the caption **"Code, data and
figures"** beside it and the repository URL beneath in text as well — so the
link survives a phone that cannot scan. A violet frame if it sits on a light
background. The layout specification places it at the foot of the poster with
the contact details.

## What the QR points at

The repository root. A visitor arriving there sees the root `README.md` and
`POSTER.md`, both of which were rewritten on this branch to be accurate as of
this state — what the project is, what it found *including the negative
results*, what it does not claim, and where the figures, analysis plan and
freeze tag are.
