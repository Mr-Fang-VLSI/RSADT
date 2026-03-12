#!/usr/bin/env python3
from PIL import Image, ImageDraw, ImageFont


W, H = 1400, 1540
BG = (255, 255, 255)
GRID = (120, 128, 140)
BLUE = (39, 84, 142)
ORANGE = (233, 145, 0)
RED = (198, 40, 40)
GRAY = (110, 110, 110)
BLACK = (35, 35, 35)
GREEN = (46, 125, 50)


def font(size: int, bold: bool = False):
    candidates = [
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf" if bold else "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/local/share/fonts/DejaVuSans-Bold.ttf" if bold else "/usr/local/share/fonts/DejaVuSans.ttf",
    ]
    for path in candidates:
        try:
            return ImageFont.truetype(path, size=size)
        except OSError:
            continue
    return ImageFont.load_default()


IMG = Image.new("RGB", (W, H), BG)
D = ImageDraw.Draw(IMG)
TITLE = font(46, True)
SUB = font(30, True)
TEXT = font(24, False)
LABEL = font(28, True)
SMALL = font(20, False)


def center_text(x0, y0, x1, text, fnt, fill=BLACK):
    bbox = D.textbbox((0, 0), text, font=fnt)
    tw = bbox[2] - bbox[0]
    th = bbox[3] - bbox[1]
    D.text(((x0 + x1 - tw) / 2, y0), text, fill=fill, font=fnt)
    return th


def draw_grid(origin, cols=5, rows=4, cell=86):
    x, y = origin
    for c in range(cols + 1):
        D.line((x + c * cell, y, x + c * cell, y + rows * cell), fill=GRID, width=3)
    for r in range(rows + 1):
        D.line((x, y + r * cell, x + cols * cell, y + r * cell), fill=GRID, width=3)
    return cell


def box(origin, rc, label_text, fill, cell=86, pad=14, text_fill=(255, 255, 255)):
    x, y = origin
    r, c = rc
    x0 = x + c * cell + pad
    y0 = y + r * cell + pad
    x1 = x + (c + 1) * cell - pad
    y1 = y + (r + 1) * cell - pad
    D.rounded_rectangle((x0, y0, x1, y1), radius=8, fill=fill)
    bbox = D.textbbox((0, 0), label_text, font=LABEL)
    tw = bbox[2] - bbox[0]
    th = bbox[3] - bbox[1]
    D.text(((x0 + x1 - tw) / 2, (y0 + y1 - th) / 2 - 2), label_text, fill=text_fill, font=LABEL)
    return ((x0 + x1) / 2, (y0 + y1) / 2)


def arrow(p0, p1, color=RED, width=6, dash=False):
    if not dash:
        D.line((p0, p1), fill=color, width=width)
    else:
        x0, y0 = p0
        x1, y1 = p1
        steps = 18
        for i in range(0, steps, 2):
            xa = x0 + (x1 - x0) * i / steps
            ya = y0 + (y1 - y0) * i / steps
            xb = x0 + (x1 - x0) * (i + 1) / steps
            yb = y0 + (y1 - y0) * (i + 1) / steps
            D.line((xa, ya, xb, yb), fill=color, width=width)
    x0, y0 = p0
    x1, y1 = p1
    vx, vy = x1 - x0, y1 - y0
    mag = (vx * vx + vy * vy) ** 0.5 or 1
    ux, uy = vx / mag, vy / mag
    px, py = -uy, ux
    tip = (x1, y1)
    left = (x1 - 18 * ux + 10 * px, y1 - 18 * uy + 10 * py)
    right = (x1 - 18 * ux - 10 * px, y1 - 18 * uy - 10 * py)
    D.polygon([tip, left, right], fill=color)


def segment(p0, p1, color=RED, width=7):
    D.line((p0, p1), fill=color, width=width)


def callout(x, y, text, color):
    D.rounded_rectangle((x, y, x + 16, y + 16), radius=3, fill=color)
    D.text((x + 24, y - 4), text, font=SMALL, fill=BLACK)


D.text((90, 42), "Orange: swapped DSPs  |  Blue: neighboring DSPs on affected nets  |  Red: nets whose span enters the swap-local HPWL comparison", font=TEXT, fill=BLACK)
D.line((80, 88, W - 80, 88), fill=(190, 190, 190), width=2)

panels = [
    ("Case 1: A and B on the same row or column", "Affected local nets: (A,a4) and (B,b4)"),
    ("Case 2: A and B on different rows and columns", "Affected local nets: (A,a1), (A,a4), (B,b3), and (B,b4)"),
    ("Case 3: Boundary or corner degeneration", "Boundary case: one side is truncated, but the same local-net comparison remains"),
]

y_offsets = [120, 560, 1000]

for idx, ((title, subtitle), y0) in enumerate(zip(panels, y_offsets)):
    D.rounded_rectangle((70, y0, W - 70, y0 + 410), radius=18, outline=(205, 205, 205), width=3, fill=(252, 252, 252))
    D.text((100, y0 + 28), title, font=SUB, fill=BLACK)
    D.text((100, y0 + 74), subtitle, font=TEXT, fill=GRAY)
    origin = (120, y0 + 130)
    cell = draw_grid(origin)

    if idx == 0:
        pb4 = box(origin, (1, 0), "b4", BLUE)
        pA = box(origin, (1, 1), "A", ORANGE)
        pB = box(origin, (1, 3), "B", ORANGE)
        pa4 = box(origin, (1, 4), "a4", BLUE)
        box(origin, (0, 1), "", BLUE)
        box(origin, (2, 1), "", BLUE)
        arrow((pA[0] + 25, pA[1]), (pB[0] - 25, pB[1]), color=GRAY, width=5, dash=True)
        segment((pb4[0] + 24, pb4[1]), (pA[0] - 24, pA[1]))
        segment((pB[0] + 24, pB[1]), (pa4[0] - 24, pa4[1]))
        D.text((640, y0 + 165), "swap(A,B)", font=TEXT, fill=GRAY)
        D.text((640, y0 + 220), "Only the two boundary nets of the local\nordered window can change their span.", font=TEXT, fill=BLACK)

    elif idx == 1:
        pa1 = box(origin, (0, 1), "a1", BLUE)
        pA = box(origin, (1, 2), "A", ORANGE)
        pa4 = box(origin, (1, 4), "a4", BLUE)
        pb4 = box(origin, (2, 0), "b4", BLUE)
        pB = box(origin, (2, 2), "B", ORANGE)
        pb3 = box(origin, (2, 4), "b3", BLUE)
        arrow((pA[0] - 12, pA[1] + 18), (pB[0] - 12, pB[1] - 18), color=GRAY, width=5, dash=True)
        segment((pA[0], pA[1] - 24), (pa1[0], pa1[1] + 24))
        segment((pA[0] + 24, pA[1]), (pa4[0] - 24, pa4[1]))
        segment((pb4[0] + 24, pb4[1]), (pB[0] - 24, pB[1]))
        segment((pB[0] + 24, pB[1]), (pb3[0] - 24, pb3[1]))
        D.text((640, y0 + 165), "swap(A,B)", font=TEXT, fill=GRAY)
        D.text((640, y0 + 220), "The cross-pattern inequality compares the\nfour local nets incident on A and B.", font=TEXT, fill=BLACK)

    else:
        pA = box(origin, (0, 0), "A", ORANGE)
        pa4 = box(origin, (0, 2), "a4", BLUE)
        pB = box(origin, (2, 0), "B", ORANGE)
        pb3 = box(origin, (2, 2), "b3", BLUE)
        box(origin, (1, 1), "", BLUE)
        arrow((pA[0], pA[1] + 24), (pB[0], pB[1] - 24), color=GRAY, width=5, dash=True)
        segment((pA[0] + 24, pA[1]), (pa4[0] - 24, pa4[1]))
        segment((pB[0] + 24, pB[1]), (pb3[0] - 24, pb3[1]))
        D.text((640, y0 + 165), "swap(A,B)", font=TEXT, fill=GRAY)
        D.text((640, y0 + 220), "At a boundary or corner, missing neighbors\nremove terms but do not change the sign logic.", font=TEXT, fill=BLACK)

callout(90, H - 95, "Affected net segment", RED)
callout(420, H - 95, "Swap direction", GRAY)
callout(700, H - 95, "Neighbor in local proof window", BLUE)
callout(1105, H - 95, "Swapped DSP", ORANGE)

IMG.save("paper/iccad2026/assets/oc_cases_precise.png")
