# Poster logos

**Status: no logo files are present in this repository.** The archive
described in the poster brief was not found anywhere on the machine this
branch was built on (`assets/` did not exist, and a filesystem search for
`*KAUST*`, `*Alagil*`, `*rzm*`, `*ntdp*` and similar returned nothing). The
files below therefore have to be supplied by the author before the poster can
be printed. Correctly sized placeholders are used in the layout specification
in the meantime, each carrying the logo name.

## The six logos

The programme brief names four required logos. Every poster in the previous
cohort's bank in fact carries six. Plan for six; the author confirms the final
set.

| # | Logo | Expected file | Present? |
|---|---|---|---|
| 1 | KAUST | `KAUST print logo-Black.png` (reversed: `KAUST print logo White-01.png`) | **no — author to supply** |
| 2 | KAUST Academy | `KAUST_Academy_logo_Full_Color.png` (black/white variants also exist) | **no — author to supply** |
| 3 | University of Cambridge | — | **no — not in the archive; author must source** |
| 4 | Department of Genetics | — | **no — not in the archive; author must source** |
| 5 | ntdp | `ntdp-logo-light.png` | **no — author to supply** |
| 6 | Alagil Foundation | `image_Alagilfoundationlogoupdated.png` (grey: `Alagail-gray.png`) | **no — author to supply** |

RZM (`rzm-invest-logo.png`) is named in the programme brief but does not appear
in the previous cohort's posters. Include it if the programme requires it; the
layout leaves room for a seventh mark in the top strip.

## Two logos are missing from the archive entirely

**University of Cambridge** and **Department of Genetics** were not in the
supplied archive and must be sourced by the author from the University's brand
resources. Do not substitute anything for them.

Two traps worth recording, both noted in the brief and worth repeating here
because both are easy to get wrong at speed:

- The archive contains logos for **other host universities** — Oxford,
  Toronto, Duke, Imperial, University College London. None of them belongs on
  this poster.
- The file named **`CREST`** is KAUST's Center for Renewable Energy and Storage
  Technologies. It is **not** the Cambridge shield. Do not use it as one.

## Rules for using them

- **Do not recolour, distort, crop, rotate or add effects to any logo.** Scale
  proportionally only.
- **Pick one treatment and use it for all six.** Either all full-colour or all
  monochrome — never a mixture. Where both a colour and a mono variant exist,
  the mono set is usually the safer choice for a single strip.
- **Never place a dark-on-transparent logo over the deep violet title band.**
  Use the white or reversed variant there, or move the strip off the band.
- Keep clear space around each mark of at least half its own height.
- Set all six to a **consistent optical height** (not a consistent bounding
  box): 34 mm in the layout specification. Wordmarks with descenders may need
  1–2 mm of adjustment by eye.

## Placement

House convention from the poster bank: a single row across the top edge,
above the title, at consistent height with clear space around each. The layout
specification in `results/reports/poster_final/POSTER_TEXT.md` places the strip
at `y 12–52 mm`, full width, 34 mm logo height. A split between a top row and
a bottom band is the accepted alternative if the top strip becomes crowded.

## When the files arrive

Drop them into this directory using the expected filenames above, then update
the "Present?" column. Nothing in the figure pipeline reads these files — they
are layout assets only — so no code needs to change.
