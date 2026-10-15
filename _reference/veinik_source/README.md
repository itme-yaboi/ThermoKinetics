# Veinik source material (reference only; only this README is tracked in git)

Raw sources behind `docs/veinik/`. **Read the curated docs first.** Come here only to check what Veinik actually wrote, or to verify an equation.

| File | What it is | Quality |
|---|---|---|
| `1964_book_front_matter_and_ch1-3_ocr.txt` | 1964 English edition: title, preface, contents, introduction, Chapters I–III (pp. iii–37, eqs. (1)–(25)) | Good prose OCR; equations sometimes garbled |
| `1964_book_ch4-end_ocr.txt` | 1964 edition, Chapter IV through the Conclusion (pp. 37–254, eqs. (26)–(522)), from the original scans | Prose fine; many equations garbled; pp. 244–245 missing here |
| `1964_book_reimaged_pages_ocr.txt` | Clearer re-images of pp. 26–31, 36–45, 50–57, 64–69, 74–79 and the missing pp. 244–247 | Better than the scans for §10–12 and §62 |
| `page_images/ch1-3/`, `page_images/reimaged/` | Screenshots of two-page spreads (file names give page numbers) | **The authority for equations.** Open the image whenever OCR looks wrong |
| `notes_by_section_ch6-18.md` | Section-by-section notes with equation numbers (Chapters VI–XVIII and the Conclusion) | Hand-checked digest |
| `1966_new_thermodynamics_RU.txt` | Veinik 1966, *New Thermodynamics* (Russian) | Only parts are sound; see `docs/veinik/05` |
| `1973_thermodynamic_couples_RU.txt` | Veinik 1973, *Thermodynamic Couples* (Russian) | Theory chapters are wrong; use only the couple descriptions listed in `docs/veinik/05` |

Rules:
- Search with grep by equation number, e.g. `grep -n "(79)"`, or by section, e.g. `grep -n "§ 46"`.
- Units in the sources are kgf, kcal, kgf·m, and γ is specific weight. Convert to SI.
- If a source equation conflicts with a non-negotiable in CLAUDE.md, the non-negotiable wins; log it in `docs/DECISIONS.md`.
- Never use Veinik's later "chronal" writings (not included here).
- This is a scanned copyrighted book for personal reference. Keep everything here except this README out of git (D-25, D-36).
