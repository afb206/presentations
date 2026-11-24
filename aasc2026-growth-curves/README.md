# Report build (AASC2026 growth-curve poster companion)

Static site: open `index.html`, or publish this folder with GitHub Pages. No build step is needed to view it.

## Rebuild
1. `cd scripts && python3 common.py <poster sources folder>`  copies the input files and checks each md5 against the poster's MD5SUMS.txt
2. `python3 figures.py`   redraws the figures (figs/*.svg) from sources/
3. `python3 assemble.py`  fills template.html -> index.html, writes ../source_list.csv

Needs Python 3 with pandas, matplotlib and pyreadr.

## How it works
- `template.html` is the layout and the text. Numbers in the text are {{tokens}}; each is read from a saved file by assemble.py and listed (file, md5, how read) in the Sources section and in source_list.csv.
- A dashed "slot" in the page means the data for it has not been supplied yet. Slots still empty are listed at the end of the assemble.py output. A slot is replaced when assemble.py is given a fragment for its token (see FILLED in assemble.py).
- Nothing here is committed or published.
