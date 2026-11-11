# ftorralba.github.io

Personal website. Plain HTML and CSS, no build step: edit a file, commit, push.

- `*.html`: one file per section. Content goes inside `<main>`.
- `js/site.js`: the sidebar list and the permanent disclaimer. Edit `PAGES` there to add or reorder sections.
- `css/style.css`: all styling.
- `files/`: PDFs to link from pages (create it when the first one is added).
- `_archive/`: the previous site, kept for reference. Folders starting with `_` are not published.

To preview locally: `python -m http.server 8000` in this folder, then open http://localhost:8000.
