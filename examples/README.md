# Fabricated event format examples

`events.csv` and `events.tsv` are hand-written format examples. They contain no participant data and do not correspond to any supplied recording.

TIME and DURATION are seconds relative to the beginning of the matching BDF. LABEL maps pass1 to active, pass2 to passive, rest1/resting to resting. Use `--two-resting` to preserve rest2 as resting2. Every requested condition must have a matching row, and all intervals must fit inside the user's recording.

Use `--events-dir`, `--bdf-pattern`, `--events-pattern`, and `--subjects` to describe your files. The script does not silently truncate out-of-bounds intervals or overwrite existing outputs.
