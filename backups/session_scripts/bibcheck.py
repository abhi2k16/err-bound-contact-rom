import re
base = r"C:\Users\abhij\Desktop\Physical_ML\proposed-research\reduced-contact-reliability"
tex = open(base + r"\manuscript_master_v2.tex", encoding="utf-8").read()
bib = open(base + r"\references.bib", encoding="utf-8").read()
cited = set()
for m in re.finditer(r"\\cite[a-z]*\{([^}]*)\}", tex):
    cited |= {k.strip() for k in m.group(1).split(",")}
entries = re.findall(r"@(\w+)\{([^,]+),(.*?)\n\}", bib, flags=re.S)
keys = {k for _, k, _ in entries}
print("cited not in bib:", sorted(cited - keys))
print("bib not cited:", sorted(keys - cited))
for typ, k, body in entries:
    def f(n):
        m = re.search(n + r"\s*=\s*[{\"](.*?)[}\"],?\s*\n", body, flags=re.S)
        return m.group(1).replace("\n", " ")[:70] if m else ""
    print(f"{k:26s}{typ:8s} yr={f('year'):5s} doi={'Y' if 'doi' in body else 'N'} | {f('author')[:38]} | {f('journal') or f('booktitle') or f('publisher') or f('howpublished') or f('note')}")
