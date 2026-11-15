import re
t = open(r"C:\Users\abhij\Desktop\Physical_ML\proposed-research\reduced-contact-reliability\manuscript_master_v2.tex", encoding="utf-8").read().replace("\r\n", "\n")
# restrict to the main body prose: drop table/algorithm/equation environments
body = t[t.index(r"\section{Introduction}"):t.index(r"\bibliographystyle")]
body = re.sub(r"\\begin\{(equation|align|aligned|algorithm|longtable|tabular|table|figure)\}.*?\\end\{\1\}", " ", body, flags=re.S)
rows = []
for para in body.split("\n\n"):
    flat = " ".join(para.split())
    for sent in re.split(r"(?<=[a-z0-9\)\}\$])\. (?=[A-Z\\])", flat):
        n = len(sent.split())
        if n > 60:
            rows.append((n, sent[:150]))
rows.sort(reverse=True)
print(len(rows), "sentences over 60 words")
for n, s in rows[:25]:
    print(n, "|", s)
