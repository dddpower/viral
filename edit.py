with open("hiv.csv", "rw") as f:
    a = f.read()
    s = "pos,   Codon"
    enumerate(a)
    for i, ch in enumerate(a):
        s += str(i)
        s += ch
        if i % 2 == 2:
            s += "\n"
    f.write(s)
