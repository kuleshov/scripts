def n50(l):
    total_length = sum(l)
    half_length = total_length / 2
    length = 0

    l.sort(reverse=True)

    while length < half_length:
        x = l.pop()
        length += x
        # print "%d/%d" % (length, total_length)

    return x

def nx(L,x):
    assert 0.0 <= x <= 1.0
    assert L
    l=list(L)
    total_length = float(sum(l))
    x_length = total_length * (1.0-x)
    length = 0

    l.sort(reverse=True)

    while length < x_length:
        t = l.pop()
        length += t
        # print "%d/%d" % (length, total_length)

    return t


