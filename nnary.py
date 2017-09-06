
def nary(alist, n):
    results = alist[len(alist)/n/2::len(alist)/n]
    return results

alist = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]+range(1,20)
print len(alist)

r = nary(alist, 3)

print r
