#!/usr/bin/python

def sc2d(n):
    print 'default ', 37.0 * n**2 / (9 * n**2)
    print 'shared X', 31.0 * n**2 / (9 * n**2)
    print 'two pass', (10.0 * (n-1)**2 + (4*n-4) + 19 * (n-1)**2) / (9.0*(n-1)*(n-1))
    print 'one pass', (9.0 * n**2 + 9 * (n-1)**2) / (9 * (n-1)**2)

def sc3d(n):
    print 'default ', 77.0 * n**2 / (19 * n**2)
    print 'two pass', (20.0*2 * (n-1)**2 + (4*n-4) + 19*(n-1)**2) / (19 * (n-1)**2)
    print 'one pass', (19.0 * n**2 + 19 * (n-1)**2) / (19 * (n-1)**2)
    print '2 one pass', (19.0*2 * n**2 + 19*2*(n-1)**2) / (19 * 2 * (n-1)**2)

print '* sc2d'
for x in [16, 24, 32, 40, 48, 56]:
    print '**', x
    sc2d(x)
    print

print '* sc3d'
for x in [16, 24, 32]:
    print '**', x
    sc3d(x)
    print
