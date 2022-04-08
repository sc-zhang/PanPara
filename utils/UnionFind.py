class UnionFind:
    def __init__(self, n):
        self._f = [_ for _ in range(n)]
    
    def find(self, x):
        t = x
        while t != self._f[t]:
            t = self._f[t]
        while x != self._f[x]:
            tmp = self._f[x]
            self._f[x] = t
            x = tmp
        return self._f[x]
    
    def union(self, x, y):
        fx = self.find(x)
        fy = self.find(y)
        if fx != fy:
            self._f[fy] = fx
    
        
