import ecc
import curves
import sys

# Certicom challenge points
QPOINTS = {
    "ecp79":ecc.ECPoint(0x0679834CEFB7215DC365, 0x4084BC50388C4E6FDFAB)
}

'''
Inversion mod p using Fermat method
'''
def invm(x, m):
    y = x % m
    return pow(y, m-2, m)

def from_montgomery(x, p):
    r = 2**160

    rinv = invm(r, p)

    return (x * rinv) % p

def to_montgomery(x, p):
    r = 2**160

    return x * r % p

def get_r_points(name):

    curve = curves.getCurveByName(name)
    if name == 'ecp79':

        a_str = [
        "c806168f60b90d7c758",
        "1984517914b97ccc1364",
        "2deed643606e2a22002c",
        "af0f8b9174d05c07ff7",
        "52248cab7d2291be77c5",
        "1eb1936a6303310fac6",
        "3306f509824ae973bf13",
        "20f62dcd7a058e7db9be",
        "120f8c2c38be74942c23",
        "338ae6f3737a41124f74",
        "187af490e26e243b3342",
        "205a50c11193447299e8",
        "5bc57ded78e9916896b9",
        "251b8826f8d9148bbda0",
        "5b90888d2b82e1721d4f",
        "333a7f456008f1b48773",
        "3e8782c66d7308c3e084",
        "54acb0faaec3c08319ad",
        "359843717ee3e17ebf0b",
        "4af891f0c1e7efbe66b8",
        "2e895547a16c1087c764",
        "41b522d8cda0ce5b924e",
        "430d751c31f1009d6faf",
        "27c9f8c726d3c43f2138",
        "5cc55299d57f2a654270",
        "49ec897ac98154b79d6c",
        "651c04fac912c07c709",
        "4487e2b97379c31598a6",
        "48a3bf746dc70b4e662a",
        "254577444d1878c51498",
        "3e0fdfc769e22a80ed7c",
        "3d06c06c898a002c96ea",
        ]

        b_str = [
        "23573b8c0594b154a7fc",
        "1211e94deb5d9f959750",
        "5cadcdce5e6cac5dd968",
        "600e7edc5eb8a2af3ce5",
        "2fc433085787f95c8ccd",
        "11e05f6e995f65e4a401",
        "318e67b52daa9b357062",
        "2b12f4fb6f1b65b77469",
        "38b864a4a65da340606a",
        "1e0b684473b6df8556b4",
        "470b7bed297003c91919",
        "226c78e7dd14d0cd3096",
        "428a4e1fa65f3e35fbc2",
        "1894e84e4e962cdfbb42",
        "1859387671918e78902a",
        "59399b32e3c2e7afda6a",
        "41df7b3657c29c5e26e3",
        "204ddc5ee834d5e81b3d",
        "59c8a416e159b19f51d7",
        "30353613c9028867a174",
        "45a1ea9f070ae65d7f9a",
        "a861f66444221e79dd4",
        "3932db73ffdff4526882",
        "ec3ac32591ae83e6558",
        "69a23a39ef0de86e928",
        "250b525ba49a1bc18810",
        "108e6248774e4309b24f",
        "2845de39f3531db60878",
        "1f0ab0eb7640e52e89df",
        "4d8d32889a552bf91898",
        "500a3f277a391ef173ed",
        "3162535bf0766f03758b"
        ]

        q = QPOINTS["ecp79"]
        rwpoints = []
        for i in range(32):
            a = int(a_str[i], 16)
            b = int(b_str[i], 16)

            ap = curve.multiply(a, curve.bp)
            bq = curve.multiply(b, q)

            r = curve.add(ap, bq)
            if not curve.verifyPoint(r):
                raise "Invalid point"

            rwpoints.append({'a':a, 'b':b, 'r':r})
    else:
        raise "Invalid curve name"
    
    return rwpoints


class RhoSolver:

    def __init__(self, curve, rwpoints, a1, a2, len1, len2, dp):
        self.curve = curve
        self.rwpoints = rwpoints

        self.dp = dp
        self.a1 = a1
        self.a2 = a2
        self.len1 = len1
        self.len2 = len2


    def iterate(self, a, b, point):

        # Convert to montgomery form so we can get the correct
        # distinguished bits
        xm = to_montgomery(point.x, self.curve.p)

        idx = xm % 32

        new_a = (a + self.rwpoints[idx]['a']) % self.curve.n
        new_b = (b + self.rwpoints[idx]['b']) % self.curve.n

        new_p = self.curve.add(self.rwpoints[idx]['r'], point)
        if not self.curve.verifyPoint(new_p):
            raise "Invalid point"
        return new_a, new_b, new_p

    def solve(self):

        # Ensure the first walk is longer
        if self.len1 < self.len2:
            self.len1, self.len2 = self.len2, self.len1
            self.a1, self.a2 = self.a2, self.a1

        # The starting points only have one exponent, a. b = 0
        p1 = self.curve.multiply(self.a1, self.curve.bp)
        b1 = 0

        p2 = self.curve.multiply(self.a2, self.curve.bp)
        b2 = 0

        print(f"x={p1.x:033x}    x={p2.x:033x}")
        print(f"y={p1.y:033x}    y={p2.y:033x}")

        # Advance the longer walk until boths walks have equal remaining steps
        diff = self.len1 - self.len2
        for _ in range(diff):
            self.a1, b1, p1 = self.iterate(self.a1, b1, p1)

        print('Starting points:')
        print(f"x={p1.x:033x}    x={p2.x:033x}")
        print(f"y={p1.y:033x}    y={p2.y:033x}")
        
        # Advance both walks together until a collision
        remaining = self.len2
        while p1.x != p2.x and remaining > 0:
            self.a1, b1, p1 = self.iterate(self.a1, b1, p1)
            self.a2, b2, p2 = self.iterate(self.a2, b2, p2)
            remaining -= 1

        print('Colliding point:')
        print(f"x={p1.x:033x}    x={p2.x:033x}")
        print(f"y={p1.y:033x}    y={p2.y:033x}")
        print(f"a={self.a1:033x}    a={self.a2:033x}")
        print(f"b={b1:033x}    b={b2:033x}")

        # Calculate the discrete logarithm
        # a1G + b1Q = a2G + b2Q
        # (a1 - a2)G = (b2 - b1)Q
        # (a1 - a2)G = (b2 - b1)kG
        #(a1 - a2) = (b2 - b1)k
        #(a1 - a2) / (b2 - b1) = k

        d1 = (self.a1 - self.a2) % self.curve.n
        d2 = (b2 - b1) % self.curve.n
        k = (d1 * invm(d2, self.curve.n)) % self.curve.n

        # Validate
        result_q = self.curve.multiply(k, self.curve.bp)
        real_q = QPOINTS["ecp79"]

        if result_q != real_q:
            print("Something went wrong")
            print('Expected {} {}'.format(hex(real_q.x), hex(real_q.y)))
            print('Got      {} {}'.format(hex(result_q.x), hex(result_q.y)))
        else:
            print('k={}'.format(hex(k)))



def usage():
    print('Usage:')
    print('solve.py [input file]')

def main():

    curve = curves.getCurveByName("ecp79")
    rwpoints = get_r_points("ecp79")

    if len(sys.argv) == 1:
        usage()
        return 1

    filepath = sys.argv[1]

    with open(filepath, 'rt') as f:
        line = f.readline()

        lines = line.split(' ')

        x = int(lines[0], 16)
        y = int(lines[1], 16)
        a1 = int(lines[2], 16)
        len1 = int(lines[3])

        a2 = int(lines[4], 16)
        len2 = int(lines[5])

        # Convert from Montgomery form
        x = from_montgomery(x, curve.p)
        y = from_montgomery(y, curve.p)

        dp = ecc.ECPoint(x, y)
        if not curve.verifyPoint(dp):
            raise "Invalid point"

    solver = RhoSolver(curve, rwpoints, a1, a2, len1, len2, dp)
    solver.solve()

if __name__ == "__main__":
    main()
 