import ecc
import curves
import sys

# Certicom challenge points
QPOINTS = {
    "ecp79":ecc.ECPoint(0x0679834CEFB7215DC365, 0x4084BC50388C4E6FDFAB),
    "ecp89":ecc.ECPoint(0x00DE1AA94FF94DB64E763E2D, 0x002A44C4C2D4EE27FA0A4BA9)
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
    elif name == "ecp89":
        a_str = [
        "1a4a0dada8ff9b3353dbf2",
        "90944f2e1223ada2fe6373",
        "df4a0580717e733f6509ee",
        "eccf0e1eecc67e7d2d51ba",
        "34a7cecd5622ec1fccb805",
        "3db9106faef4f211b86f46",
        "206d628fd82f3f8cf7c08c",
        "2b721e761ae96179ffcc52",
        "10273929637074fab430b93",
        "82157440614fb7efc1f7a8",
        "70cbe069aed013279d9ec6",
        "3e584b772af8b8cb4e643e",
        "437535a44320d27861f68f",
        "18d37c5dc2837b351f98b6",
        "12f22ef037079792dd1a615",
        "e0257fe8b74b2faf11bea7",
        "10c012a09941aabf6c08359",
        "3720c8660fd83687484cbb",
        "69e5c7153e1d5c25bfb59",
        "124a06e492da5faf7693ee9",
        "1b922194e8550be06bf5b",
        "dd47cf5735c318aee15b40",
        "10701daa485386eadb2a892",
        "3a602074819aa0f2104bb7",
        "11c6ab1b0958bf9c38180f0",
        "32f38636e9e9ab27fb67f6",
        "10d0ec7513fc759ead97aa7",
        "711464d1088c61f8f18b52",
        "a1cd764b84bb55914fe2c2",
        "4658eb940363aa67c2b211",
        "cd1e7e6cda0eadc0a053e5",
        "ddc4e005bd535fb2f36283",
        ]

        b_str = [
        "a7f6a7f2e76f256c095115",
        "152691186eec9a4d1c336",
        "9ffe439631edbc7564e0b8",
        "11667595e7ed7b81ec2e68d",
        "b211f196f744455e8e1386",
        "e3d20c11bf0de746353e0c",
        "bd8ab6cc0f2bad3fd554ad",
        "14bde4d32d0e83e3b9d7c58",
        "784d8489b3d7c6303faa4f",
        "bb37cc89fbae126905ddc9",
        "98859bdd38b1b418dd12e4",
        "a90491a0351165cf876c5c",
        "fce0bece5e75d66e80f8c3",
        "ff39d4d8fae4ebbf6f7733",
        "b127c7dad6245f6d9f4e2",
        "1527ebed54230c960c860d3",
        "1541e14ca34b027b0ca76b4",
        "4b6c4ba15a6024c7190cd0",
        "b97c64b5c0afcb08d99385",
        "a0f3786fe87887935f9405",
        "836263f1f57926211949d1",
        "2d1d6f6c49a3065f9d58f8",
        "d6195039176c2f580fd817",
        "95135bf1863af1d5061bce",
        "ebe66e88f9bc0611bae27",
        "14d07436aab9d12a8e4a552",
        "118e2ede252b7f4aec373aa",
        "50144ab259be479b8a2ab3",
        "e3f222243072f9dbca29a4",
        "9900a0e7aef2842f6af0fb",
        "783c50fcbe09e5274b1ff8",
        "846c000769a327ee944dc6",
        ]
        q = QPOINTS["ecp89"]
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

    #curve = curves.getCurveByName("ecp89")
    #rwpoints = get_r_points("ecp89")

    if len(sys.argv) == 1:
        usage()
        return 1

    filepath = sys.argv[1]

    with open(filepath, 'rt') as f:

        # First line is curve name
        curve_name = f.readline().strip()
        curve = curves.getCurveByName(curve_name)
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

    rwpoints = get_r_points(curve_name)

    solver = RhoSolver(curve, rwpoints, a1, a2, len1, len2, dp)
    solver.solve()

if __name__ == "__main__":
    main()
 