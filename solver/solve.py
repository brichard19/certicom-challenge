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
        "1614aeb8e92db814df57",
        "4027a11e368e956d92af",
        "5167d26cf01aca8bf67c",
        "732e3e470ae7638ed05",
        "5a1b3cb7c1c5d0631ab6",
        "23bb3fd83dbc9352d5e1",
        "3b499686e5018ffe3e2b",
        "2de349f117f79bf77dbc",
        "494bed2ef651d0fe4b6e",
        "3a3caadb82edb51b377b",
        "391d345ae163b3edb2a5",
        "374a4f58940039749c8f",
        "60ddf2e5033edc753ef9",
        "3663cf9cc6d126c292fb",
        "12c308adf7f92a59fe65",
        "79084e620ce46337b79",
        "5334c520d01bcca5f438",
        "271dfc84af3c4191bd53",
        "3bbb41fcba037324ef30",
        "30a6d7ea510a40bdbd36",
        "17fc76572f0e53b9b860",
        "627d40ad93ece4b4187",
        "291d2f8d459f90aaa910",
        "24a5d671bb553afde95a",
        "606447988cc6778d59d0",
        "11b85dc7f5d8d8fe0db9",
        "2e3ddc6562ff7f880357",
        "29072a23fad221d674b1",
        "3e6099679a21ab541c45",
        "430932e0c44731c97bc3",
        "30de088da3da65a8da14",
        "33c89818c30aa03c38bd",
        ]

        b_str = [
        "26ad00a3569c3785004e",
        "2320f5a64d245a6e48e1",
        "25318ac49ecc98ccb8e8",
        "4c9af6b038efb323b91e",
        "12bff8347f3f2844d9a9",
        "381b2920dfe5941c94c0",
        "415ab7e9fe8ba6914814",
        "94ee5c57b318e7e416c",
        "13dd9410cec72986530b",
        "5ef10216635657136ff9",
        "ddbf26d833d9c4879bd",
        "c9c6249ea691d445b12",
        "15557da3f36a91ebfc0d",
        "27f8e3c12af04271c534",
        "49d1535780bceb4fd6bf",
        "1e784e1c9252ac1e790c",
        "21215b6fef72e6328835",
        "32d94875ecb56b920dcf",
        "4ac05cb562e67f279fec",
        "43c291909a8afd751912",
        "17b829d768e142338e05",
        "677d94839904aad3720",
        "4fabf6b53722e7c46de3",
        "2d4404483feffaddc45c",
        "f3031b646c1245e52cc",
        "97d2be05cadebd9cb1d",
        "5e7195eefb021fbba33b",
        "1197fb72ac16016f55e0",
        "c14c2e53d858a989c0a",
        "169c545411f32d9b958",
        "2a0a3c3ec79b5120a83f",
        "1110024670d75ffcf4e6",
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

        # The walk was originally done in Montgomery form, convert so we can
        # get the correct bits
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

        p1 = self.curve.multiply(self.a1, self.curve.bp)
        b1 = 0

        p2 = self.curve.multiply(self.a2, self.curve.bp)
        b2 = 0

        # Advance the longer walk so that both walks have the same remaining iterations
        # left. 
        diff = self.len1 - self.len2
        for _ in range(diff):
            self.a1, b1, p1 = self.iterate(self.a1, b1, p1)
        
        # Advance both walks together until a collision
        while p1.x != p2.x:
            self.a1, b1, p1 = self.iterate(self.a1, b1, p1)
            self.a2, b2, p2 = self.iterate(self.a2, b2, p2)

        print('Colliding point:')
        print('{} {}'.format(hex(p1.x), hex(p1.y)))
        print('a1 = {} b1 = {}'.format(hex(self.a1), b1))
        print('a2 = {} b2 = {}'.format(hex(self.a2), b2))

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
 